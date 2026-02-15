import streamlit as st
import io
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio import Restriction  # インポートエラー回避のため一括インポート
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from dna_features_viewer import BiopythonTranslator

# 自作モジュールのインポート
from data.constants import CULTURE_DB
from modules.analyzer import analyze_vector_resistance
from modules.designer import design_gibson_primers, check_dimer
from modules.exporter import create_pdf

# --- UI 基本設定 ---
st.set_page_config(page_title="BioDesigner v1.3", layout="wide", initial_sidebar_state="expanded")
st.title("🧬 Cloning & Expression Designer")
st.markdown("ベクター解析からプロトコル作成までを自動化します。")

# --- 1. サイドバー: 宿主の設定 ---
st.sidebar.header("1. Host Selection")
host_choice = st.sidebar.selectbox("宿主を選択", list(CULTURE_DB.keys()))

# --- 2. メインパネル: ファイルアップロード ---
st.header("📂 Step 1: Vector Upload")
uploaded_file = st.file_uploader("ベクターファイル (.gb, .fasta) をアップロード", type=["gb", "fasta"])

if uploaded_file:
    # 配列データの読み込み
    content = uploaded_file.getvalue().decode("utf-8")
    handle = io.StringIO(content)
    file_format = "genbank" if uploaded_file.name.endswith(".gb") else "fasta"
    record = SeqIO.read(handle, file_format)
    
    st.success(f"解析成功: {record.id} ({len(record.seq)} bp)")
    st.divider()

    # --- 3. セクション1: クローニング設計 ---
    st.header("✂️ Step 2: Cloning Design")
    col1, col2 = st.columns([1, 1])
    
    with col1:
        st.subheader("制限酵素サイトの選択")
        
        # 安定した制限酵素検索 (Commandoセットを使用)
        try:
            search_results = Restriction.Commando.search(record.seq)
        except AttributeError:
            # 万が一のフォールバック
            search_results = Restriction.AllEnzymes.search(record.seq)
            
        # 1箇所だけ切る(Unique Cutter)酵素のリストを作成
        unique_cutters = sorted([str(enz) for enz, sites in search_results.items() if len(sites) == 1])
        
        if unique_cutters:
            selected_ez_name = st.selectbox("Unique Cutter を選択", unique_cutters)
            
            # 選択された酵素のオブジェクトを取得して切断位置(0-based)を特定
            ez_obj = getattr(Restriction, selected_ez_name)
            # Biopythonの結果(1-based)を、スライス用の0-basedに変換
            cut_pos = search_results[ez_obj][0] - 1 
            st.info(f"📍 {selected_ez_name} の切断位置: {cut_pos} bp 目")
        else:
            st.error("利用可能なUnique Cutterが見つかりませんでした。")
            st.stop()

    with col2:
        st.subheader("インサート配列の入力")
        ins_raw = st.text_area("挿入する遺伝子配列(ATGC)を入力", placeholder="ATGGT...", height=150)
    
    if ins_raw and len(ins_raw) >= 40:
        # 空白や改行を除去
        ins_seq = Seq(ins_raw.strip().replace("\n", "").replace(" ", "").upper())
        
        # プライマー設計
        primers = design_gibson_primers(record.seq, ins_seq, cut_pos)
        
        st.subheader("✅ 設計されたプライマー")
        p_col1, p_col2 = st.columns(2)
        with p_col1:
            st.write("**Forward (5'->3')**")
            st.code(str(primers["Forward"]["seq"]))
            st.caption(f"Tm: {primers['Forward']['Tm']}°C / {check_dimer(primers['Forward']['seq'], primers['Forward']['seq'])[1]}")
        with p_col2:
            st.write("**Reverse (5'->3')**")
            st.code(str(primers["Reverse"]["seq"]))
            st.caption(f"Tm: {primers['Reverse']['Tm']}°C / {check_dimer(primers['Reverse']['seq'], primers['Reverse']['seq'])[1]}")

        # --- 4. セクション2: 仮想プラスミド図 ---
        st.divider()
        st.header("🗺️ Step 3: Final Construct Map")
        # ベクターを分割してインサートを挿入
        final_seq = record.seq[:cut_pos] + ins_seq + record.seq[cut_pos:]
        final_rec = SeqRecord(final_seq, id="Construct")
        
        # フィーチャー（視覚パーツ）の追加
        final_rec.features.append(SeqFeature(FeatureLocation(0, cut_pos), type="misc_feature", qualifiers={"label": ["Vector_Up"]}))
        final_rec.features.append(SeqFeature(FeatureLocation(cut_pos, cut_pos + len(ins_seq)), type="CDS", qualifiers={"label": ["INSERT"], "color": ["#ff4b4b"]}))
        final_rec.features.append(SeqFeature(FeatureLocation(cut_pos + len(ins_seq), len(final_seq)), type="misc_feature", qualifiers={"label": ["Vector_Down"]}))
        
        fig, ax = plt.subplots(figsize=(10, 2))
        translator = BiopythonTranslator()
        graphic_record = translator.translate_record(final_rec)
        graphic_record.plot(ax=ax, with_ruler=True)
        st.pyplot(fig)

    # --- 5. セクション3: 実験プロトコル ---
    st.divider()
    st.header("📋 Step 4: Experimental Protocol")
    
    # 耐性遺伝子の特定
    res_genes = analyze_vector_resistance(record)
    h_data = CULTURE_DB.get(host_choice, {})
    
    # KeyErrorを防ぐために .get() で安全にデータを取得
    protocol_disp = {
        "Transformation": h_data.get("trans_method", "データなし"),
        "Media": h_data.get("media", "データなし"),
        "Incubation": h_data.get("incubation", "データなし"),
        "Selection": ", ".join([h_data.get("antibiotics", {}).get(g, "不明") for g in res_genes]) if res_genes else "None detected"
    }
    
    st.json(protocol_disp)
    
    # インサート入力がある場合のみ、PDFダウンロードボタンを表示
    if ins_raw and 'primers' in locals():
        pdf_data = create_pdf(
            host_choice, 
            selected_ez_name, 
            str(primers["Forward"]["seq"]), 
            str(primers["Reverse"]["seq"]), 
            protocol_disp
        )
        st.download_button(
            label="📄 プロトコルをPDFでダウンロード", 
            data=bytes(pdf_bytes if 'pdf_bytes' in locals() else pdf_data), 
            file_name="cloning_protocol.pdf",
            mime="application/pdf"
        )