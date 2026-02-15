import streamlit as st
import io
import matplotlib.pyplot as plt
from Bio import SeqIO, Restriction
from Bio.Restriction import Commando
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
st.set_page_config(page_title="BioDesigner v1.1", layout="wide", initial_sidebar_state="expanded")
st.title("🧬 Cloning & Expression Designer")
st.markdown("ベクター解析からプロトコル作成までを自動化します。")

# --- 1. サイドバー: 宿主の設定 ---
st.sidebar.header("1. 宿主の選択")
host_choice = st.sidebar.selectbox("対象となる宿主細胞", list(CULTURE_DB.keys()))

if host_choice == "哺乳類細胞 (Mammalian)":
    cell_line = st.sidebar.selectbox("細胞株 (推奨条件表示用)", ["HEK293T", "CHO", "HeLa", "Other"])

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
        # インポートエラー回避のための安定した解析ロジック
        all_results = Commando.search(record.seq)
        unique_cutters = sorted([str(enz) for enz, sites in all_results.items() if len(sites) == 1])
        
        if unique_cutters:
            selected_ez_name = st.selectbox("Unique Cutter (1箇所のみ切断) を選択", unique_cutters)
            
            # 選択された酵素オブジェクトを取得して切断位置を特定
            ez_obj = getattr(Restriction, selected_ez_name)
            cut_pos = all_results[ez_obj][0] # 1-based index
            st.info(f"📍 {selected_ez_name} の切断位置: {cut_pos} bp")
        else:
            st.error("利用可能なUnique Cutterが見つかりませんでした。")
            st.stop()

    with col2:
        st.subheader("インサート配列の入力")
        ins_raw = st.text_area("挿入する遺伝子配列を入力 (ATGC)", placeholder="ATGGTGAGCA...", height=150)
    
    if ins_raw and len(ins_raw) >= 40:
        # 配列のクリーニング
        ins_seq = Seq(ins_raw.strip().replace("\n", "").replace(" ", "").upper())
        
        # プライマー設計 (15bp Overhang)
        primers = design_gibson_primers(record.seq, ins_seq, cut_pos)
        
        st.subheader("✅ 設計されたプライマー (Gibson/In-Fusion用)")
        p_col1, p_col2 = st.columns(2)
        
        with p_col1:
            st.write("**Forward Primer (5'->3')**")
            st.code(str(primers["Forward"]["seq"]))
            is_dimer, msg = check_dimer(primers["Forward"]["seq"], primers["Forward"]["seq"])
            st.caption(f"Tm: {primers['Forward']['Tm']}°C / {msg}")
            
        with p_col2:
            st.write("**Reverse Primer (5'->3')**")
            st.code(str(primers["Reverse"]["seq"]))
            is_dimer, msg = check_dimer(primers["Reverse"]["seq"], primers["Reverse"]["seq"])
            st.caption(f"Tm: {primers['Reverse']['Tm']}°C / {msg}")

        # --- 4. セクション2: 仮想プラスミド図 ---
        st.divider()
        st.header("🗺️ Step 3: Final Construct Map")
        
        # 配列の結合
        final_seq = record.seq[:cut_pos] + ins_seq + record.seq[cut_pos:]
        final_rec = SeqRecord(final_seq, id="Construct", name="Construct")
        
        # 可視化用フィーチャーの追加
        final_rec.features.append(SeqFeature(FeatureLocation(0, cut_pos), type="misc_feature", qualifiers={"label": ["Vector_Upstream"], "color": ["#cff0ff"]}))
        final_rec.features.append(SeqFeature(FeatureLocation(cut_pos, cut_pos + len(ins_seq)), type="CDS", qualifiers={"label": ["INSERT_GENE"], "color": ["#ff4b4b"]}))
        final_rec.features.append(SeqFeature(FeatureLocation(cut_pos + len(ins_seq), len(final_seq)), type="misc_feature", qualifiers={"label": ["Vector_Downstream"], "color": ["#cff0ff"]}))
        
        # 図の描画
        fig, ax = plt.subplots(figsize=(12, 2))
        translator = BiopythonTranslator()
        graphic_record = translator.translate_record(final_rec)
        graphic_record.plot(ax=ax, with_ruler=True)
        st.pyplot(fig)
        st.caption("※リニア（線形）表示での構築確認用マップです。")

        # --- 5. セクション3: 実験プロトコル ---
        st.divider()
        st.header("📋 Step 4: Experimental Protocol")
        
        res_genes = analyze_vector_resistance(record)
        h_data = CULTURE_DB[host_choice]
        
        # 薬剤情報の紐付け
        selected_abs = []
        for g in res_genes:
            ab_name = h_data["antibiotics"].get(g)
            if ab_name:
                selected_abs.append(f"{g.upper()}: {ab_name}")

        c1, c2, c3 = st.columns(3)
        with c1:
            st.info("**💡 導入条件**")
            st.write(f"- 手法: {h_data['trans_method']}")
        with c2:
            st.success("**🧫 培養条件**")
            st.write(f"- 培地: {h_data['media']}")
            st.write(f"- 期間: {h_data['incubation']}")
        with c3:
            st.warning("**💊 選択薬剤**")
            if selected_abs:
                for ab in selected_abs: st.write(f"- {ab}")
            else:
                st.write("- 配列から耐性を特定できませんでした")

        # PDF出力
        protocol_dict = {
            "Transformation": h_data["trans_method"],
            "Media": h_data["media"],
            "Antibiotics": ", ".join(selected_abs) if selected_abs else "None"
        }
        
        pdf_bytes = create_pdf(
            host_choice, 
            selected_ez_name, 
            str(primers["Forward"]["seq"]), 
            str(primers["Reverse"]["seq"]), 
            protocol_dict
        )
        
        st.download_button(
            label="📄 プロトコルをPDFでダウンロード",
            data=bytes(pdf_bytes),
            file_name="cloning_protocol.pdf",
            mime="application/pdf"
        )

    elif ins_raw:
        st.warning("インサート配列が短すぎます (40bp以上必要です)")