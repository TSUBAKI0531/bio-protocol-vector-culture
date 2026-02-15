import streamlit as st
import io
import matplotlib.pyplot as plt
from Bio import SeqIO, Restriction
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from dna_features_viewer import BiopythonTranslator

# 自作モジュールのインポート
from data.constants import CULTURE_DB
from modules.analyzer import analyze_vector_resistance
from modules.designer import design_gibson_primers, check_dimer
from modules.exporter import create_pdf

# --- UI 設定 ---
st.set_page_config(page_title="BioDesigner v1.0", layout="wide")
st.title("🧬 Cloning & Expression Designer")

# --- サイドバー: 宿主設定 ---
st.sidebar.header("1. Host Selection")
host_choice = st.sidebar.selectbox("宿主を選択", list(CULTURE_DB.keys()))

# --- メイン: ファイルアップロード ---
uploaded_file = st.file_uploader("ベクターファイルをアップロード (GenBank推奨)", type=["gb", "fasta"])

if uploaded_file:
    # 配列読み込み
    content = uploaded_file.getvalue().decode("utf-8")
    handle = io.StringIO(content)
    file_format = "genbank" if uploaded_file.name.endswith(".gb") else "fasta"
    record = SeqIO.read(handle, file_format)
    
    st.divider()

    # --- セクション1: 制限酵素 & クローニング設計 ---
    col1, col2 = st.columns([1, 1])
    
    with col1:
        st.subheader("✂️ Cloning Site")
        analysis = Restriction.Analysis(Restriction.Commando, record.seq)
        unique_cutters = sorted([str(e) for e in analysis.unique_cutters()])
        
        if unique_cutters:
            selected_ez = st.selectbox("Unique Cutter を選択", unique_cutters)
            rb = Restriction.RestrictionBatch([selected_ez])
            # 切断箇所の取得（1-basedを0-basedへ）
            cut_pos = rb.search(record.seq)[getattr(Restriction, selected_ez)][0]
            st.info(f"切断位置: {cut_pos} bp")
        else:
            st.error("利用可能なUnique Cutterが見つかりません。")
            st.stop()

    with col2:
        st.subheader("🧬 Insert Sequence")
        ins_raw = st.text_area("インサート配列を入力", placeholder="ATGC...", height=100)
    
    if ins_raw and len(ins_raw) >= 40:
        ins_seq = Seq(ins_raw.strip().upper())
        primers = design_gibson_primers(record.seq, ins_seq, cut_pos)
        
        # プライマー表示 & チェック
        st.subheader("✅ Designed Primers")
        p1, p2 = st.columns(2)
        with p1:
            st.write("**Forward (5'->3')**")
            st.code(str(primers["Forward"]["seq"]))
            st.caption(f"Tm: {primers['Forward']['Tm']}°C / {check_dimer(primers['Forward']['seq'], primers['Forward']['seq'])[1]}")
        with p2:
            st.write("**Reverse (5'->3')**")
            st.code(str(primers["Reverse"]["seq"]))
            st.caption(f"Tm: {primers['Reverse']['Tm']}°C / {check_dimer(primers['Reverse']['seq'], primers['Reverse']['seq'])[1]}")

        # --- 仮想マップ表示 ---
        st.divider()
        st.subheader("🗺️ Virtual Construct Map")
        final_seq = record.seq[:cut_pos] + ins_seq + record.seq[cut_pos:]
        final_rec = SeqRecord(final_seq, name="Construct")
        # フィーチャー追加
        final_rec.features.append(SeqFeature(FeatureLocation(0, cut_pos), type="misc_feature", qualifiers={"label": ["Vector_Part1"]}))
        final_rec.features.append(SeqFeature(FeatureLocation(cut_pos, cut_pos+len(ins_seq)), type="CDS", qualifiers={"label": ["INSERT"], "color": ["#ff4b4b"]}))
        final_rec.features.append(SeqFeature(FeatureLocation(cut_pos+len(ins_seq), len(final_seq)), type="misc_feature", qualifiers={"label": ["Vector_Part2"]}))
        
        fig, ax = plt.subplots(figsize=(10, 1.5))
        graphic_record = BiopythonTranslator().translate_record(final_rec)
        graphic_record.plot(ax=ax, with_ruler=True)
        st.pyplot(fig)

    # --- セクション2: プロトコル出力 ---
    st.divider()
    st.subheader("📋 Experimental Protocol")
    
    res_genes = analyze_vector_resistance(record)
    h_data = CULTURE_DB[host_choice]
    
    # UI表示用データの整理
    protocol_disp = {
        "Transformation": h_data["trans_method"],
        "Media": h_data["media"],
        "Incubation": h_data["incubation"],
        "Selection": ", ".join([h_data["antibiotics"].get(g, "Unknown") for g in res_genes]) if res_genes else "None detected"
    }
    
    st.json(protocol_disp)

    # PDFダウンロード
    if ins_raw:
        pdf_data = create_pdf(host_choice, selected_ez, str(primers["Forward"]["seq"]), str(primers["Reverse"]["seq"]), protocol_disp)
        st.download_button("📄 PDF プロトコルをダウンロード", data=bytes(pdf_data), file_name="protocol.pdf")