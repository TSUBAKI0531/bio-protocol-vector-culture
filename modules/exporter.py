from fpdf import FPDF

def create_pdf(host, enzyme, f_primer, r_primer, protocol_dict):
    """実験プロトコルのPDF生成 (Unicodeエラー回避・英語変換版)"""
    pdf = FPDF()
    pdf.add_page()
    
    # 翻訳マップ（PDFの文字化けを防ぐため、出力時のみ英語に変換）
    trans = {
        "大腸菌 (E. coli)": "E. coli",
        "酵母 (S. cerevisiae)": "S. cerevisiae",
        "哺乳類細胞 (Mammalian)": "Mammalian Cells",
        "ヒートショック法 (42°C, 30-45秒)": "Heat Shock (42C, 30-45s)",
        "酢酸リチウム法": "Lithium Acetate Method",
        "リポフェクション法 または 電解穿孔法": "Lipofection / Electroporation",
        "LB培地 / 37°C": "LB Media / 37C",
        "YPD または SD培地 / 30°C": "YPD or SD Media / 30C",
        "DMEM+10%FBS / 37°C (5% CO2)": "DMEM+10%FBS / 37C (5% CO2)",
        "16-18時間 (振とう培養)": "16-18h (Shaking Culture)",
        "2-3日間": "2-3 Days",
        "一過性: 2-3日 / 安定株: 1-2週間": "Transient: 2-3d / Stable: 1-2w",
        "Ampicillin (100 µg/mL)": "Ampicillin (100 ug/mL)",
        "Kanamycin (50 µg/mL)": "Kanamycin (50 ug/mL)",
        "Puromycin (1-10 µg/mL)": "Puromycin (1-10 ug/mL)"
    }

    # タイトル
    pdf.set_font("Arial", "B", 16)
    pdf.cell(0, 10, "Molecular Cloning Protocol Summary", ln=True, align="C")
    pdf.ln(10)
    
    # 基本情報
    pdf.set_font("Arial", "B", 12)
    pdf.cell(0, 10, f"Target Host: {trans.get(host, host)}", ln=True)
    pdf.cell(0, 10, f"Cloning Site: {enzyme}", ln=True)
    pdf.ln(5)
    
    # 1. プライマー情報
    pdf.cell(0, 10, "1. Primer Sequences (5' -> 3')", ln=True)
    pdf.set_font("Courier", "", 10)
    pdf.multi_cell(0, 8, f"Forward: {f_primer}\nReverse: {r_primer}")
    pdf.ln(5)
    
    # 2. 培養・導入条件
    pdf.set_font("Arial", "B", 12)
    pdf.cell(0, 10, "2. Experimental Conditions", ln=True)
    pdf.set_font("Arial", "", 11)
    for key, value in protocol_dict.items():
        # 英語に変換。変換できない場合は非ASCII文字を?に置換してクラッシュを防ぐ
        val_en = trans.get(value, value)
        val_safe = "".join([c if ord(c) < 256 else "?" for c in str(val_en)])
        pdf.cell(0, 8, f"- {key}: {val_safe}", ln=True)
    
    return pdf.output()