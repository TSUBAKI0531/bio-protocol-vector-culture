from fpdf import FPDF

def create_pdf(host, enzyme, f_primer, r_primer, protocol_text):
    pdf = FPDF()
    pdf.add_page()
    pdf.set_font("Arial", "B", 16)
    pdf.cell(0, 10, "Bio-Protocol Summary", ln=True, align="C")
    pdf.ln(10)
    
    pdf.set_font("Arial", "B", 12)
    pdf.cell(0, 10, f"Host: {host}", ln=True)
    pdf.cell(0, 10, f"Enzyme: {enzyme}", ln=True)
    
    pdf.ln(5)
    pdf.cell(0, 10, "1. Primers", ln=True)
    pdf.set_font("Courier", "", 10)
    pdf.multi_cell(0, 10, f"Fwd: {f_primer}\nRev: {r_primer}")
    
    pdf.ln(5)
    pdf.set_font("Arial", "B", 12)
    pdf.cell(0, 10, "2. Conditions", ln=True)
    pdf.set_font("Arial", "", 11)
    pdf.multi_cell(0, 10, protocol_text)
    
    return pdf.output()