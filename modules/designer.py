from Bio.SeqUtils import MeltingTemp as mt

def design_gibson_primers(vector_seq, insert_seq, cut_site_index, overhang_len=15):
    upstream_vec = vector_seq[cut_site_index - overhang_len : cut_site_index]
    downstream_vec = vector_seq[cut_site_index : cut_site_index + overhang_len]
    ins_start = insert_seq[:20]
    ins_end = insert_seq[-20:]
    
    fwd = upstream_vec + ins_start
    rev = (ins_end + downstream_vec).reverse_complement()
    
    return {
        "Forward": {"seq": fwd, "Tm": round(mt.Tm_NN(fwd), 1)},
        "Reverse": {"seq": rev, "Tm": round(mt.Tm_NN(rev), 1)}
    }

def check_dimer(seq1, seq2, threshold=5):
    s1_str = str(seq1).upper()
    s2_rc_str = str(seq2.reverse_complement()).upper()
    three_prime_end = s1_str[-threshold:]
    if three_prime_end in s2_rc_str:
        return True, f"⚠️ 3'末端({three_prime_end})が結合リスクあり"
    return False, "✅ OK"