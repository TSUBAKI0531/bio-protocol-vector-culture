CULTURE_DB = {
    "大腸菌 (E. coli)": {
        "trans_method": "ヒートショック法 (42°C, 30-45秒)",
        "media": "LB培地 / 37°C",
        "incubation": "16-18時間 (振とう培養)",
        "antibiotics": {
            "amp": "Ampicillin (100 µg/mL)", 
            "kan": "Kanamycin (50 µg/mL)", 
            "puro": "Puromycin (N/A)", 
            "bla": "Blasticidin (50 µg/mL)",
            "tet": "Tetracycline (12.5 µg/mL)"
        }
    },
    "酵母 (S. cerevisiae)": {
        "trans_method": "酢酸リチウム法",
        "media": "YPD または SD培地 / 30°C",
        "incubation": "2-3日間",
        "antibiotics": {
            "kan": "G418 (200 µg/mL)", 
            "neo": "G418 (200 µg/mL)", 
            "hyg": "Hygromycin B (200 µg/mL)",
            "zeo": "Zeocin (100 µg/mL)"
        }
    },
    "哺乳類細胞 (Mammalian)": {
        "trans_method": "リポフェクション法 または 電解穿孔法",
        "media": "DMEM+10%FBS / 37°C (5% CO2)",
        "incubation": "一過性: 2-3日 / 安定株: 1-2週間",
        "antibiotics": {
            "puro": "Puromycin (1-10 µg/mL)", 
            "neo": "G418 (500 µg/mL)", 
            "bla": "Blasticidin (10 µg/mL)",
            "hyg": "Hygromycin B (200 µg/mL)"
        }
    }
}