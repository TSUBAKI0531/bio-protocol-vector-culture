🧬 Bio-Protocol Designer
Cloning Design & Experimental Protocol Generator

クローニング（Gibson/In-Fusion）の設計、プライマーのチェック、そして宿主に応じた培養条件の出力を一気通貫で行うWebアプリケーションです。研究者が実験設計にかける時間を短縮し、計算ミスを未然に防ぐことを目的に開発しました。

🌟 主な機能
ベクター解析: GenBank/FASTAファイルを読み込み、耐性遺伝子や制限酵素サイトを自動抽出。

シームレスクローニング設計: Gibson AssemblyやIn-Fusionに最適な15bpオーバーハング付きプライマーを自動生成。

プライマーチェック: 3'末端の相補性（ダイマーリスク）を簡易スクリーニング。

仮想マップ表示: 挿入後のプラスミド図（リニアマップ）をリアルタイムに視覚化。

プロトコル生成: 宿主（大腸菌・酵母・哺乳類）に合わせた導入手法や培地組成をデータベースから提案。

PDF出力: 設計したプライマーと実験条件を、そのまま実験ノートに貼れるPDF形式でエクスポート。

🚀 デプロイ済みURL
[ここにStreamlit CloudのURLを貼り付けてください]

例: https://bio-protocol-designer.streamlit.app

🛠️ 技術スタック
Language: Python 3.10+

Frontend: Streamlit

Bioinformatics: Biopython

Visualization: dna_features_viewer, Matplotlib

Document Export: fpdf2

📂 ディレクトリ構造
Plaintext
bio-protocol-designer/
├── app.py                  # メインアプリケーション
├── modules/                # ロジックモジュール
│   ├── analyzer.py         # 配列解析
│   ├── designer.py         # プライマー設計
│   └── exporter.py         # PDF出力
├── data/                   # データベース
│   └── constants.py        # 培養条件DB
├── requirements.txt        # 依存ライブラリ
└── LICENSE                 # MIT License
💻 ローカルでの実行方法
リポジトリをクローン

Bash
git clone https://github.com/あなたのユーザー名/bio-protocol-designer.git
cd bio-protocol-designer
仮想環境の構築とライブラリのインストール

Bash
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
pip install -r requirements.txt
アプリの起動

Bash
streamlit run app.py
📝 ライセンス
このプロジェクトは MIT License の下で公開されています。

👤 著者
[Sukeda Masaki/TSUBAKI0531]

GitHub: @TSUBAKI0531