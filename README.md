---

### 1. 程式原始碼 (`pubmed_scanner.py`)

請將以下程式碼存為專案根目錄下的 `pubmed_scanner.py`。

```python
import os
import sqlite3
import datetime
from Bio import Entrez
import smtplib
from email.mime.text import MIMEText
from email.header import Header
import google.generativeai as genai # 範例使用 Gemini API，也可換成 OpenAI

# 配置區 (透過環境變數讀取)
EMAIL_SENDER = os.getenv('EMAIL_SENDER')
EMAIL_PASSWORD = os.getenv('EMAIL_PASSWORD')
EMAIL_RECEIVER = os.getenv('EMAIL_RECEIVER')
NCBI_EMAIL = os.getenv('NCBI_EMAIL')
NCBI_API_KEY = os.getenv('NCBI_API_KEY')
GEMINI_API_KEY = os.getenv('GEMINI_API_KEY')
KEYWORDS = "CRISPR AND Alzheimer" # 您可以修改此關鍵字

# 初始化 API
Entrez.email = NCBI_EMAIL
Entrez.api_key = NCBI_API_KEY
genai.configure(api_key=GEMINI_API_KEY)
model = genai.GenerativeModel('gemini-1.5-flash')

def init_db():
    conn = sqlite3.connect('papers.db')
    c = conn.cursor()
    c.execute('''CREATE TABLE IF NOT EXISTS papers (pmid TEXT PRIMARY KEY)''')
    conn.commit()
    return conn

def search_pubmed(keywords):
    # 搜尋過去 1 天內的論文
    handle = Entrez.esearch(db="pubmed", term=keywords, reldate=1, datetype="pdat")
    record = Entrez.read(handle)
    handle.close()
    return record["IdList"]

def fetch_details(pmid):
    handle = Entrez.efetch(db="pubmed", id=pmid, retmode="xml")
    records = Entrez.read(handle)
    handle.close()
    
    article = records['PubmedArticle'][0]['MedlineCitation']['Article']
    title = article['ArticleTitle']
    abstract = "".join(article.get('Abstract', {}).get('AbstractText', ["無摘要"]))
    doi = ""
    for id_obj in records['PubmedArticle'][0]['PubmedData']['ArticleIdList']:
        if id_obj.attributes.get('IdType') == 'doi':
            doi = str(id_obj)
    return title, abstract, doi

def summarize_ai(title, abstract):
    prompt = f"請將以下學術論文摘要轉化為淺顯易懂的科普內容（繁體中文），讓非專業人士也能理解核心貢獻。不要使用艱澀術語。\n標題：{title}\n摘要：{abstract}"
    response = model.generate_content(prompt)
    return response.text

def send_email(content):
    msg = MIMEText(content, 'plain', 'utf-8')
    msg['Subject'] = Header(f"今日 PubMed 學術快報 - {datetime.date.today()}", 'utf-8')
    msg['From'] = EMAIL_SENDER
    msg['To'] = EMAIL_RECEIVER

    with smtplib.SMTP_SSL("smtp.gmail.com", 465) as server:
        server.login(EMAIL_SENDER, EMAIL_PASSWORD)
        server.sendmail(EMAIL_SENDER, [EMAIL_RECEIVER], msg.as_string())

def main():
    conn = init_db()
    c = conn.cursor()
    
    pmids = search_pubmed(KEYWORDS)
    new_papers = []
    
    for pmid in pmids:
        c.execute("SELECT pmid FROM papers WHERE pmid=?", (pmid,))
        if not c.fetchone():
            title, abstract, doi = fetch_details(pmid)
            summary = summarize_ai(title, abstract)
            new_papers.append(f"標題：{title}\n連結：https://doi.org/{doi}\n【科普摘要】：\n{summary}\n" + "="*30)
            c.execute("INSERT INTO papers VALUES (?)", (pmid,))
    
    if new_papers:
        full_report = "\n\n".join(new_papers)
        send_email(full_report)
        print(f"已寄送 {len(new_papers)} 篇新論文摘要。")
    else:
        print("今日無新論文。")
    
    conn.commit()
    conn.close()

if __name__ == "__main__":
    main()

```

---

### 2. GitHub Actions 設定 (`.github/workflows/daily_check.yml`)

此檔案負責每天定時執行程式，並將更新後的資料庫（`papers.db`）存回 GitHub 倉庫。

```yaml
name: Daily PubMed Scanner

on:
  schedule:
    - cron: '0 0 * * *' # 每天 UTC 00:00 (台灣時間早上 08:00) 執行
  workflow_dispatch: # 允許手動觸發測試

jobs:
  build:
    runs-on: ubuntu-latest
    steps:
      - name: Checkout code
        uses: actions/checkout@v3

      - name: Set up Python
        uses: actions/setup-python@v4
        with:
          python-version: '3.9'

      - name: Install dependencies
        run: |
          pip install biopython google-generativeai

      - name: Run scanner
        env:
          EMAIL_SENDER: ${{ secrets.EMAIL_SENDER }}
          EMAIL_PASSWORD: ${{ secrets.EMAIL_PASSWORD }}
          EMAIL_RECEIVER: ${{ secrets.EMAIL_RECEIVER }}
          NCBI_EMAIL: ${{ secrets.NCBI_EMAIL }}
          NCBI_API_KEY: ${{ secrets.NCBI_API_KEY }}
          GEMINI_API_KEY: ${{ secrets.GEMINI_API_KEY }}
        run: python pubmed_scanner.py

      - name: Commit and push if database updated
        run: |
          git config --global user.name "github-actions[bot]"
          git config --global user.email "github-actions[bot]@users.noreply.github.com"
          git add papers.db
          git commit -m "Update papers database [skip ci]" || echo "No changes to commit"
          git push

```

---

### 3. 設定步驟指南

#### 第一步：建立 GitHub 倉庫

1. 在 GitHub 建立一個私有倉庫 (Private Repo)。
2. 將 `pubmed_scanner.py` 和上述的 `.yml` 檔案上傳。

#### 第二步：設定 Secrets (金鑰管理)

到 GitHub 倉庫的 `Settings` > `Secrets and variables` > `Actions` > `New repository secret`，新增以下變數：

* `EMAIL_SENDER`: 你的發信 Gmail 地址。
* `EMAIL_PASSWORD`: **Gmail 應用程式專用密碼** (非登入密碼，需在 Google 帳號安全性設定中產生)。
* `EMAIL_RECEIVER`: 你的收信地址。
* `NCBI_EMAIL`: 你的電子郵件（用於 NCBI API 識別）。
* `NCBI_API_KEY`: 從 NCBI 帳戶申請的 API Key。
* `GEMINI_API_KEY`: 從 [Google AI Studio](https://aistudio.google.com/) 申請的免費金鑰。

#### 第三步：關於 SQLite 的持續性

程式碼結尾有一段 `git push`，這非常重要。GitHub Actions 的運作環境是暫時性的，如果不把 `papers.db` 重新 commit 回去，隔天系統就會忘記已經讀過哪些論文。透過這種方式，您的 GitHub 倉庫本身就是您的資料庫存儲點。

---

### 關鍵優點說明

1. **無成本**：GitHub Actions 對私有倉庫有充足的每月免費額度，Gemini API (Flash 版本) 也有極高的免費呼叫額度。
2. **自動去重**：利用 SQLite 檔案儲存 PMID，確保每篇論文只會寄送一次摘要。
3. **靈活調整**：您隨時可以修改 `KEYWORDS` 或調整 `summarize_ai` 函式中的 Prompt，來改變科普摘要的口吻。

您需要我針對其中某一部分（例如如何取得 Gmail 應用程式密碼）做更詳細的說明嗎？
