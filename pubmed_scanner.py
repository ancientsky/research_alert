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
