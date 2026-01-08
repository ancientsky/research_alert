import os
import sqlite3
import datetime
import smtplib
from email.mime.text import MIMEText
from email.header import Header

# 第三方套件
from Bio import Entrez
import google.generativeai as genai

# --- 設定區 ---

# 1. 讀取環境變數
EMAIL_SENDER = os.getenv('EMAIL_SENDER')
EMAIL_PASSWORD = os.getenv('EMAIL_PASSWORD')
EMAIL_RECEIVER = os.getenv('EMAIL_RECEIVER')
NCBI_EMAIL = os.getenv('NCBI_EMAIL')
NCBI_API_KEY = os.getenv('NCBI_API_KEY')
GEMINI_API_KEY = os.getenv('GEMINI_API_KEY')

# 2. 搜尋關鍵字與參數
KEYWORDS = "artificial intelligence AND infectious disease"  # 請修改為您感興趣的關鍵字

# 3. LLM 模型設定
# 如果 Google AI Studio 中已有 'gemini-3.0-flash'，請直接將下方字串改為 'gemini-3.0-flash'
MODEL_NAME = 'gemini-3-flash-preview' 

# --- 初始化 ---
Entrez.email = NCBI_EMAIL
Entrez.api_key = NCBI_API_KEY
genai.configure(api_key=GEMINI_API_KEY)
model = genai.GenerativeModel(MODEL_NAME)

def init_db():
    """初始化 SQLite 資料庫，若不存在則建立"""
    conn = sqlite3.connect('papers.db')
    c = conn.cursor()
    c.execute('''CREATE TABLE IF NOT EXISTS papers (pmid TEXT PRIMARY KEY)''')
    conn.commit()
    return conn

def search_pubmed(keywords):
    """搜尋過去 1 天內的論文 PMID"""
    try:
        handle = Entrez.esearch(db="pubmed", term=keywords, reldate=1, datetype="pdat")
        record = Entrez.read(handle)
        handle.close()
        return record["IdList"]
    except Exception as e:
        print(f"PubMed 搜尋失敗: {e}")
        return []

def fetch_details(pmid):
    """根據 PMID 獲取標題、摘要與 DOI"""
    try:
        handle = Entrez.efetch(db="pubmed", id=pmid, retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        
        if not records['PubmedArticle']:
            return None, None, None

        article = records['PubmedArticle'][0]['MedlineCitation']['Article']
        title = article['ArticleTitle']
        
        # 處理摘要可能為列表的情況
        abstract_list = article.get('Abstract', {}).get('AbstractText', ["無摘要"])
        abstract = "".join([str(x) for x in abstract_list]) # 確保轉為字串
        
        doi = ""
        for id_obj in records['PubmedArticle'][0]['PubmedData']['ArticleIdList']:
            if id_obj.attributes.get('IdType') == 'doi':
                doi = str(id_obj)
                break
                
        return title, abstract, doi
    except Exception as e:
        print(f"獲取論文詳情失敗 (PMID: {pmid}): {e}")
        return None, None, None

def summarize_ai(title, abstract):
    """使用 Gemini Flash 進行科普摘要"""
    prompt = (
        f"你是一位專業的科普作家。請閱讀以下醫學論文摘要，"
        f"用繁體中文寫一段約 150 字的「科普摘要」。"
        f"重點：1. 這項研究解決了什麼問題？ 2. 有什麼新發現？ 3. 對未來有什麼影響？"
        f"請避免艱澀術語，讓一般大眾也能看懂。\n\n"
        f"標題：{title}\n"
        f"原始摘要：{abstract}"
    )
    try:
        response = model.generate_content(prompt)
        return response.text
    except Exception as e:
        return f"摘要生成失敗: {e}"

def send_email(content):
    """發送匯總郵件"""
    if not content: return

    msg = MIMEText(content, 'plain', 'utf-8')
    msg['Subject'] = Header(f"【每日論文速遞】{datetime.date.today()}", 'utf-8')
    msg['From'] = EMAIL_SENDER
    msg['To'] = EMAIL_RECEIVER

    try:
        with smtplib.SMTP_SSL("smtp.gmail.com", 465) as server:
            server.login(EMAIL_SENDER, EMAIL_PASSWORD)
            server.sendmail(EMAIL_SENDER, [EMAIL_RECEIVER], msg.as_string())
        print("郵件發送成功！")
    except Exception as e:
        print(f"郵件發送失敗: {e}")

def main():
    print(f"開始執行 - 模型: {MODEL_NAME} - Python 3.12")
    conn = init_db()
    c = conn.cursor()
    
    pmids = search_pubmed(KEYWORDS)
    print(f"找到 {len(pmids)} 篇相關論文 (含舊資料)")
    
    new_papers_content = []
    
    for pmid in pmids:
        # 檢查是否已存在資料庫
        c.execute("SELECT pmid FROM papers WHERE pmid=?", (pmid,))
        if c.fetchone():
            continue # 已處理過，跳過
            
        print(f"處理新論文: {pmid}")
        title, abstract, doi = fetch_details(pmid)
        
        if title and abstract:
            summary = summarize_ai(title, abstract)
            link = f"https://doi.org/{doi}" if doi else f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
            
            entry = (
                f"📄 標題：{title}\n"
                f"🔗 連結：{link}\n"
                f"💡 科普摘要：\n{summary}\n"
                f"{'-'*40}"
            )
            new_papers_content.append(entry)
            
            # 成功處理後才寫入 DB
            c.execute("INSERT INTO papers VALUES (?)", (pmid,))
            conn.commit() # 每次成功都存檔，避免中斷遺失
    
    if new_papers_content:
        full_report = "以下是您訂閱的最新論文科普摘要：\n\n" + "\n\n".join(new_papers_content)
        send_email(full_report)
    else:
        print("今日無新論文，未發送郵件。")
    
    conn.close()

if __name__ == "__main__":
    main()
