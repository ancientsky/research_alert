import os
import sqlite3
import datetime
import time
import requests
from Bio import Entrez
import google.generativeai as genai

# --- 設定區 ---

WEBHOOK_URL = os.getenv('WEBHOOK_URL')
NCBI_EMAIL = os.getenv('NCBI_EMAIL')
NCBI_API_KEY = os.getenv('NCBI_API_KEY')
GEMINI_API_KEY = os.getenv('GEMINI_API_KEY')

# 搜尋關鍵字
KEYWORDS = "Artificial Intelligence AND Epidemics"

# LLM 模型
MODEL_NAME = 'gemini-3-flash-preview' 

# --- 初始化 ---
Entrez.email = NCBI_EMAIL
Entrez.api_key = NCBI_API_KEY
genai.configure(api_key=GEMINI_API_KEY)
model = genai.GenerativeModel(MODEL_NAME)

def init_db():
    """
    初始化 SQLite 資料庫
    修改點：增加欄位檢測，若使用舊版 DB 會自動升級 Schema
    """
    conn = sqlite3.connect('papers.db')
    c = conn.cursor()
    
    # 建立表格 (如果完全不存在)
    c.execute('''CREATE TABLE IF NOT EXISTS papers 
                 (pmid TEXT PRIMARY KEY, 
                  title TEXT, 
                  abstract TEXT, 
                  summary TEXT, 
                  processed_date TEXT)''')
    
    # --- 自動遷移邏輯 (針對舊版資料庫) ---
    # 檢查目前有哪些欄位
    c.execute("PRAGMA table_info(papers)")
    existing_columns = [info[1] for info in c.fetchall()]
    
    # 如果缺欄位，就動態補上 (Migration)
    new_columns = {
        'title': 'TEXT',
        'abstract': 'TEXT',
        'summary': 'TEXT',
        'processed_date': 'TEXT'
    }
    
    for col_name, col_type in new_columns.items():
        if col_name not in existing_columns:
            print(f"資料庫升級：新增欄位 {col_name}")
            c.execute(f"ALTER TABLE papers ADD COLUMN {col_name} {col_type}")
            
    conn.commit()
    return conn

def search_pubmed(keywords):
    """搜尋過去 1 天內的論文 (限制 10 篇最新)"""
    try:
        handle = Entrez.esearch(
            db="pubmed", 
            term=keywords, 
            reldate=1, 
            datetype="pdat", 
            retmax=10, 
            sort='date'
        )
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
        
        abstract_list = article.get('Abstract', {}).get('AbstractText', ["無摘要"])
        abstract = "".join([str(x) for x in abstract_list])
        
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
        f"用繁體中文寫一段約 100-150 字的「科普摘要」。"
        f"結構請包含：1. 背景與問題 2. 核心發現 3. 意義。"
        f"請使用條列式或分段，使其在聊天軟體中易於閱讀。\n\n"
        f"標題：{title}\n"
        f"原始摘要：{abstract}"
    )
    try:
        response = model.generate_content(prompt)
        return response.text
    except Exception as e:
        return f"摘要生成失敗: {e}"

def send_chat_message(text):
    if not WEBHOOK_URL:
        return

    headers = {'Content-Type': 'application/json; charset=UTF-8'}
    data = {"text": text}
    
    try:
        requests.post(WEBHOOK_URL, json=data, headers=headers)
    except Exception as e:
        print(f"Webhook 連線錯誤: {e}")

def main():
    print(f"開始執行 - 模型: {MODEL_NAME} (詳細存檔版)")
    conn = init_db()
    c = conn.cursor()
    
    pmids = search_pubmed(KEYWORDS)
    print(f"找到 {len(pmids)} 篇相關論文")
    
    new_count = 0
    today_str = datetime.date.today().isoformat() # 格式：YYYY-MM-DD
    
    for pmid in pmids:
        c.execute("SELECT pmid FROM papers WHERE pmid=?", (pmid,))
        if c.fetchone():
            continue 
            
        print(f"處理新論文: {pmid}")
        title, abstract, doi = fetch_details(pmid)
        
        if title and abstract:
            summary = summarize_ai(title, abstract)
            link = f"https://doi.org/{doi}" if doi else f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
            
            # 1. 發送通知
            message = (
                f"📄 *{title}*\n"
                f"{'-'*20}\n"
                f"{summary}\n\n"
                f"🔗 <{link}|點擊閱讀原文>" 
            )
            send_chat_message(message)
            new_count += 1
            
            # 2. 存入詳細資料 (修改點)
            # 資料結構：(pmid, title, abstract, summary, processed_date)
            try:
                c.execute(
                    "INSERT INTO papers (pmid, title, abstract, summary, processed_date) VALUES (?, ?, ?, ?, ?)", 
                    (pmid, title, abstract, summary, today_str)
                )
                conn.commit()
            except sqlite3.OperationalError as e:
                # 預防性的錯誤捕捉，雖然 init_db 已經處理過遷移
                print(f"資料庫寫入錯誤: {e}")

            time.sleep(1) 
    
    if new_count > 0:
        send_chat_message(f"✅ 今日更新完畢，共推送 {new_count} 篇新論文。")
    else:
        print("今日無新論文，未發送訊息。")
    
    conn.close()

if __name__ == "__main__":
    main()
