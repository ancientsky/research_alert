import os
import sqlite3
import datetime
import time
import requests # 新增：用於呼叫 Webhook
from Bio import Entrez
import google.generativeai as genai

# --- 設定區 ---

# 1. 讀取環境變數 (移除 Email 相關，新增 Webhook)
WEBHOOK_URL = os.getenv('WEBHOOK_URL')
NCBI_EMAIL = os.getenv('NCBI_EMAIL')
NCBI_API_KEY = os.getenv('NCBI_API_KEY')
GEMINI_API_KEY = os.getenv('GEMINI_API_KEY')

# 2. 搜尋關鍵字
KEYWORDS = "artificial intelligence AND infecious disease"

# 3. LLM 模型設定
MODEL_NAME = 'gemini-3-flash-preview' 

# --- 初始化 ---
Entrez.email = NCBI_EMAIL
Entrez.api_key = NCBI_API_KEY
genai.configure(api_key=GEMINI_API_KEY)
model = genai.GenerativeModel(MODEL_NAME)

def init_db():
    """初始化 SQLite 資料庫"""
    conn = sqlite3.connect('papers.db')
    c = conn.cursor()
    c.execute('''CREATE TABLE IF NOT EXISTS papers (pmid TEXT PRIMARY KEY)''')
    conn.commit()
    return conn

def search_pubmed(keywords):
    """
    搜尋過去 1 天內的論文 PMID
    修改點：限制只回傳最新的 10 篇 (retmax=10, sort='date')
    """
    try:
        handle = Entrez.esearch(
            db="pubmed", 
            term=keywords, 
            reldate=1, 
            datetype="pdat", 
            retmax=10,    # <--- 限制回傳最大數量為 10
            sort='date'   # <--- 強制按日期排序，確保是「最新」的 10 篇
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
    """發送訊息到 Google Chat Webhook"""
    if not WEBHOOK_URL:
        print("錯誤：未設定 WEBHOOK_URL")
        return

    headers = {'Content-Type': 'application/json; charset=UTF-8'}
    
    # 建立 payload
    data = {"text": text}
    
    try:
        response = requests.post(WEBHOOK_URL, json=data, headers=headers)
        if response.status_code != 200:
            print(f"Webhook 發送失敗: {response.text}")
    except Exception as e:
        print(f"Webhook 連線錯誤: {e}")

def main():
    print(f"開始執行 - 模型: {MODEL_NAME}")
    conn = init_db()
    c = conn.cursor()
    
    pmids = search_pubmed(KEYWORDS)
    print(f"找到 {len(pmids)} 篇相關論文")
    
    new_count = 0
    
    # 如果有新論文，先發送一個開頭訊息
    # (為了避免洗版，這裡我們先計算未讀數量，若要即時發送可省略此步驟，直接進入迴圈)
    
    for pmid in pmids:
        c.execute("SELECT pmid FROM papers WHERE pmid=?", (pmid,))
        if c.fetchone():
            continue 
            
        print(f"處理新論文: {pmid}")
        title, abstract, doi = fetch_details(pmid)
        
        if title and abstract:
            summary = summarize_ai(title, abstract)
            link = f"https://doi.org/{doi}" if doi else f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
            
            # 組合單篇論文訊息
            message = (
                f"📄 *{title}*\n"
                f"{'-'*20}\n"
                f"{summary}\n\n"
                f"🔗 <{link}|點擊閱讀原文>" 
            )
            
            send_chat_message(message)
            new_count += 1
            
            # 寫入 DB
            c.execute("INSERT INTO papers VALUES (?)", (pmid,))
            conn.commit()
            
            # 避免觸發 API 速率限制，稍作停頓
            time.sleep(1) 
    
    if new_count > 0:
        send_chat_message(f"✅ 今日更新完畢，共推送 {new_count} 篇新論文。")
    else:
        print("今日無新論文，未發送訊息。")
    
    conn.close()

if __name__ == "__main__":
    main()
