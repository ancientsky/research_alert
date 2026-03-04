import os
import sqlite3
import datetime
import time
import requests
from Bio import Entrez

# --- 新版 SDK Import ---
from google import genai
from google.genai import errors
from google.genai import types
# --- 設定區 ---

WEBHOOK_URL = os.getenv('WEBHOOK_URL')
NCBI_EMAIL = os.getenv('NCBI_EMAIL')
NCBI_API_KEY = os.getenv('NCBI_API_KEY')
GEMINI_API_KEY = os.getenv('GEMINI_API_KEY')

# 搜尋關鍵字
#KEYWORDS = "Artificial Intelligence"
KEYWORDS = '("artificial intelligence"[Title/Abstract] OR "LLM"[Title/Abstract] OR "large language model"[Title/Abstract]) AND hasabstract[text]'

# LLM 模型名稱 (新版 SDK 通用)
#MODEL_NAME = 'gemini-3-flash-preview'
MODEL_NAME = 'gemini-3.1-flash-lite-preview'

# --- 安全限制設定 (針對 Free Tier) ---
API_DELAY_SECONDS = 12  # 每次呼叫 AI 後休息 12 秒 (確保低於 5 RPM)
MAX_DAILY_PAPERS = 20    # 每次執行最多處理 20 篇 (確保低於 20 RPD)

# --- 初始化 ---
Entrez.email = NCBI_EMAIL
Entrez.api_key = NCBI_API_KEY

# 新版 Client 初始化
client = genai.Client(api_key=GEMINI_API_KEY)

def init_db():
    """初始化 SQLite 資料庫與自動遷移"""
    conn = sqlite3.connect('papers.db')
    c = conn.cursor()
    
    c.execute('''CREATE TABLE IF NOT EXISTS papers 
                 (pmid TEXT PRIMARY KEY, 
                  title TEXT, 
                  abstract TEXT, 
                  summary TEXT, 
                  processed_date TEXT)''')
    
    # 自動檢查並新增欄位 (Migration)
    c.execute("PRAGMA table_info(papers)")
    existing_columns = [info[1] for info in c.fetchall()]
    new_columns = {
        'title': 'TEXT', 'abstract': 'TEXT', 
        'summary': 'TEXT', 'processed_date': 'TEXT'
    }
    for col_name, col_type in new_columns.items():
        if col_name not in existing_columns:
            try:
                c.execute(f"ALTER TABLE papers ADD COLUMN {col_name} {col_type}")
            except sqlite3.OperationalError:
                pass 
            
    conn.commit()
    return conn

def search_pubmed(keywords):
    """搜尋過去 1 天內的論文 (限制數量)"""
    try:
        handle = Entrez.esearch(
            db="pubmed", 
            term=keywords, 
            reldate=1, 
            datetype="pdat", 
            retmax=20, 
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
        time.sleep(1) # 禮貌性延遲
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


####old summarize_ai
# def summarize_ai(title, abstract):
#     """使用新版 Google Gen AI SDK 進行科普摘要"""
#     prompt = (
#         f"你是一位專業的科普作家。請閱讀以下醫學論文摘要，"
#         f"用繁體中文寫一段約 100-150 字的「科普摘要」。"
#         f"結構請包含：1. 背景與問題 2. 核心發現 3. 意義。"
#         f"請使用條列式或分段，使其在聊天軟體中易於閱讀。\n\n"
#         f"標題：{title}\n"
#         f"原始摘要：{abstract}"
#     )
    
#     # 新版呼叫方式
#     response = client.models.generate_content(
#         model=MODEL_NAME,
#         contents=prompt
#     )
#     return response.text
####end_old summarize_ai

def summarize_ai(title, abstract):
    """使用新版 Google Gen AI SDK 進行科普摘要"""
    
    prompt = f"""[角色任務] 擔任頂尖科技與醫學科普作家，專為AI推動辦公室團隊解析前沿論文。

[背景資訊] 團隊同仁需要快速掌握 Pubmed 上最新 AI 應用的期刊發展，以便在日常通訊軟體中迅速吸收新知，並評估技術落地的潛在價值。

[具體指令] 閱讀提供的論文標題與原始摘要，將核心資訊提煉為 150 至 200 字的繁體中文重點摘要。請依序產出三個區塊：1. 💡 痛點與背景、2. 🔍 核心AI發現、3. 🚀 應用與意義。

標題：{title}
原始摘要：{abstract}

[約束條件] 全面使用精簡的條列式或短分段排版，段落間保留空白行。每段開頭結合適當的 Emoji 提升視覺引導與手機閱讀體驗。語氣保持專業、易懂且充滿活力，字數嚴格控制在 150 至 200 字之間。"""
    
    # 建立生成設定，啟用 MINIMAL 思考層級
    generate_content_config = types.GenerateContentConfig(
        thinking_config=types.ThinkingConfig(
            thinking_level="MINIMAL",
        ),
    )
    
    # 呼叫模型 (維持非串流以配合 Webhook 傳送，並帶入 config)
    response = client.models.generate_content(
        model=MODEL_NAME,
        contents=prompt,
        config=generate_content_config
    )
    return response.text




def send_chat_message(text):
    if not WEBHOOK_URL: return
    headers = {'Content-Type': 'application/json; charset=UTF-8'}
    try:
        requests.post(WEBHOOK_URL, json={"text": text}, headers=headers)
    except Exception as e:
        print(f"Webhook 連線錯誤: {e}")

def main():
    print(f"開始執行 - SDK: google-genai - 模型: {MODEL_NAME}")
    print(f"限制模式: 每次最多 {MAX_DAILY_PAPERS} 篇，間隔 {API_DELAY_SECONDS} 秒")
    
    conn = init_db()
    c = conn.cursor()
    
    pmids = search_pubmed(KEYWORDS)
    print(f"PubMed 找到 {len(pmids)} 篇候選論文")
    
    new_count = 0
    today_str = datetime.date.today().isoformat()
    
    for pmid in pmids:
        if new_count >= MAX_DAILY_PAPERS:
            print(f"⚠️ 已達到單次執行上限 ({MAX_DAILY_PAPERS} 篇)。")
            send_chat_message(f"⚠️ 今日論文較多，僅推送前 {MAX_DAILY_PAPERS} 篇以節省額度。")
            break

        c.execute("SELECT pmid FROM papers WHERE pmid=?", (pmid,))
        if c.fetchone():
            continue 
            
        print(f"處理新論文: {pmid}")
        title, abstract, doi = fetch_details(pmid)
        
        if title and abstract:
            try:
                # 3. 呼叫 AI (新版錯誤處理)
                summary = summarize_ai(title, abstract)
                # 👇 新增這行：將標準 Markdown 的雙星號替換為 Google Chat 的單星號
                summary = summary.replace('**', '*')
                
                link = f"https://doi.org/{doi}" if doi else f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
                message = (
                    f"📄 *{title}*\n"
                    f"{'-'*20}\n"
                    f"{summary}\n\n"
                    f"🔗 <{link}|點擊閱讀原文>" 
                )
                send_chat_message(message)
                
                c.execute(
                    "INSERT INTO papers (pmid, title, abstract, summary, processed_date) VALUES (?, ?, ?, ?, ?)", 
                    (pmid, title, abstract, summary, today_str)
                )
                conn.commit()
                new_count += 1
                
                print(f"✅ 處理成功，休息 {API_DELAY_SECONDS} 秒...")
                time.sleep(API_DELAY_SECONDS)

            except errors.ClientError as e:
                # 新版 SDK 的錯誤處理，通常 429 會包含在 ClientError 中
                if e.code == 429:
                    print("❌ API 額度已用盡 (429 Resource Exhausted)。")
                    send_chat_message("❌ 今日 AI 額度已用盡，停止後續任務。")
                    break
                else:
                    print(f"⚠️ API ClientError: {e}")
                    time.sleep(5)
            except Exception as e:
                print(f"⚠️ 未預期錯誤: {e}")
                time.sleep(5)
    
    if new_count > 0:
        send_chat_message(f"✅ 今日更新完畢，共推送 {new_count} 篇新論文。")
    else:
        send_chat_message("今日無新增論文。")
    
    conn.close()

if __name__ == "__main__":
    main()
