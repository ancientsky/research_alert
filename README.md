# 🧬 PubMed AI 每日學術論文速遞 (Google Chat 版)

![Python](https://img.shields.io/badge/Python-3.12-blue?logo=python)
![GitHub Actions](https://img.shields.io/badge/Automated-GitHub%20Actions-green?logo=github-actions)
![Gemini AI](https://img.shields.io/badge/AI-Gemini%203%20Flash%20Preview-orange?logo=google)

這是一個全自動化的學術追蹤工具，專為忙碌的研究人員與工程師設計。它會每天自動搜尋 PubMed 上最新的醫學論文，利用 Google Gemini AI 將艱澀的摘要轉化為易懂的中文科普短文，並透過 Webhook 推送到您的 Google Chat 群組。

## ✨ 核心功能

* **自動化搜尋**：每天定時（預設 UTC 00:00）搜尋指定關鍵字（如 Dengue, COVID-19 等）。
* **AI 智慧摘要**：使用 Google 最新 `gemini-3-flash-preview` 模型，生成高品質中文科普摘要。
* **智慧去重**：內建 SQLite 資料庫，記錄已讀論文，確保不會收到重複通知。
* **額度保護**：針對免費用戶優化，包含 API 速率限制（Rate Limiting）與每日推播上限（預設 20 篇），避免帳號被鎖。
* **完全免費**：利用 GitHub Actions 與 Google AI Free Tier，零成本運作。

---

## 🚀 如何部署 (給同事的快速指南)

如果你想複製這份專案並追蹤自己感興趣的領域，請跟著以下步驟操作：

### 第一步：Fork 專案
1. 在本頁面右上方點擊 **「Fork」** 按鈕。
2. 將專案複製到你自己的 GitHub 帳號下。

### 第二步：申請必要的 API Key
你需要準備以下三組金鑰（詳細申請教學請見文件下方）：
1.  **Google Chat Webhook URL**: 用於接收通知。
2.  **NCBI API Key**: 用於存取 PubMed 資料庫。
3.  **Google Gemini API Key**: 用於 AI 摘要生成。

### 第三步：設定 GitHub Secrets
為了安全起見，不要將金鑰直接寫在程式碼裡。請在你的 GitHub 倉庫中設定環境變數：

1. 進入你的倉庫頁面，點選上方 **「Settings」**。
2. 左側選單找到 **「Secrets and variables」** > **「Actions」**。
3. 點擊 **「New repository secret」**，依序新增以下 4 個變數：

| Name | Value (填入你的內容) | 說明 |
| :--- | :--- | :--- |
| `WEBHOOK_URL` | `https://chat.googleapis.com/...` | Google Chat 的 Webhook 連結 |
| `NCBI_EMAIL` | `your_email@example.com` | 你的 Email (NCBI 要求聯絡用) |
| `NCBI_API_KEY` | `sk-xxxx...` | NCBI PubMed API Key |
| `GEMINI_API_KEY` | `AIzaSy...` | Google AI Studio API Key |

### 第四步：啟動自動化
1. 點選上方 **「Actions」** 分頁。
2. 左側點選 **「Daily PubMed Scanner」**。
3. 如果看到 "Workflows aren't being run on this forked repository"，請點擊綠色按鈕啟用。
4. 你可以手動測試：點擊右側的 **「Run workflow」** 按鈕，系統會立即執行一次檢查。

---

## ⚙️ 如何修改關鍵字與設定

所有設定都位於 `pubmed_scanner.py` 檔案的前幾行。你可以直接在 GitHub 網頁上編輯此檔案：

1. 開啟 `pubmed_scanner.py`。
2. 點擊右上角的 ✏️ (Edit) 圖示。
3. 修改 `KEYWORDS` 變數：


# 修改前
KEYWORDS = "Artificial Intelligence"

# 修改後 (例如你想追蹤結核病和HIV)
KEYWORDS = "Tuberculosis AND HIV"

# 進階搜尋 (例如限定 2024 年之後)
KEYWORDS = '"Machine Learning"[Title] AND 2024[Date - Publication]'



4. **提交變更**：捲動到頁面最下方，點擊 **「Commit changes」**。

---

## 🛠️ 各API Key及Webhook詳細申請教學

### 1. 申請 Google AI Studio API Key (Gemini)

這是免費的，但在 Free Tier 下有速率限制（每分鐘 5 次請求，每天 50 次請求）。本程式已針對此限制優化。

1. 前往 [Google AI Studio](https://aistudio.google.com/)。
2. 登入 Google 帳號。
3. 點擊左上角的 **"Get API key"**。
4. 點擊 **"Create API key"**。
5. 複製生成的 `AIzaSy...` 開頭的字串，這就是 `GEMINI_API_KEY`。

### 2. 申請 NCBI PubMed API Key

雖然不加 Key 也能搜尋，但加上 Key 可以提高穩定性與存取頻率（從 3 次/秒提升至 10 次/秒）。

1. 註冊/登入 [NCBI 帳戶](https://www.ncbi.nlm.nih.gov/account/)。
2. 點擊右上角的帳號名稱，進入 **"Dashboard"** (或 Account Settings)。
3. 在 "API Key Management" 區塊，點擊 **"Create an API Key"**。
4. 複製生成的長字串，這就是 `NCBI_API_KEY`。

### 3. 建立 Google Chat Webhook

這是讓機器人把訊息傳到你聊天室的管道。

1. 開啟 Google Chat，進入（或建立）一個 **Space (聊天室)**。
2. 點擊上方的聊天室名稱 > 選擇 **「應用程式與整合功能」(Apps & integrations)**。
3. 點擊 **「管理 Webhook」(Manage webhooks)**。
4. 輸入名稱（例如：`PubMed Bot`）並貼上一個可愛的頭像 URL（選填）。
5. 點擊儲存，複製產生的連結（`https://chat.googleapis.com/...`），這就是 `WEBHOOK_URL`。

---

## 💻 本地開發 (Local Development)

如果你想在自己的電腦上開發或測試：

**1. 安裝環境**
確保你有安裝 Python 3.9+。

```bash
git clone [https://github.com/你的帳號/你的專案名.git](https://github.com/你的帳號/你的專案名.git)
cd 你的專案名
pip install -r requirements.txt

```

**2. 設定環境變數**
在 Windows PowerShell:

```powershell
$env:WEBHOOK_URL="你的連結"
$env:GEMINI_API_KEY="你的Key"
$env:NCBI_API_KEY="你的Key"
$env:NCBI_EMAIL="你的Email"
python pubmed_scanner.py

```

在 Mac/Linux:

```bash
export WEBHOOK_URL="你的連結"
export GEMINI_API_KEY="你的Key"
export NCBI_API_KEY="你的Key"
export NCBI_EMAIL="你的Email"
python pubmed_scanner.py

```

---

## 📂 檔案結構說明

* `pubmed_scanner.py`: 主程式核心邏輯。
* `.github/workflows/daily_check.yml`: 設定 GitHub Actions 定時任務的腳本。
* `requirements.txt`: Python 套件清單 (biopython, google-genai 等)。
* `papers.db`: SQLite 資料庫，儲存已處理過的論文 ID 與摘要（自動生成，勿手動刪除）。

## ⚠️ 注意事項

* **API 額度**：由於使用 Gemini 免費版，程式預設每天最多處理 **20 篇** 新論文，且每篇之間會強制休息 12 秒，以避免觸發 "429 Quota Exceeded" 錯誤。
* **資料庫更新**：GitHub Actions 執行完畢後會自動將更新後的 `papers.db` 推送回倉庫，請確保 Actions 有寫入權限（已在 yml 中設定）。


## ⏱️ 進階設定：調整執行時間 (YAML 詳解)

本專案的自動化排程是由 `.github/workflows/daily_check.yml` 控制的。如果你想要更改機器人運作的時間（例如改成下午執行），請依照以下說明修改。

### 1. 修改執行時間 (Cron Schedule)
GitHub Actions 使用 **UTC 時間**，台灣時間 (CST) 是 **UTC+8**。
請編輯 `.github/workflows/daily_check.yml` 檔案，找到 `on: schedule: - cron:` 這一行。

```yaml
on:
  schedule:
    # 格式為: '分 時 日 月 星期' (使用 UTC 時間)
    # 範例：UTC 00:00 = 台灣時間 08:00
    - cron: '0 0 * * *'

```

**常見時間設定範例：**

| 台灣時間 (UTC+8) | UTC 設定值 | YAML 寫法 (`cron`) |
| --- | --- | --- |
| **早上 08:00** (預設) | 00:00 | `'0 0 * * *'` |
| **早上 09:00** | 01:00 | `'0 1 * * *'` |
| **中午 12:00** | 04:00 | `'0 4 * * *'` |
| **晚上 18:00** | 10:00 | `'0 10 * * *'` |

> 💡 **小工具**：不確定 Cron 怎麼寫？可以使用 [Crontab.guru](https://crontab.guru/) 網站進行測試。

### 2. 理解權限設定 (Permissions)

在 YAML 檔案中，你會看到這段設定：

```yaml
permissions:
  contents: write

```

**這是做什麼的？**
這是為了賦予 GitHub Actions **「寫入權限」**。因為程式執行完畢後，需要將更新後的資料庫 (`papers.db`) 推送 (git push) 回倉庫保存。如果沒有這一行，程式會出現 `403 Forbidden` 錯誤，導致機器人會一直重複推送相同的舊論文。

### 3. 環境變數對應 (Environment Variables)

在 `jobs:` 區塊下方，定義了 Python 程式如何讀取 GitHub Secrets：

```yaml
      - name: Run PubMed Scanner
        env:
          # 左邊是 Python 程式碼用的變數名稱
          # 右邊是 GitHub Settings 裡面的 Secret 名稱
          WEBHOOK_URL: ${{ secrets.WEBHOOK_URL }}
          NCBI_EMAIL: ${{ secrets.NCBI_EMAIL }}
          NCBI_API_KEY: ${{ secrets.NCBI_API_KEY }}
          GEMINI_API_KEY: ${{ secrets.GEMINI_API_KEY }}
        run: python pubmed_scanner.py

```

如果你在 Python 程式碼中新增了其他的 API Key 需求（例如增加了 Slack 通知），記得也要在這裡補上對應的設定，否則程式會讀不到金鑰。


## License

MIT License

