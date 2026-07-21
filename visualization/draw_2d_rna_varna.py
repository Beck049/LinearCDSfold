import os
import subprocess

# 1. 測試資料
# seq = "CAUUGGCACAGGCAAUGACUGUGUCUCAGACGACACUAGCAGUCGGUGUCAGUCGACACUAGCAGUCAGUGUCAGAUUGCAUGACCGGGUUCCCUU"
# ss = "..(((((...(((..((...(((((......))))).....(((((......))))).....)))))...))))).((....))..((....)).."

seq = "AUGGGGCUGCGGUACAGCAAGGAGGUGAGGGAUCGUCAUGGGGAUAAGGAUCCUGAGGGCAGGAUCCCUAUCACCCAGACGAUGCCUCAGACACUGUACGGCAGAUAUAAUUGUAAGAGCUGUUGGUUCGCCAACAAGGGGCUGAUCAAGUGCUCUAACCAUUAUAUCUGCCUGCGGUGUCUGACCUCCAUGCUGAACCGCACUGACUACUGUGAGAUCUGUGGGGAGGUGCUGCCCAAGCGGCUCACCUUUGAGACCACUCCCACAGCACCUCCCUACACCCCG"
ss  = "..((((.((((((.(((((.((((((((((.((((((.((((((((.((((((((....)))))))).)))).)))))))))).))))((((((((((.((((((((((((.((.(((((..((((((.(((.......))).))))))..))))).)).)))))))))))))))))))))).)))))).))))).))))))((.((....)).))...(((((((((((((((....((.((((((....)))).)).))....)))))))))))))))))))."

# 2. 取得 VARNAv3-93.jar 的完整絕對路徑 (避免 Java Classpath 找不到檔案)
jar_name = "VARNAv3-93.jar"
jar_path = os.path.abspath(jar_name)

if not os.path.exists(jar_path):
    print(f"❌ 錯誤：在目錄中找不到 {jar_name}！")
    print(f"請確認檔案名稱是否為 {jar_name}")
    exit(1)

output_image = "rna_structure.png"

# 3. 組裝 Java 執行指令 (針對 Java 26 與 Classpath 修正)
cmd = [
    "java",
    "-Djava.awt.headless=true",  # Headless 繪圖模式
    "--add-opens",
    "java.desktop/java.awt=ALL-UNNAMED",  # 解除 Java 26 的 AWT 模組限制
    "--add-opens",
    "java.desktop/javax.swing=ALL-UNNAMED",
    "-cp",
    jar_path,  # 使用完整絕對路徑
    "fr.orsay.lri.varna.applications.VARNAcmd",  # 正確的 Entry Class
    "-sequenceDBN",
    seq,
    "-structureDBN",
    ss,
    "-algorithm",
    "naview",  # 放射狀折疊演算法
    "-bpStyle",
    "line",  # 配對線條樣式
    "-resolution",
    "2.0",  # 高畫質輸出
    "-o",
    output_image,
]

try:
    print(f"🚀 正在讀取 `{jar_path}` 並呼叫 VARNA 繪製...")
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    print(f"✨ 繪製成功！圖片已儲存至：{os.path.abspath(output_image)}")
except subprocess.CalledProcessError as e:
    print("❌ 繪製失敗，錯誤訊息如下：")
    print("STDOUT:", e.stdout)
    print("STDERR:", e.stderr)