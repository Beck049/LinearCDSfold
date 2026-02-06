import pandas as pd
import sys
from typing import Any
from collections.abc import Iterator

def read_protein_pandas(file_path: str) -> pd.DataFrame:
    # 自動偵測標頭並讀取
    df = pd.read_csv(file_path, skipinitialspace=True)
    return df

def yield_protein_blocks(file_path: str) -> Iterator[dict[str, Any]]:
    with open(file_path, 'r', encoding='utf-8') as f:
        block = {}
        capture_next_lines = 0  # 用來追蹤是否正在讀取 CDS 與 Structure
        
        for line in f:
            line = line.strip()
            
            # 1. 分隔線處理：吐出資料並重置
            if line.startswith('====='):
                if block:
                    yield block
                    block = {}
                capture_next_lines = 0
                continue
            
            # 2. 處理「標籤在下一行」的特殊欄位
            if line.startswith("Coding sequence and its secondary structure:"):
                capture_next_lines = 2 # 接下來兩行分別是 CDS 和 Structure
                continue
            
            if capture_next_lines == 2:
                block["CDS"] = line
                capture_next_lines -= 1
                continue
            elif capture_next_lines == 1:
                block["Structure"] = line
                capture_next_lines -= 1
                continue
            
            # 3. 解析冒號分隔的一般欄位
            if ':' in line:
                key, value = line.split(':', 1)
                key, value = key.strip(), value.strip()
                
                if key == "Protein ID":
                    block["UniPort ID"] = value
                elif key == "Protein length":
                    block["Protein length"] = int(value)
                elif key == "Lambda":
                    block["Lambda"] = float(value)
                elif key == "Folding free energy":
                    block["MFE"] = float(value.split()[0])
                elif key == "CAI":
                    block["CAI"] = float(value)
                elif key == "Total runtime":
                    block["Run time"] = float(value.split()[0])
        
        # 處理檔案末尾
        if block:
            yield block

def check_structure(structure: str) -> bool:
    stack = 0
    for char in structure:
        if char == '(':
            stack += 1
        elif char == ')':
            stack -= 1
            if stack < 0:
                return False
    return stack == 0

# main function
if __name__ == "__main__":
    if len(sys.argv) >= 3:
        txt_file = sys.argv[1]
        csv_file = sys.argv[2]
    else:
        txt_file = 'multiple_seqs.txt'
        csv_file = 'result.csv'

    df = read_protein_pandas(csv_file)
    txt = yield_protein_blocks(txt_file)

    counter = 0
    for block in txt:
        # block 是一個 dict
        print(f"------------- 正在處理 {block['UniPort ID']} -------------")
        df_row = df[counter:counter+1]
        if not df_row.empty and df_row.iloc[0]['UniPort ID'] == block['UniPort ID']:
            # 比對 DataFrame 中的資料
            is_match = True
            is_match &= (df_row.iloc[0]['Protein length'] == block['Protein length'])
            is_match &= (abs(df_row.iloc[0]['Lambda'] - block['Lambda']) < 1e-6)
            is_match &= (round(df_row.iloc[0]['MFE'], 3) == block['MFE'])
            is_match &= (round(df_row.iloc[0]['CAI'], 3) == block['CAI'])
            if not is_match:
                print("DataFrame 中的資料與預期不符：")
                print(f"received: {block['Lambda']}\t{block['MFE']}\t{block['CAI']}")
                print(f"expected: {df_row.iloc[0]['Lambda']}\t{df_row.iloc[0]['MFE']}\t{df_row.iloc[0]['CAI']}")
            # 比對 2nd Structure
            structure = block['Structure']
            is_valid_structure = check_structure(structure)
            if not is_valid_structure:
                print("2nd Structure 格式不正確")
        else:
            print("DataFrame 中找不到對應的 UniPort ID")
            sys.exit(1)

        counter += 1