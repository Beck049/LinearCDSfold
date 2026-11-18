import csv
import math
import sys

# 指定要比對的目標欄位（順序會依照此列表輸出）
TARGET_COLUMNS = ["MFE", "CAI", "crossTalkCnt", "cMFE", "cCAI", "cCrossTalkCnt"]

# 作為主鍵比對的 ID 欄位名稱
ID_COLUMN = "UniPort ID"

# 浮點數比對容許誤差 (可依需求調整，如 1e-4)
TOLERANCE = 1e-4


def normalize_value(val: str):
    """清理字串並嘗試轉為 float 比對數值；若非數值則回傳修剪後的字串。"""
    val = val.strip()
    try:
        return float(val)
    except ValueError:
        return val


def is_equal(v1, v2) -> bool:
    """比較兩值是否相同（支援浮點數誤差容許）。"""
    if isinstance(v1, float) and isinstance(v2, float):
        return math.isclose(v1, v2, abs_tol=TOLERANCE)
    return v1 == v2


def load_csv_data(filepath: str):
    """讀取 CSV 並以 UniPort ID 為 key 建立字典，同時回傳欄位名稱列表。"""
    data = {}
    with open(filepath, mode="r", encoding="utf-8-sig") as f:
        reader = csv.reader(f)
        try:
            raw_headers = next(reader)
        except StopIteration:
            return {}, []

        headers = [h.strip() for h in raw_headers]

        if ID_COLUMN not in headers:
            print(f"[錯誤] 檔案 {filepath} 找不到主鍵欄位 '{ID_COLUMN}'！現有欄位: {headers}")
            sys.exit(1)

        for line_num, row in enumerate(reader, start=2):
            if not row or not any(row):
                continue
            row_dict = {
                h: normalize_value(v)
                for h, v in zip(headers, row)
            }
            row_id = row_dict.get(ID_COLUMN)

            if row_id in data:
                print(f"[警告] 檔案 {filepath} 第 {line_num} 行出現重複 {ID_COLUMN}: {row_id}")
            data[row_id] = row_dict

    return data, headers


def main():
    if len(sys.argv) < 3:
        print("Usage: python check.py <test.csv> <answer.csv>")
        sys.exit(1)

    test_file = sys.argv[1]
    ans_file = sys.argv[2]

    test_data, test_headers = load_csv_data(test_file)
    ans_data, ans_headers = load_csv_data(ans_file)

    test_ids = set(test_data.keys())
    ans_ids = set(ans_data.keys())

    has_error = False

    # 1. 檢查缺少或多出的 UniPort ID
    missing_ids = ans_ids - test_ids
    extra_ids = test_ids - ans_ids

    if missing_ids:
        has_error = True
        print(f"[錯誤] 比對序列缺少以下 {ID_COLUMN} ({len(missing_ids)} 個):")
        for mid in sorted(missing_ids):
            print(f"  - 缺少: {mid}")

    if extra_ids:
        has_error = True
        print(f"[錯誤] 比對序列多出以下 {ID_COLUMN} ({len(extra_ids)} 個):")
        for eid in sorted(extra_ids):
            print(f"  - 多出: {eid}")

    # 2. 篩選出 test.csv 與 answer.csv 皆存在的目標欄位（test.csv 無此欄位則跳過）
    test_header_set = set(test_headers)
    ans_header_set = set(ans_headers)

    compare_columns = []
    for col in TARGET_COLUMNS:
        if col in test_header_set and col in ans_header_set:
            compare_columns.append(col)

    # 3. 針對共同存在的 UniPort ID 比對指定欄位數值
    common_ids = test_ids & ans_ids

    diff_count = 0
    for row_id in sorted(common_ids):
        row_test = test_data[row_id]
        row_ans = ans_data[row_id]

        row_errors = []
        for col in compare_columns:
            val_test = row_test.get(col)
            val_ans = row_ans.get(col)

            if not is_equal(val_test, val_ans):
                row_errors.append(f"  欄位 [{col}]: 輸出={val_test} | 預期={val_ans}")

        if row_errors:
            has_error = True
            diff_count += 1
            print(f"[數值不符] {ID_COLUMN}: {row_id}")
            for err in row_errors:
                print(err)

    if not has_error:
        print(f"比對完全正確")
        sys.exit(0)
    else:
        print(f"\n比對結束：發現錯誤，共有 {diff_count} 筆資料欄位不符。")
        sys.exit(1)


if __name__ == "__main__":
    main()