import argparse
import pandas as pd
import re
import os

def extract_uniprot_id(header):
    # 使用正則表達式提取 >sp|ID| 或 >tr|ID| 中的 ID 部分
    match = re.search(r'>[a-z]{2}\|([^|]+)\|', header)
    if match:
        return match.group(1)
    return None

def disassemble_fasta(filename):
    if not os.path.exists(filename):
        print(f"Error: File {filename} not found.")
        return

    current_id = None
    current_seq = []

    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # 如果之前有抓到序列，先寫入前一個檔案
                if current_id and current_seq:
                    with open(f"{current_id}.fasta", 'w') as out_f:
                        out_f.write(f">{current_id}\n" + "".join(current_seq) + "\n")
                
                # 提取新的 ID
                current_id = extract_uniprot_id(line)
                current_seq = []
                if not current_id:
                    print(f"Warning: Could not parse ID from header: {line}")
            else:
                current_seq.append(line)

        # 寫入最後一條序列
        if current_id and current_seq:
            with open(f"{current_id}.fasta", 'w') as out_f:
                out_f.write(f">{current_id}\n" + "".join(current_seq) + "\n")
    
    print(f"Disassemble complete. Individual FASTA files created.")

def assemble_fasta(csv_filename, output_name):
    if not os.path.exists(csv_filename):
        print(f"Error: CSV file {csv_filename} not found.")
        return

    df = pd.read_csv(csv_filename)
    
    # 確保輸出檔名有 .fasta 副檔名
    if not output_name.endswith('.fasta'):
        output_name += '.fasta'

    with open(output_name, 'w') as out_f:
        for uniprot_id in df['UniProt ID']:
            file_name = f"{uniprot_id}.fasta"
            if os.path.exists(file_name):
                with open(file_name, 'r') as f:
                    out_f.write(f.read().strip() + "\n")
            else:
                print(f"Warning: File {file_name} not found, skipping...")

    print(f"Assemble complete. Combined file saved as: {output_name}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="FASTA Disassemble/Assemble Tool")
    parser.add_argument("filename", nargs='?', help="Input FASTA or CSV filename")
    parser.add_argument("--disassemble", action="store_true", help="Split FASTA into individual files based on ID")
    parser.add_argument("--assemble", action="store_true", help="Merge files based on CSV list")
    parser.add_argument("--name", help="Output filename for assembly")

    args = parser.parse_args()

    if args.disassemble and args.filename:
        disassemble_fasta(args.filename)
    elif args.assemble and args.filename and args.name:
        assemble_fasta(args.filename, args.name)
    else:
        parser.print_help()