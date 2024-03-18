import pandas as pd



with open('place_seq.txt', 'r') as f:
    data = f.read().split("//")
    
    result_list = []
    for item in data:
        item_list = []
        lines = item.split("\n")
        description = ''
        for line in lines:
            if line.startswith("ID"):
                item_list.append(line.split()[1].strip())
            elif line.startswith("AC"):
                item_list.append(line.split()[1].strip())
            elif line.startswith("DE") or line.startswith("KW") or line.startswith("R") or line.startswith("OS"):
                description = description + line[5:].strip() + ' '
            elif line.startswith("SQ"):
                item_list.append(lines[-2].strip())
                item_list.append(description)
                
            elif line.startswith("XX"):
                continue
            else:
                continue
        result_list.append(item_list)

df = pd.DataFrame(result_list, columns=["ID", "AC", "Sequence", "Description"])
df.to_csv("parse_place_seq.txt", sep="\t", index=False)