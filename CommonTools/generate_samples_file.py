import os
import re
import argparse
import pandas as pd


parser = argparse.ArgumentParser()
parser.add_argument('i', help='输入目录')
args = parser.parse_args()

# 指定文件夹路径
folder_path = args.i

# 获取所有文件名
file_names = os.listdir(folder_path)

# 初始化字典来存储样本名和文件名对
samples = {}

# 正则表达式匹配数字
pattern = re.compile(r'\d+')

for file_name in file_names:
    # 从文件名中提取除去数字部分的文件名
    base_name = re.sub(pattern, '', file_name)

    # 提取文件名中的数字部分
    numbers = pattern.findall(file_name)

    # 以除去数字部分的文件名为键，将文件名添加到对应样本的列表中
    key = base_name
    if key not in samples:
        samples[key] = []
    samples[key].append({'file_name': file_name, 'numbers': numbers})

# 将字典转换为 DataFrame
data = {'Sample Name': [], 'File 1': [], 'File 2': []}

for base_name, files in samples.items():
    # 如果有两个以上文件名以相同的文件名开头，则将其视为一对
    if len(files) >= 2:
        for i in range(len(files)):
            for j in range(i + 1, len(files)):
                # 检查数字是否只有一个不同
                diff_count = sum(a != b for a, b in zip(files[i]['numbers'], files[j]['numbers']))
                if diff_count == 1:
                    data['Sample Name'].append(base_name)
                    data['File 1'].append(files[i]['file_name'])
                    data['File 2'].append(files[j]['file_name'])

df = pd.DataFrame(data)

# 保存 DataFrame 为 CSV 文件
df.to_csv('sample_file_pairs.txt', index=False, sep='\t')

print(df)
