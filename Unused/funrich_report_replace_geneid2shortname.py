#!/usr/bin/env python
import pandas as pd
import os


def replace_shortname(funrich_file, shortname):
    def my_rep(s):
        s = s.split(';')
        result = []
        for i in s:
            i = i.strip()
            if i in shortname.keys():
                result.append(shortname[i])
            else:
                result.append(i)
        result_str = ';'.join(result)
        return result_str

    file_xlsx = pd.ExcelFile(funrich_file, engine='openpyxl')
    sheet_names = file_xlsx.sheet_names
    with pd.ExcelWriter("result/" + funrich_file, engine='xlsxwriter') as writer:
        for each_sheet in sheet_names:
            df = pd.read_excel(file_xlsx, each_sheet, skiprows=8, engine='openpyxl')
            df['gene/proteins mapped (from input data set)'] = df['gene/proteins mapped (from input data set)'].apply(lambda x: my_rep(x))
            df.to_excel(excel_writer=writer, sheet_name=each_sheet, index=False) # type: ignore


if __name__ == '__main__':
    shortname = pd.read_csv('Sus_scrofa_shortname.txt', sep='\t', index_col=0)
    shortname.dropna(subset='Gene_shortname', inplace=True)
    shortname_dict = shortname.to_dict()['Gene_shortname']

    for each_file in os.listdir():
        if not each_file.endswith('xlsx'):
            continue
        replace_shortname(each_file, shortname_dict)

