import os
import pandas as pd
cur_path=os.getcwd()
dir_list=os.listdir(cur_path)
for each in dir_list:
    if '_p_filtered.txt' in each:
        sample=each.split('_p_filtered.txt')[0]
        dic1={}
        ls=[]
        with open(each,'r') as f1:
            num=0
            for each_line1 in f1:
                num+=1
                if num==1:
                    continue
                gene=each_line1.split('\t')[0]
                ls.append(gene)
        fpkm=sample+'.xlsx'
        df=pd.read_excel(fpkm)
        df_fltered = df[df['GeneID'].isin(ls)]
        df_fltered.reset_index(drop=True, inplace=True)
        sample_group=sample+'_group.txt'
        with open(sample_group,'r') as f2:
            dic_group={}
            num2=0
            for each_line2 in f2:
                num2+=1
                if num2==1:
                    continue
                group=each_line2.split('\t')[0]
                sample2=each_line2.split('\t')[-1].strip()
                if group in dic_group:
                    dic_group[group]+='\t'+sample2
                else:
                    dic_group[group]=sample2
        df_result=pd.DataFrame()
        df_result['GeneID']=df_fltered['GeneID']
        for each_key in dic_group:


            sample_list = dic_group[each_key].strip().split('\t')
            sum = df.loc[:, sample_list[0]]
            if sample_list[0]!=sample_list[-1]:
                for each_sample in sample_list[1:]:
                    sum += df.loc[:, each_sample]
            ave = round(sum / len(sample_list), 2)
            df_result[each_key] = round(sum / len(sample_list), 2)
        out_doc=sample+'_annovafiltered_fpkm.txt'
        df_result.to_csv(out_doc, sep='\t',index=None)




