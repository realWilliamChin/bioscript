import os

cur_path=os.getcwd()
dir_list=os.listdir(cur_path)

f=open('cluster_gene_id_count.txt','w')
f.write('group'+'\t'+'total'+'\t'+'count'+'\n')
for each in dir_list:
    if ('_cluster' in each)&('txt' in each):
        dic1 = {}
        total=0
        sample=each.split('_cluster')[0]
        num = 0
        f1 = open(each, 'r')
        #with open(each,'r') as f1:
        for each_line1 in f1:
            num+=1
            if num==1:
                title=each_line1.strip()
            else:
                id=each_line1.split('\t')[0]
                cluster=each_line1.split('\t')[-1].strip()
                if cluster in dic1:
                    dic1[cluster]+='\n'+id
                else:
                    dic1[cluster]=id
        f1.close()

        for each_key in dic1:
            new_doc=sample+'_group'+each_key+'.txt'

            f2=open(new_doc,'w')
            f2.write(dic1[each_key])
            count=len(dic1[each_key].split('\n'))
            f2.close()
            sampel_name=sample+'_group'+each_key
            f.write(sampel_name.replace('gene_id_','')+'\t'+str(count)+'\t'+str(round(count/num,2))+'\n')









