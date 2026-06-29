import os
cur_path=os.getcwd()
dir_list=os.listdir(cur_path)
f3=open('summary.txt','w')
f3.write('sample'+'\t'+'count'+'\t'+'total'+'\t'+'percentage'+'\n')
for each in dir_list:
    if '_p.txt' in each:
        dic1={}
        out_doc=each.replace('_p','_p_filtered')
        f2=open(out_doc,'w')
        with open(each,'r') as f1:
            num=0
            count=0
            for each_line in f1:
                num+=1

                if num==1:
                    title=each_line.strip()
                    f2.write(each_line)
                else:
                    flag=each_line.strip().split('\t')[-2]
                    if float(flag)<0.05:
                        f2.write(each_line)
                        count+=1
        f3.write(out_doc.replace('.txt','')+'\t'+str(count)+'\t'+str(num)+'\t'+str(round(count/num,4))+'\n')
        f2.close()






