


f2=open('gene_id_seed_vs_germinated1.txt','w')
dic = {}
with open('gene_id_seed_vs_germinated.txt','r') as f1:

    for each_line in f1:
        if each_line not in dic:
            dic[each_line]=''
            f2.write(each_line)
        else:
            print(each_line)

f2.close()


















# #################################################################################33
# def rpsblast_annotation(input_file):
#     dic_cdd = {}
#     # cdd  numeber对应的content
#     with open('/opt/biosoft/RpsbProc-x64-linux/data/cddid.tbl', 'r') as f1:
#         for each_line in f1:
#             cdd_number = each_line.split('\t')[0]
#             content = each_line.strip().strip('\t', 1)[-1]
#             dic_cdd[cdd_number] = content
#
#     parser = argparse.ArgumentParser(description="把当前文件夹中的Gene_def文件夹中的所有注释文件，匹配_ID.txt文件里的geneid，输出到_ID_all_def.txt中 ")
#
#     parser.add_argument('--f', dest='input',
#                         help='cog/smart/tifr/...的rpsblast结果')
#
#     with open (input_file) as f2:
#         new_doc=input_file.replace('.out','_annotation.blast')
#         f3=open(new_doc,'w')
#         for each_line2 in f2:
#             cdd=each_line2.strip().split('\t')[1].replace('CDD:','')
#             if cdd in dic_cdd:
#                 f3.write(each_line2.strip()+'\t'+dic_cdd[cdd]+'\n')
#         f3.close()




