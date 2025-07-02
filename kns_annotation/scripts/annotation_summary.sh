#!/bin/bash

# 检查参数数量
if [ $# -ne 7 ]; then
    echo "Usage: $0 <anno_fasta> <swiss_goid> <swiss_gene_def> <kegg_gene_def> <nr_gene_def> <output_pic> <output_summary_file>"
    exit 1
fi

# 接收参数
ANNO_FASTA=$1
SWISS_GOID=$2
SWISS_GENE_DEF=$3
KEGG_GENE_DEF=$4
NR_GENE_DEF=$5
OUTPUT_PIC=$6
OUTPUT_SUMMARY_FILE=$7

# 创建临时文件
GO_ID="GO.txt"
SWISS_ID="Swiss.txt"
KEGG_ID="KEGG.txt"
NR_ID="NR.txt"

# 处理输入文件
tail -n +2 ${SWISS_GOID} | cut -f 1 > ${GO_ID}
tail -n +2 ${SWISS_GENE_DEF} | cut -f 1 > ${SWISS_ID}
tail -n +2 ${KEGG_GENE_DEF} | cut -f 1 > ${KEGG_ID}
tail -n +2 ${NR_GENE_DEF} | cut -f 1 > ${NR_ID}

# 绘制Venn图
python /home/colddata/qinqiang/script/Plot/Venn/venn_plot.py \
    -i ${GO_ID} ${SWISS_ID} ${KEGG_ID} ${NR_ID} \
    --pic-name ${OUTPUT_PIC}

# 计算统计信息
total_id_count=$(grep -c '>' ${ANNO_FASTA})
kegg_count=$(wc -l ${KEGG_ID} | awk '{print $1}')
go_count=$(wc -l ${GO_ID} | awk '{print $1}')
nr_count=$(wc -l ${NR_ID} | awk '{print $1}')
swiss_count=$(wc -l ${SWISS_ID} | awk '{print $1}')

# 生成汇总文件
echo -e "Database\tTotal\tAnnotated\tPercentage" > ${OUTPUT_SUMMARY_FILE}
echo -e "KEGG\t$total_id_count\t$((kegg_count - 1))\t$(( (kegg_count - 1) * 100 / total_id_count ))%" >> ${OUTPUT_SUMMARY_FILE}
echo -e "GO\t$total_id_count\t$((go_count - 1))\t$(( (go_count - 1) * 100 / total_id_count ))%" >> ${OUTPUT_SUMMARY_FILE}
echo -e "NR\t$total_id_count\t$((nr_count - 1))\t$(( (nr_count - 1) * 100 / total_id_count ))%" >> ${OUTPUT_SUMMARY_FILE}
echo -e "Swiss Protein\t$total_id_count\t$((swiss_count - 1))\t$(( (swiss_count - 1) * 100 / total_id_count ))%" >> ${OUTPUT_SUMMARY_FILE}

# 清理临时文件
rm ${GO_ID} ${SWISS_ID} ${KEGG_ID} ${NR_ID} 
