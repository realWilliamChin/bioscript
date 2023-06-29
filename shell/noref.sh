#!/bin/bash
# set -e -o pipefail

check_conda_env() {
    conda_env=$(conda env list | grep "*" | cut -d " " -f 1)
    echo $conda_env
}

check_compareinfo_and_samplesdescribed() {
    # 如果 compare info 和 samples described 文件不存在则跳过检查
    if [[ ! -f ${work_dir}/compare_info.txt ]] || [[ ! -f ${work_dir}/samples_described.txt ]]; then
        echo "compare_info.txt 或 samples_described.txt 文件不存在"
        return 0
    fi
    echo "检查 compare_info.txt 和 samples_described.txt 文件是否有不对应的地方"
    python ${python_script}/check_SampDesAndCompInfo.py
}

### fastqc
exec_fastqc() {
    mkdir ${pinjiedata}/fastqc
    echo "正在后台生成 fastqc 报告"
    # 检查 pinjiedata 目录下都是 fq 文件还是文件夹，如果是文件夹，则循环文件夹下执行 fastqc
    # 如果是文件，则直接执行 fastqc
    if [[ $(ls ${pinjiedata} | grep ".fq" | wc -l) -gt 1 ]]; then
        echo "nohup fastqc ${pinjiedata}/*.fq \
        -o ${pinjiedata}/fastqc > ${log}/fastqc.log 2>&1 && echo "fastqc 已完成" & " | bash
    elif [[ $(ls ${pinjiedata} | grep ".fq" | wc -l) -lt 1 ]]; then
        for i in $(ls ${pinjiedata}); do
            echo "nohup fastqc ${pinjiedata}/${i}/*.fq \
            --o ${pinjiedata}/fastqc > ${log}/fastqc.log 2>&1 && echo "fastqc 已完成" & " | bash
        done
    fi
}

### pinjie 步骤
exec_pinjie() {
    # 生成 sample trinity 文件, 可能需要修改一下 sample trinity，太多了拼接太慢了
    source /home/train/miniconda3/bin/activate base
    python /home/colddata/qinqiang/script/noreference/trinity_samples_file.py
    cat ${work_dir}/samples_trinity.txt
    echo "##############################################\n"
    read -p "是否更改 samples_trinity.txt，回车继续"
    mv ${work_dir}/samples_trinity.txt ${work_dir}/${specie}_samples_trinity.txt

    mkdir ${assemble_trinity}
    # 生成 trinity.fasta 文件
    echo "正在执行pinjie流程"
    Trinity --seqType fq \
        --max_memory ${max_memory}G \
        --no_salmon \
        --no_version_check \
        --samples_file ${specie}_samples_trinity.txt \
        --output ${assemble_trinity} \
        --CPU ${num_threads} --SS_lib_type RF --normalize_reads \
        --min_contig_length 500

    perl /home/train/trainingStuff/bin/extract_longest_isoforms_from_TrinityFasta.pl \
        ${assemble_trinity}/Trinity.fasta >${assemble_trinity}/unigene_longest.fasta

    cd-hit-est -i ${assemble_trinity}/unigene_longest.fasta \
        -o ${unigene_fasta} -c 0.95 && echo "cd-hit-est 已完成，已生成 ${specie}_unigene.fasta"

    /opt/biosoft/Trinity-v2.8.5/util/TrinityStats.pl \
        ${unigene_fasta} >${assemble_trinity}/assemble_stat.txt
    seqkit stats ${unigene_fasta} >>${assemble_trinity}/assemble_stat.txt
    grep '>' ${unigene_fasta} | cut -d ' ' -f 1 | tr -d '>' > ${specie}_all_gene_id.txt
    echo "pinjie 流程已完成"
    assemble_report
}
# assemble_stat.txt for report
assemble_report() {
    Total_sequence_num=$(tail -n 1 assemble_stat.txt | awk -F' ' '{print $4}' | tr -d ',')
    Total_sequence_bases=$(grep 'Total assembled bases' assemble_stat.txt | head -n 1 | awk -F':' '{print $2}' | tr -d ' ')
    Percent_GC=$(grep 'Percent GC' assemble_stat.txt | awk -F':' '{print $2}' | tr -d ' ')
    Largest_transcript=$(tail -n 1 assemble_stat.txt | awk -F' ' '{print $8}' | tr -d ',')
    Smallest_transcript=$(tail -n 1 assemble_stat.txt | awk -F' ' '{print $6}' | tr -d ',')
    Average_length=$(tail -n 1 assemble_stat.txt | awk -F' ' '{print $7}' | tr -d ',')
    N50=$(grep N50 assemble_stat.txt | head -n 1 | awk -F':' '{print $2}' | tr -d ' ')
    echo -e "Total_sequence_num\t${Total_sequence_num}
    Total_sequence_bases\t${Total_sequence_bases}
    Percent_GC\t${Percent_GC}
    Largest_transcript\t${Largest_transcript}
    Smallest_transcript\t${Smallest_transcript}
    Average_length\t${Average_length}
    N50\t${N50}" > assemble_stat_report.txt
}

### Annotation
## Swiss
py_swiss() {
    cd ${annotation} || exit
    python ${python_script}/annotation/swiss.py
    cd ${work_dir} || exit
}
exec_swiss() {
    mkdir ${annotation}
    echo "执行 Annotation - Swiss 步骤"
    /opt/biosoft/ncbi-blast-2.9.0+/bin/blastx \
        -db /home/data/ref_data/Linux_centos_databases/2019_Unprot_databases/swissprot \
        -query ${unigene_fasta} \
        -out ${annotation}/${specie}_unigene_swiss.blast \
        -max_target_seqs 20 \
        -evalue 1e-5 \
        -num_threads ${num_threads} \
        -outfmt "6 qacc sacc pident qcovs qcovhsp ppos length mismatch gapopen qstart qend sstart send evalue bitscore stitle" 
    py_swiss && echo "swiss 注释完成"
}

## nr
py_nr() {
    cd ${annotation} || exit
    python ${python_script}/annotation/nr.py
    cd ${work_dir} || exit
}
exec_nr() {
    echo "执行 Annotation - nr 步骤"
    mkdir -p ${annotation}/temp
    diamond blastx --db /home/data/ref_data/db/diamond_nr/diamond_nr \
        --query ${unigene_fasta} \
        --out ${annotation}/${specie}_unigene_nr_diamond.blast \
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle \
        --sensitive \
        --max-target-seqs 1 \
        --evalue 1e-5 \
        --id 30 \
        --block-size 20.0 \
        --tmpdir ${annotation}/temp \
        --index-chunks 1 && echo 'nr 注释完成'
}

### cog
# emappey.py youhuma 45分钟，运行完需要 conda deactivate
exec_cog() {
    cd ${annotation} || exit
    echo "正在切换到 python27 conda 环境"
    source /home/train/miniconda3/bin/activate python27
    check_conda_env
    emapper.py -i ${specie}_unigene.fasta \
        --output ${specie}_unigene \
        --cpu ${num_threads} \
        --data_dir /home/data/xuezhen/download/eggNOG/eggnog.db \
        -m diamond --translate
    echo "正在切换到 base conda 环境"
    source /home/train/miniconda3/bin/activate base
    check_conda_env
    cd ${work_dir} || exit
}

### kegg
kegg() {
    # TODO: 后面继续写爬虫程序
    # 拿着 $specie_unigene.fasta 去 kegg 网站生成 keg 文件，再做下面的东西 -t plant/animal 动物或植物
    python ${python_script}/annotation/kegg.py -t ${specie_type} -i ${all_gene_id}
}

### transdecoder
transdecoder() {
    # 也可以自动生成到报告里
    cd ${annotation} || exit
    mkdir transdecoder
    cd transdecoder || exit
    TransDecoder.LongOrfs -t ${unigene_fasta}
    # 先运行完上面那句，才能运行下面这句
    TransDecoder.Predict -t ${unigene_fasta}
    cd ${work_dir} || exit
}
exec_annotation() {
    mkdir ${annotation}
    cp ${unigene_fasta} ${annotation}
    exec_swiss
    exec_nr
    # exec_cog
    transdecoder
}

### rsem
rsem() {
    mkdir ${reference}
    # 建库
    cd ${reference} || exit
    rsem-prepare-reference --bowtie2 ${unigene_fasta} ${specie}
    cd ${work_dir} || exit
    mkdir ${mapping}
    # 读取 samples_described.txt 循环比对
    cat ${work_dir}/samples_described.txt | grep -v sample | while read line; do
        sample=$(echo $line | awk '{print $2}')
        echo "正在执行 $sample 的 rsem 流程"
        rsem-calculate-expression -p ${num_threads} \
            --bowtie2 --strandedness reverse --paired-end \
            ${pinjiedata}/${sample}_1.clean.fq.gz ${pinjiedata}/${sample}_2.clean.fq.gz \
            ${reference}/${specie} ${mapping}/${sample} \
            1>${mapping}/${sample}_mapping.txt \
            2>${mapping}/${sample}.log
    done
    cd ${mapping} || exit
    # 提出 rawreads
    python /home/colddata/chen/03_transcript/No_Reference_transcriptome/all_sample_raw_reads.py
    # 提出 fpkm
    python /home/colddata/chen/03_transcript/No_Reference_transcriptome/all_sample_fpkm.py
    # 过滤掉 rawreads 小于 50 的 gene
    python /home/colddata/chen/03_transcript/No_Reference_transcriptome/count_filtered.py
    # 提出过滤后的基因的 fpkm 值
    python /home/colddata/chen/03_transcript/No_Reference_transcriptome/fpkm_filtered.py
    cd ${work_dir} || exit

    # 列名按照 samples_described.txt 重新排序
    python ${python_script}/realignment_fpkm_columns.py -r \
    -s ${work_dir}/samples_described.txt -f ${mapping}/reads_matrix_filtered.txt
    python ${python_script}/realignment_fpkm_columns.py -r \
    -s ${work_dir}/samples_described.txt -f ${mapping}/fpkm_matrix_filtered.txt
    python ${python_script}/realignment_fpkm_columns.py -r \
    -s ${work_dir}/samples_described.txt -f ${mapping}/gene_count_matrix.txt
    python ${python_script}/realignment_fpkm_columns.py -r \
    -s ${work_dir}/samples_described.txt -f ${mapping}/gene_fpkm_matrix.txt

    # 合并 fpkm 和 reads 矩阵
    cd ${mapping} || exit
    echo -e "正在合并 fpkm 和 reads 矩阵，Def 使用\n"
    echo -e "${kegg_gene_def}"
    echo -e "${nr_gene_def}"
    echo -e "${swiss_gene_def}\n"
    python "${python_script}"/merge_fpkm_reads_matrix.py \
        -f fpkm_matrix_filtered.txt \
        -r reads_matrix_filtered.txt \
        -k "$kegg_gene_def" \
        -n "$nr_gene_def" \
        -s "$swiss_gene_def" \
        -o fpkm_and_reads_matrix_filtered_data_def.txt
    python "${python_script}"/merge_fpkm_reads_matrix.py \
        -f gene_fpkm_matrix.txt \
        -r gene_count_matrix.txt \
        -k "$kegg_gene_def" \
        -n "$nr_gene_def" \
        -s "$swiss_gene_def" \
        -o fpkm_and_reads_data_def.txt
    cd ${work_dir} || exit
}

### multi deseq
### 需要创建一个 compare.txt 和 samples_described.txt
multi_deseq() {
    mkdir ${multideseq}
    cd ${multideseq} || exit
    cp ${mapping}/*matrix_filtered.txt ${multideseq}/
    cp ${work_dir}/samples_described.txt ${multideseq}/
    cp ${work_dir}/compare_info.txt ${multideseq}/
    Rscript /home/data/command/Reference_transcriptome/multiple_samples_DESeq2_V7.r ${rlog_number}

    cd ${multideseq} || exit
    python ${python_script}/de_results_add_def.py \
        -k "$kegg_gene_def" \
        -n "$nr_gene_def" \
        -s "$swiss_gene_def"
    cd ${work_dir} || exit
}

### jiaofu
jiaofu() {
    # 设置目录变量
    Background_materials_dir=${jiaofu}/00_Background_materials
    Funrich_software_def_files_dir=${jiaofu}/00_Funrich_software_def_files
    Reference_genome_annotation_files_dir=${jiaofu}/00_Reference_genome_annotation_files
    Sequencing_quality_dir=${jiaofu}/00_测序质控文件
    Original_expression_data_dir=${jiaofu}/01_Original_expression_data
    DEG_analysis_dir=${jiaofu}/02_DEG_analysis
    PPI_analysis_KEGG_pathways_dir=${jiaofu}/03_PPI_analysis_KEGG_pathways
    Softwares_for_analysis_dir=${jiaofu}/04_Softwares_for_analysis
    Prep_files_dir=${jiaofu}/Prep_files_dir

    transdecoder_file_dir=${Reference_genome_annotation_files_dir}/拼接转录本序列/
    DEG_analysis_Analysis_dir=${DEG_analysis_dir}/Analysis
    DEG_data_dir=${DEG_analysis_Analysis_dir}/Expression_data
    DEG_data_graphs_dir=${DEG_analysis_Analysis_dir}/Expression_data_graphs
    DEG_analysis_Funrich_report_dir=${DEG_analysis_dir}/Analysis_Funrich_report
    Funrich_Barplot_dir=${DEG_analysis_Funrich_report_dir}/DEG_Barplot_graphs
    Funrich_Bubble_dir=${DEG_analysis_Funrich_report_dir}/DEG_Bubble_graphs
    Funrich_report_dir=${DEG_analysis_Funrich_report_dir}/Funrich_DEG_analysis_report

    wego=${Prep_files_dir}/wego

    # 创建目录
    mkdir -p "$Background_materials_dir" \
        "$Funrich_software_def_files_dir" \
        "$Reference_genome_annotation_files_dir" \
        "$Sequencing_quality_dir" \
        "$Original_expression_data_dir" \
        "$DEG_analysis_dir" \
        "$PPI_analysis_KEGG_pathways_dir" \
        "$Softwares_for_analysis_dir" \
        "$Prep_files_dir" \
        "$transdecoder_file_dir" \
        "$DEG_analysis_Analysis_dir" \
        "$DEG_data_dir" \
        "$DEG_data_graphs_dir" \
        "$DEG_analysis_Funrich_report_dir" \
        "$Funrich_Barplot_dir" \
        "$Funrich_Bubble_dir" \
        "$Funrich_report_dir"

    echo -e "\ncopy files to 00_Funrich_software_def_files"
    cp "$annotation"/*GO*ID.txt $Funrich_software_def_files_dir
    cp "$annotation"/*KEGG.txt "$Funrich_software_def_files_dir"
    cp "$annotation"/*shortname.txt "$Funrich_software_def_files_dir"
    cp "$assemble_trinity/${specie}_all_gene_id.txt" "$Funrich_software_def_files_dir"

    echo -e "\ncopy files to 00_Reference_genome_annotation_files"
    cp "$annotation"/*_gene_def.txt "$Reference_genome_annotation_files_dir"

    echo -e "\ncopy files to transdecoder_file_dir"
    cp "$annotation"/transdecoder/*.pep "$transdecoder_file_dir"
    cp "$annotation"/transdecoder/*.cds "$transdecoder_file_dir"

    echo -e "\ncopy files to 00_测序质控文件"
    cp "$pinjiedata"/fastqc/*.zip "$Sequencing_quality_dir"

    echo -e "\ncopy files to 01_Original_expression_data"
    cp "$mapping"/*_data_def.txt "$Original_expression_data_dir"

    echo -e "\ncopy files to 02_DEG_analysis"
    cp "$multideseq"/*Down_ID.txt "$DEG_analysis_Analysis_dir"
    cp "$multideseq"/*Up_ID.txt "$DEG_analysis_Analysis_dir"
    cp "$multideseq"/*_DEG_data.txt "$DEG_data_dir"

    # Prep_files
    echo -e "\ncopy files to Prep_files"
    cp "$assemble_trinity"/assemble_stat_report.txt "$Prep_files_dir"

    # wego
    echo -e "\ncopy files to wego"
    cp "$annotation"/*GO*ID.txt "$wego"
    cp "$annotation"/*idNo_def.txt "$wego"/all_geneid_GO.txt
    cp "$multideseq"/*Down_ID.txt "$wego"
    cp "$multideseq"/*Up_ID.txt "$wego"
    cd "$wego" || exit
    # python ${python_script}/wego.py
    # rm *Down_ID.txt *Up_ID.txt *GO*ID.txt
    cd "$word_dir" || exit
}

### 变量
# 目录
python_script=/home/colddata/qinqiang/script
work_dir=$(pwd)
log=${work_dir}/log
rawdata=${work_dir}/01_Rawdata
pinjiedata=${work_dir}/02_Pinjiedata
assemble_trinity=${work_dir}/03_Assemble_trinity
annotation=${work_dir}/04_Annotation
reference=${work_dir}/05_Reference
mapping=${work_dir}/06_Mapping_mapping
multideseq=${work_dir}/07_Multi_deseq
jiaofu=${work_dir}/jiaofu_perpare

# 文件
unigene_fasta=${assemble_trinity}/${specie}_unigene.fasta
kegg_gene_def=$(realpath -s  ${annotation}/*KEGG_gene_def.txt)
nr_gene_def=$(realpath -s ${annotation}/*nr_gene_def.txt)
swiss_gene_def=$(realpath -s ${annotation}/*swiss_gene_def.txt)

### 读取输入参数
# 默认参数
max_memory=500
num_threads=60
rlog_number=1

while getopts ":h-:" opt; do
    case ${opt} in
    -)
        case "${OPTARG}" in
        specie)
            specie="${!OPTIND}"
            OPTIND=$((OPTIND + 1))
            ;;
        threads)
            num_threads="${!OPTIND}"
            OPTIND=$((OPTIND + 1))
            ;;
        max_memory)
            max_memory="${!OPTIND}"
            OPTIND=$((OPTIND + 1))
            ;;
        rlog_number)
            rlog_number="${!OPTIND}"
            OPTIND=$((OPTIND + 1))
            ;;
        help)
            echo "使用 -h 查看帮助信息"
            exit 0
            ;;
        *)
            echo "无效的选项：--$OPTARG"
            exit 1
            ;;
        esac
        ;;
    h)
        echo "帮助信息："
        echo "使用方法：./noref.sh [选项]"
        echo "选项："
        echo "  --specie <specie> 指定物种"
        echo "  --threads <threads> 指定线程数"
        echo "  --max_memory <max_memory> 指定最大内存"
        echo "  --rlog_number <rlog_number> multi deseq 的倍数"
        echo "  -h 显示帮助信息"
        exit 0
        ;;
    \?)
        echo "无效的选项：-$OPTARG"
        exit 1
        ;;
    esac
done

show_menu() {
    echo "
    无参执行程序
    物种名: $specie

    0 执行所有流程
    ————————————————
    1 执行 merge
    2 执行 fastqc
    3 执行 pinjie
    4 执行 pinjie result
    5 执行 Annotation
        5.1 执行 swiss
        5.2 执行 nr
        5.3 执行 cog
        5.4 执行 kegg
        5.5 执行 transdecoder
    6 执行 rsem
    7 执行 multi_deseq
    8 整理交付目录
 "
    # show_status
    read -p "请输入选择(其他任意退出): " select

    case "${select}" in
    0)
        echo "程序开始执行"
        # check_source_data
        exec_fastqc
        exec_pinjie
        pinjie_result
        exec_annotation
        rsem
        multi_deseq
        jiaofu
        ;;
    1)
        merge
        ;;
    2)
        exec_fastqc
        ;;
    3)
        exec_pinjie
        ;;
    4)
        pinjie_result
        ;;
    5)
        exec_annotation
        ;;
    5.1)
        exec_swiss
        ;;
    5.2)
        exec_nr
        ;;
    5.3)
        exec_cog
        ;;
    5.4)
        exec_kegg
        ;;
    5.5)
        transdecoder
        ;;
    6)
        rsem
        ;;
    7)
        multi_deseq
        ;;
    8)
        jiaofu
        ;;
    *)
        # LOGE "请输入正确的数字: "
        exit 0
        ;;
    esac
}
show_menu
