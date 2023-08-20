#!/bin/bash
# set -e -o pipefail

check_conda_env() {
    conda_env=$(conda env list | grep "*" | cut -d " " -f 1)
    echo $conda_env
}

check_compareinfo_and_samplesdescribed() {
    # 如果 compare info 和 samples described 文件不存在则跳过检查
    echo "检查 compare_info.txt 和 samples_described.txt 中..."
    if [[ ! -f ${work_dir}/compare_info.txt ]] || [[ ! -f ${work_dir}/samples_described.txt ]]; then
        echo "compare_info.txt 或 samples_described.txt 文件不存在或没有放在工作目录"
        return 0
    fi
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
    # python ${python_script}/noreference/trinity_samples_file.py
    cat ${work_dir}/samples_trinity.txt
    echo -e "\n##############################################\n"
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
        --CPU ${num_threads} \
        --SS_lib_type RF \
        --normalize_reads \
        --min_contig_length 500 \
        && echo "已生成 Trinity.fasta"

    echo "正在执行 extract_longest_isoforms_from_TrinityFasta.pl"
    perl /home/train/trainingStuff/bin/extract_longest_isoforms_from_TrinityFasta.pl \
        ${assemble_trinity}/Trinity.fasta \
        > ${assemble_trinity}/unigene_longest.fasta \
        && echo "已生成 unigene_longest.fasta"

    echo "cd-hit-est 正在执行"
    cd-hit-est -i ${assemble_trinity}/unigene_longest.fasta \
        -o ${assemble_trinity}/${specie}_unigene.fasta \
        -c 0.95 \
        && echo "cd-hit-est 已完成，已生成 ${specie}_unigene.fasta"

    /opt/biosoft/Trinity-v2.8.5/util/TrinityStats.pl \
        ${assemble_trinity}/${specie}_unigene.fasta \
        > ${assemble_trinity}/assemble_stat.txt
    seqkit stats ${assemble_trinity}/${specie}_unigene.fasta >> ${assemble_trinity}/assemble_stat.txt
    grep '>' ${assemble_trinity}/${specie}_unigene.fasta | cut -d ' ' -f 1 | tr -d '>' > ${specie}_all_gene_id.txt
    echo "pinjie 流程已完成"
    cd ${assemble_trinity} || exit
    assemble_report
    cd ${work_dir} || exit
}
# assemble_stat.txt for report
assemble_report() {
    # 针对 TrinityStats.pl 和 seqkit stats 生成的 assemble_stat.txt 进行处理，可直接插入报告使用
    Total_sequence_num=$(tail -n 1 assemble_stat.txt | awk -F' ' '{print $4}' | tr -d ',')
    Total_sequence_bases=$(grep 'Total assembled bases' assemble_stat.txt | head -n 1 | awk -F':' '{print $2}' | tr -d ' ')
    Percent_GC=$(grep 'Percent GC' assemble_stat.txt | awk -F':' '{print $2}' | tr -d ' ')
    Largest_transcript=$(tail -n 1 assemble_stat.txt | awk -F' ' '{print $8}' | tr -d ',')
    Smallest_transcript=$(tail -n 1 assemble_stat.txt | awk -F' ' '{print $6}' | tr -d ',')
    Average_length=$(tail -n 1 assemble_stat.txt | awk -F' ' '{print $7}' | tr -d ',')
    N50=$(grep N50 assemble_stat.txt | head -n 1 | awk -F':' '{print $2}' | tr -d ' ')
    echo -en "Total_sequence_num\t${Total_sequence_num}\n
    Total_sequence_bases\t${Total_sequence_bases}\n
    Percent_GC\t${Percent_GC}\n
    Largest_transcript\t${Largest_transcript}\n
    Smallest_transcript\t${Smallest_transcript}\n
    Average_length\t${Average_length}\n
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
        -query ${assemble_trinity}/${specie}_unigene.fasta \
        -out ${annotation}/${specie}_unigene_swiss.blast \
        -max_target_seqs 20 \
        -evalue 1e-5 \
        -num_threads ${num_threads} \
        -outfmt "6 qacc sacc pident qcovs qcovhsp ppos length mismatch gapopen qstart qend sstart send evalue bitscore stitle" 
    # 如果 swiss 注释成功，则对 swiss blast 结果进行处理
    if [ -f ${annotation}/${specie}_unigene_swiss.blast ]; then
        py_swiss && echo "swiss 注释完成"
    else
        echo "swiss blast 失败"
    fi
    
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
        --query ${assemble_trinity}/${specie}_unigene.fasta \
        --out ${annotation}/${specie}_unigene_nr_diamond.blast \
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle \
        --sensitive \
        --max-target-seqs 1 \
        --evalue 1e-5 \
        --id 30 \
        --block-size 20.0 \
        --tmpdir ${annotation}/temp \
        --index-chunks 1 && echo 'nr 注释完成'
    # 如果 nr 执行成功，则对 nr 结果进行处理
    if [ -f ${annotation}/${specie}_unigene_nr_diamond.blast ]; then
        py_nr && echo "nr 注释完成"
    else
        echo "nr 注释失败"
    fi
}

### cog
# emappey.py youhuma 45分钟，运行完需要 conda deactivate
exec_cog() {
    mkdir ${annotation}
    cd ${annotation} || exit
    echo "正在切换到 python27 conda 环境"
    source /home/train/miniconda3/bin/activate python27
    check_conda_env
    emapper.py -i ${assemble_trinity}/${specie}_unigene.fasta \
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
exec_kegg() {
    # 判断 $specie_unigene.fasta 是否大于 50000 条，如果大于 50000 条，则需要分批生成 keg 文件
    if [[ $(grep -c '>' ${assemble_trinity}/${specie}_unigene.fasta) -gt 50000 ]]; then
        echo "unigene.fasta 文件大于 50000 条"
        split -l 100000 ${assemble_trinity}/${specie}_unigene.fasta ${annotation}/${specie}_unigene.fasta_
        for i in $(ls ${annotation} | grep "${specie}_unigene.fasta_"); do
            echo "正在生成 ${i} 的 keg 文件"
            python ${python_script}/annotation/kegg_annotation.py \
            -f ${annotation}/${i} \
            -o ${annotation}/${i}_keg \
            -l "${kegg_org}" && echo "${i} 的 keg 文件已生成"
            # 检查文件是否生成，如果生成则把临时分割的文件删除
            if [[ -f ${annotation}/${i}_keg ]]; then
                cat ${annotation}/${i}_keg >> ${annotation}/${specie}_unigene.keg
                rm ${annotation}/${i}_keg
                rm ${annotation}/${i}
            elif [[ ! -f ${annotation}/${i}_keg ]]; then
                echo "${i} 的 keg 文件生成失败"
            fi
        done
    else
        echo "unigene.fasta 文件小于 50000 条"
        python ${python_script}/annotation/kegg_annotation.py \
            -f ${assemble_trinity}/${specie}_unigene.fasta \
            -o ${annotation}/${specie}_unigene.keg \
            -l "${kegg_org}" && echo "${specie}_unigene.fasta 的 keg 文件已生成"
    fi

    # 拿着 $specie_unigene.fasta 去 kegg 网站生成 keg 文件，再做下面的东西 -t plant/animal 动物或植物
    cd ${annotation} || exit
    python ${python_script}/annotation/kegg.py -t ${specie_type} -i ${word_dir}/${specie}_all_gene_id.txt
    cd ${work_dir} || exit
}

### transdecoder
transdecoder() {
    # 也可以自动生成到报告里
    cd ${annotation} || exit
    mkdir transdecoder
    cd transdecoder || exit
    TransDecoder.LongOrfs -t ${assemble_trinity}/${specie}_unigene.fasta
    # 先运行完上面那句，才能运行下面这句
    TransDecoder.Predict -t ${assemble_trinity}/${specie}_unigene.fasta
    cd ${work_dir} || exit
}
exec_annotation() {
    mkdir ${annotation}
    cp ${assemble_trinity}/${specie}_unigene.fasta ${annotation}
    exec_swiss
    exec_nr
    exec_cog
    transdecoder
}

### rsem
rsem() {
    mkdir ${reference}
    # 建库
    cd ${reference} || exit
    echo "正在建库 ${assemble_trinity}/${specie}_unigene.fasta ${specie}"
    sleep 60
    rsem-prepare-reference --bowtie2 ${assemble_trinity}/${specie}_unigene.fasta ${specie}
    cd ${work_dir} || exit
    mkdir ${mapping}
    # 读取 samples_described.txt 循环比对
    tail -n +2 ${work_dir}/samples_described.txt | grep -v '^$' | while read line; do
        sample=$(echo $line | awk '{print $2}')
        sample1=$sample$(ls ${pinjiedata} | grep ${sample} | awk -F"${sample}" '{print $2}' | grep '1')
        sample2=$sample$(ls ${pinjiedata} | grep ${sample} | awk -F"${sample}" '{print $2}' | grep '2')
        echo "正在执行 $sample 的 rsem 流程, $sample1 $sample2 正在比对"
        rsem-calculate-expression -p ${num_threads} \
            --bowtie2 \
            --strandedness reverse \
            --paired-end \
            ${pinjiedata}/${sample1} ${pinjiedata}/${sample2} \
            ${reference}/${specie} ${mapping}/${sample} \
            1>${mapping}/${sample}_mapping.txt \
            2>${mapping}/${sample}.log
        echo "$sample 的 rsem 比对流程已完成"
    done

    # 使用 base 的 python3 环境运行 python 脚本
    source /home/train/miniconda3/bin/activate base
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
        -s ${work_dir}/samples_described.txt \
        -f ${mapping}/reads_matrix_filtered.txt
    python ${python_script}/realignment_fpkm_columns.py -r \
        -s ${work_dir}/samples_described.txt \
        -f ${mapping}/fpkm_matrix_filtered.txt
    python ${python_script}/realignment_fpkm_columns.py -r \
        -s ${work_dir}/samples_described.txt \
        -f ${mapping}/gene_count_matrix.txt
    python ${python_script}/realignment_fpkm_columns.py -r \
        -s ${work_dir}/samples_described.txt \
        -f ${mapping}/gene_fpkm_matrix.txt
    # 合并 fpkm 和 reads 矩阵
    kegg_gene_def=$(realpath -s  ${annotation}/*KEGG_gene_def.txt)
    nr_gene_def=$(realpath -s ${annotation}/*nr_gene_def.txt)
    swiss_gene_def=$(realpath -s ${annotation}/*swiss_gene_def.txt)
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
    cp ${mapping}/reads_matrix_filtered.txt ${multideseq}/
    cp ${mapping}/fpkm_matrix_filtered.txt ${multideseq}/
    cp ${work_dir}/samples_described.txt ${multideseq}/
    cp ${work_dir}/compare_info.txt ${multideseq}/
    echo "正在执行 multi_deseq 流程，Rlog number 为 ${rlog_number}"
    Rscript /home/data/command/Reference_transcriptome/multiple_samples_DESeq2_V7.r ${rlog_number}

    kegg_gene_def=$(realpath -s  ${annotation}/*KEGG_gene_def.txt)
    nr_gene_def=$(realpath -s ${annotation}/*nr_gene_def.txt)
    swiss_gene_def=$(realpath -s ${annotation}/*swiss_gene_def.txt)
    cd ${multideseq} || exit
    python ${python_script}/de_results_add_def.py \
        -k "$kegg_gene_def" \
        -n "$nr_gene_def" \
        -s "$swiss_gene_def"
    cat DEG_summary.txt
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
# unigene_fasta=

### 读取输入参数
# 默认参数
max_memory=500
num_threads=60
rlog_number=1

run_program() {
    case "$1" in
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
        echo "无效的选项：$1"
        exit 0
        ;;
    esac
}

help_info() {
    echo "帮助信息：
    ————————————————————————————————————————————————————————
    0 执行所有流程

    2 执行 fastqc
    3 执行 pinjie（需指定 --specie, --threads, --max-memory）
    5 执行 Annotation（需指定 --specie, --threads）
        5.1 执行 swiss（需指定 --specie, --threads)
        5.2 执行 nr（需指定 --specie)
        5.3 执行 cog（需指定 --specie，--threads）
        5.4 执行 kegg（需指定 --specie, --specie-type [plant, animal], --kegg-org）
        5.5 执行 transdecoder（需指定 --specie）
    6 执行 rsem（需指定 --specie, --threads）
    7 执行 multi_deseq（需指定 --rlog-number）
    8 整理交付目录
    ————————————————————————————————————————————————————————
    使用方法：./noref.sh [参数]
    参数：
        -r, --run <run> 指定运行的流程 (指上面的数字，例如 -r 5)
        -s, --specie <specie> 指定物种
        -c, --threads <threads> 指定线程数
        -m, --max-memory <max_memory> 指定最大内存
        --rlog-number <rlog_number> multi deseq 的倍数
        --trinity-type <trinity_type> 设置 trinity 的类型 all/planA/custom
            all：全部拼接
            planA：只拼接每组中最长的
            custom：自定义，生成文件后手动修改
        -h 显示帮助信息
"
}
while [[ $# -gt 0 ]]; do
    case "$1" in
        -r | --run)
            run="$2"
            echo "run: $2"
            shift 2
            ;;
        -s | --specie)
            specie="$2"
            echo "specie: $2"
            shift 2
            ;;
        -c | --threads)
            num_threads="$2"
            echo "num_threads: $2"
            shift 2
            ;;
        -m | --max-memory)
            max_memory="$2"
            echo "max_memory: $2"
            shift 2
            ;;
        --rlog-number)
            rlog_number="$2"
            echo "rlog_number: $2"
            shift 2
            ;;
        --trinity-type)
            trinity_type="$2"
            echo "trinity_type: $2"
            shift 2
            ;;
        --kegg-org)
            kegg_org="$2"
            echo "kegg_org: $2"
            shift 2
            ;;
        --specie-type)
            specie_type="$2"
            echo "specie_type: $2"
            shift 2
            ;;
        -h)
            help_info
            exit 0
            ;;
        *)
            echo "无效的选项：$1"
            exit 1
            ;;
    esac
done
if [[ -z $run ]]; then
    echo "未指定运行流程，使用 -h 查看帮助"
    exit 0
fi
run_program $run
