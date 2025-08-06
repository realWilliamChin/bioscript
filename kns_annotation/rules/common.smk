import os
import pandas as pd


# 定义工作目录
WORK_DIR = os.getcwd()

# 加载配置文件
configfile: os.path.join(WORK_DIR, 'config.yaml')

# 定义各目录路径
KEGG_DIR = os.path.join(WORK_DIR, config.get('kegg_dir', ''))
NR_DIR = os.path.join(WORK_DIR, config.get('nr_dir', ''))
SWISS_DIR = os.path.join(WORK_DIR, config.get('swiss_dir', ''))
COG_DIR = os.path.join(WORK_DIR, config.get('cog_dir', ''))
BIOGRID_DIR = os.path.join(WORK_DIR, config.get('biogrid_dir', ''))
GSEA_DIR = os.path.join(WORK_DIR, config.get('gsea_dir', ''))

SPECIE_NAME = config.get('specie', '')
SPECIE_TYPE = config.get('specie_type', '')

ANNO_FASTA = config.get('anno_fasta', '')
ALL_GENEID = config.get('all_geneid_file', '').format(specie=SPECIE_NAME)
BASICINFO_FILE = config.get('basicinfo_file', '')

# 条件函数
def get_kegg_outputs():
    if config.get('kegg_dir', '') != '':
        return [
            f'{KEGG_DIR}/{SPECIE_NAME}_KEGG_clean.txt',
            f'{KEGG_DIR}/{SPECIE_NAME}_KEGG_gene_def.txt'
        ]
    return []

def get_swiss_outputs():
    if config.get('swiss_dir', '') != '':
        return [
            f'{SWISS_DIR}/{SPECIE_NAME}_GO_BP_ID.txt',
            f'{SWISS_DIR}/{SPECIE_NAME}_GO_CC_ID.txt',
            f'{SWISS_DIR}/{SPECIE_NAME}_GO_MF_ID.txt',
            f'{SWISS_DIR}/{SPECIE_NAME}_swiss.blast',
            f'{SWISS_DIR}/{SPECIE_NAME}_swiss_gene_def.txt'
        ]
    return []

def get_nr_outputs():
    if config.get('nr_dir', '') != '':
        return [
            f'{NR_DIR}/{SPECIE_NAME}_nr.blast',
            f'{NR_DIR}/{SPECIE_NAME}_nr_gene_def.txt'
        ]
    return []

def get_cog_outputs():
    if config.get('cog_dir', '') != '':
        return [
            f'{COG_DIR}/{SPECIE_NAME}.emapper.seed_orthologs'
        ]
    return []

def get_biogrid_outputs():
    if config.get('biogrid_dir', '') != '':
        return [
            f'{BIOGRID_DIR}/Biogrid_PPI_relation_from_{SPECIE_NAME}.txt'
        ]
    return []

def get_gsea_outputs():
    if config.get('gsea_dir', '') != '':
        return [
            f'{GSEA_DIR}/{SPECIE_NAME}_GO_BP_ID.GMT',
            f'{GSEA_DIR}/{SPECIE_NAME}_GO_CC_ID.GMT',
            f'{GSEA_DIR}/{SPECIE_NAME}_GO_MF_ID.GMT',
            f'{GSEA_DIR}/{SPECIE_NAME}_KEGG_tier2.GMT',
            f'{GSEA_DIR}/{SPECIE_NAME}_KEGG_tier3.GMT'
        ]
    return []

def get_kns_outputs():
    if all(not config.get(dir_key, '') for dir_key in ['kegg_dir', 'swiss_dir', 'nr_dir']):
        return []
    return [f'{WORK_DIR}/{SPECIE_NAME}_kns_gene_def.txt']

def get_annotation_summary_outputs():
    if all(not config.get(dir_key, '') for dir_key in ['kegg_dir', 'swiss_dir', 'nr_dir', 'cog_dir']):
        return []
    return [
        f'{WORK_DIR}/{SPECIE_NAME}_annotation_summary_venn.jpeg',
        f'{WORK_DIR}/{SPECIE_NAME}_annotation_summary.txt'
    ]

# 定义最终输出文件
rule all_outputs:
    input:
        get_kegg_outputs(),
        get_nr_outputs(),
        get_swiss_outputs(),
        get_kns_outputs(),
        get_cog_outputs(),
        get_biogrid_outputs(),
        # get_transdecoder_outputs(),
        get_gsea_outputs(),
        get_annotation_summary_outputs()

rule gsea:
    input:
        go_bp = f'{SWISS_DIR}/{SPECIE_NAME}_GO_BP_ID.txt',
        go_cc = f'{SWISS_DIR}/{SPECIE_NAME}_GO_CC_ID.txt',
        go_mf = f'{SWISS_DIR}/{SPECIE_NAME}_GO_MF_ID.txt',
        kegg_tier2 = f'{KEGG_DIR}/{SPECIE_NAME}_KEGG_tier2.txt',
        kegg_tier3 = f'{KEGG_DIR}/{SPECIE_NAME}_KEGG.txt'
    output:
        go_bp = f'{GSEA_DIR}/{SPECIE_NAME}_GO_BP_ID.GMT',
        go_cc = f'{GSEA_DIR}/{SPECIE_NAME}_GO_CC_ID.GMT',
        go_mf = f'{GSEA_DIR}/{SPECIE_NAME}_GO_MF_ID.GMT',
        kegg_tier2 = f'{GSEA_DIR}/{SPECIE_NAME}_KEGG_tier2.GMT',
        kegg_tier3 = f'{GSEA_DIR}/{SPECIE_NAME}_KEGG_tier3.GMT'
    conda: 'python310'
    shell:
        """
        mkdir -p {GSEA_DIR}
        python /home/colddata/qinqiang/script/kns_annotation/scripts/gsea.py \
            --gobp {input.go_bp} \
            --gocc {input.go_cc} \
            --gomf {input.go_mf} \
            --kegg-tier2 {input.kegg_tier2} \
            --kegg-tier3 {input.kegg_tier3} \
            --output-dir {GSEA_DIR}
        """
    

rule annotation_summary_venn_plot:
    input:
        anno_fasta = f'{ANNO_FASTA}',
        swiss_goid = f'{SWISS_DIR}/{SPECIE_NAME}_swiss_idNo_def.txt',
        swiss_gene_def = f'{SWISS_DIR}/{SPECIE_NAME}_swiss_gene_def.txt',
        kegg_gene_def = f'{KEGG_DIR}/{SPECIE_NAME}_KEGG_gene_def.txt',
        nr_gene_def = f'{NR_DIR}/{SPECIE_NAME}_nr_gene_def.txt'
    output:
        output_pic = f'{WORK_DIR}/{SPECIE_NAME}_annotation_summary_venn.jpeg',
        output_summary_file = f'{WORK_DIR}/{SPECIE_NAME}_annotation_summary.txt'
    conda: 'python310'
    shell:
        """
        bash /home/colddata/qinqiang/script/kns_annotation/scripts/annotation_summary.sh \
            {input.anno_fasta} \
            {input.swiss_goid} \
            {input.swiss_gene_def} \
            {input.kegg_gene_def} \
            {input.nr_gene_def} \
            {output.output_pic} \
            {output.output_summary_file}
        """


rule kns_def_merge:
    input:
        all_geneid = f'{ALL_GENEID}',
        basicinfo_file = f'{BASICINFO_FILE}',
        swiss_gene_def = f'{SWISS_DIR}/{SPECIE_NAME}_swiss_gene_def.txt',
        kegg_gene_def = f'{KEGG_DIR}/{SPECIE_NAME}_KEGG_gene_def.txt',
        nr_gene_def = f'{NR_DIR}/{SPECIE_NAME}_nr_gene_def.txt'
    output:
        all_geneid_def = f'{WORK_DIR}/{SPECIE_NAME}_kns_gene_def.txt'
    conda: 'python310'
    shell:
        """
        if [ -s {input.basicinfo_file} ]; then
            input_file={input.basicinfo_file}
            header_param=""
        else
            input_file={input.all_geneid}
            header_param="--input-header 0"
        fi
        
        python /home/colddata/qinqiang/script/transcriptome/genedf_add_expression_and_def.py \
            -i $input_file \
            $header_param \
            -s {input.swiss_gene_def} \
            -k {input.kegg_gene_def} \
            -n {input.nr_gene_def} \
            -o {output.all_geneid_def}
        """
        
