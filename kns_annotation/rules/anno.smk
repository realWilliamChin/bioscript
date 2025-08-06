rule kegg_anno:
    input:
        all_geneid_file = f'{ALL_GENEID}',
        anno_fasta = f'{ANNO_FASTA}'
    output:
        ko_file = f'{KEGG_DIR}/{SPECIE_NAME}_ko.txt',
        clean_file = f'{KEGG_DIR}/{SPECIE_NAME}_KEGG_clean.txt',
        gene_def = f'{KEGG_DIR}/{SPECIE_NAME}_KEGG_gene_def.txt',
        kegg_tier2 = f'{KEGG_DIR}/{SPECIE_NAME}_KEGG_tier2.txt',
        kegg_tier3 = f'{KEGG_DIR}/{SPECIE_NAME}_KEGG.txt'
    priority: 10
    conda: "python310"
    params:
        kegg_org_list = config.get('kegg_anno_org_list', ''),
        split_anno = config.get('split_anno', 1),
        fasta_type = SPECIE_TYPE
    shell:
        """
        python /home/colddata/qinqiang/script/kns_annotation/scripts/kegg_annotation.py \
            -f {input.anno_fasta} \
            -k {output.ko_file} \
            -t {params.fasta_type} \
            -l {params.kegg_org_list} \
            -s {params.split_anno} \
            -o {KEGG_DIR}/{SPECIE_NAME}
        """
    

# nr annotation
rule swiss_anno:
    input:
        f'{ANNO_FASTA}'
    output:
        blast_file = f'{SWISS_DIR}/{SPECIE_NAME}_swiss.blast',
        gene_def = f'{SWISS_DIR}/{SPECIE_NAME}_swiss_gene_def.txt',
        go_bp = f'{SWISS_DIR}/{SPECIE_NAME}_GO_BP_ID.txt',
        go_cc = f'{SWISS_DIR}/{SPECIE_NAME}_GO_CC_ID.txt',
        go_mf = f'{SWISS_DIR}/{SPECIE_NAME}_GO_MF_ID.txt',
        swiss_goid = f'{SWISS_DIR}/{SPECIE_NAME}_swiss_idNo_def.txt',
    conda: "python310"
    params:
        output_prefix = f'{SWISS_DIR}/{SPECIE_NAME}'
    threads: lambda wildcards: max(1, workflow.cores // 3)
    shell:
        """
        python /home/colddata/qinqiang/script/kns_annotation/scripts/swiss.py \
            --fasta {input} \
            -b {output.blast_file} \
            -o {params.output_prefix} \
            -t {threads}
        """


rule nr_anno:
    input:
        f'{ANNO_FASTA}'
    output:
        blast_file = f'{NR_DIR}/{SPECIE_NAME}_nr.blast',
        gene_def = f'{NR_DIR}/{SPECIE_NAME}_nr_gene_def.txt'
    conda: "python310"
    params:
        output_prefix = f'{NR_DIR}/{SPECIE_NAME}'
    shell:
        """
        python /home/colddata/qinqiang/script/kns_annotation/scripts/nr.py \
            --fasta {input} \
            -b {output.blast_file} \
            -o {params.output_prefix}
        """


rule cog_anno:
    input:
        f'{ANNO_FASTA}'
    output:
        annotations = f'{COG_DIR}/{SPECIE_NAME}.emapper.annotations',
        seed_orthologs = f'{COG_DIR}/{SPECIE_NAME}.emapper.seed_orthologs'
    conda: "eggnog2"
    threads: lambda wildcards: max(1, workflow.cores // 3)
    shell:
        """
        python /home/train/miniconda3/envs/eggnog2/bin/emapper.py \
            -i {input} \
            -o {COG_DIR}/{SPECIE_NAME} \
            --cpu {threads} \
            --data_dir /opt/biosoft/eggnog5.0.0 \
            --dmnd_db /opt/biosoft/eggnog5.0.0/eggnog_proteins.dmnd \
            -m diamond \
            --translate \
            --itype genome \
            --genepred search
        """


rule biogrid:
    input:
        f'{ANNO_FASTA}'
    output:
        f'{BIOGRID_DIR}/Biogrid_PPI_relation_from_{SPECIE_NAME}.txt'
    conda: 'python310'
    threads: 5
    params:
        genome_type = SPECIE_TYPE
    shell:
        """
        python /home/colddata/qinqiang/script/kns_annotation/scripts/biogrid.py \
            -f {input} \
            -d {params.genome_type} \
            -c {threads} \
            -p {BIOGRID_DIR}/{SPECIE_NAME}
        """


# rule transdecoder:
#     input:
#         f'{ANNO_FASTA}'
#     output:
#         orfs = f'{TRANSDECODER_DIR}/{SPECIE_NAME}.transdecoder_dir/longest_orfs.pep',
#         cds = f'{TRANSDECODER_DIR}/{SPECIE_NAME}.transdecoder_dir/longest_orfs.cds',
#         gff3 = f'{TRANSDECODER_DIR}/{SPECIE_NAME}.transdecoder_dir/longest_orfs.gff3',
#         bed = f'{TRANSDECODER_DIR}/{SPECIE_NAME}.transdecoder_dir/longest_orfs.bed'
#     conda: 'base'
#     threads: lambda wildcards: max(1, workflow.cores // 2)
#     shell:
#         """
#         TransDecoder.LongOrfs -t {input} -O {TRANSDECODER_DIR}/{SPECIE_NAME}.transdecoder_dir --gene_trans_map {input}.trans_map
#         TransDecoder.Predict -t {input} -O {TRANSDECODER_DIR}/{SPECIE_NAME}.transdecoder_dir \
#             --cpu {threads} \
#             --single_best_only \
#             --retain_long_orfs 1000 \
#             --retain_pfam_hits
#         """