#!/home/train/miniconda3/bin/python
# -*- coding: utf-8 -*-
# Created Time  : 2025/07/28 13:25
# Author        : William GoGo
import os

sample = config.get("sample")
input_file = config.get("input_file")
output_dir = config.get("output_dir", os.getcwd())
humandb_path = config.get("humandb_path", "/home/data/opt/biosoft/annovar/humandb")
buildver = config.get("buildver", "hg38")
protocols = config.get('protocols', 'refGene,gnomad41_exome,clinvar_20240917')
operations = config.get('operations', 'g,f,f')
arguments = config.get('arguments', ',,-infoasscore')

rule all:
    input:
        expand("{output_dir}/{sample}.hg38_multianno.txt", 
               output_dir=output_dir, 
               sample=sample)

rule convert_to_annovar:
    input:
        vcf = input_file
    output:
        avinput = os.path.join(output_dir, '{sample}.avinput')
    shell:
        """
        perl /home/data/opt/biosoft/annovar/convert2annovar.pl \
            -format vcf4 \
            -allsample \
            -includeinfo \
            -withfreq \
            {input.vcf} \
            -outfile {output.avinput}
        """

rule annotate_variants:
    input:
        avinput = rules.convert_to_annovar.output.avinput
    output:
        annotated = os.path.join(output_dir, '{sample}.hg38_multianno.txt')
    shell:
        """
        perl /home/data/opt/biosoft/annovar/table_annovar.pl \
            {input.avinput} \
            {humandb_path} \
            --buildver {buildver} \
            --out {output_dir}/{wildcards.sample} \
            --protocol {protocols} \
            --operation {operations} \
            --argument {arguments} \
            --polish \
            --intronhgvs 100 \
            --otherinfo \
            --remove
        """