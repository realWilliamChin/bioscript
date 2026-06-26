#!/bin/bash
# 【项目模板脚本】本脚本中的 VCF_DIR 是按具体项目硬编码的客户数据路径，
# 每次新项目使用时需要修改 VCF_DIR 与样本目录匹配模式(IonXpress_*)等，
# 不是一个通用 CLI 工具。

# 处理目录
VCF_DIR="/home/colddata/qinqiang/Project/2026_04_15_安贞医院_外显子/01_数据分析/IonTorrent/01_Variant_caller_result"
BCFTOOLS="/home/data/opt/biosoft/bcftools-1.9/bin/bcftools"
TABIX="/home/data/opt/biosoft/htslib-1.9/bin/tabix"

# 遍历所有IonXpress目录
for sample_dir in ${VCF_DIR}/IonXpress_*; do
    if [ -d "${sample_dir}" ]; then
        sample_name=$(basename ${sample_dir})
        vcf_file="${sample_dir}/lifted_over.vcf"
        
        if [ -f "${vcf_file}" ]; then
            echo "正在处理样本: ${sample_name}"
            
            # 备份原始文件
            if [ ! -f "${vcf_file}.bak" ]; then
                cp ${vcf_file} ${vcf_file}.bak
                echo "已备份原始文件到 ${vcf_file}.bak"
            fi
            
            # 修改样本名
            ${BCFTOOLS} reheader -s <(echo ${sample_name}) ${vcf_file}.bak > ${vcf_file}
            
            # 压缩VCF
            ${BCFTOOLS} view -Oz -o ${vcf_file}.gz ${vcf_file}
            
            # 生成索引
            ${TABIX} -p vcf ${vcf_file}.gz
            
            echo "样本 ${sample_name} 处理完成，生成文件: ${vcf_file}.gz ${vcf_file}.gz.tbi"
            echo "----------------------------------------"
        else
            echo "警告: 样本 ${sample_name} 下不存在 lifted_over.vcf 文件，跳过"
        fi
    fi
done

echo "所有样本处理完成！"