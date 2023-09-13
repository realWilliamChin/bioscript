#!/bin/bash
# 0 执行所有流程
# 
# 1 执行 fastqc
# 2 执行 pinjie
# 3 执行 Annotation
# 3.1 执行 swiss
# 3.2 执行 nr
# 3.3 执行 cog
# 3.4 执行 kegg
# 3.5 执行 transdecoder
# 3.6 执行 biogrid 程序
# 3.7 执行 annotation report
# 4 执行 rsem
# 5 执行 multi_deseq
# 7 整理交付目录

# ！！！ 此脚本需在项目根目录运行
# 选择运行流程
run=0
# 指定物种，例如 Medicago_truncatula
specie=zihuamuxu
# 指定线程数
threads=60
# 指定最大内存
max_memory=500
# 指定 Trinity 的类型 [all(所有), max(每组中最大的)]
trinity_type=max
# 指定物种类型，例如 [plant, animal]
specie_type=plant
# kegg 注释的物种，例如 ath, osa
kegg_org="ath, osa"
# multi deseq 的倍数
rlog_number=1
