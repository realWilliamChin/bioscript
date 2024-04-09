#!/bin/bash
find ./ -name "*1.fastq.gz" >filelist.txt
find ./ -name "*1.fq.gz" >> filelist.txt

cleandata_dir="./02_Cleandata"
report_dir="./fastp_report"
mkdir $cleandata_dir
mkdir $report_dir

while IFS= read -r line; do
    f1=$line
    f2=$(echo $line | sed 's/1.fastq.gz/2.fastq.gz/')
    key_name=$(echo ${f1##*/})
    key_name=$(echo ${key_name} | sed 's/_1.fq.gz//')
    echo "fastp -G -q 5 -u 50 -A -n 15 -l 150 \
        --overlap_diff_limit 1 \
        --overlap_diff_percent_limit 10 \
        -i ${f1} -I ${f2} \
        -o ${cleandata_dir}/${f1##*/} \
        -O ${cleandata_dir}/${f2##*/} \
        -j ${report_dir}/${f1##*/}_fastp.json &" | /bin/bash
done <filelist.txt

cd ${report_dir}
rawdata2clean_json_result.py
