#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2023/07/14 16:50
# Author        : William GoGo
import time
import re
import os, sys
import argparse
import requests
import poplib
import pandas as pd
from Bio import SeqIO
from lxml import etree
from fake_useragent import UserAgent
from email import encoders
from email.header import Header, decode_header
from email.mime.text import MIMEText
from email.parser import Parser
from email.utils import parseaddr, formataddr

import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger

sys.path.append('/home/colddata/qinqiang/script/CommonTools/')
sys.path.append('/home/colddata/qinqiang/script/transcriptome/')
from merge_fpkm_reads_matrix import merge_fpkm_reads
from load_input import load_table, write_output_df

class EmailReceiver:
    def __init__(self, server, port):
        self.server = server
        self.port = port
    
    
    def guess_charset(self, msg):
        charset = msg.get_charset()
        if charset is None:
            content_type = msg.get('Content-Type', '').lower()
            pos = content_type.find('charset=')
            if pos >= 0:
                charset = content_type[pos + 8:].strip()
        return charset


    def decode_str(self, s):
        value, charset = decode_header(s)[0]
        if charset:
            value = value.decode(charset)
        return value


    def receive_mail(self, user, pswd):
        server = poplib.POP3(self.server, self.port)
        max_retry = 0
        while True:
            try:
                server.user(user)
                server.pass_(pswd)
                resp, mails, octets = server.list()
                if len(mails) == 0:
                    raise RuntimeError('No mail found from eamil address!')
                break
            except Exception as e:
                if max_retry <= 3:
                    logger.warning(f'获取邮件错误，10s 后尝试重新获取')
                    time.sleep(10)
                    max_retry += 1
                else:
                    logger.critical(f'尝试 3 次失败')
                    raise RuntimeError(f'Error message: {e}')
        
        # 从最近 3 封或 5 封邮件选取 @genome 的邮件，以防收取到其他邮件
        msg = None
        for i in [0, -1, -2]:
            resp, lines, octets = server.retr(len(mails) - i)
            msg = Parser().parsestr(b'\r\n'.join(lines).decode('utf-8'))
            if "genome.jp" in msg.get('From'):
                break
        # header_from = msg.get('From', '')
        # mail_text = self.decode_str(msg.get('Subject', ''))
        msg_date = msg.get('Date', '')
        result = self.get_info(msg)
        result['date'] = msg_date
        server.quit()
        
        return result


    def get_info(self, msg):
        result = {}
        for header in ['From', 'To', 'Subject']:
            value = msg.get(header)
            if not value:
                continue
            
            if header == 'Subject':
                value = self.decode_str(value)
            else:
                hdr, addr = parseaddr(value)
                name = self.decode_str(hdr)
                value = f'{name} <{addr}>'
            # print(f'{header}: {value}')
            result[header.lower()] = value.strip()
        # print(msg.is_multipart())
        if msg.is_multipart():
            parts = msg.get_payload()
            for part in parts:
                if part.get_content_type() == 'text/plain':
                    # self.print_info(part)
                    content = part.get_payload(decode=True)
                    charset = self.guess_charset(part)
                    if charset:
                        content = content.decode(charset)
                    # print(f'Text: {content}...')
                    result['text'] = content
         
        else:
            content_type = msg.get_content_type()
            # print(content_type)
            if content_type == 'text/plain' or content_type == 'text/html':
                content = msg.get_payload(decode=True)
                charset = self.guess_charset(msg)
                if charset:
                    content = content.decode(charset)
                    # print(f'Text: {content}...')
                    result['text'] = content.strip()
        return result
    
    
class KeggUploadError(Exception):
    """KEGG 上传相关的异常"""
    pass

def upload_fasta_to_kegg(fasta_file, org_lst, eamil_address):
    fasta_file_name = fasta_file.split('/')[-1]
    # step.1 读取文件为上传的格式
    with open(fasta_file, 'r') as fasta_file:
        fasta_content = '\r\n'.join(fasta_file.read().strip().split('\n'))
    upload_key_name = fasta_file_name.split('.')[0]
    upload_file_url = 'https://www.genome.jp/kaas-bin/kaas_main'
    boundary = '----WebKitFormBoundaryu7jixhGFRtH1qvAA'
    headers = {
        'Host': 'www.genome.jp',
        'Connection': 'keep-alive',
        'Cache-Control': 'max-age=0',
        'sec-ch-ua': '"Not.A/Brand";v="8", "Chromium";v="114", "Microsoft Edge";v="114"',
        'sec-ch-ua-mobile': '?0',
        'Origin': 'https://www.genome.jp',
        'Content-Type': f'multipart/form-data; boundary={boundary}',
        'User-Agent': "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36 Edg/120.0.0.0",
        'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3;q=0.7',
        'Sec-Fetch-Site': 'same-origin',
        'Sec-Fetch-Mode': 'navigate',
        'Sec-Fetch-Dest': 'document',
        'Referer': 'https://www.genome.jp/kaas-bin/kaas_main'
    }
    boundary = '--' + boundary
    payload = '\r\n'.join([boundary,
                           'Content-Disposition: form-data; name="continue"',
                           '',
                           '1',
                            boundary,
                            'Content-Disposition: form-data; name="prog"',
                            '',
                            'BLAST',
                            boundary,
                            'Content-Disposition: form-data; name="text"',
                            '',
                            '',
                            boundary,
                            'Content-Disposition: form-data; name="uptype"',
                            '',
                            'q_file',
                            boundary,
                            'Content-Disposition: form-data; name="peptide2"',
                            '',
                            'n',
                            boundary,
                            f'Content-Disposition: form-data; name="file"; filename="{fasta_file_name}"',
                            'Content-Type: application/octet-stream',
                            '',
                            f'{fasta_content}',
                            '',
                            boundary,
                            'Content-Disposition: form-data; name="qname"',
                            '',
                            f'{upload_key_name}',
                            boundary,
                            'Content-Disposition: form-data; name="mail"',
                            '',
                            f'{eamil_address}',
                            boundary,
                            'Content-Disposition: form-data; name="dbmode"',
                            '',
                            'manual',
                            boundary,
                            'Content-Disposition: form-data; name="org_list"',
                            '',
                            f'{org_lst}',
                            boundary,
                            'Content-Disposition: form-data; name="way"',
                            '',
                            'b',
                            boundary,
                            'Content-Disposition: form-data; name="mode"',
                            '',
                            'compute',
                            boundary,
                            ''
                           ]).strip("'")
    
    # upload file
    logger.info(f'uploading {fasta_file_name}...')
    response = requests.post(upload_file_url, headers=headers, data=payload)
    if response.status_code != 200:
        raise KeggUploadError(f'上传失败 {fasta_file_name}! 状态码: {response.status_code}\n{response.text}')
    
    html = etree.HTML(response.text)
    
    # check upload status
    sub_element = html.xpath('//div[@id="main"]/p[1]/text()')[0]
    right_sub = 'Job Request'
    
    if sub_element == right_sub:
        logger.success('upload success!')
    elif 'Sorry' in sub_element:
        raise KeggUploadError(f'上传失败!\n错误信息: {sub_element}')
    else:
        raise KeggUploadError(f'上传失败!\n{response.text}')
    

def email_link_click(link):
    headers = {
        'Host': 'www.genome.jp',
        'User-Agent': "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36 Edg/120.0.0.0",
    }
    max_retry = 0
    while max_retry <= 5:
        resp = requests.get(link, headers=headers)
        if resp.status_code == 200:
            return resp.text
        else:
            max_retry += 1
            time.sleep(3)
    sys.exit(1)


def kegg_anno(mail_type, username, password, fasta_file, org_lst: str, output_file):
    logger.info(f"fasta file: {fasta_file}\norg_list：{org_lst}")
    mail = EmailReceiver(f'imap.{mail_type}.com', 110)
    # step.1 submit fasta file to kegg and get request_ID and job status
    mail_content = mail.receive_mail(username, password)
    # 用来判断是否是同一个邮件
    mail_date = mail_content['date']

    try:
        upload_fasta_to_kegg(fasta_file, org_lst, username)
    except KeggUploadError as e:
        logger.critical(f"上传失败: {str(e)}")
        raise

    time.sleep(5)
    job_ID = None
    job_key = None
    for t in range(0, 30):
        mail_content = mail.receive_mail(username, password)
        if 'Job request accepted' in mail_content['subject'] and mail_content['from'] == '<noreply@genome.jp>':
            submit_link = re.findall(r"https://.*submit.*", mail_content['text'])[0]
            job_ID = submit_link.split('&')[-2].split('=')[-1]
            job_key = submit_link.split('&')[-1].split('=')[-1]
            submit_resp = email_link_click(submit_link)
            submit_status = etree.HTML(submit_resp).xpath('//h1/text()')[0]
            if submit_status == 'Submitted':
                logger.info(f'Job requesting ID:{job_ID} ...')
                mail_date = mail_content['date']
                break
            else:
                logger.critical(f'Job request Failed! not Submitted ID:{job_ID}!')
                raise KeggUploadError(f'任务提交失败! ID:{job_ID}')
        elif mail_content['date'] == mail_date:
            logger.info(f'receiving job request email {t+1} times ...')
            time.sleep(5)
            t += 1
        else:
            logger.critical(f'Job request Failed!\n{mail_content}')
            raise KeggUploadError('任务请求失败!')
            
    # step.2 获取任务提交状态
    status_link = None
    status_ID = None
    time.sleep(30)
    for t in range(0, 12):
        mail_content = mail.receive_mail(username, password)
        if 'KAAS - Accepted' == mail_content['subject'] and mail_content['from'] == '<kaas@genome.jp>':
            logger.info(f'Job Accepted ID:{job_ID}!')
            status_link = re.findall(r"https://.*user.*", mail_content['text'])[0]
            status_ID = status_link.split('&')[-2].split('=')[-1]
            mail_date = mail_content['date']
            break
        elif mail_content['date'] == mail_date:
            logger.info(f'get job status {t+1} times')
            time.sleep(30)
            t += 1
        else:
            logger.critical(f'Job Accept Failed!\n{mail_content}')
            raise KeggUploadError('任务接受失败!')
            
    # step.3 获取任务完成状态
    # 校验是否同一个任务
    if job_ID != status_ID:
        logger.error('request_ID and status_ID are not the same!')
        raise KeggUploadError('任务ID不匹配!')
    
    # 根据序列数量动态计算等待时间
    fasta_seq_num = len([record for record in SeqIO.parse(fasta_file, "fasta")])
    base_wait_time = 60  # 基础等待时间（分钟）
    wait_times = max(int(fasta_seq_num/1000)*base_wait_time, base_wait_time)  # 至少等待基础时间
    logger.info(f'Job is running wait for {wait_times/60} minutes ...')
    
    # 使用指数退避策略进行等待
    max_retries = 60
    base_sleep_time = 300  # 5分钟
    max_sleep_time = 1800  # 30分钟
    current_sleep_time = base_sleep_time
    
    for t in range(1, max_retries):
        mail_content = mail.receive_mail(username, password)
        
        if 'Annotation was completed' in mail_content['subject'] and 'kaas@genome.jp' in mail_content['from']:
            logger.success('Job Complete!')
            
            # you should click this link then you can download the file
            logger.info(f'downloading {output_file} ...')
            visit_link = f'https://www.genome.jp/kaas-bin/kaas_main?mode=brite&id={job_ID}&key={job_key}'
            email_link_click(visit_link)
            time.sleep(5)
            
            # from download link get file
            # download_link = f'https://www.genome.jp/kegg-bin/download_htext?htext=q00001.keg&format=htext&filedir=/tools/kaas/files/log/result/{job_ID}'
            download_link = f'https://www.genome.jp/tools/kaas/files/dl/{job_ID}/query.ko'
            resp = email_link_click(download_link)
            
            # check file line > 1
            logger.info('checking file ...')
            if len(resp.split('\n')) > 1:
                logger.success('checking file done ...')
            else:
                logger.critical(f'get file failed!, file is not correct ID{job_ID},key{job_key}!\n{resp}')
                raise KeggUploadError('获取结果文件失败!')
                
            # output file
            with open(output_file, 'w') as kegg_annotation:
                kegg_annotation.write(resp)
            logger.success(f'download success! {output_file} is ready!')
            
            break
        elif mail_content['date'] == mail_date:
            logger.info(f'Job is running, checked email {t} times({current_sleep_time/60}min) ...')
            if t == max_retries:
                logger.error(f'Job maybe error, Have been waiting {max_retries*base_sleep_time/60} minutes ...')
                raise KeggUploadError('任务执行超时!')
            time.sleep(current_sleep_time)
            # 使用指数退避策略增加等待时间
            current_sleep_time = min(current_sleep_time * 1.5, max_sleep_time)
            t += 1
        else:
            logger.critical(f'Job Failed!\n{mail_content}')
            raise KeggUploadError('任务执行失败!')


def ko03000(kegg_gene_df, kegg_tier3_df):
    ko03000_def_df_file = '/home/colddata/qinqiang/script/kns_annotation/scripts/ko03000_def.txt'
    ko03000_def_df = load_table(ko03000_def_df_file)
    
    # 处理 NA/NaN 值
    kegg_tier3_df['Pathway'] = kegg_tier3_df['Pathway'].fillna('N/A')
    
    ko03000_geneid_list = kegg_tier3_df[kegg_tier3_df['Pathway'].str.contains("ko03000", na=False)]['GeneID'].tolist()
    ko03022_geneid_list = kegg_tier3_df[kegg_tier3_df['Pathway'].str.contains("ko03022", na=False)]['GeneID'].tolist()
    
    ko03000_df = kegg_gene_df[kegg_gene_df['GeneID'].isin(ko03000_geneid_list)][['GeneID', 'KEGG_ID']].copy()
    ko03022_df = kegg_gene_df[kegg_gene_df['GeneID'].isin(ko03022_geneid_list)][['GeneID', 'KEGG_ID']].copy()
    ko03000_df = pd.merge(left=ko03000_df, right=ko03000_def_df, on='KEGG_ID', how='left')
    
    return ko03000_df, ko03022_df


def filter_ko_df(ko_df, specie_type):
    # 植物物种过滤动物的一些注释
    if specie_type == 'plant':
        # 定义需要过滤的类别
        filter_categories = {
            'A09160': 'Human Diseases',
            'A09190': 'Organismal Systems',
            'A09150': {
                'keep_patterns': [
                    '09158:Development and regeneration',
                    '09159:Environmental adaptatqion'
                ]
            }
        }
        
        # 过滤掉完全不需要的类别
        for category in ['A09160', 'A09190']:
            ko_df = ko_df[~ko_df['Category'].str.contains(category, na=False)]
        
        # 处理需要部分保留的类别
        if 'A09150' in filter_categories:
            mask = ko_df['Category'].str.contains('A09150:Organismal Systems', na=False)
            if mask.any():
                keep_pattern = '|'.join(filter_categories['A09150']['keep_patterns'])
                ko_df = ko_df[~mask | (mask & ko_df['Subcategory'].str.contains(keep_pattern, na=False))]
        
        # 删除所有空值行
        ko_df = ko_df.dropna(subset=['Category'])
        
    return ko_df


def parse_keg(ko_file, specie_type, all_id_file, fpkm, reads, output_prefix):
    kegg_db = load_table('/home/colddata/qinqiang/script/kns_annotation/scripts/kegg_db.txt')
    ko_df = load_table(ko_file, header=None, names=['GeneID', 'KO_ID'], dtype=str)
    ko_df = ko_df.dropna(subset=['KO_ID'])
    ko_df = pd.merge(left=ko_df, right=kegg_db, on='KO_ID', how='left')
    ko_df = ko_df[['GeneID', 'KEGG_Pathway', 'Subcategory', 'Category', 'KO_ID', 'Gene_symbol', 'Enzyme_Description']]
    
    # 过滤，保存 KEGG_clean.txt 文件
    kegg_clean_file = output_prefix + 'KEGG_clean.txt'
    ko_df = filter_ko_df(ko_df, specie_type)
    write_output_df(ko_df, kegg_clean_file, header=False, index=False)
                
    # 生成 KEGG_gene_def 文件
    kegg_gene_def_titles = ['GeneID', 'Pathway', 'Level2', 'Level1', 'KEGG_ID', 'Gene_shortname', 'Description EC_number']
    kegg_clean_df = load_table(kegg_clean_file, header=None, names=kegg_gene_def_titles, dtype={"GeneID": str})
    gene_def_df = kegg_clean_df.drop(columns=['Pathway', 'Level2', 'Level1']).copy()
    gene_def_df.drop_duplicates(subset=['GeneID', 'KEGG_ID'], keep='first', inplace=True)
    logger.info('输出 KEGG_gene_def 文件，去重前的数量:{}，去重后的数量:{}'.format(kegg_clean_df.shape[0]-1, gene_def_df.shape[0]-1))
    
    # 提取 EC 编号
    def extract_ec_number(desc):
        try:
            ec_match = re.search(r'\[EC[^]]*\]', desc)
            if ec_match:
                return ec_match.group(0).strip('[]')
            return 'N/A'
        except (AttributeError, TypeError):
            return 'N/A'
            
    gene_def_df['EC_number'] = gene_def_df['Description EC_number'].apply(extract_ec_number)
    gene_def_df['KEGG_def'] = gene_def_df['Description EC_number'].str.split('[', expand=True)[0].str.strip()
    gene_def_df.fillna(value='N/A', inplace=True)
    
    # Gene_shortname 不能设置空为 NA，好像是某个软件识别 NA 会有问题（张老师说的）
    gene_def_df['Gene_shortname'] = gene_def_df['Gene_shortname'].str.split(',', expand=True)[0].fillna('')
    gene_def_df.drop(columns='Description EC_number', inplace=True)
    write_output_df(gene_def_df, output_prefix + 'KEGG_gene_def.txt', index=False)
    
    # 如果指定了 -i allgeneid 文件，则生成 all_gene_id + KEGG gene_short_name 新文件
    if all_id_file:
        all_gene_df = load_table(all_id_file, header=None, names=['GeneID'], dtype={"GeneID": str})
        all_gene_df = pd.merge(left=all_gene_df, right=gene_def_df[['GeneID', 'Gene_shortname']], on='GeneID', how='left')
        all_gene_df.fillna(value='', inplace=True)
        write_output_df(all_gene_df, output_prefix + 'shortname.txt', header=False, index=False)
    
    # 生成 tier2 文件
    tier2_name = output_prefix + 'KEGG_tier2.txt'
    tier2_df = kegg_clean_df.copy()
    tier2_df.drop(columns=['Pathway', 'Level1', 'KEGG_ID', 'Gene_shortname', 'Description EC_number'], inplace=True)
    tier2_df['Level2'] = tier2_df['Level2'].str.replace('\\','').str.replace(' / ', '_').str.replace('/', '_').str.replace(', ', '').str.replace(',', '_')
    tier2_df.drop_duplicates(keep='first', inplace=True)
    write_output_df(tier2_df, tier2_name, header=False, index=False)
    
    # tier3 文件，后来改名 KEGG.txt 了
    # geneid \t ko pathway
    # no header
    tier3_name = output_prefix + 'KEGG.txt'
    tier3_df = kegg_clean_df.copy()
    tier3_df.drop(columns=['Level2', 'Level1', 'KEGG_ID', 'Gene_shortname', 'Description EC_number'], inplace=True)
    tier3_df['Pathway'] = tier3_df['Pathway'].str.replace('\\','').str.replace(' / ', '_').str.replace('/', '_').str.replace(', ', '').str.replace(',', '_')
    tier3_df.drop_duplicates(keep='first', inplace=True)
    write_output_df(tier3_df, tier3_name, header=False, index=False)
    
    # ko03022_basal_transcription_factor.txt (子集) from tier3
    # 比对 ko03000_def.txt 加上定义
    ko03000_df, ko03022_df = ko03000(gene_def_df, tier3_df)
    write_output_df(ko03000_df, output_prefix + 'ko03000_transcription_factors.txt', index=False)
    write_output_df(ko03022_df, output_prefix + 'ko03022_basal_transcription_factor.txt', index=False)
    
    # ko03000 需要增加表达量，生成新文件 (2024_02_22)
    if fpkm and reads:
        ko03000_fpkm_reads_df = merge_fpkm_reads(fpkm, reads)
        ko03000_fpkm_reads_df = pd.merge(ko03000_fpkm_reads_df, ko03000_df, on='GeneID', how='inner')
        # ko03000_fpkm_reads_df = ko03000_fpkm_reads_df[ko03000_fpkm_reads_df['GeneID'].isin(ko03000_df['GeneID'])]
        write_output_df(ko03000_fpkm_reads_df, output_prefix + 'ko03000_expression_data_def.txt', index=False)


def kegg_levelb_count_barplot(keggclean_file, output_file):
    data_kegg = pd.read_csv(keggclean_file, sep="\t", header=None)

    # 统计 LevelB 和 LevelA 的组合频率
    count_kegg = data_kegg.iloc[:, 2:4].value_counts().reset_index()
    count_kegg = count_kegg[count_kegg.iloc[:, 2] > 0]

    # 清理 LevelA 和 LevelB 列
    count_kegg.iloc[:, 1] = count_kegg.iloc[:, 1].str.replace(r"A\d+:", "", regex=True)
    count_kegg.iloc[:, 0] = count_kegg.iloc[:, 0].str.replace(r"\d+:", "", regex=True)

    # 重命名列
    count_kegg.columns = ["LevelB", "LevelA", "Count"]
    count_kegg = count_kegg.sort_values(by=["LevelA", "LevelB"])

    # 计算 LevelB 的唯一值数量，用于动态设置图像高度
    n_levels = count_kegg["LevelB"].nunique()
    height = n_levels * 0.2  # 动态计算高度

    # 绘制图表
    plt.figure(figsize=(12, height))
    sns.barplot(data=count_kegg, y="LevelB", x="Count", hue="LevelA", dodge=False)
    plt.title("KEGG Count")
    plt.xlabel("Count")
    plt.ylabel("LevelB")
    plt.tight_layout()

    # 保存图表
    plt.savefig(output_file, dpi=300, bbox_inches="tight")


def parse_input():
    parser = argparse.ArgumentParser(description='kegg annotation，一次只能注释一个任务，除非添加其他邮箱')
    parser.add_argument('-f', '--fasta', help='输入 fasta file')
    parser.add_argument('-k', '--ko-file', dest='ko_file', type=str, required=True, help='指定 ko-file 文件')
    parser.add_argument('-l', '--org_lst', help='指定物种列表，以 ", " 分隔')
    parser.add_argument('-t', '--type', type=str, required=True, help='指定物种类型，植物=plant, 动物=animal')
    parser.add_argument('--allid', type=str, help='指定 all_id 文件，用来生成 shortname.txt')
    parser.add_argument('-o', '--output-prefix', dest='output_prefix', type=str, required=True, help='指定输出 prefix 例如 kegg/caomei')
    
    # 这两个参数用于生成 ko03000_expression_data.txt 文件
    parser.add_argument('--fpkm', type=str,
                        help='指定 fpkm 文件, 用于对 ko03000 添加表达量生成新文件，如果不指定，则不生成 ko03000_expression_data 文件')
    parser.add_argument('--reads', type=str,
                        help='指定 reads 文件, 用于对 ko03000 添加表达量生成新文件, 如果指定需要 fpkm 和 reads 都存在')
    
    account_parser = parser.add_argument_group("其他邮箱账号参数")
    account_parser.add_argument('-u', '--username', default='dagongrenyyds@163.com',
                        help='【默认就行】输入163邮箱用户名')
    account_parser.add_argument('-p', '--password', default='WJk6Bdg4y7tjz82W',  # IFBWCQOBTEOXODTK
                        help='【默认就行】输入163邮箱密码，密码非常规密码，查看获取方法 https://m.jqjq.net/jiqiao/108147.html')
    account_parser.add_argument('-m', '--mail_type', default='163', choices=['163', 'qq'],
                        help='【默认就行】输入邮箱类型 163 or qq，目前只支持 163 邮箱')
    
    args = parser.parse_args()
    
    
    # 检测 fpkm 和 reads 文件是否有效
    if args.fpkm or args.reads:
        if os.path.exists(args.fpkm) is False:
            raise Exception('fpkm 文件不存在')
        if os.path.exists(args.reads) is False:
            raise Exception('reads 文件不存在')
    else:
        args.fpkm, args.reads = False, False
    
    if args.output_prefix:
        if not args.output_prefix.endswith('_'):
            args.output_prefix = args.output_prefix + '_'
    else:
        args.output_prefix = ''
    
    return args


def main():
    args = parse_input()
    ko_file = args.ko_file
    # 直接运行注释
    if args.org_lst and args.fasta:
        kegg_anno(args.mail_type, args.username, args.password, args.fasta, args.org_lst, ko_file)

    # 解析 keg 文件
    parse_keg(ko_file, args.type, args.allid, args.fpkm, args.reads, args.output_prefix)
    
    # 画图
    kegg_levelb_count_barplot(f'{args.output_prefix}KEGG_clean.txt', f'{args.output_prefix}KEGG_levelB_count.jpeg')
    
    logger.success('Done!')


if __name__ == '__main__':
    main()