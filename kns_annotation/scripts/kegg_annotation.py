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
from loguru import logger

sys.path.append(os.path.abspath('/home/colddata/qinqiang/script/'))
sys.path.append(os.path.abspath('/home/colddata/qinqiang/script/transcriptome/'))
from merge_fpkm_reads_matrix import merge_fpkm_reads


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
    if response.status_code == 200:
        logger.info(f'uploaded {fasta_file_name}...')
    else:
        logger.critical(f'upload failed {fasta_file_name}!')
        logger.critical(response.text)
        exit(1)
    html = etree.HTML(response.text)
    
    # check upload status
    sub_element = html.xpath('//div[@id="main"]/p[1]/text()')[0]
    right_sub = 'Job Request'
    # An eamil has been sent to email_address for confirmation.
    # res_element = html.xpath('//p[@class="res"]/text()')[0]
    # right_res = f'An email has been sent to {eamil_address} for confirmation.'
    if sub_element == right_sub:
        logger.success('upload success!')
    elif 'Sorry' in sub_element:
        logger.critical(f'upload failed!\nerror message: {sub_element}')
        sys.exit(1)
    else:
        logger.critical(f'upload failed!\n{response.text}')
        sys.exit(1)
    

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


def kegg_anno(mail_type, username, password, fasta_file, org_lst: str, output_filename):
    logger.info(f"fasta file: {fasta_file}\norg_list：{org_lst}")
    mail = EmailReceiver(f'imap.{mail_type}.com', 110)
    # step.1 submit fasta file to kegg and get request_ID and job status
    mail_content = mail.receive_mail(username, password)
    # 用来判断是否是同一个邮件
    mail_date = mail_content['date']

    upload_fasta_to_kegg(fasta_file, org_lst, username)
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
                sys.exit(1)
        elif mail_content['date'] == mail_date:
            logger.info(f'receiving job request email {t+1} times ...')
            time.sleep(5)
            t += 1
        else:
            logger.critical(f'Job request Failed!\n{mail_content}')
            sys.exit(1)
            
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
            sys.exit(1)
            
    # step.3 获取任务完成状态
    # 校验是否同一个任务
    if job_ID != status_ID:
        logger.error('request_ID and status_ID are not the same!')
    
    fasta_seq_num = len([record for record in SeqIO.parse(fasta_file, "fasta")])
    wait_times = int(fasta_seq_num/1000)*60
    logger.info(f'Job is running wait for {wait_times/60} minutes ...')
    time.sleep(wait_times)
    
    for t in range(1, 60):
        mail_content = mail.receive_mail(username, password)
        
        if 'Annotation was completed' in mail_content['subject'] and 'kaas@genome.jp' in mail_content['from']:
            logger.success('Job Complete!')
            
            # you should click this link then you can download the file
            logger.info(f'downloading {output_filename} ...')
            visit_link = f'https://www.genome.jp/kaas-bin/kaas_main?mode=brite&id={job_ID}&key={job_key}'
            email_link_click(visit_link)
            time.sleep(5)
            
            # from download link get file
            download_link = f'https://www.genome.jp/kegg-bin/download_htext?htext=q00001.keg&format=htext&filedir=/tools/kaas/files/log/result/{job_ID}'
            resp = email_link_click(download_link)
            
            # check file if correct
            logger.info('checking file ...')
            if resp.split('\n')[0] == '+D\tKO':
                logger.info('checking file done ...')
            else:
                resp_file_header = resp.split('\n')[0:5]
                logger.critical(f'get file failed!, file is not correct ID{job_ID},key{job_key}!\n{resp_file_header}')
                sys.exit(1)
                
            # output file
            with open(output_filename, 'w') as kegg_annotation:
                kegg_annotation.write(resp)
            logger.success(f'download success! {output_filename} is ready!')
            
            break
        elif mail_content['date'] == mail_date:
            logger.info(f'Job is running, checked email {t} times(5min) ...')
            if t == 60:
                logger.error(f'Job maybe error, Have been waiting 360 minutes ...')
                sys.exit(1)
            time.sleep(300)
            t += 1
        else:
            logger.critical(f'Job Failed!\n{mail_content}')
            sys.exit(1)


def ko03000(kegg_gene_df, kegg_tier3_df):
    ko03000_def_df_file = '/home/colddata/qinqiang/script/transcriptome/annotation/ko03000_def.txt'
    ko03000_def_df = pd.read_csv(ko03000_def_df_file, sep='\t')
    ko03000_geneid_list = kegg_tier3_df[kegg_tier3_df['Pathway'].str.contains("ko03000")]['GeneID'].tolist()
    ko03022_geneid_list = kegg_tier3_df[kegg_tier3_df['Pathway'].str.contains("ko03022")]['GeneID'].tolist()
    
    ko03000_df = kegg_gene_df[kegg_gene_df['GeneID'].isin(ko03000_geneid_list)][['GeneID', 'KEGG_ID']].copy()
    ko03022_df = kegg_gene_df[kegg_gene_df['GeneID'].isin(ko03022_geneid_list)][['GeneID', 'KEGG_ID']].copy()
    # ko03000_df = ko03000_df
    ko03000_df = pd.merge(left=ko03000_df, right=ko03000_def_df, on='KEGG_ID', how='left')
    
    return ko03000_df, ko03022_df


def parse_keg(keg_file, specie_type, all_id_file, fpkm, reads):
    key_name = keg_file.replace('.keg', '')
    kegg_file = keg_file.replace('.keg','_KEGG_original.txt')
    # 初始处理 keg 文件
    command = 'perl /home/colddata/chen/03_transcript/annotation/kegg/kaas_parse_keg_file.pl -i {} -o {}'.format(keg_file, kegg_file)
    os.system(command)
    
    # 使用 python 处理 keg 文件
    # process_keg(keg_file, kegg_file)
    
    with open(kegg_file,'r') as f1:
        kegg_clean_file = key_name + '_KEGG_clean.txt'
        f2 = open(kegg_clean_file, 'w')
        
        # 植物物种过滤动物的一些注释
        if specie_type == 'plant':
            for each_line in f1:
                if 'A09160' in each_line:
                    continue
                if 'A09190' in each_line:
                    continue
                if 'A09150' in each_line:
                    if '09158:Development and regeneration' not in each_line and '09159:Environmental adaptation' not in each_line:
                        continue
                f2.write(each_line)
                
        # 其他物种暂时不用过滤，直接写入
        else:
            for each_line in f1:
                f2.write(each_line)
        f2.close()
    
    # 生成 KEGG_gene_def 文件
    kegg_gene_def_titles = ['GeneID', 'Pathway', 'Level2', 'Level1', 'KEGG_ID', 'Gene_shortname', 'Description EC_number']
    kegg_clean_df = pd.read_csv(kegg_clean_file, sep='\t', header=None, names=kegg_gene_def_titles, dtype={"GeneID": str})
    gene_def_df = kegg_clean_df.drop(columns=['Pathway', 'Level2', 'Level1']).copy()
    gene_def_df.drop_duplicates(subset=['GeneID', 'KEGG_ID'], keep='first', inplace=True)
    logger.info('生成 KEGG_gene_def 文件，去重前的数量:{}，去重后的数量:{}'.format(kegg_clean_df.shape[0]-1, gene_def_df.shape[0]-1))
    gene_def_df['EC_number'] = gene_def_df['Description EC_number'].str.split('[', expand=True)[1].str.replace(']', '').str.split(' ', expand=True)[0]
    gene_def_df['KEGG_def'] = gene_def_df['Description EC_number'].str.split('[', expand=True)[0].str.strip()
    gene_def_df.fillna(value='NA', inplace=True)
    # Gene_shortname 不能设置空为 NA，设置为空（张老师说的）好像是某个软件识别 NA 会有问题
    gene_def_df['Gene_shortname'] = gene_def_df['Gene_shortname'].str.split(',', expand=True)[0].fillna(value='')
    gene_def_df.drop(columns='Description EC_number', inplace=True)
    gene_def_df.to_csv(key_name + '_KEGG_gene_def.txt', sep='\t', index=False)
    
    # 如果指定了 -i allgeneid 文件，则生成 all_gene_id + KEGG gene_short_name 新文件
    if all_id_file:
        all_gene_df = pd.read_csv(all_id_file, sep='\t', names=['GeneID'], dtype={"GeneID": str})
        all_gene_df = pd.merge(left=all_gene_df, right=gene_def_df[['GeneID', 'Gene_shortname']], on='GeneID', how='left')
        all_gene_df.fillna(value='', inplace=True)
        all_gene_df.to_csv(key_name + '_shortname.txt', sep='\t', index=False)
    
    # 生成 tier2 文件
    tier2_name = key_name + '_KEGG_tier2.txt'
    tier2_df = kegg_clean_df.copy()
    tier2_df.drop(columns=['Pathway', 'Level1', 'KEGG_ID', 'Gene_shortname', 'Description EC_number'], inplace=True)
    tier2_df['Level2'] = tier2_df['Level2'].str.replace('\\','').str.replace(' / ', '_').str.replace('/', '_').str.replace(', ', '').str.replace(',', '_')
    tier2_df.drop_duplicates(keep='first', inplace=True)
    tier2_df.to_csv(tier2_name, sep='\t', index=False, header=None)
    
    # tier3 文件，后来改名 KEGG.txt 了
    # geneid \t ko pathway
    # no header
    tier3_name = key_name + '_KEGG.txt'
    tier3_df = kegg_clean_df.copy()
    tier3_df.drop(columns=['Level2', 'Level1', 'KEGG_ID', 'Gene_shortname', 'Description EC_number'], inplace=True)
    tier3_df['Pathway'] = tier3_df['Pathway'].str.replace('\\','').str.replace(' / ', '_').str.replace('/', '_').str.replace(', ', '').str.replace(',', '_')
    tier3_df.drop_duplicates(keep='first', inplace=True)
    tier3_df.to_csv(tier3_name, sep='\t', index=False, header=None)
    
    # ko03022_basal_transcription_factor.txt (子集) from tier3
    # 比对 ko03000_def.txt 加上定义
    ko03000_df, ko03022_df = ko03000(gene_def_df, tier3_df)
    ko03000_df.to_csv(key_name + '_ko03000_transcription_factors.txt', sep='\t', index=False)
    ko03022_df.to_csv(key_name + '_ko03022_basal_transcription_factor.txt', sep='\t', index=False)
    
    # ko03000 需要增加表达量，生成新文件 (2024_02_22)
    if fpkm and reads:
        ko03000_fpkm_reads_df = merge_fpkm_reads(fpkm, reads)
        ko03000_fpkm_reads_df = pd.merge(ko03000_fpkm_reads_df, ko03000_df, on='GeneID', how='inner')
        # ko03000_fpkm_reads_df = ko03000_fpkm_reads_df[ko03000_fpkm_reads_df['GeneID'].isin(ko03000_df['GeneID'])]
        ko03000_fpkm_reads_df.to_csv(key_name + '_ko03000_expression_data_def.txt', sep='\t', index=False)


def split_fasta(fasta_file, split_parts: int):
    split_cmd = f'seqkit split -p {split_parts} {fasta_file}'
    logger.info(f'运行命令 {split_cmd}')
    os.system(split_cmd)


def parse_input():
    parser = argparse.ArgumentParser(description='kegg annotation，一次只能注释一个任务，除非添加其他邮箱')
    parser.add_argument('-f', '--fasta', help='输入 fasta file')
    parser.add_argument('-k', '--keg', type=str, help='指定 keg 文件')
    parser.add_argument('-l', '--org_lst', help='指定物种列表，以 ", " 分隔')
    parser.add_argument('-s', '--split', type=int, help='切割 fasta 文件，如果太大的话')
    parser.add_argument('-t', '--type', type=str, help='指定物种类型，植物=plant, 动物=animal')
    parser.add_argument('--allid', type=str, help='指定 all_id 文件，用来生成 shortname.txt')
    
    # 这两个参数用于生成 ko03000_expression_data.txt 文件
    parser.add_argument('--fpkm', type=str,
                        help='指定 fpkm 文件, 用于对 ko03000 添加表达量生成新文件，如果不指定，则不生成 ko03000_expression_data 文件')
    parser.add_argument('--reads', type=str,
                        help='指定 reads 文件, 用于对 ko03000 添加表达量生成新文件, 如果指定需要 fpkm 和 reads 都存在')
    
    account_parser = parser.add_argument_group("其他邮箱账号参数")
    account_parser.add_argument('-u', '--username', default='dagongrenyyds@163.com',
                        help='【默认就行】输入163邮箱用户名')
    account_parser.add_argument('-p', '--password', default='IFBWCQOBTEOXODTK',
                        help='【默认就行】输入163邮箱密码，密码非常规密码，查看获取方法 https://m.jqjq.net/jiqiao/108147.html')
    account_parser.add_argument('-m', '--mail_type', default='163', choices=['163', 'qq'],
                        help='【默认就行】输入邮箱类型 163 or qq，目前只支持 163 邮箱')
    
    args = parser.parse_args()
    
    if args.fasta and not args.keg:
        logger.critical('未输入 keg 文件名')
        sys.exit(1)
    
    if args.keg or args.fasta:
        pass
    else:
        try:
            args.keg = [x for x in os.listdir() if x.endswith('.keg')][0]
        except Exception:
            logger.critical('未输入 keg 文件名，尝试使用 keg 结尾文件未找到')
            sys.exit(1)
    
    if args.split == None:
        pass
    elif args.split == 1:
        delattr(args, "split")
    elif args.split >= 5:
        logger.warning(f'fasta 文件分割大于 5 次，将会延长注释时间')
    
    # 检测 fpkm 和 reads 文件是否有效
    if args.fpkm or args.reads:
        if os.path.exists(args.fpkm) is False:
            raise Exception('fpkm 文件不存在')
        if os.path.exists(args.reads) is False:
            raise Exception('reads 文件不存在')
    else:
        args.fpkm, args.reads = False, False
    
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_input()
    keg_file = args.keg
    if args.split:
        split_fasta(args.fasta, args.split)
        split_dir = f'{args.fasta}.split'
        fasta_file_list = [f'{os.path.join(split_dir, x)}' for x in os.listdir(split_dir)]
        for fasta_file in fasta_file_list:
            keg_sp_file = f'{fasta_file.split(".")[-2]}_{keg_file}'
            kegg_anno(args.mail_type, args.username, args.password, fasta_file, args.org_lst, keg_sp_file)
        os.system(f'cat *.keg > {keg_file}')
    elif args.org_lst and args.fasta:
        kegg_anno(args.mail_type, args.username, args.password, args.fasta, args.org_lst, keg_file)

    parse_keg(keg_file, args.type, args.allid, args.fpkm, args.reads)
    
    logger.success('Done!')



    