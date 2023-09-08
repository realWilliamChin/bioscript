#!/home/train/miniconda3/bin/python
# -*- coding: utf-8 -*-
# Created Time  : 2023/07/14 16:50
# Author        : William GoGo
import time
import re
import os, sys
import argparse
import requests
import poplib
from lxml import etree
from fake_useragent import UserAgent
from email import encoders
from email.header import Header, decode_header
from email.mime.text import MIMEText
from email.parser import Parser
from email.utils import parseaddr, formataddr


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
        try:
            server.user(user)
            server.pass_(pswd)
            resp, mails, octets = server.list()
            if len(mails) == 0:
                raise RuntimeError('No mail found from eamil address!')
        except Exception as e:
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
        'User-Agent': UserAgent().random,
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
    print(f'uploading {fasta_file_name}...')
    response = requests.post(upload_file_url, headers=headers, data=payload)
    if response.status_code == 200:
        print(f'uploaded {fasta_file_name}...')
    else:
        print(f'upload failed {fasta_file_name}!')
        print(response.text)
        exit(1)
    html = etree.HTML(response.text)
    
    # check upload status
    sub_element = html.xpath('//div[@id="main"]/p[1]/text()')[0]
    right_sub = 'Job Request'
    # An eamil has been sent to email_address for confirmation.
    # res_element = html.xpath('//p[@class="res"]/text()')[0]
    # right_res = f'An email has been sent to {eamil_address} for confirmation.'
    if sub_element == right_sub:
        print('upload success!')
    elif 'Sorry' in sub_element:
        print(f'upload failed!\nerror message: {sub_element}')
        sys.exit(1)
    else:
        print(f'upload failed!\n{response.text}')
        sys.exit(1)
    

def email_link_click(link):
    headers = {
        'Host': 'www.genome.jp',
        'User-Agent': UserAgent().random,
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
    

def parse_input():
    argparser = argparse.ArgumentParser(description='kegg annotation，一次只能注释一个任务，除非添加其他邮箱')
    argparser.add_argument('-f', '--fasta', required=True,
                           help='输入 fasta file')
    argparser.add_argument('-o', '--output', default='kegg_annotation.keg',
                           help='指定输出文件名')
    argparser.add_argument('-l', '--org_lst', required=True,
                           help='指定物种列表，以 ", " 分隔')
    argparser.add_argument('-u', '--username', default='dagongrenyyds@163.com',
                           help='【默认就行】输入163邮箱用户名')
    argparser.add_argument('-p', '--password', default='IFBWCQOBTEOXODTK',
                           help='【默认就行】输入163邮箱密码，密码非常规密码，查看获取方法 https://m.jqjq.net/jiqiao/108147.html')
    argparser.add_argument('-m', '--mail_type', default='163', choices=['163', 'qq'],
                           help='【默认就行】输入邮箱类型 163 or qq，目前只支持 163 邮箱')
    return argparser.parse_args()


def kegg_anno(mail_type, username, password, fasta_file, org_lst, output_filename):
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
                print(f'Job requesting ID:{job_ID} ...')
                mail_date = mail_content['date']
                break
            else:
                print(f'Job request Failed! not Submitted ID:{job_ID}!')
                sys.exit(1)
        elif mail_content['date'] == mail_date:
            print(f'receiving job request email {t+1} times ...')
            time.sleep(5)
            t += 1
        else:
            print(f'Job request Failed!\n{mail_content}')
            sys.exit(1)
            
    # step.2 获取任务提交状态
    status_link = None
    status_ID = None
    time.sleep(30)
    for t in range(0, 12):
        mail_content = mail.receive_mail(username, password)
        if 'KAAS - Accepted' == mail_content['subject'] and mail_content['from'] == '<kaas@genome.jp>':
            print(f'Job Accepted ID:{job_ID}!')
            status_link = re.findall(r"https://.*user.*", mail_content['text'])[0]
            status_ID = status_link.split('&')[-2].split('=')[-1]
            mail_date = mail_content['date']
            break
        elif mail_content['date'] == mail_date:
            print(f'get job status {t+1} times')
            time.sleep(30)
            t += 1
        else:
            print(f'Job Accept Failed!\n{mail_content}')
            sys.exit(1)
            
    # step.3 获取任务完成状态
    # 校验是否同一个任务
    if job_ID != status_ID:
        print('request_ID and status_ID are not the same!')
    print('Job is running wait for 5 minutes ...')
    time.sleep(5*60)
    
    for t in range(1, 360):
        mail_content = mail.receive_mail(username, password)
        
        if 'Annotation was completed' in mail_content['subject'] and 'kaas@genome.jp' in mail_content['from']:
            print('Job Complete!')
            
            # you should click this link then you can download the file
            print(f'downloading {output_filename} ...')
            visit_link = f'https://www.genome.jp/kaas-bin/kaas_main?mode=brite&id={job_ID}&key={job_key}'
            email_link_click(visit_link)
            time.sleep(5)
            
            # from download link get file
            download_link = f'https://www.genome.jp/kegg-bin/download_htext?htext=q00001.keg&format=htext&filedir=/tools/kaas/files/log/result/{job_ID}'
            resp = email_link_click(download_link)
            
            # check file if correct
            print('checking file ...')
            if resp.split('\n')[0] == '+D\tKO':
                print('checking file done ...')
            else:
                resp_file_header = resp.split('\n')[0:5]
                print(f'get file failed!, file is not correct ID{job_ID},key{job_key}!\n{resp_file_header}')
                sys.exit(1)
                
            # output file
            with open(output_filename, 'w') as kegg_annotation:
                kegg_annotation.write(resp)
            print(f'download success! {output_filename} is ready!')
            
            break
        elif mail_content['date'] == mail_date:
            if str(t / 5).split('.')[-1] == '0':
                print(f'Job is running, Have been waiting {t+5} minutes ...')
            elif t == 360:
                print(f'Job maybe error, Have been waiting 360 minutes ...')
                sys.exit(1)
            time.sleep(60)
            t += 1
        else:
            print(f'Job Failed!\n{mail_content}')
            sys.exit(1)


def main():
    args = parse_input()
    print(f"fasta file: {args.fasta}\norg_list：{args.org_lst}")
    kegg_anno(args.mail_type, args.username, args.password, args.fasta, args.org_lst, args.output)


if __name__ == '__main__':
    main()
    