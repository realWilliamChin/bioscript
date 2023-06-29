#!/home/train/miniconda3/bin/python
# -*- coding: utf-8 -*-
# Created Time  : 2023/06/27 11:00
# Author        : William GoGo
import os
import pandas as pd
import argparse
import requests
import json
from urllib.parse import quote
from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.support.ui import Select
from selenium.webdriver.common.by import By
from fake_useragent import UserAgent
import time


def parse_input():
    argparser = argparse.ArgumentParser(description='')
    argparser.add_argument('-a', '--all', default='all_geneid_Go.txt',
                           help='指定所有基因的 GO 文件, 默认当前文件夹下的 _all_geneid_Go 文件')
    args = argparser.parse_args()
    return args


def geneid_goid(all_geneid_goid_file, updownid_file):
    """
    生成每个样本的上下调基因的 GO 文件
    """
    updown_go_filename = updownid_file.replace('_ID.txt', '_Go.txt')
    # 读取 updownid_file 为 list
    with open(updownid_file, 'r') as file:
        id_lst = list(file.read().splitlines())
    # 循环读取 all_geneid_goid_file，如果 GeneID 在 id_lst 中，就写入 updown_go_filename
    with open(all_geneid_goid_file, 'r') as file:
        for line in file:
            geneid = line.split('\t')[0]
            if geneid in id_lst:
                with open(updown_go_filename, 'a') as updown_go_file:
                    updown_go_file.write(line)


def get_wego_pic(geneid_goid_file):
    # step.1 上传文件到 wego
    with open(geneid_goid_file, 'r') as file:
        lines = file.read().splitlines()
        geneid_goid_content = '\r\n'.join(lines)
    
    upload_file_url = 'https://wego.genomics.cn/file/upload'
    headers = {
        'Host': 'wego.genomics.cn',
        'Connection': 'keep-alive',
        'sec-ch-ua': '"Not.A/Brand";v="8", "Chromium";v="114", "Microsoft Edge";v="114"',
        'sec-ch-ua-mobile': '?0',
        'User-Agent': UserAgent().random,
        'Content-Type': 'multipart/form-data; boundary=----WebKitFormBoundarynPWPEF3MomzY7ppI',
        'Accept': 'application/json',
        'Cache-Control': 'no-cache',
        'X-Requested-With': 'XMLHttpRequest',
        'sec-ch-ua-platform': '"Windows"',
        'Origin': 'https://wego.genomics.cn',
        'Sec-Fetch-Site': 'same-origin',
        'Sec-Fetch-Mode': 'cors',
        'Sec-Fetch-Dest': 'empty',
        'Referer': 'https://wego.genomics.cn/',
        'Accept-Encoding': 'gzip, deflate, br',
        'Accept-Language': 'en-US,en;q=0.9'
    }
    payload = f'------WebKitFormBoundarynPWPEF3MomzY7ppI\r\nContent-Disposition: form-data; name=\"file[0]\"; filename=\"{geneid_goid_file.replace("_Go.txt", "")}\"\r\nContent-Type: text/plain\r\n\r\n{geneid_goid_content}\r\n\r\n------WebKitFormBoundarynPWPEF3MomzY7ppI--\r\n'
    
    response = requests.post(upload_file_url, headers=headers, data=payload)
    upload_code = response.json()['code']
    
    # step.2 获取 wego 网页
    get_view_url = 'https://wego.genomics.cn/savejob'
    payload = quote(f'bundle[dir]={dir}&bundle[archive]=2018-11-01&bundle[format]=native', safe='=&')
    headers['Content-Type'] = 'application/x-www-form-urlencoded; charset=UTF-8'
    headers['Accept'] = '*/*'
    response = requests.post(get_view_url, headers=headers, data=payload)
    getpic_code = response.json()['code']
    print(f'{geneid_goid_file}：\t\tupload:{upload_code}\tgetview:{getpic_code}')
    view_ID = response.json()['jid']
    
    # step.3 加载网页
    chrome_options = Options()
    chrome_options.add_argument('--headless')
    chrome_options.add_argument('--disable-gpu')
    chrome_options.add_argument("--no-sandbox")
    chrome_options.add_argument("--disable-dev-shm-usage")
    driver = webdriver.Chrome(service=Service('./lib/chromedriver'), options=chrome_options)
    # driver.implicitly_wait(10)
    
    driver.set_window_size(1920, 1440)
    wego_view_link = f'https://wego.genomics.cn/view/{view_ID}'
    driver.get(wego_view_link)
    
    # step.4 get job status
    headers['Referer'] = wego_view_link
    get_status_url = 'https://wego.genomics.cn/getJobDetail'
    max_retry = 0
    while max_retry < 500:
        response = requests.post(get_status_url, headers=headers, data=f'jobname={view_ID}&level=3')
        if response.json()['code'] == 'success':
            break
        time.sleep(1)
        max_retry += 1
    driver.find_element(By.ID, "graph_tab").click()
    # 设置大小
    picsize_element = driver.find_element(By.ID, "adj_size_1")
    picsize_element.clear()
    picsize_element.send_keys("1500*1000")
    # 设置标题（改上面文件名就好了，不用改这个了
    # driver.find_element(By.CLASS_NAME, 'legend_name'):
    # pictitle_element.clear()
    # pictitle_element.send_keys(geneid_goid_file.replace('_Go.txt', ''))
    # 设置图片类型
    pictype_element = driver.find_element(By.ID, "g1_type_sel")
    Select(pictype_element).select_by_visible_text("image/jpeg")
    time.sleep(1)
    
    download_max_retry = 0
    driver.save_screenshot(f'{geneid_goid_file}.png')
    while download_max_retry < 5:
        driver.find_element(By.ID, "download_btn1").click()
        time.sleep(3)
        for jepg_file in os.listdir():
            if jepg_file == 'graph1.jpeg':
                os.rename(jepg_file, geneid_goid_file.replace('.txt', '.jpeg'))
                driver.quit()
                return 0
        download_max_retry += 1
    driver.quit()
    print(f'download wego picture failed {geneid_goid_file}')


def main():
    args = parse_input()
    for each_file in os.listdir():
        if 'ID.txt' not in each_file:
            continue
        if 'Up' or 'Down' in each_file:
            geneid_goid(args.all, each_file)
    for each_file in os.listdir():
        if each_file.endswith('Go.txt'):
            get_wego_pic(each_file)
    

if __name__ == '__main__':
    main()
    