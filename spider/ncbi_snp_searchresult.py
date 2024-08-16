#!/home/train/miniconda3/bin/python
# -*- coding: utf-8 -*-
# Created Time  : 2024/08/12 11:26
# Author        : William GoGo
import os
import argparse
import requests
import json
from urllib.parse import quote
from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.support.ui import Select
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from fake_useragent import UserAgent
import time
from loguru import logger


def ncbi_snp_search(search_item, output_dir='./'):
    
    chrome_options = Options()
    chrome_options.add_argument('--headless')
    chrome_options.add_argument('--disable-gpu')
    chrome_options.add_argument("--no-sandbox")
    chrome_options.add_argument("--disable-dev-shm-usage")
    chrome_options.add_argument("--ignore-certificate-errors")
    chrome_options.add_experimental_option("prefs", {"download.default_directory": output_dir})
    # driver = webdriver.Chrome(service=Service('C:\\Users\\PC\Downloads\\chromedriver-win64\\chromedriver-win64\\chromedriver.exe'), options=chrome_options)
    driver = webdriver.Chrome(service=Service('/home/colddata/qinqiang/script/lib/chromedriver'), options=chrome_options)
    # driver.implicitly_wait(10)
    driver.set_window_size(1920, 1440)
    ncbi_link = f'https://www.ncbi.nlm.nih.gov/'
    driver.get(ncbi_link)

    select_element = driver.find_element(By.ID, "database")
    sellect_object = Select(select_element)
    sellect_object.select_by_value("snp")
    
    input_element = driver.find_element(By.ID, "term")
    input_element.clear()
    input_element.send_keys(search_item)
    input_element.send_keys(Keys.RETURN)

    chromosome_element = driver.find_element(By.XPATH, "//dt[text()='Chromosome: ']/following-sibling::dd")
    chromosome_text = chromosome_element.text
    chromosome_1_values = chromosome_text.split('\n')[0].split("(")[0].split(':')[1]
    chromosome_2_values = chromosome_text.split('\n')[1].split("(")[0].split(':')[1]
    print(chromosome_1_values, chromosome_2_values)

    varview_element = driver.find_element(By.XPATH, "//span[@class='snpsum_hgvs']/a")
    varview_element.click()
    
    def find_and_click(driver, element_name, xpath):
        try:
            logger.info(f'定位 {element_name}')
            element = WebDriverWait(driver, 300).until(EC.presence_of_element_located((By.XPATH, xpath)))
            element.click()
        except TimeoutError:
            print(f"{element_name} 元素未在规定时间内出现，程序退出")
        return element

    for chromosome in [chromosome_1_values, chromosome_2_values]:
        chromosome_values = chromosome
        
        find_and_click(driver, "tools 按钮", '//span[text()="Tools"]')
        find_and_click(driver, 'search 按钮', '//span[text()="Search"]')
        search_box = find_and_click(driver, 'search 输入框', '//input[@role="combobox"]')
        search_box.send_keys(f'{int(chromosome_values)-500}-{int(chromosome_values)+500}')
        find_and_click(driver, 'search ok 按钮', '//span[text()="OK"]')
        time.sleep(15)
        find_and_click(driver, 'download 按钮', '//span[text()="Download"]')
        find_and_click(driver, 'download fasta 按钮', '//span[text()="Download FASTA"]')
        time.sleep(3)
        find_and_click(driver, 'download visible range 按钮', '//span[text()="FASTA (Visible Range)"]')
        time.sleep(5)
    
    all_download_files = [x for x in os.listdir(output_dir) if x.endswith('fa')]
    for f in all_download_files:
        data = open(f, 'r').read()
        open(f, 'w').write('>' + search_item + '_' + data.replace('>', ''))


if __name__ == '__main__':
    ncbi_snp_search("rs1131454")