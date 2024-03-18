#!/home/train/miniconda3/bin/python
# -*- coding: utf-8 -*-
# Created Time  : 2023/09/11 10:07
# Author        : William GoGo
import requests
import math
import time
import pandas as pd
from lxml import etree
from fake_useragent import UserAgent


def search_pathway_object(search_list, org_type):
    
    search_list = '\r\n'.join(search_list)
    # POST
    boundary = '----WebKitFormBoundarygRLF3yT3KUXk1BB5'
    url = 'https://www.kegg.jp/kegg-bin/search_pathway_object'
    header = {
        'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3;q=0.7',
        'Accept-Encoding': 'gzip, deflate, br',
        'Accept-Language': 'en-US,en;q=0.9',
        'Cache-Control': 'max-age=0',
        'Content-Type': f'multipart/form-data; boundary={boundary}',
        'Content-Length': '49276',
        'Host': 'www.kegg.jp',
        'Origin': 'https://www.kegg.jp',
        'Referer': 'https://www.kegg.jp/kegg/mapper/search.html',
        'Upgrade-Insecure-Requests': '1',
        'User-Agent': UserAgent().random,
        'sec-ch-ua': '"Chromium";v="116", "Not)A;Brand";v="24", "Microsoft Edge";v="116"',
    }
    boundary = '--' + boundary
    payload = '\r\n'.join([
        boundary,
        'Content-Disposition: form-data; name="org_type"',
        '',
        f'{org_type}',
        boundary,
        'Content-Disposition: form-data; name="org"',
        '',
        '',
        boundary,
        'Content-Disposition: form-data; name="unclassified"',
        '',
        f'{search_list}',
        '',
        boundary,
        'Content-Disposition: form-data; name="color_list"; filename=""',
        'Content-Type: application/octet-stream',
        '',
        '',
        boundary,
        'Content-Disposition: form-data; name="default"',
        '',
        'pink',
        boundary,
        'Content-Disposition: form-data; name="org_name"',
        '',
        f'{org_type}',
        boundary,
        'Content-Disposition: form-data; name="s_sample"',
        '',
        '',
        boundary
    ]).strip("'")
    
    max_retry = 0
    while max_retry <= 5:
        resp = requests.post(url, headers=header, data=payload)
        if resp.status_code != 200:
            max_retry += 1
            time.sleep(5)
            continue
        else:
            break
    # open('save.txt', 'w').write(resp.text)
    html = etree.HTML(resp.text)
    pathway_lst = html.xpath('//div[@class="box2"]//li/a[1]/text()')
    print(f'搜索到 pathway list 有 {len(pathway_lst)} 个')
    
    return pathway_lst,


def search_compound_path(kegg_pathway, page_number=1):
    
    url = f'https://www.genome.jp/dbget-bin/get_linkdb?-t+compound+-p+{page_number}+path:map{kegg_pathway}'
    header = {
        "Accept": 'text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3;q=0.7',
        "Accept-Encoding": 'gzip, deflate, br',
        "Host": "www.genome.jp",
        "User-Agent": UserAgent().random,
    }
    
    max_retry = 0
    while max_retry <= 5:
        resp = requests.get(url, headers=header)
        if resp.status_code != 200:
            max_retry += 1
            time.sleep(5)
            continue
        else:
            break
    
    html = etree.HTML(resp.text)
    
    return html


def process_search_compound_path(kegg_pathway, org_type):
    if kegg_pathway.startswith(org_type):
        kegg_pathway = kegg_pathway.replace(org_type, '')
    html = search_compound_path(kegg_pathway)
   
    # 检查总数，超过 1000 会分页 
    all_compound_displayed = int([x for x in html.xpath('/html/body/text()') if 'Hits:' in x][0].split(':')[-1].split('from')[0].strip())
    # 总数为 0，返回 None，跳过此页
    if all_compound_displayed == 0:
        print("检查页面总数为 0，跳过")
        return None
    page_count = math.ceil(all_compound_displayed / 1000)

    all_kegg_compound_definition_lst = []
    all_kegg_compound_id_lst = []
    if page_count > 1:
        for i in range(1, page_count+1):
            html = search_compound_path(kegg_pathway, i) 

            kegg_compound_id_lst = html.xpath('//pre/a/text()')
            kegg_compound_definition_lst = html.xpath('//pre/text()[position() >= 2]')

            all_kegg_compound_id_lst = all_kegg_compound_id_lst + kegg_compound_id_lst
            all_kegg_compound_definition_lst = all_kegg_compound_definition_lst + kegg_compound_definition_lst
    else:
        html = search_compound_path(kegg_pathway) 

        all_kegg_compound_id_lst = html.xpath('//pre/a/text()')
        all_kegg_compound_definition_lst = html.xpath('//pre/text()[position() >= 2]')

    all_kegg_compound_definition_lst = [x.strip() for x in all_kegg_compound_definition_lst]

    # 判断
    if all_compound_displayed == len(all_kegg_compound_definition_lst) == len(all_kegg_compound_id_lst):
        print(f'检查页面显示总数:{all_compound_displayed}, 写入总数:{len(all_kegg_compound_id_lst)}')
    else:
        print(all_compound_displayed, len(all_kegg_compound_id_lst), len(all_kegg_compound_definition_lst))
        print("数量不匹配，程序退出")
        exit(1)
    
    df = pd.DataFrame({
        "kegg_compound_id": all_kegg_compound_id_lst,
        f"{org_type}": f"{org_type}{kegg_pathway}",
        "kegg_compound_definition": all_kegg_compound_definition_lst
    })

    return df


def main():
    # 2023_09_11
    # org_type = 'hsa'
    # file = '/home/colddata/qinqiang/work/2023_09_11_crawl_kegg_pathway/2023_03_06_KEGG_ID_class_def.txt'
    # search_lst = list(pd.read_csv(file, sep='\t', usecols=[0], header=None)[0])
    # print(f'input count {len(search_lst)}')
    # pathway_list = search_pathway_object(search_lst, org_type)[0]
    # df_lst = []

    # 2023_09_13
    org_type = 'hsa'
    filee = '/home/colddata/qinqiang/work/2023_09_11_crawl_kegg_pathway/Hsa_KEGG_Pathway_ID.txt'
    dff = pd.read_csv(filee, sep='\t', header=None, names=[org_type, org_type + '_def'])
    pathway_list = list(dff[org_type].values)
    df_lst = [] 

    for i, each_pathway in enumerate(pathway_list):
        print(f'processing {i}-{each_pathway} ...')
        df = process_search_compound_path(each_pathway, org_type)
        # 如果没有会返回 None
        # if df != None:
        df_lst.append(df)
        
    
    result_df = pd.concat(df_lst, axis=0)

    # 2023_09_13 替换列
    result_df = pd.merge(left=result_df, right=dff, left_on=org_type, right_on=org_type, how='left')
    result_df.drop(columns=[org_type], inplace=True)
    result_df = result_df.reindex(columns=['kegg_compound_id', org_type+'_def', 'kegg_compound_definition'])


    result_df.to_csv(f'/home/colddata/qinqiang/work/2023_09_11_crawl_kegg_pathway/Hsa_KEGG_Pathway_ID_result.txt', sep='\t', index=False)
    print(result_df)


if __name__ == '__main__':
    main()