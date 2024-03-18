import requests
import pandas as pd
import time
from lxml import etree
import re


def get_all_module(ko_number):
    url = f'https://www.kegg.jp/entry/pathway+map{ko_number}'
    
    header = {
        "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/121.0.0.0 Safari/537.36 Edg/121.0.0.0"
    }
    
    resp = requests.get(url, headers=header)
    html = etree.HTML(resp.text)
    # print(resp.text)
    module_list_element = html.xpath('//td[@class="fr3 w1"]//tr[5]/td')[0]
    module_list = module_list_element.xpath('.//td[@class="vtop pd0"]//a/text()')
    print(module_list)
    
    return module_list


def get_module_info(module):
    url = f'https://www.kegg.jp/module/{module}'
    
    header = {
        "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/121.0.0.0 Safari/537.36 Edg/121.0.0.0"
    }
    
    resp = requests.get(url, headers=header)
    html = etree.HTML(resp.text)
    entry = html.xpath('//div[@class="definition"]//tr[1]/td[2]/text()')[0].strip()
    name = html.xpath('//div[@class="definition"]//tr[2]/td[2]/text()')[0].strip()
    definition_list = html.xpath('//div[@class="definition"]//tr[3]/td/a/text()')
    rc_number_list = html.xpath('//div[@class="definition"]//tr[6]/td[2]/a/text()')
    c_number_list = [x for x in rc_number_list if x.startswith('C')]
    r_number_list = [x for x in rc_number_list if x.startswith('R')]
    c_number_list = list(set(c_number_list))
    r_number_list = list(set(r_number_list))
    definition = ";".join(definition_list)
    
    return entry, name, definition, c_number_list, r_number_list
    

def main():
    ko_number = '01100'
    module_list = get_all_module(ko_number)
    
    df = pd.DataFrame(columns=['module', 'entry', 'name', 'definition', 'c_number', 'r_number'])
    for module in module_list:
        print(f'正在处理 {module}')
        entry, name, definition_list, c_number_list, r_number_list = get_module_info(module)
        name = name.replace(' ', '_')
        name = re.sub(r'[^\w\s]', '_', name)
        # 把 name 字符串所有的连续的下划线改为一个下划线
        while True:
            if '__' in name:
                name = name.replace('__', '_')
            else:
                break
        name = name.strip('_')
        c_number_list = ";".join(c_number_list)
        r_number_list = ";".join(r_number_list)
        print(f'entry: {entry}, name: {name}\ndefinition: {definition_list}, \nc_number: {c_number_list}\nr_number_list: {r_number_list}')
        df = df.append({
            'module': module,
            'entry': entry,
            'name': name,
            'definition': definition_list,
            'c_number': c_number_list,
            'r_number': r_number_list
            }, ignore_index=True)
        time.sleep(1)

    df.to_csv(f'ko01100.txt', sep="\t", index=False)


if __name__ == '__main__':
    main()


