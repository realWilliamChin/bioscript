import os, sys
import requests
import pandas as pd
import time
from lxml import etree
import re


def get_tcmh_page(tcmh_number):

    url = f"https://www.bidd.group/TCMID/herb.php?herb=TCMH{tcmh_number}"

    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/122.0.0.0 Safari/537.36 Edg/122.0.0.0'
    }
    try:
        resp = requests.request("GET", url, headers=headers)
        if resp.status_code == 200:
            html = etree.HTML(resp.text)
        else:
            print(f"获取网页错误, code={resp.status_code}")
            time.sleep(10)
            get_tcmh_page(tcmh_number)
    except Exception as e:
        print(f"获取网页错误, 30 秒后重试{e}")
        time.sleep(30)
        get_tcmh_page(tcmh_number)
    
    return html


def get_all_tcmh_desc():
    open('tcmh.txt', 'w').write('Latin Name\tEnglish Name\t中文名\tTCMH_ID(Component ID)\n')
    tcmh = {"Latin Name": [], "English Name": [], "中文名": [], "TCMH_ID(Component ID)": []}
    for tcmh_number in range(1,2752):
        each_tcmh_info = []
        print(f"获取 {tcmh_number} 的详细信息中")
        tcmh_html_page = get_tcmh_page(tcmh_number)
        tcmh_info_section_dict = {
            1: "TCMH_ID(Component ID)",
            2: "Latin Name",
            3: "English Name",
            5: "中文名"
        }
        for each_section in tcmh_info_section_dict.keys():
            try:
                section_value = tcmh_html_page.xpath(f'//section[@class="section-padding"]/div//table/tr[{each_section}]/td[2]/text()')[0]
                each_tcmh_info.append(section_value)
                #tcmh[tcmh_info_section_dict[each_section]].append(section_value)
                #print(section_value)
            except Exception as e:
                print(f'获取 {tcmh_number} 的 {tcmh_info_section_dict[each_section]} 错误，{e}，设置为 NA')
        open('tcmh.txt', 'a').write('\t'.join(each_tcmh_info) + '\n')
        
    #tcmh_df = pd.DataFrame(tcmh)
    #tcmh_df.to_csv("tcmh.txt", sep='\t', index=False)

def download_img(src, output_name):
    url = f"https://www.bidd.group/TCMID/{src}"
    
    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/122.0.0.0 Safari/537.36 Edg/122.0.0.0'
    }
    try:
        resp = requests.request("GET", url, headers=headers)
        if resp.status_code == 200:
            open(output_name, 'wb').write(resp.content)
            return True
        elif resp.status_code == 404:
            print(f"没有这个图片, {output_name}，跳过")
            return False
        else:
            print(f"获取错误， code = {resp.status_code}")
            time.sleep(10)
            download_img(src, output_name)
    except Exception as e:
        print(f"获取网页错误, {e}")
        time.sleep(10)
        download_img(src, output_name)


def get_Molecular_Weight(src):
    print(f'正在获取 molecular_weight {src}')
    url = f"https://www.ebi.ac.uk/chembl/interface_api/es_proxy/es_data/get_es_document/chembl_molecule/{src.split('/')[-1]}"
    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/122.0.0.0 Safari/537.36 Edg/122.0.0.0'
    }
    try:
        resp = requests.request("GET", url, headers=headers)
        if resp.status_code == 200:
            data = resp.json()
            molecular_weight = data["_source"]['molecule_properties']['full_mwt']
            return molecular_weight
        elif resp.status_code == 404:
            print(f"The ID has not been included in the latest release.")
            return "INACTIVE"
        else:
            print(f"获取错误， code = {resp.status_code}")
            time.sleep(10)
            get_Molecular_Weight(src)
    except Exception as e:
        print(f"获取错误，{e}")
        time.sleep(30)
        get_Molecular_Weight(src)


def get_tcmh_chemical_ingredients(tcmh_number: int):
    # chemical ingredients 图片输出目录
    img_output_dir = f"TCMH{tcmh_number}"
    os.mkdir(img_output_dir)
    
    tcmh_html_page = get_tcmh_page(tcmh_number)
    
    # 卡片上的信息
    chemical_ingredients_ele_list = tcmh_html_page.xpath('//div[@id="structures"]/div[@class="container"]/div[@class="row"]/div[@class="col-sm-4 col-md-3 col-lg-2"]')
    # detail 里面的信息
    chemical_ingredients_info_ele_list = tcmh_html_page.xpath('//div[@id="structures"]/div[@class="container"]/div[@class="row"]/div[@class="modal"]')
    
    if len(chemical_ingredients_info_ele_list) != len(chemical_ingredients_ele_list):
        print(f'获取到的卡片和 detail 信息 list 不相等，请检查一下，tcmh_number = {tcmh_number}')
        return None
    
    ingredient_id_list = []
    formula_list = []
    chemical_structure_list = []
    common_name_list = []
    iupac_name_list = []
    molecular_weight_list = []
    
    # 循环卡片上的信息
    for each_ingredient in chemical_ingredients_ele_list:
        ingredient_ID_text_list = each_ingredient.xpath('.//div[@class="card-block"]/p/text()')
        ingredient_ID_text_list = [x.strip() for x in ingredient_ID_text_list if x.strip() != '']
        if 'Ingredient' in ingredient_ID_text_list[0] and 'Formula' in ingredient_ID_text_list[1]:
            ingredient_ID = ingredient_ID_text_list[0].split(':')[-1]
            formula_ID = ingredient_ID_text_list[1].split(':')[-1]
            img_src = each_ingredient.xpath('.//img[@class="card-img-top"]/@src')
            if len(img_src) != 1:
                print(f'当前 tcmh_number {tcmh_number} 的 ingredient ID {ingredient_ID} 获取图片 src 出现错误，不等于 1，跳过')
                print(img_src)
                continue
            else:
                img_src = img_src[0]
            output_img_name = os.path.join(img_output_dir, img_src.split('/')[-1])
            if download_img(img_src, output_img_name):
                chemical_structure_list.append("True")
            else:
                chemical_structure_list.append("False")
            ingredient_id_list.append(ingredient_ID)
            formula_list.append(formula_ID)
    
    # 循环 detail 里面的信息
    for each_ingredient in chemical_ingredients_info_ele_list:
        ingredient_common_name = each_ingredient.xpath('.//div[@class="modal-body"]/table[1]/tr[1]/td[2]//text()')
        ingredient_iupac_name = each_ingredient.xpath('.//div[@class="modal-body"]/table[1]/tr[2]/td[2]//text()')
        external_identifiers = each_ingredient.xpath('.//div[@class="modal-body"]/table[1]/tr[6]/td[2]/a[1]/@href')
        if len(ingredient_common_name) == 0:
            print(f'没有获取到 common_name')
            ingredient_common_name = "NA"
        elif len(ingredient_common_name) > 1:
            print(f'获取到的 common_name list 不是 1 {ingredient_common_name}')
            print(f'尝试选取第一个作为 common_name')
            ingredient_common_name = ingredient_common_name[0].strip()
        else:
            ingredient_common_name = ingredient_common_name[0].strip()
            
        if len(ingredient_iupac_name) == 0:
            print(f'没有获取到 iupac_name')
        elif len(ingredient_iupac_name) > 1:
            print(f'获取到的 iupac_name list 不是 1 {ingredient_iupac_name}')
            print(f'尝试选取第一个作为 iupac_name')
            ingredient_iupac_name = ingredient_iupac_name[0].strip()
        else:
            ingredient_iupac_name = ingredient_iupac_name[0].strip()
        
        if len(external_identifiers) == 0:
            print(f"没有获取到 Molecular Weight 的链接")
            external_identifiers = "NA"
        elif len(external_identifiers) > 1:
            print(f'获取到的 Molecular Weight 链接不是一个，{external_identifiers}')
            print(f'尝试选取第一个链接去获取')
            external_identifiers = external_identifiers[0].strip()
            external_identifiers = get_Molecular_Weight(external_identifiers)
        elif external_identifiers[0].strip().endswith("N/A"):
            external_identifiers = "NA"
        else:
            external_identifiers = external_identifiers[0].strip()
            external_identifiers = get_Molecular_Weight(external_identifiers)
        
        common_name_list.append(ingredient_common_name)
        iupac_name_list.append(ingredient_iupac_name)
        molecular_weight_list.append(external_identifiers)

    with open(f"TCMH{tcmh_number}_info.txt", 'w') as f:
        f.write(f"Common Name\tIUPAC Name\tChemical Structure\tIngredient ID\tFormula ID\tMolecular Weight\n")
        for i in range(0, len(ingredient_id_list)):
            f.write(f"{common_name_list[i]}\t{iupac_name_list[i]}\t{chemical_structure_list[i]}\t{ingredient_id_list[i]}\t{formula_list[i]}\t{molecular_weight_list[i]}\n")


def main():
    #get_all_tcmh_desc()
    tcmh_num_list = open('tcmh_number.txt', 'r').readlines()
    for tcmh_number in tcmh_num_list:
        tcmh_number = tcmh_number.strip()
        print(f'\n\n正在获取 {tcmh_number} 的详细信息')
        get_tcmh_chemical_ingredients(tcmh_number)


if __name__ == '__main__':
    main()