import requests
import pandas as pd
import time


def crawler(product_number):
    search_url = f'http://www.coolaber.com/Json/GetSearchProduct.asp?rows=20&page=1&DisplyObj=field1%2Cfield8%2Cfield15%2Cprice_spec_count%2Cprice_spec&keyword={product_number}'
    
    header = {
        "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/121.0.0.0 Safari/537.36 Edg/121.0.0.0"
    }
    
    search_result = requests.get(search_url, headers=header)
    search_result_json = search_result.json()
    product_id = search_result_json['rows'][0]['P_ID']
    column_id = search_result_json['rows'][0]['P_CurrColumn']
    
    # product_info_url = f'http://www.coolaber.com/Column.asp?Model=Product_Detail&Column_ID={column_id}&P_ID={product_id}'
    product_info_url = f"http://www.coolaber.com/Json/GetSearchProduct.asp?action=GetShow&P_ID={product_id}&DisplyObj=field%2Cpicture%2Cvideo%2Cfile%2Cdescription%2Cprice_spec_count%2Cprice_spec"
    
    # 储存条件 OVF_Field8
    # 级别 OVF_Field4
    # 纯度 OVF_Field7
    
    product_info_resp = requests.get(product_info_url, headers=header)
    product_info_json = product_info_resp.json()
    chucun_tiaojian = product_info_json['OVF_Field8']
    jibie = product_info_json['OVF_Field4']
    chundu = product_info_json['OVF_Field7']
    
    return jibie, chundu, chucun_tiaojian
    

def main():
    product_list_file = 'product.txt'
    product_list = pd.read_csv(product_list_file, sep='\t', header=None, encoding='iso-8859-1')[0].values.tolist()
    
    with open('result.csv', 'w') as f:
        f.write('型号\t级别\t纯度\t储存条件\n')
        for product in product_list:
            try:
                jibie, chundu, chucun_tiaojian = crawler(product)
            except Exception:
                jibie, chundu, chucun_tiaojian = '未知', '未知', '未知'
            print(f"{product}\t{jibie}\t{chundu}\t{chucun_tiaojian}")
            f.write(f"{product}\t{jibie}\t{chundu}\t{chucun_tiaojian}\n")
            time.sleep(1)
            



if __name__ == '__main__':
    main()