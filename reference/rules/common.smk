import os
import pandas as pd


# 加载配置文件
configfile: 'config.yaml'

# 定义工作目录，运行目录
WORK_DIR = os.getcwd()
LOG_DIR = config['paths']['log_dir'].format(work_dir=WORK_DIR)
DATA_DIR = config['paths']['input_data'].format(work_dir=WORK_DIR)
MAPPING_DIR = config['paths']['mapping_output'].format(work_dir=WORK_DIR)
EXPRESSION_DIR = config['paths']['expression_output'].format(work_dir=WORK_DIR)
SAMPLES_FILE = config['paths']['samples_file'].format(work_dir=WORK_DIR)

SPECIE_NAME = config['specie']

samples_df = pd.read_csv(f'{SAMPLES_FILE}', sep='\t')
samples_list = samples_df['sample'].tolist()


# --------------------------
# run prep
# --------------------------
rule check_samples_data:
    input:
        samples_file = f'{SAMPLES_FILE}',
        data_dir = f'{DATA_DIR}'
    priority: 100
    output:
        touch(f'{WORK_DIR}/logs/samples_check.done')
    shell:
        """
        python /home/colddata/qinqiang/script/CommonTools/check_cleandata_samples.py \
            -s {input.samples_file} \
            --cd {input.data_dir}
        """