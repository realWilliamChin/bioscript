import os
import pandas as pd

# 加载配置文件
configfile: 'config.yaml'

# 定义工作目录，运行目录
WORK_DIR = os.getcwd()
LOG_DIR = os.path.join(WORK_DIR, config['paths']['log_dir'])
DATA_DIR = os.path.join(WORK_DIR, config['paths']['input_data'])
TRINITY_DIR = os.path.join(WORK_DIR, config['paths']['trinity_output'])
REFERENCE_DIR = os.path.join(WORK_DIR, config['paths']['reference_dir'])
MAPPING_DIR = os.path.join(WORK_DIR, config['paths']['mapping_output'])
EXPRESSION_DIR = os.path.join(WORK_DIR, config['paths']['expression_output'])
SAMPLES_FILE = os.path.join(WORK_DIR, config['paths']['samples_file'])

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
    output:
        touch(f'{WORK_DIR}/logs/samples_check.done')
    shell:
        """
        python /home/colddata/qinqiang/script/CommonTools/check_SampDesAndCompInfo.py \
            -s {input.samples_file} \
            -d {input.data_dir}
        """