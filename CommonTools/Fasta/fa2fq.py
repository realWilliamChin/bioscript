import argparse
import random
from Bio import SeqIO

def fasta_to_fastq(input_fasta, output_fastq, quality_type='uniform', 
                   uniform_score=33, min_random_score=20, max_random_score=40, 
                   offset=33, id_suffix="/1"):
    """
    将FASTA文件转换为FASTQ文件
    
    参数:
    input_fasta: 输入FASTA文件路径
    output_fastq: 输出FASTQ文件路径
    quality_type: 质量分数类型 ('uniform'或'random')
    uniform_score: 统一质量分数值 (default=33)
    min_random_score: 随机质量分数最小值 (default=20)
    max_random_score: 随机质量分数最大值 (default=40)
    offset: FASTQ质量分数偏移量 (33=Phred+33, 64=Phred+64)
    id_suffix: 为每个read ID添加的后缀 (可选)
    """
    # 生成质量分数的函数
    def generate_quality_sequence(length, qual_type, uniform_score, min_score, max_score):
        if qual_type == 'uniform':
            # 创建统一质量分数的字符串
            return chr(uniform_score + offset) * length
        
        elif qual_type == 'random':
            # 创建随机质量分数 (默认20-40之间)
            return ''.join([chr(random.randint(min_score, max_score) + offset) 
                            for _ in range(length)])
        
        else:
            raise ValueError(f"不支持的质量分数类型: {qual_type}")
    
    # 记录转换统计
    converted_count = 0
    total_length = 0
    min_length = float('inf')
    max_length = 0
    
    # 创建输出文件
    with open(output_fastq, 'w') as output_handle:
        # 遍历输入FASTA文件
        for record in SeqIO.parse(input_fasta, "fasta"):
            seq = str(record.seq).upper()
            seq_len = len(seq)
            
            # 更新长度统计
            min_length = min(min_length, seq_len)
            max_length = max(max_length, seq_len)
            total_length += seq_len
            
            # 创建read ID (可选添加后缀)
            read_id = record.id
            if id_suffix and not read_id.endswith(id_suffix):
                read_id += id_suffix
            
            # 生成质量分数序列
            quality_string = generate_quality_sequence(
                seq_len, quality_type, uniform_score, min_random_score, max_random_score
            )
            
            # 写入FASTQ格式记录
            output_handle.write(f"@{read_id}\n")
            output_sequence = '\n'.join([seq[i:i+100] for i in range(0, seq_len, 100)])
            output_handle.write(f"{output_sequence}\n")
            output_handle.write(f"+\n")
            
            # 每行质量分数不超过100个字符
            for i in range(0, len(quality_string), 100):
                output_handle.write(f"{quality_string[i:i+100]}\n")
            
            converted_count += 1
    
    # 打印转换结果摘要
    avg_length = total_length / converted_count if converted_count > 0 else 0
    
    print(f"\n转换完成:")
    print(f"  FASTA文件: {input_fasta}")
    print(f"  FASTQ文件: {output_fastq}")
    print(f"  转换记录数: {converted_count}")
    print(f"  序列长度: 最小={min_length}bp, 最大={max_length}bp, 平均={avg_length:.1f}bp")
    print(f"  质量分数: 类型={quality_type}")
    
    if quality_type == 'uniform':
        print(f"  统一质量分数: {uniform_score} (Phred值), 字符={chr(uniform_score + offset)}")
    else:
        print(f"  随机质量范围: Phred {min_random_score}-{max_random_score}")
    
    print(f"  质量分数偏移: Phred+{offset}")
    if id_suffix:
        print(f"  Read ID后缀: '{id_suffix}'")
    
    print(f"\n提示: 该转换生成了模拟的质量分数, 并非原始测序质量!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='将FASTA文件转换为FASTQ格式')
    parser.add_argument('-i', '--input', required=True, help='输入FASTA文件路径')
    parser.add_argument('-o', '--output', required=True, help='输出FASTQ文件路径')
    parser.add_argument('-q', '--quality', choices=['uniform', 'random'], default='uniform', 
                       help='质量分数类型: uniform(默认)或random')
    parser.add_argument('--uniform_score', type=int, default=33, 
                       help='统一质量分数的数值(0-41之间, 默认=33)')
    parser.add_argument('--min_score', type=int, default=20, 
                       help='随机质量分数最小值(默认=20)')
    parser.add_argument('--max_score', type=int, default=40, 
                       help='随机质量分数最大值(默认=40)')
    parser.add_argument('--offset', type=int, choices=[33, 64], default=33, 
                       help='质量分数偏移: 33(Phred+33,默认)或64(Phred+64)')
    parser.add_argument('--suffix', type=str, default="/1", 
                       help='添加到read ID后的后缀(默认="/1")')
    
    args = parser.parse_args()
    
    # 输入验证
    if args.uniform_score < 0 or args.uniform_score > 41:
        parser.error("uniform_score 必须在0-41之间")
    if args.min_score < 0 or args.max_score > 41:
        parser.error("质量分数必须在0-41之间")
    if args.min_score > args.max_score:
        parser.error("min_score 不能大于 max_score")
    
    # 执行转换
    fasta_to_fastq(
        args.input, 
        args.output,
        quality_type=args.quality,
        uniform_score=args.uniform_score,
        min_random_score=args.min_score,
        max_random_score=args.max_score,
        offset=args.offset,
        id_suffix=args.suffix
    )