#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw(sum);
#example:perl summary.pl Phylum.txt Class.txt Order.txt  Family.txt Genus.txt  Species.txt
my $output_file = "column_sums.txt";
my @input_files = @ARGV or die "Usage: $0 <file1> [file2 ...]\n";




# 写入输出文件
open(my $out, '>', $output_file) or die "无法创建输出文件 $output_file: $!";
#close $out;
# 写入输出文件
#open(my $out, '>>', $output_file) or die "无法创建输出文件 $output_file: $!";
        
# 用于存储列总和的哈希
my %column_sums;
my @headers;
my $first_file = 1;
my $fileName="";

foreach my $file (@input_files) {
    open(my $fh, '<', $file) or die "无法打开文件 $file: $!";
    if ($file =~ /(.*?)\.txt/) {$fileName=$1;};
    
    my $line = <$fh>;  # 读取表头
    chomp $line;
    my @current_headers = split(/\t/, $line);
    
    if ($first_file) {
        @headers = @current_headers;
        $first_file = 0;
        
        print $out join("\t",@headers), "\n";
        #close $out;       
        
        
    } else {
        # 检查表头是否一致
        unless (@headers ~~ @current_headers) {
            die "错误: 文件 $file 的表头与其他文件不匹配\n";
        }
    }
    
    # 初始化列总和哈希
    unless (keys %column_sums) {
        foreach my $h (@headers) {
            $column_sums{$h} = 0;
        }
    }
    
    # 处理数据行
    while ($line = <$fh>) {
        chomp $line;
        my @fields = split(/\t/, $line);
        next unless @fields;  # 跳过空行
        
        # 第一列是样品名称，跳过不统计
        my $sample_name = shift @fields;
        
        # 累加各列值
        for my $i (0..$#fields) {
            my $header = $headers[$i+1];  # +1因为第一列是样品名
            my $value = $fields[$i] || 0;  # 空值视为0
            $column_sums{$header} += $value;
        }
    }
     print $out  "$fileName\t";
     #print $out  "ggg\t";
    close $fh;
   
    # 写入各列总和
    foreach my $header2 (@headers[1..$#headers]) {  # 跳过第一列(样品名)
     #print $out join("\t", $header, $column_sums{$header}), "\t";
     print $out "$column_sums{$header2}\t";
    }
    print $out "\n";
    foreach my $header3 (@headers[1..$#headers]) {  # 跳过第一列(样品名)
     #print $out join("\t", $header, $column_sums{$header}), "\t";
     $column_sums{$header3}=0;
    }
}



# 写入表头
#print $out join("\t", "Column", "Total_Sum"), "\n";

# 写入各列总和
foreach my $header (@headers[1..$#headers]) {  # 跳过第一列(样品名)
    #print $out join("\n", $column_sums{$header}), "\t";
}

close $out;

print "处理完成! 结果已保存到 $output_file\n";