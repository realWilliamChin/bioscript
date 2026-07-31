#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created Time  : 2026/07/30 17:41
# Author        : William GoGo
import argparse
import csv
import sys
from pathlib import Path
from loguru import logger


def parse():
    p = argparse.ArgumentParser()
    p.add_argument('--in',  dest='infile',  required=True,
                   help='输入的原始 regulation 文件(两列:SYMBOL<TAB>value,无表头)')
    p.add_argument('--out', dest='outfile', required=True,
                   help='输出 regulation.txt(带表头:KEGG_ID<TAB>regulation,K 号)')
    p.add_argument('--kns', required=True,
                   help='<Species>_kns_gene_basicinfo_def.txt(需含 GeneSymbol、KEGG_ID 列)')
    return p.parse_args()


def main():
    args = parse()
    sym2k = {}
    sym2k_ci = {}
    with open(args.kns, encoding='utf-8') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            sym = (row.get('GeneSymbol') or '').strip()
            k   = (row.get('KEGG_ID')    or '').strip()
            if not sym or sym == 'N/A' or not k or k == 'N/A':
                continue
            for s in (x.strip() for x in sym.split(',') if x.strip()):
                sym2k.setdefault(s, k)
                sym2k_ci.setdefault(s.lower(), k)
    logger.info(f'加载 SYMBOL→K 映射:{len(sym2k)} 条 (from {args.kns})')

    hits, miss = [], []
    with open(args.infile, encoding='utf-8') as f:
        for line in f:
            line = line.rstrip('\n')
            if not line:
                continue
            parts = line.split('\t')
            sym = parts[0]
            val = parts[1] if len(parts) > 1 else ''
            k = sym2k.get(sym) or sym2k_ci.get(sym.lower())
            if k:
                hits.append((k, val))
            else:
                miss.append(sym)

    # K 号去重
    seen = set()
    uniq = []
    for k, v in hits:
        if k in seen:
            continue
        seen.add(k)
        uniq.append((k, v))

    Path(args.outfile).parent.mkdir(parents=True, exist_ok=True)
    with open(args.outfile, 'w', encoding='utf-8') as f:
        f.write('KEGG_ID\tregulation\n')
        for k, v in uniq:
            f.write(f'{k}\t{v}\n')

    logger.success(f'命中 {len(hits)} / 总 {len(hits) + len(miss)},去重后 {len(uniq)} 行 → {args.outfile}')
    if miss:
        logger.warning(f'未映射到 K 号的 SYMBOL({len(miss)}):{", ".join(miss)}')


if __name__ == '__main__':
    main()
