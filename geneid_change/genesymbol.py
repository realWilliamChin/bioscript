#!/usr/bin/env python
import pandas as pd


gene_info = 'geneid_def/gene_info'
gene2ensembl = 'geneid_def/gene2ensembl'
species = 'species.txt'

species_df = pd.read_csv(species, sep='\t', header=None, names=['Specie', '#tax_id'])
gene_info_df = pd.read_csv(gene_info, sep='\t', usecols=[0, 1, 2, 3, 8, 9], low_memory=False)
gene2ensembl_df = pd.read_csv(gene2ensembl, sep='\t', usecols=[0, 1, 2])


for each_specie in species_df.itertuples():
    print(f'processing{each_specie[1]}')
    each_specie_df = gene_info_df[gene_info_df['#tax_id'] == each_specie[2]].copy()
    each_specie_df.drop(columns=['#tax_id'], inplace=True)
    each_specie_df.to_csv('02_NCBI2GeneSymbol/' + each_specie[1] + '_Ncbi2genesymbol.txt', sep='\t', index=False)
    
    each_specie_df2 = gene2ensembl_df[gene2ensembl_df['#tax_id'] == each_specie[2]].copy()
    each_specie_df2.drop(columns=['#tax_id'], inplace=True)
    each_specie_df2.to_csv('01_NCBI2Embl/' + each_specie[1] + '_Ncbi2embl.txt', sep='\t', index=False)



for each_specie in species_df.itertuples():
#     NCBI2EMbl_filelist = [x for x in os.listdir('01_NCBI2Embl')]
#     NCBI2GeneSymbol_filelist = [x for x in os.listdir('02_NCBI2GeneSymbol')]
#     print(f'{each_specie[1]}_Ncbi2embl.txt', f'{each_specie[1]}_Ncbi2genesymbo.txt')
#     if f'{each_specie[1]}_Ncbi2embl.txt' in NCBI2EMbl_filelist and f'{each_specie[1]}_Ncbi2genesymbo.txt' in NCBI2GeneSymbol_filelist:
#         print(f'processing {each_specie}')
    embl_df = pd.read_csv(f'01_NCBI2Embl/{each_specie[1]}_Ncbi2embl.txt', sep='\t', dtype={"GeneID": str}, names=['GeneID', 'Embl_ID'], skiprows=1)
    symbol_df = pd.read_csv(f'02_NCBI2GeneSymbol/{each_specie[1]}_Ncbi2genesymbol.txt', usecols=[0, 1], names=['GeneID', 'GeneSymbol'], sep='\t', dtype={"GeneID": str})

    result_df = pd.merge(left=embl_df, right=symbol_df, how='inner', on='GeneID')
    result_df.drop(columns=['GeneID'], inplace=True)
    result_df.drop_duplicates(subset=['Embl_ID', 'GeneSymbol'], inplace=True)
    result_df.rename(columns={'Embl_ID': 'GeneID'}, inplace=True)
    result_df.to_csv(f'03_Result/{each_specie[1]}_Embl2GeneSymbol.txt', sep='\t', index=False)
    
    
# for i in os.listdir('03_Result'):
#     df = pd.read_csv('03_Result/' + i, sep='\t', names=['GeneID', 'GeneSymbol'], skiprows=0)
#     df.to_csv('03_Result/' + i, sep='\t', index=False)