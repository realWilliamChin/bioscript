20260608
#服务器：
cd /home/colddata/qinqiang/Project/Test/singleCellPipeline_test/pipelineTest1/

conda activate test_scanpy_python311
python

import os
import scanpy as sc
import anndata as ad
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import shutil
import pickle
import scanpy.external as sce
from itertools import combinations
from adjustText import adjust_text  # 关键：防基因名重叠
import sc_utils as scu
from harmonypy import run_harmony

import celltypist
from celltypist import models

from scipy.stats import kruskal, mannwhitneyu
import scikit_posthocs as sp
from statsmodels.stats.multitest import multipletests
from scipy import stats

rawData_dir='/NAS/CurrentData/2026_data/2026_05_08_单细胞数据分析_狼疮性肾炎_杭州/Data5/GSE285773_RAW'

sample2mtxDir_dict={e.split('_')[0]:os.path.join(rawData_dir, e.split('_')[0]) for e in os.listdir(rawData_dir)}

targetSampleIDs_lst=['GSM8708836', 'GSM8708847', 'GSM8708839', 'GSM8708850', 'GSM8708840', 'GSM8708843', 'GSM8708841', 'GSM8708837', 'GSM8708838', 'GSM8708842', 'GSM8708828', 'GSM8708830', 'GSM8708833']

#如每个样本已有下面三个文件（barcodes.tsv.gz、features.tsv.gz和matrix.mtx.gz）下面生成sample2mtxFiles_dict的代码不用运行：
sample2mtxFiles_dict={}
for e in os.listdir(rawData_dir):
	sampleID=e.split('_')[0]
	mtxFile=os.path.join(os.getcwd(), os.path.join(rawData_dir, e))
	if sampleID in sample2mtxFiles_dict: 
		sample2mtxFiles_dict[sampleID].append(mtxFile)
	else:
		sample2mtxFiles_dict[sampleID]=[mtxFile]
		
sample2adata_dict={}
for k, v in sample2mtxDir_dict.items():
	if k in targetSampleIDs_lst:
		#adata_doubletsRmv(os.path.join(rawData_dir, k)) 
		#for e in sample2mtxFiles_dict[k]: #如每个样本已有下面三个文件（barcodes.tsv.gz、features.tsv.gz和matrix.mtx.gz）此行及下面2行不用运行
		#	shutil.copy(e, os.path.join(rawData_dir, k, e.split('_')[-1]))
		sample2adata_dict[k]=sc.read_10x_mtx(os.path.join(rawData_dir, k))

#读取样品相关信息：sampleID、groupID和dataID。
samplesInfos_df=pd.read_excel('/NAS/CurrentData/2026_data/2026_05_08_单细胞数据分析_狼疮性肾炎_杭州/Data5/Data5_record.xlsx', sheet_name='分析数据整理', header=1)
samplesMeta_df = samplesInfos_df[['Sample_ID', 'Group_ID', 'Data_ID']]
samplesMeta_df=samplesMeta_df.dropna()

samplesMeta_df['dataID']=samplesMeta_df['Data_ID']
samplesMeta_df['Group_ID']=samplesMeta_df['Group_ID'].apply(lambda x:x.replace(' ', ''))

samplesMeta_df.set_index('Data_ID', inplace=True)
sample_meta={}
for k, v in sample2adata_dict.items():
	sample_meta[k]=samplesMeta_df.loc[k].to_dict()
	sample_meta[k].update({'path':v})

adatas = []
for sample_id, info in sample_meta.items():
    # 1. 读取数据
    adata = sample2adata_dict[sample_id]
    
    # 2. 强制修复细胞ID重复（核心步骤！）
    # 生成绝对唯一的细胞ID：sample_id + 序号（避免原始ID格式问题导致重复）
    new_obs_names = [f"{sample_id}_cell_{i}" for i in range(adata.n_obs)]
    
    # 3. 修复基因名重复（统一基因名格式，避免大小写/版本号导致的重复）
    adata.var_names = adata.var_names.str.upper()  # 全部转为大写
    adata.var_names = adata.var_names.str.replace(r'\.\d+$', '', regex=True)  # 去除基因名后的版本号（如CD3D.1 → CD3D）
    # 去除重复基因（保留第一个）
    adata = adata[:, ~adata.var_names.duplicated()]
    
    # 4. 添加元信息
    adata.obs['sample'] = sample_id
    adata.obs['sampleID'] = info['Sample_ID']
    adata.obs['group'] = info['Group_ID']
    adata.obs['cellID'] = new_obs_names  # 覆盖原始细胞ID
    
    # 5. 检查当前样品的标签唯一性（调试用，可删除）
    assert len(adata.obs_names.unique()) == adata.n_obs, f"{sample_id}细胞ID重复！"
    assert len(adata.var_names.unique()) == adata.n_vars, f"{sample_id}基因名重复！"
    
    adatas.append(adata)

adata_merged = ad.concat(
     adatas,
     label='sample',
     keys=list(sample_meta.keys()), join='outer', fill_value=0
 )

# 标记线粒体与核糖体基因（人类）
adata_merged.var['mt'] = adata_merged.var_names.str.startswith('MT-')
adata_merged.var['ribo'] = adata_merged.var_names.str.startswith(('RPS', 'RPL'))

hemoglobin_genes = ['HBA1', 'HBA2', 'HBB', 'HBD', 'HBE1', 'HBG1', 'HBG2', 'HBZ', 'HBQ1']
# 仅保留数据中真实存在的血红蛋白基因
hb_genes_found = [g for g in hemoglobin_genes if g in adata_merged.var_names]
adata_merged.var['hb'] = False
if hb_genes_found:
    adata_merged.var.loc[hb_genes_found, 'hb'] = True
    print(f"检测到血红蛋白基因: {hb_genes_found}")
else:
    print("警告：未检测到任何血红蛋白基因，请检查基因注释！")

# 计算QC指标，包含线粒体、核糖体和血红蛋白
sc.pp.calculate_qc_metrics(adata_merged, qc_vars=['mt', 'ribo', 'hb'], inplace=True)

os.mkdir('pipelineTest1')
os.chdir('pipelineTest1')
os.mkdir(os.path.join(os.getcwd(), '00_Raw_data_after_quality_control'))

# 分样本QC小提琴图（核心诊断步骤）
fig1, axes1 = plt.subplots(2, 2, figsize=(18, 10))
for i, metric in enumerate(['n_genes_by_counts', 'total_counts']):
    sc.pl.violin(adata_merged, metric, groupby='sample', rotation=45, stripplot=False, ax=axes1[0, i], show=False, legend=None)
    # 不做全局统一阈值，而是观察每个样本的尾巴

for i, metric in enumerate(['pct_counts_mt', 'pct_counts_hb']):
    sc.pl.violin(adata_merged, metric, groupby='sample', rotation=45, stripplot=False, ax=axes1[1, i], show=False, legend=None)

fig1.tight_layout()
fig1.savefig(os.path.join(os.getcwd(), '00_Raw_data_after_quality_control', 'QC_metrics_beforeQC_violin.png'), bbox_inches='tight')

fig2, axes2 = plt.subplots(1, 3, figsize=(21, 5))   # 一行两列
sc.pl.scatter(adata_merged, x='n_genes_by_counts', y='pct_counts_mt',
              color='sampleID', ax=axes2[0], show=False)
axes2[0].set_title('Genes vs MT%')

sc.pl.scatter(adata_merged, x='n_genes_by_counts', y='pct_counts_hb',
              color='sampleID', ax=axes2[1], show=False)
axes2[1].set_title('Genes vs HB%')

sc.pl.scatter(adata_merged, x='total_counts', y='n_genes_by_counts',
              color='sampleID', ax=axes2[2], show=False)
axes2[2].set_title('UMI vs Genes')
fig2.tight_layout()
fig2.savefig(os.path.join(os.getcwd(), '00_Raw_data_after_quality_control', 'QC_metrics_beforeQC_scatter.png'), bbox_inches='tight')

min_genes=200,
max_genes=7500,
min_counts=500,
max_counts=50000,
max_mt_percent=20,
# 创建筛选掩码
cell_mask = (
        (adata_merged.obs.n_genes_by_counts >= min_genes) &
        (adata_merged.obs.n_genes_by_counts <= max_genes) &
        (adata_merged.obs.total_counts >= min_counts) &
        (adata_merged.obs.total_counts <= max_counts) &
        (adata_merged.obs.pct_counts_mt <= max_mt_percent) &
	(adata_merged.obs.pct_counts_hb <= 0.5)
    )

adata_filtered = adata_merged[cell_mask, :].copy()

adata_filtered.raw=adata_filtered
sc.pp.normalize_total(adata_filtered, target_sum=1e4)
sc.pp.log1p(adata_filtered)
sc.pp.highly_variable_genes(adata_filtered, n_top_genes=3000)
adata_hvg = adata_filtered[:, adata_filtered.var.highly_variable].copy()

adata_hvg.obs['doublet_score'] = np.nan
adata_hvg.obs['predicted_doublet'] = False
for s in adata_hvg.obs['sample'].unique():
    idx = adata_hvg.obs['sample'] == s
    tmp = adata_hvg[idx].copy()
    
    sc.pp.scrublet(tmp, expected_doublet_rate=0.06)
    
    adata_hvg.obs.loc[idx, 'doublet_score'] = tmp.obs['doublet_score']
    adata_hvg.obs.loc[idx, 'predicted_doublet'] = tmp.obs['predicted_doublet']

# 把结果回填
adata_filtered.obs['doublet_score'] = adata_hvg.obs['doublet_score']
adata_filtered.obs['predicted_doublet'] = adata_hvg.obs['predicted_doublet']

adata_doubletsRmv = adata_filtered[~adata_filtered.obs['predicted_doublet'].fillna(False).astype(bool).values].copy()

fig1, axes1 = plt.subplots(2, 2, figsize=(18, 10))
for i, metric in enumerate(['n_genes_by_counts', 'total_counts']):
    sc.pl.violin(adata_doubletsRmv , metric, groupby='sample', rotation=45, stripplot=False, ax=axes1[0, i], show=False, legend=None)
    # 不做全局统一阈值，而是观察每个样本的尾巴

for i, metric in enumerate(['pct_counts_mt', 'pct_counts_hb']):
    sc.pl.violin(adata_doubletsRmv , metric, groupby='sample', rotation=45, stripplot=False, ax=axes1[1, i], show=False, legend=None)

axes1[1,0].set_ylim(0, 100)
axes1[1,1].set_ylim(0, 100)
fig1.tight_layout()
fig1.savefig(os.path.join(os.getcwd(), '00_Raw_data_after_quality_control', 'QC_metrics_afterQC_violin.png'), bbox_inches='tight')

fig2, axes2 = plt.subplots(1, 3, figsize=(21, 5))   # 一行两列
sc.pl.scatter(adata_doubletsRmv, x='n_genes_by_counts', y='pct_counts_mt',
              color='sampleID', ax=axes2[0], show=False)

axes2[0].set_title('Genes vs MT%')

sc.pl.scatter(adata_doubletsRmv, x='n_genes_by_counts', y='pct_counts_hb',
              color='sampleID', ax=axes2[1], show=False)
axes2[1].set_title('Genes vs HB%')

sc.pl.scatter(adata_doubletsRmv, x='total_counts', y='n_genes_by_counts',
              color='sampleID', ax=axes2[2], show=False)
axes2[2].set_title('UMI vs Genes')
fig2.tight_layout()
fig2.savefig(os.path.join(os.getcwd(), '00_Raw_data_after_quality_control', 'QC_metrics_afterQC_scatte.png'), bbox_inches='tight')

os.mkdir(os.path.join(os.getcwd(), 'adataRDS'))
adata_doubletsRmv.write_h5ad(os.path.join(os.getcwd(), 'adataRDS', 'adataFiltrd.h5ad'))

QC_bySample_stat_df=pd.DataFrame(adata_merged.obs.groupby('sampleID')['sampleID'].count())
QC_bySample_stat_df.columns=['rawCellsCount']
QC_bySample_stat_df['pct_counts_mt>=20']=adata_merged.obs.query('pct_counts_mt>=20').groupby('sampleID')['sampleID'].count()
QC_bySample_stat_df['doublets']=adata_filtered.obs.query('predicted_doublet=="True"').groupby('sampleID')['sampleID'].count()
QC_bySample_stat_df['filterdCellsCount']=adata_doubletsRmv.obs.groupby('sampleID')['sampleID'].count()
QC_bySample_stat_df['pct_counts_hb>=0.5']=adata_merged.obs.query('pct_counts_hb>=0.5').groupby('sampleID')['sampleID'].count()
QC_bySample_stat_df['total_counts < 500 or >50000']=adata_merged.obs.query('total_counts <500 | total_counts >50000').groupby('sampleID')['sampleID'].count()
QC_bySample_stat_df['n_genes_by_counts < 200 or >7500']=adata_merged.obs.query('n_genes_by_counts <200 | total_counts >7500').groupby('sampleID')['sampleID'].count()
QC_bySample_stat_df['filterdCellsFraction(%)']=QC_bySample_stat_df[['rawCellsCount', 'filterdCellsCount']].apply(lambda x:round(100*x['filterdCellsCount']/x['rawCellsCount'], 2), axis=1)

QC_bySample_stat_df = QC_bySample_stat_df[['rawCellsCount', 'pct_counts_mt>=20', 'pct_counts_hb>=0.5', 'doublets', 'total_counts < 500 or >50000', 'n_genes_by_counts < 200 or >7500', 'filterdCellsCount', 'filterdCellsFraction(%)']]
QC_bySample_stat_df.to_excel(os.path.join(os.getcwd(), '00_Raw_data_after_quality_control', 'QC_bySample_stat_df.xlsx'))

sampleIDLst=adata_doubletsRmv.obs['sample'].unique().to_list()
for e in sampleIDLst:
	os.mkdir(os.path.join(os.getcwd(), '00_Raw_data_after_quality_control', e))	

for e in sampleIDLst:
	e
	cellIDs=adata_doubletsRmv.obs[adata_doubletsRmv.obs['sample']==e]['cellID'].to_list()
	filterdAdata = adata_doubletsRmv[adata_doubletsRmv.obs['cellID'].isin(cellIDs), :].copy()
	scu.write_mtx(filterdAdata, os.path.join(os.getcwd(), '00_Raw_data_after_quality_control', e))

new_adata = sc.pp.filter_genes(adata, min_cells=3, copy=True)
adata_doubletsRmv.n_vars
33282
new_adata.n_vars
18149


#鉴定高变基因
sc.pp.normalize_total(adata_doubletsRmv, target_sum=1e4)
sc.pp.log1p(adata_doubletsRmv)
sc.pp.highly_variable_genes(
    adata_doubletsRmv,
    n_top_genes=2000,
    min_mean=0.01,    # 必须放宽！拯救 TBX21/BCL6/RORC
    max_mean=5,
    #batch_key='group', 
	min_disp=0.5,
    flavor="seurat"
)

#文献中的marker genes！！！
cd4_markers = {
    "T cell":       ["TRAC", "TRBC1", "TRBC2", "CD4", "CD8A"],
    "Naive T":	    ["TCF7", "LEF1", "CD28", "BCL2"],
    "Tcm":          ["CCR7", "SELL", "CD27", "IL7R", "CD44"],
    "Tem":          ["CD44", "GZMK", "CXCR3", "CCR6"],
    "Effector_CD4": ["GZMK", "GZMB", "NKG7", "CCL5", "PRF1"],
    "Temra":        ["GZMB", "PRF1", "KLRG1", "GNLY"],
    "TRM":          ["ITGAE", "CXCR6", "ITGA1", "ZNF683"],
    "Th1":          ["TBX21", "CXCR3", "IFNG"],
    "Th2":          ["GATA3", "IL4", "IL5", "IL13"],
    "Th17":         ["RORC", "CCR6", "IL17A"],
    "Treg":         ["FOXP3", "IL2RA", "CTLA4", "IKZF2"],
    "cTfh":         ["BCL6", "CXCR5", "ICOS", "IL21", "SH2D1A"],
}

markersLst=[]
for k, v in cd4_markers.items():
	markersLst.extend(v)

for gene in markersLst:
    if gene in adata_doubletsRmv.var_names:
        adata_doubletsRmv.var.loc[gene, 'highly_variable'] = True


#批次矫正（需选择）
os.makedirs('01_Primary_clustering/01_批次校正')
os.chdir('01_Primary_clustering/01_批次校正')

sc.pp.scale(adata_doubletsRmv, max_value=10)
sc.tl.pca(adata_doubletsRmv, n_comps=50, svd_solver="arpack", use_highly_variable=True)
sc.pp.neighbors(adata_doubletsRmv, key_added='noBC')
sc.tl.leiden(adata_doubletsRmv, neighbors_key='noBC', key_added='noBC')
sc.tl.umap(
	adata_doubletsRmv, neighbors_key='noBC', key_added='X_noBC'
)

sc.pl.embedding(
    adata_doubletsRmv,
    basis="noBC",   # 对应 obsm['X_noBC']
    color=["noBC", "group"],
    title="UMAP (no Batch Correction)",
    save="_umap.png"
)

#批次矫正
pca_matrix = adata_doubletsRmv.obsm["X_pca"]
meta_data = adata_doubletsRmv.obs
# 运行 Harmony
harmony_obj = run_harmony(
    data_mat=pca_matrix.T,  # 必须转置
    meta_data=meta_data,
    vars_use=["group"],     # 你的批次列
    nclust=10,
    max_iter_harmony=10,
    verbose=True,
)
adata_doubletsRmv.obsm["X_pca_harmony"] = harmony_obj.Z_corr

sc.pp.neighbors(adata_doubletsRmv, key_added='BC', use_rep="X_pca_harmony")
sc.tl.leiden(adata_doubletsRmv, neighbors_key='BC', key_added='BC')
sc.tl.umap(
	adata_doubletsRmv, neighbors_key='BC', key_added='BC'
)

sc.pl.embedding(
    adata_doubletsRmv,
    basis="BC",   # 对应 obsm['X_noBC']
    color=["BC", "group"],
    title="UMAP (Batch Correction)",
    save="_umap.png"
)


#细胞聚类，！
os.makedirs('../02_细胞聚类/count')
os.chdir('../')
sc.pp.neighbors(adata_doubletsRmv, use_rep="X_pca")

resLst=[5*(e+3)/100 for e in range(8)] #分辨率范围
for r in resLst:
	'leiden_r0'+ str(int(r*100))
	sc.tl.leiden(adata_doubletsRmv, resolution=r, key_added='leiden_r0'+ str(int(r*100)))

adata_doubletsRmv_r2counts_df_dict={}
for r in resLst:
	res='leiden_r0'+ str(int(r*100))
	adata_doubletsRmv_r2counts_df=pd.DataFrame(adata_doubletsRmv.obs.groupby(res)[res].count())
	adata_doubletsRmv_r2counts_df.columns=[res]
	adata_doubletsRmv_r2counts_df_dict[res]=adata_doubletsRmv_r2counts_df

for k, v in adata_doubletsRmv_r2counts_df_dict.items():
	fileName=k.replace('leiden_r', 'resolution')+'_clusterCellCount_re_20260606.csv'
	fileName
	#v.shape
	v.columns=['cellsCount']
	v.reset_index(inplace=True)
	v.rename(columns={v.columns[0]:'clusterID'}, inplace=True)
	v.to_csv(os.path.join(os.getcwd(), '02_细胞聚类/count', fileName), index=False)

clusterDir=os.path.join(os.getcwd(), '02_细胞聚类')
adata_doubletsRmv_r2clustersCount_dict={format(int(k.replace('leiden_r', ''))/100, '.2f'):v.shape[0] for k, v in adata_doubletsRmv_r2counts_df_dict.items()}
adata_doubletsRmv_r2clustersCount_df = pd.DataFrame.from_dict(adata_doubletsRmv_r2clustersCount_dict, orient='index')
adata_doubletsRmv_r2clustersCount_df.columns=['clusterCount']
adata_doubletsRmv_r2clustersCount_df.reset_index(inplace=True)
adata_doubletsRmv_r2clustersCount_df.rename(columns={'index':'resolution'}, inplace=True)
adata_doubletsRmv_r2clustersCount_df.to_csv(os.path.join(clusterDir, 'count', 'resolution_allCount.csv'), index=False)

os.chdir(clusterDir)
sc.pp.neighbors(adata_doubletsRmv, use_rep="X_pca")
sc.tl.umap(adata_doubletsRmv)
for k in adata_doubletsRmv_r2counts_df_dict:	
	sc.pl.umap(
		adata_doubletsRmv,
		color=[k],
		title=['UMAP (Clusters)'],
		save='_' + k + '.png'
	)

for k in adata_doubletsRmv_r2counts_df_dict:
	sc.pl.dotplot(adata_doubletsRmv, groupby=k, standard_scale="var", var_names=cd4_markers, save=k+'_allTmarkers.png')

for k in adata_doubletsRmv_r2counts_df_dict:
	# 差异基因分析
	sc.tl.rank_genes_groups(adata_doubletsRmv, groupby=k, method="wilcoxon")
	# 查看结果
	sc.pl.rank_genes_groups_dotplot(adata_doubletsRmv, groupby=k, standard_scale="var", n_genes=5, save=f'{k}_top5.png')
	

#选定分辨率
os.makedirs(os.path.join(clusterDir, 'marker_gene'))
sc.tl.rank_genes_groups(adata_doubletsRmv, groupby='leiden_r025', use_raw=False, pts=True, method="wilcoxon")
r025_markers_df=sc.get.rank_genes_groups_df(adata_doubletsRmv, group=None)
r025_markers_df.rename(columns={'group':'clusterID'}, inplace=True)
r025_markers_df.to_excel(os.path.join(clusterDir, 'marker_gene/leiden_r025_all.xlsx'), index=False)

r025_markers_df.rename(columns={'names':'genes'}, inplace=True)
r025_clustersLst=r025_markers_df['clusterID'].unique().tolist()
cluster2marker_top5_dict={}
for e in r025_clustersLst:
	top5_df=r025_markers_df[r025_markers_df['clusterID']==e].head(5)
	top200_df=r025_markers_df[r025_markers_df['clusterID']==e].head(200)
	cluster2marker_top5_dict[e]=top5_df['genes'].tolist()
	top5_df.to_csv(os.path.join(clusterDir, f'marker_gene/r025_Cluster{e}_top5.csv'), index=False)
	top200_df['genes'].to_csv(os.path.join(clusterDir, f'marker_gene/r025_Cluster{e}_top200.csv'), index=False)

#手工选定聚类的分辨率后：
for cluster in r025_clustersLst:
	f'cluster{cluster}'
	top5=cluster2marker_top5_dict[cluster]
	sc.pl.umap(
		adata_doubletsRmv,
		color=top5,
		ncols=3,
		save='_r025_cluster' +cluster + '_top5.png'
	)

correlation_matrix, cluster_expr_df = calculate_cluster_correlation_expression(adata_doubletsRmv, 'leiden_r025')

os.chdir('../..')
plot_correlation_heatmap(correlation_matrix, out_file=os.path.join(os.getcwd(), "01_Primary_clustering/02_细胞聚类/figures/leidenR025_cellClusters_correlationHeatmap_20260605.png"), title="Cluster Correlation Heatmap", annot_kws={'size':8})


#03_细胞注释: 需手动选择模型！！！
os.makedirs('01_Primary_clustering/03_细胞注释')
os.makedirs('01_Primary_clustering/03_细胞注释/stat')
os.chdir('01_Primary_clustering/03_细胞注释')

adata_doubletsRmv.obs.reset_index(inplace=True)
adata_doubletsRmv.obs.rename(columns={'index':'barcode'}, inplace=True)
adata_doubletsRmv.obs.set_index('cellID', inplace=True)

immuneLow_mod=models.Model.load(model='Immune_All_Low.pkl') 
immuneLowPred_mv = celltypist.annotate(adata_doubletsRmv, model=immuneLow_mod, majority_voting=True, mode='prob match')
adata_doubletsRmv.obs[['mv_predLables', 'mv_over_clustering', 'majority_voting']]=immuneLowPred_mv.predicted_labels[['predicted_labels', 'over_clustering', 'majority_voting']]

mvPred_old2new_dict={e:f'mv_{e}'  for e in ['predicted_labels', 'over_clustering', 'majority_voting', 'conf_score']}
adata_doubletsRmv.obs.rename(columns=mvPred_old2new_dict, inplace=True)

bmPreds = celltypist.annotate(adata_doubletsRmv, model=immuneLow_mod, mode='best match')
bmPred_adata = bmPreds.to_adata()
adata_doubletsRmv.obs[['bm_predLables', 'bm_conf_score']]=bmPred_adata.obs[['predicted_labels', 'conf_score']]

os.chdir('../../')
adata_doubletsRmv.write_h5ad(os.path.join(os.getcwd(), 'adataRDS', 'adataFiltrd_celltypist.h5ad'))

# 先查看每组细胞数量（可选）
print(adata_doubletsRmv.obs['bm_predLables'].value_counts())

# 保留至少 10 个细胞的群（数字自己改）
min_cells_per_cluster = 10
clusters_to_keep = adata_doubletsRmv.obs['bm_predLables'].value_counts()
clusters_to_keep = clusters_to_keep[clusters_to_keep >= min_cells_per_cluster].index.tolist()

# 过滤adata
adata_filtered = adata_doubletsRmv[adata_doubletsRmv.obs['bm_predLables'].isin(clusters_to_keep)].copy()

# 再跑差异分析
sc.tl.rank_genes_groups(
    adata_filtered,
    groupby="bm_predLables",
    use_raw=True
)

sc.pl.rank_genes_groups_dotplot(adata_filtered, groupby="bm_predLables", standard_scale="var", n_genes=5, dendrogram=False, save='leidenR025_bm_noDendrogram.png')

correlation_matrix, cluster_expr_df = calculate_cluster_correlation_expression(adata_doubletsRmv, 'bm_predLables')
plot_correlation_heatmap(correlation_matrix, out_file=os.path.join(os.getcwd(), "figures/celltypistBM_cellClusters_correlationHeatmap.png"), title="Cluster Correlation Heatmap", annot_kws={'size':8})

#
adata_doubletsRmv.obs.groupby('bm_predLables')['bm_predLables'].count().sort_values(ascending=False).to_excel('stat/celltypist_cellAnnotStat.xlsx')
correlation_matrix.to_excel('stat/celltypist_cellGrps_correlation_matrix.xlsx')


#04_样本或组汇总
grpComp_dir=os.path.join(os.getcwd(), '01_Primary_clustering/04_样本或组汇总')
os.mkdir(grpComp_dir)

os.chdir(grpComp_dir)
sc.pl.umap(
	adata_doubletsRmv,
	color='sampleID',
	title=['UMAP (Clusters)'],
	save='_r025_sampleID.png'
	)

sc.pl.umap(
	adata_doubletsRmv,
	color='group',
	title=['UMAP (Clusters)'],
	save='_r025_group.png'
	)

os.mkdir(os.path.join(os.getcwd(), 'count'))

group_r025_clusterCellsCounts_df = pd.crosstab(adata_doubletsRmv.obs['group'], adata_doubletsRmv.obs['leiden_r025'])
group_r025_clusterCellsProp_df = group_r025_clusterCellsCounts_df.T.div(group_r025_clusterCellsCounts_df.sum(axis=1), axis=1)
group_r025_clusterCellsProp_df.reset_index(inplace=True)
group_r025_clusterCellsProp_df.rename(columns={'leiden_r025':'clusterID'}, inplace=True)
group_r025_clusterCellsProp_df.shape
(5, 6)

group_r025_clusterCellsCounts_df=group_r025_clusterCellsCounts_df.T
group_r025_clusterCellsCounts_df.shape
(5, 5)
group_r025_clusterCellsCounts_df.reset_index(inplace=True)
group_r025_clusterCellsCounts_df.rename(columns={'leiden_r025':'clusterID'}, inplace=True)

list(cluster2marker_top5_dict.items())[:2]
[('0', ['RPS4Y1', 'RPS8', 'RPL6', 'RPL17', 'RPL27A']), ('1', ['S100A4', 'AHNAK', 'SYNE2', 'KLF6', 'ITGB1'])]

group_r025_clusterCellsCounts_df['Core_Markers']=group_r025_clusterCellsCounts_df['clusterID'].apply(lambda x:','.join(cluster2marker_top5_dict[x]))
group_r025_clusterCellsCounts_df.set_index('clusterID', inplace=True)
group_r025_clusterCellsCounts_df['Cell_Count']=group_r025_clusterCellsCounts_df[group_r025_clusterCellsCounts_df.columns[0:5]].sum(axis=1)

group_r025_clusterCellsCounts_df['Cell_Percentage']=group_r025_clusterCellsCounts_df['Cell_Count'].apply(lambda x:round(100*x/138683, 2))
group_r025_clusterCellsCounts_df.columns
Index(['Class_III', 'Class_III-IV', 'Class_IV-V', 'HD', 'NoLN', 'Core_Markers',
       'Cell_Count', 'Cell_Percentage'],
      dtype='object', name='group')

group_r025_clusterCellsCounts1_df = pd.crosstab(adata_doubletsRmv.obs['group'], adata_doubletsRmv.obs['leiden_r025'])
group_r025_clusterCellsProp_df = group_r025_clusterCellsCounts1_df.T.div(group_r025_clusterCellsCounts1_df.sum(axis=1), axis=1)
group_r025_clusterCellsProp_df.reset_index(inplace=True)
group_r025_clusterCellsProp_df.rename(columns={'leiden_r025':'clusterID'}, inplace=True)

grpOld2new_dict={e:f'{e.rstrip()}(%)' for e in group_r025_clusterCellsProp_df.columns[1:]}
grpOld2new_dict
{'Class_III': 'Class_III(%)', 'Class_III-IV': 'Class_III-IV(%)', 'Class_IV-V': 'Class_IV-V(%)', 'HD ': 'HD(%)', 'No LN': 'No LN(%)'}
group_r025_clusterCellsProp_df.rename(columns=grpOld2new_dict, inplace=True)
group_r025_clusterCellsProp_df.set_index('clusterID', inplace=True)
group_r025_clusterCellsProp_df=group_r025_clusterCellsProp_df.map(lambda x:round(100*x, 2))

group_r025_clusterCellsStat_df=pd.concat([group_r025_clusterCellsCounts_df, group_r025_clusterCellsProp_df], axis=1)
group_r025_clusterCellsStat_df.reset_index(inplace=True)
group_r025_clusterCellsStat_df.rename(columns={'clusterID':'Cluster'}, inplace=True)

new_order = ['Cluster', 'Cell_Count', 'Cell_Percentage', 'Core_Markers', 'HD', 'NoLN', 'Class_III', 'Class_III-IV', 'Class_IV-V', 'HD(%)', 'NoLN(%)', 'Class_III(%)', 'Class_III-IV(%)', 'Class_IV-V(%)']
group_r025_clusterCellsStat_df = group_r025_clusterCellsStat_df[new_order]
group_r025_clusterCellsStat_df.columns
Index(['Cluster', 'Cell_Count', 'Cell_Percentage', 'Core_Markers', 'HD',
       'NoLN', 'Class_III', 'Class_III-IV', 'Class_IV-V', 'HD(%)', 'NoLN(%)',
       'Class_III(%)', 'Class_III-IV(%)', 'Class_IV-V(%)'],
      dtype='object', name='group')

group_r025_clusterCellsCounts_df.to_excel('count/group_r025_clusterCellsCountAnno_df.xlsx')

grpsLst=['HD', 'NoLN', 'Class_III', 'Class_III-IV', 'Class_IV-V']
grpPairs = list(combinations(grpsLst, 2))  # 2表示两两组合

os.mkdir(os.path.join(os.getcwd(), 'diffExpre'))
for e in grpPairs:
	e
	sc.tl.rank_genes_groups(adata_doubletsRmv, groupby='group', reference=e[0], groups=[e[1]], pts=True, method="wilcoxon")
	grpComp_markers_df=sc.get.rank_genes_groups_df(adata_doubletsRmv, group=None)
	#grpComp_markers_df.to_excel(os.path.join(os.getcwd(), 'diffExpre', f'{e[0].rstrip()}_vs_{e[1]}_all.xlsx'))


sample_r025_clusterCellsCounts_df = pd.crosstab(adata_doubletsRmv.obs['sampleID'], adata_doubletsRmv.obs['leiden_r025'])
sample_r025_clusterCellsProp_df = sample_r025_clusterCellsCounts_df.T.div(sample_r025_clusterCellsCounts_df.sum(axis=1), axis=1)
sample_r025_clusterCellsProp_df.reset_index(inplace=True)
sample_r025_clusterCellsProp_df.rename(columns={'leiden_r025':'clusterID'}, inplace=True)

sample_r025_clusterCellsProp_df.set_index('clusterID', inplace=True)
sample_r025_clusterCellsProp_df.T.plot(kind='bar', stacked=True, colormap='tab20')
plt.title('Cluster composition of samples', fontweight='bold')
plt.xlabel('sample')
plt.ylabel('Proportion')
plt.legend(title='Cluster', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tick_params(axis='x', rotation=45)
plt.savefig("figures/cellsamples_bars.png", bbox_inches="tight")

sample_r025_clusterCellsCounts_df=sample_r025_clusterCellsCounts_df.T
sample_r025_clusterCellsCounts_df['Cell_Count']=sample_r025_clusterCellsCounts_df[sample_r025_clusterCellsCounts_df.columns].sum(axis=1)
sample_r025_clusterCellsCounts_df['Cell_Percentage']=sample_r025_clusterCellsCounts_df['Cell_Count'].apply(lambda x:round(100*x/138683, 2))

sample_r025_clusterCellsCounts_df.reset_index(inplace=True)
sample_r025_clusterCellsCounts_df['Core_Markers']=sample_r025_clusterCellsCounts_df['leiden_r025'].apply(lambda x:','.join(cluster2marker_top5_dict[x]))

colsSorted=['leiden_r025', 'Cell_Count', 'Cell_Percentage', 'Core_Markers']
colsSorted.extend(sample_r025_clusterCellsCounts_df.columns[1:14])

sample_r025_clusterCellsCounts_df = sample_r025_clusterCellsCounts_df.loc[:, colsSorted]
sample_r025_clusterCellsCounts_df['Cluster']=sample_r025_clusterCellsCounts_df['leiden_r025'].apply(lambda x:x)
sample_r025_clusterCellsCounts_df.set_index('leiden_r025', inplace=True)

sample_r025_clusterCellsProp_df.set_index('clusterID', inplace=True)
sample_r025_clusterCellsProp_df1=sample_r025_clusterCellsProp_df.map(lambda x:round(100*x, 2))
sample_r025_clusterCellsProp_df1.rename(columns={e:f'{e}(%)' for e in sample_r025_clusterCellsProp_df1.columns}, inplace=True)

sample_r025_clusterCellsStat_df=pd.concat([sample_r025_clusterCellsCounts_df, sample_r025_clusterCellsProp_df1], axis=1)
sample_r025_clusterCellsStat_df.shape
(5, 30)
sample_r025_clusterCellsStat_df.columns
Index(['Cell_Count', 'Cell_Percentage', 'Core_Markers', 'HD3', 'HD5', 'HD8',
       'SLE1', 'SLE2', 'SLE3', 'SLE4', 'SLE5', 'SLE6', 'SLE7', 'SLE8', 'SLE12',
       'SLE15', 'Cluster', 'HD3(%)', 'HD5(%)', 'HD8(%)', 'SLE1(%)', 'SLE2(%)',
       'SLE3(%)', 'SLE4(%)', 'SLE5(%)', 'SLE6(%)', 'SLE7(%)', 'SLE8(%)',
       'SLE12(%)', 'SLE15(%)'],
      dtype='object', name='sampleID')

sample_r025_clusterCellsStat_df.reset_index(inplace=True)
sample_r025_clusterCellsStat_df.rename(columns={'index':'clusterID'}, inplace=True)
sample_r025_clusterCellsStat_df.to_csv('count/sample_r025_clusterCellsStatAnno_df.csv', index=False)

padj_threshold = 0.05
logfc_threshold = 1
plt.figure(figsize=(10, 8))

grpComp2markers_dict={}
for e in grpPairs:
	e
	f'04_样本或组汇总/diffExpre/{e[0].rstrip().replace(" ", "-")}_vs_{e[1].rstrip().replace(" ", "-")}_all_20260522.xlsx'
	
	sc.tl.rank_genes_groups(adata_doubletsRmv, groupby='group', reference=e[0], groups=[e[1]], use_raw=True, pts=True, rankby_abs=True, method="wilcoxon")
	grpComp_markers_df=sc.get.rank_genes_groups_df(adata_doubletsRmv, group=None)
	
	grpsComp_cells=(
	(adata_doubletsRmv.obs.group ==e[0]) |
        (adata_doubletsRmv.obs.group ==e[1])
	)
	adata_grpsComp=adata_doubletsRmv[grpsComp_cells, :].copy()
	
	grpComp_markers_df['abs_logFC']=grpComp_markers_df['logfoldchanges'].apply(lambda x:abs(x))
	top50=grpComp_markers_df.query('pvals_adj<0.05 & logfoldchanges.abs() >=1').sort_values('abs_logFC', ascending=False)['names'].head(50).tolist()
	
	#sc.pl.heatmap(adata_grpsComp, top50, use_raw=False, groupby='group', swap_axes=True, vmin=-5, vmax=5, save=f'_{e[0].rstrip().replace(" ", "-")}_vs_{e[1].rstrip().replace(" ", "-")}_top50.png')
	grpComp2markers_dict[e]=grpComp_markers_df

samplesLst=sample_r025_clusterCellsCounts_df.columns[3:16].tolist()
for e in samplesLst:
	adataGrp=adata_doubletsRmv[adata_doubletsRmv.obs.sampleID ==e, :].copy()
	e
	sc.pl.umap(
		adataGrp,
		color='leiden_r025',
		title=['UMAP (Clusters)'],
		save=f'_r025_{e}.png'
	)

os.chdir('04_样本或组汇总')

sample_r025_clusterCellsPro_df=sample_r025_clusterCellsStat_df[[e for e in sample_r025_clusterCellsStat_df.columns if '%' in e]]
sample_r025_clusterCellsPro_df['clusterID']=sample_r025_clusterCellsStat_df['clusterID']
sample_r025_clusterCellsPro_df.set_index('clusterID', inplace=True)

clusterProp_amongConditions_stat_dict={}
for cluster in sample_r025_clusterCellsPro_df.index:
	clusterProp_amongConditions_stat_dict[cluster]=celltypeProp_amongConditions_KWTesting_20260609(sample_r025_clusterCellsPro_df, cluster, sample2grp_dict)

sampleGrp_df=adata_doubletsRmv.obs[['sampleID', 'group']]
sampleGrp_df=sampleGrp_df.drop_duplicates()
sampleGrp_df.set_index('sampleID', inplace=True)
sample2grp_dict=sampleGrp_df.to_dict()

clusterProp_amongConditions_stat_dict={}
for cluster in sample_r025_clusterCellsPro_df.index:
	clusterProp_amongConditions_stat_dict[cluster]=celltypeProp_amongConditions_KWTesting_20260605(sample_r025_clusterCellsPro_df, cluster, sample2grp_dict)
	
conditions_clusterProp_kw_dict={k:v[:2] for k, v in clusterProp_amongConditions_stat_dict.items()}
conditions_clusterProp_kw_df=pd.DataFrame.from_dict(conditions_clusterProp_kw_dict, orient='index')
conditions_clusterProp_kw_df.columns=['stat', 'p_value']

os.mkdir('stat')
conditions_clusterProp_kw_df['clusterID']=conditions_clusterProp_kw_df.index
conditions_clusterProp_kw_df=conditions_clusterProp_kw_df[['clusterID', 'stat', 'p_value']]
conditions_clusterProp_kw_df.to_excel('stat/conditions_clusterProp_kw_df.xlsx', index=False)

clusterProp_amongConditions_sigDunn_dict={k:v[2] for k, v in clusterProp_amongConditions_stat_dict.items() if v[1] <0.05}
for k, v in clusterProp_amongConditions_stat_dict.items():
	v[2].to_excel(f'stat/Cluster{k}_proportion_conditionsDunn.xlsx')


def calculate_cluster_correlation_expression(adata, cluster_key, use_raw=True):
    """
    基于基因表达均值计算细胞亚群相关性
    
    参数:
    adata: AnnData对象
    cluster_key: 细胞亚群注释的obs列名
    use_raw: 是否使用raw属性中的表达数据
    
    返回:
    correlation_matrix: 相关性矩阵
    """
    # 确保数据已归一化
    if use_raw and adata.raw is not None:
        expr_matrix = adata.raw.X.T
        var_names = adata.raw.var_names
    else:
        expr_matrix = adata.X.T
        var_names = adata.var_names
    
    # 将表达矩阵转换为DataFrame
    expr_df = pd.DataFrame(expr_matrix.toarray() if hasattr(expr_matrix, 'toarray') 
                          else expr_matrix, 
                          index=var_names, 
                          columns=adata.obs_names)
    
    # 按细胞亚群分组计算均值表达
    cluster_means = []
    clusters = adata.obs[cluster_key].unique()
    
    for cluster in clusters:
        cluster_cells = adata.obs[adata.obs[cluster_key] == cluster].index
        cluster_mean = expr_df[cluster_cells].mean(axis=1)
        cluster_means.append(cluster_mean)
    
    # 构建亚群表达矩阵
    cluster_expr_df = pd.concat(cluster_means, axis=1)
    cluster_expr_df.columns = clusters
    
    # 计算相关性矩阵
    correlation_matrix = cluster_expr_df.corr(method='pearson')
    
    return correlation_matrix, cluster_expr_df 

def plot_correlation_heatmap(correlation_matrix, out_file, annot_kws, title="Cluster Correlation Heatmap", figsize=(16, 12), cmap='coolwarm', annot=True):
    """
    绘制相关性热图
    """
    plt.figure(figsize=figsize)
    
    # 创建掩码矩阵，隐藏上三角（可选）
    mask = np.triu(np.ones_like(correlation_matrix, dtype=bool))
    
    sns.heatmap(correlation_matrix, 
                #mask=mask if annot else None,  # 只在需要标注时使用掩码
                cmap=cmap,
                annot=annot,
		annot_kws=annot_kws,
                fmt='.2f',
                square=True,
                linewidths=0.5,
                cbar_kws={'label': 'Pearson correlation coefficient​'})
    
    plt.title(title, fontsize=16, pad=20)
    plt.tight_layout()
    plt.savefig(out_file, bbox_inches="tight")

def celltypeProp_amongConditions_KWTesting_20260609(samples_celltypesProp_df, cellType, sample2condition_dict):
	sample2celltypeProp_dict=samples_celltypesProp_df.loc[cellType].to_dict()
	
	grp2celltypeProp_dict={}
	for k, v in sample2celltypeProp_dict.items():
		grp=sample2condition_dict[k.split('(')[0]]
		grp
		if grp in grp2celltypeProp_dict:
			grp2celltypeProp_dict[grp].append(v)
		else:
			grp2celltypeProp_dict[grp]=[v]
	
	kw_stat, kw_p = stats.kruskal(*grp2celltypeProp_dict.values())
	
	dfs = []
	for group_name, values in grp2celltypeProp_dict.items():
		df_temp = pd.DataFrame({
			"value": values, "group": group_name
			})
		dfs.append(df_temp)
	condition2celltypeProp_df = pd.concat(dfs, ignore_index=True)
	condition2celltypeProp_dunn = sp.posthoc_dunn(condition2celltypeProp_df, val_col='value', group_col='group', p_adjust='holm')
	return kw_stat, kw_p, condition2celltypeProp_dunn

