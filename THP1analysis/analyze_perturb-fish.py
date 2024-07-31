import numpy as np
import pandas as pd
import anndata as ad
from matplotlib import pyplot as plt
from scipy.spatial import distance
import os
from collections import defaultdict
import seaborn as sns
from scipy.stats import ranksums
from scipy.cluster import hierarchy 

#this script explores the data to find the best parameters to use when running FR-perturb
#then it runs it on each dataset to also check correlation. #also looks at cell density effects
def print_correlations(cov_string,dir='zombie_sample2',return_values=False,q_thresh=0.1):
	LFC = pd.read_csv('%s/%s_LFCs.txt' % (dir,cov_string),sep='\t',index_col=0).drop(columns=['Control'])
	Q = pd.read_csv('%s/%s_qvals.txt' % (dir,cov_string),sep='\t',index_col=0).drop(columns=['Control'])
	x = LFC.values.flatten()
	q = Q.values.flatten()
	for zz in [(q_thresh,1),(1,q_thresh),(q_thresh,q_thresh)]:
		idx = np.where((q<zz[0])*(qd<zz[1]))[0]
		corr = 1-distance.correlation(x[idx],y[idx])
		print(zz, len(idx),corr)
	if return_values:
		return LFC,Q,x,q


def write_data_subset(idx,path):
	adata = ad.AnnData(M.iloc[idx])
	adata.obs_names = [f"Cell_{i+1:d}" for i in idx]
	adata.obs['area'] = C.iloc[idx]['area'].values
	adata.obs['eccentricity'] = C.iloc[idx]['eccentricity'].values
	adata.obs['circularity'] = C.iloc[idx]['circularity'].values
	adata.obs['extent'] = C.iloc[idx]['extent'].values
	adata.obs['num_neighbors'] = N.values[idx,0]
	adata.obs['total_counts'] = adata.X.sum(1)
	adata.write('%scount_data.h5ad' % path)
	P.iloc[:,idx].to_csv('%sperturbation_design.csv' % path,sep=' ')


def write_data_subset_nocovs(idx,path):
	adata = ad.AnnData(M.iloc[idx])
	adata.obs_names = [f"Cell_{i+1:d}" for i in idx]
	adata.obs['num_neighbors'] = N.values[idx,0]
	adata.obs['total_counts'] = adata.X.sum(1)
	adata.write('%scount_data.h5ad' % path)
	P.iloc[:,idx].to_csv('%sperturbation_design.csv' % path,sep=' ')


def check_consistency(dir,q_thresh=0.1):
	LFC1 = pd.read_csv('%s/split1_LFCs.txt' % (dir),sep='\t',index_col=0).drop(columns=['Control'])
	Q1 = pd.read_csv('%s/split1_qvals.txt' % (dir),sep='\t',index_col=0).drop(columns=['Control'])
	LFC2 = pd.read_csv('%s/split2_LFCs.txt' % (dir),sep='\t',index_col=0).drop(columns=['Control'])
	Q2 = pd.read_csv('%s/split2_qvals.txt' % (dir),sep='\t',index_col=0).drop(columns=['Control'])
	x1 = LFC1.values.flatten()
	q1 = Q1.values.flatten()
	x2 = LFC2.values.flatten()
	q2 = Q2.values.flatten()
	idx_1 = np.where((q1<q_thresh)*(q2<1))[0]
	corr_1 = 1-distance.correlation(x1[idx_1],x2[idx_1])
	idx_both = np.where((q1<q_thresh)*(q2<q_thresh))[0]
	corr_both = 1-distance.correlation(x1[idx_both],x2[idx_both])
	return corr_1,corr_both,LFC1,Q1,LFC2,Q2


def run_random_split(dir,num_permute=10000,covariates='total_counts',rank=30):
	num_cells = M.shape[0]
	idx_sort = np.argsort(np.random.random(num_cells))
	idx1 = idx_sort[:int(num_cells/2)]
	idx2 = idx_sort[len(idx1):]
	write_data_subset(idx1,'%s/split1_' % dir)
	write_data_subset(idx2,'%s/split2_' % dir)
	_=os.system('python ./run_FR_Perturb.py --input-h5ad %s/split1_count_data.h5ad --input-perturbation-matrix %s/split1_perturbation_design.csv --control-perturbation-name Control --out %s/split1 --covariates %s --compute-pval --num-perms %d --rank %d --lambda1 0.0 --lambda2 0.0' % (dir,dir,dir,covariates,num_permute,rank))
	_=os.system('python ./run_FR_Perturb.py --input-h5ad %s/split2_count_data.h5ad --input-perturbation-matrix %s/split2_perturbation_design.csv --control-perturbation-name Control --out %s/split2 --covariates %s --compute-pval --num-perms %d --rank %d --lambda1 0.0 --lambda2 0.0' % (dir,dir,dir,covariates,num_permute,rank))

def run_random_split_withoutnorm(dir,num_permute=10000,covariates='total_counts',rank=30):
	num_cells = M.shape[0]
	idx_sort = np.argsort(np.random.random(num_cells))
	idx1 = idx_sort[:int(num_cells/2)]
	idx2 = idx_sort[len(idx1):]
	write_data_subset_nocovs(idx1,'%s/split1_' % dir)
	write_data_subset_nocovs(idx2,'%s/split2_' % dir)
	_=os.system('python ./run_FR_PerturbnoNorm.py --input-h5ad %s/split1_count_data.h5ad --input-perturbation-matrix %s/split1_perturbation_design.csv --control-perturbation-name Control --out %s/split1 --covariates %s --compute-pval --num-perms %d --rank %d --lambda1 0.0 --lambda2 0.0' % (dir,dir,dir,covariates,num_permute,rank))
	_=os.system('python ./run_FR_PerturbnoNorm.py --input-h5ad %s/split2_count_data.h5ad --input-perturbation-matrix %s/split2_perturbation_design.csv --control-perturbation-name Control --out %s/split2 --covariates %s --compute-pval --num-perms %d --rank %d --lambda1 0.0 --lambda2 0.0' % (dir,dir,dir,covariates,num_permute,rank))

def perturbation_violins(gene,target):
	Z = M[[gene]] + 1
	p = []
	for c in P:
		if P.loc['Control',c]:
			p.append('Control')
		elif P.loc[target,c]:
			p.append(target)
		else:
			p.append('Other')
	Z['perturbation'] = p
	Z = Z.drop(Z[Z.perturbation == 'Other'].index)
	w = ranksums(Z[(Z['perturbation'] == 'Control')][gene], Z[(Z['perturbation'] == target)][gene])
	print('wilcoxon: %f' % w[1], 'FR-Perturb qvalue: %f' % res[1].loc[gene,target])
	_=sns.violinplot(data=Z, x='perturbation', y=gene, log_scale=True)
	plt.show()


def density_violins(gene,target,show=False,savepath=None):
	Z = M[[gene]] + 1
	Z['density'] = ['low' if n<3 else 'high' for n in N.values]
	p = []
	for c in P:
		if P.loc['Control',c]:
			p.append('Control')
		elif P.loc[target,c]:
			p.append(target)
		else:
			p.append('Other')
	Z['perturbation'] = p
	Z = Z.drop(Z[Z.perturbation == 'Other'].index)
	w_low = ranksums(Z[(Z['density'] == 'low') & (Z['perturbation'] == 'Control')][gene], Z[(Z['density'] == 'low') & (Z['perturbation'] == target)][gene])
	w_high = ranksums(Z[(Z['density'] == 'high') & (Z['perturbation'] == 'Control')][gene], Z[(Z['density'] == 'high') & (Z['perturbation'] == target)][gene])
	print('low density wilcoxon: %f' % w_low[1], 'FR-Perturb qvalue: %f' % Q_low.loc[gene,target])
	print('high density wilcoxon: %f' % w_high[1], 'FR-Perturb qvalue: %f' % Q_high.loc[gene,target])
	_=sns.violinplot(data=Z, x="density", y=gene, hue="perturbation", log_scale=True, split=True, gap=0.1)
	if savepath is not None:
		plt.savefig(savepath)
		plt.close()
	elif show:
		plt.show()

def density_violins_eps(gene,target,show=False,savepath=None):
	Z = M[[gene]] + 1
	Z['density'] = ['low' if n<3 else 'high' for n in N.values]
	p = []
	for c in P:
		if P.loc['Control',c]:
			p.append('Control')
		elif P.loc[target,c]:
			p.append(target)
		else:
			p.append('Other')
	Z['perturbation'] = p
	Z = Z.drop(Z[Z.perturbation == 'Other'].index)
	w_low = ranksums(Z[(Z['density'] == 'low') & (Z['perturbation'] == 'Control')][gene], Z[(Z['density'] == 'low') & (Z['perturbation'] == target)][gene])
	w_high = ranksums(Z[(Z['density'] == 'high') & (Z['perturbation'] == 'Control')][gene], Z[(Z['density'] == 'high') & (Z['perturbation'] == target)][gene])
	print('low density wilcoxon: %f' % w_low[1], 'FR-Perturb qvalue: %f' % Q_low.loc[gene,target])
	print('high density wilcoxon: %f' % w_high[1], 'FR-Perturb qvalue: %f' % Q_high.loc[gene,target])
	_=sns.violinplot(data=Z, x="density", y=gene, hue="perturbation", log_scale=True, split=True, gap=0.1)
	if savepath is not None:
		plt.savefig(savepath, format='eps')
		plt.close()
	elif show:
		plt.show()



M = pd.read_csv('zombie_sample1/merfish.csv') # cells x genes
P = pd.read_csv('zombie_sample1/intrinsicdesign.csv', index_col=0) # target x cells
C = pd.read_csv('zombie_sample1/covariates.csv') # cells x covariates
N = pd.read_csv('zombie_sample1/covnumberofneighbors.csv', header=None)

# 'HN1' is not in doug data
M = M.drop(columns=['HN1'])

adata = ad.AnnData(M)
adata.obs_names = [f"Cell_{i+1:d}" for i in range(adata.n_obs)]
adata.obs['area'] = C['area'].values
adata.obs['eccentricity'] = C['eccentricity'].values
adata.obs['circularity'] = C['circularity'].values
adata.obs['extent'] = C['extent'].values
adata.obs['num_neighbors'] = N.values[:,0]
adata.obs['total_counts'] = adata.X.sum(1)
adata.write('zombie_sample1/count_data.h5ad')

P.to_csv('zombie_sample1/perturbation_design.csv',sep=' ')


LFC_doug = pd.read_csv('perturbseq/neweffectsfromDoug.csv',index_col=0)
Q_doug = pd.read_csv('perturbseq/dougsQvalues.csv',index_col=0)
y = LFC_doug.values.flatten()
qd = Q_doug.values.flatten()


# evaluate self-consistency using different covariates, ranks, lambdas
Corr_1 = defaultdict(list)
Corr_both = defaultdict(list)
for cov in adata.obs.columns:
	for _ in range(4):
		run_random_split('zombie_sample1',covariates=cov)
		res = check_consistency('zombie_sample1')
		Corr_1[cov].append(res[0])
		Corr_both[cov].append(res[1])

for cov in ['total_counts,circularity','total_counts,num_neighbors','area,circularity']:
	for _ in range(4):
		run_random_split('zombie_sample1',covariates=cov)
		res = check_consistency('zombie_sample1')
		Corr_1[cov].append(res[0])
		Corr_both[cov].append(res[1])


for rank in np.arange(8,36,2):
	for _ in range(4):
		run_random_split('zombie_sample1',rank=rank)
		res = check_consistency('zombie_sample1')
		Corr_1[rank].append(res[0])
		Corr_both[rank].append(res[1])

for rank in np.arange(8,36,2):
	for _ in range(4):
		run_random_split_withoutnorm('zombie_sample1',rank=rank)
		res = check_consistency('zombie_sample1')
		Corr_1[rank].append(res[0])
		Corr_both[rank].append(res[1])

for key in Corr_1.keys():
	print(key,np.average(Corr_1[key]),np.average(Corr_both[key]))


# conclusions:
# 	not a huge difference among covariates, but total_counts or area seem fine
# 	lambda1, lambda2 = 0 seems best
# 	rank 32 seems best


###
### compare with perturb-seq
###

dir = 'zombie_sample1'
covariates,num_permute,rank = 'total_counts',10000,34
_=os.system('python ./run_FR_PerturbnoNorm.py --input-h5ad %s/count_data.h5ad --input-perturbation-matrix %s/perturbation_design.csv --control-perturbation-name Control --out %s/all_data --covariates %s --compute-pval --num-perms %d --rank %d --lambda1 0.0 --lambda2 0.0' % (dir,dir,dir,covariates,num_permute,rank))
res = print_correlations('all_data',dir='zombie_sample2',return_values=True)
row_order = hierarchy.leaves_list(hierarchy.linkage(LFC_doug.values,optimal_ordering=True))
col_order = hierarchy.leaves_list(hierarchy.linkage(LFC_doug.values.T,optimal_ordering=True))
_=sns.heatmap(LFC_doug.iloc[row_order,col_order], cmap='vlag',center=0,vmin=-0.4,vmax=0.4)
_=sns.heatmap(res[0].iloc[row_order,col_order], cmap='vlag',center=0,vmin=-0.4,vmax=0.4)


LFC_melt = LFC_doug.stack().reset_index()
LFC_melt.columns = ['Gene','Target','Perturb-Seq effect']
LFC_melt['Perturb-FISH effect'] = res[0].stack().reset_index()[0]
idx_q = (res[3] < 0.1)*(qd < 0.1)
_=sns.scatterplot(data=LFC_melt[idx_q],x='Perturb-Seq effect',y='Perturb-FISH effect')
plt.savefig('%s/LFC_scatter_both_significant.pdf' % dir)
plt.close()

Z = []
for target in LFC_doug.columns:
	z = LFC_melt[idx_q][LFC_melt[idx_q].Target == target]
	if len(z) > 3:
		c = 1-distance.correlation(z['Perturb-Seq effect'], z['Perturb-FISH effect'])
		print(target,z.shape[0],c)
		Z.append(z)


Z = pd.concat(Z, ignore_index=True)
g = sns.FacetGrid(Z, col="Target",  col_wrap=3)
_=g.map(sns.scatterplot, "Perturb-Seq effect", "Perturb-FISH effect")
plt.savefig('%s/LFC_scatter_both_significant_individual_genes.pdf' % dir)
plt.close()


###
### density effects
###
n = N.values
n[n>7] = 7
Z1 = M[['TNF']] + 1
Z1['Neighbors'] = N
Z1.columns = ['Expression','Neighbors']
Z1['Gene'] = 'TNF'
Z1=Z1[Z1.Expression<20]
Z2 = M[['IL1A']] + 1
Z2['Neighbors'] = N
Z2.columns = ['Expression','Neighbors']
Z2['Gene'] = 'IL1A'
Z2=Z2[Z2.Expression<20]
Z = pd.concat([Z1,Z2], ignore_index=True)
g = sns.FacetGrid(Z, col='Gene')
_=g.map(sns.violinplot, "Neighbors", "Expression", log_scale=True)
plt.savefig('%s/density_violins_TNF_IL1A.eps' % dir, format='eps')
plt.close()



idx1 = np.where(N[0].values < 3)[0]
idx2 = np.where(N[0].values >= 3)[0]
dir = 'zombie_sample1'
covariates,num_permute,rank = 'area',10000,34
write_data_subset(idx1,'%s/split1_' % dir)
write_data_subset(idx2,'%s/split2_' % dir)
_=os.system('python ~/Documents/Code/FR-Perturb/run_FR_Perturb.py --input-h5ad %s/split1_count_data.h5ad --input-perturbation-matrix %s/split1_perturbation_design.csv --control-perturbation-name Control --out %s/split1 --covariates %s --compute-pval --num-perms %d --rank %d --lambda1 0.0 --lambda2 0.0' % (dir,dir,dir,covariates,num_permute,rank))
_=os.system('python ~/Documents/Code/FR-Perturb/run_FR_Perturb.py --input-h5ad %s/split2_count_data.h5ad --input-perturbation-matrix %s/split2_perturbation_design.csv --control-perturbation-name Control --out %s/split2 --covariates %s --compute-pval --num-perms %d --rank %d --lambda1 0.0 --lambda2 0.0' % (dir,dir,dir,covariates,num_permute,rank))
res = check_consistency(dir)
LFC_low, Q_low = res[2:4]
LFC_high,Q_high = res[4:6]
Q_sum = (Q_low < 0.2) + (Q_high < 0.2)
idx_target = np.where(Q_sum.values.sum(0) > 0)[0]
idx_exp = np.where(Q_sum.values.sum(1) > 0)[0]
LFC_diff = LFC_low.iloc[idx_exp,idx_target] - LFC_high.iloc[idx_exp,idx_target]

LFC_diff.stack().index[np.argsort(abs(LFC_diff.values.flatten()))[-20:]]

density_violins('EBI3','MAP3K7',savepath='%s/density_violins_EBI3_MAP3K7.pdf' % dir)
density_violins('TNF','JUN',savepath='%s/density_violins_TNF_JUN.pdf' % dir)

_=sns.clustermap(LFC_diff, cmap="vlag", center=0)





#
# sample 2
#
M = pd.read_csv('zombie_sample2/finalmerfish.csv') # cells x genes
P = pd.read_csv('zombie_sample2/finalzombie.csv') # needs to be target x cells
C = pd.read_csv('zombie_sample2/covariates.csv') # cells x covariates
N = pd.read_csv('zombie_sample2/filterednumberofneighbor.csv', header=None)

# 'HN1' is not in doug data
M = M.drop(columns=['HN1'])

P = P.rename(columns={'control': 'Control'}).T
P = P.rename(columns= dict([(i,'Cell_%d' % (i+1)) for i in P.columns.values]))

adata = ad.AnnData(M)
adata.obs_names = [f"Cell_{i+1:d}" for i in range(adata.n_obs)]
adata.obs['area'] = C['area'].values
adata.obs['eccentricity'] = C['eccentricity'].values
adata.obs['circularity'] = C['circularity'].values
adata.obs['num_neighbors'] = N.values[:,0]
adata.obs['total_counts'] = adata.X.sum(1)
adata.write('zombie_sample2/count_data.h5ad')

P.to_csv('zombie_sample2/perturbation_design.csv',sep=' ')

# counts
print_correlations('counts',dir='zombie_sample2')

# (0.1, 1) 355 0.14940818438317516
# (1, 0.1) 720 0.2285378178095273
# (0.1, 0.1) 84 0.332672986605987


Corr_1 = []
Corr_both = []
for _ in range(10):
	run_random_split('zombie_sample2')
	res = check_consistency('zombie_sample2')
	Corr_1.append(res[0])
	Corr_both.append(res[1])


#sample2 me

dir = 'zombie_sample2'
covariates,num_permute,rank = 'total_counts',10000,34
_=os.system('python ./run_FR_PerturbnoNorm.py --input-h5ad %s/count_data.h5ad --input-perturbation-matrix %s/perturbation_design.csv --control-perturbation-name Control --out %s/all_data --covariates %s --compute-pval --num-perms %d --rank %d --lambda1 0.0 --lambda2 0.0' % (dir,dir,dir,covariates,num_permute,rank))
res = print_correlations('all_data',dir='zombie_sample2',return_values=True)
row_order = hierarchy.leaves_list(hierarchy.linkage(LFC_doug.values,optimal_ordering=True))
col_order = hierarchy.leaves_list(hierarchy.linkage(LFC_doug.values.T,optimal_ordering=True))
_=sns.heatmap(LFC_doug.iloc[row_order,col_order], cmap='vlag',center=0,vmin=-0.4,vmax=0.4)
_=sns.heatmap(res[0].iloc[row_order,col_order], cmap='vlag',center=0,vmin=-0.4,vmax=0.4)


LFC_melt = LFC_doug.stack().reset_index()
LFC_melt.columns = ['Gene','Target','Perturb-Seq effect']
LFC_melt['Perturb-FISH effect'] = res[0].stack().reset_index()[0]
idx_q = (res[3] < 0.1)*(qd < 0.1)
_=sns.scatterplot(data=LFC_melt[idx_q],x='Perturb-Seq effect',y='Perturb-FISH effect')
plt.savefig('%s/LFC_scatter_both_significant.pdf' % dir)
plt.close()

Z = []
for target in LFC_doug.columns:
	z = LFC_melt[idx_q][LFC_melt[idx_q].Target == target]
	if len(z) > 3:
		c = 1-distance.correlation(z['Perturb-Seq effect'], z['Perturb-FISH effect'])
		print(target,z.shape[0],c)
		Z.append(z)


Z = pd.concat(Z, ignore_index=True)
g = sns.FacetGrid(Z, col="Target",  col_wrap=3)
_=g.map(sns.scatterplot, "Perturb-Seq effect", "Perturb-FISH effect")
plt.savefig('%s/LFC_scatter_both_significant_individual_genes.pdf' % dir)
plt.close()

#pooled data me

dir = 'zombie_bothsamples'
M = pd.read_csv('zombie_bothsamples/merfishpooleddata.csv') # cells x genes
P = pd.read_csv('zombie_bothsamples/intrinsicdesignpooleddata.csv', index_col=0) # target x cells
N = pd.read_csv('zombie_bothsamples/covnumberofneighborspooled.csv', header=None)

# 'HN1' is not in doug data
M = M.drop(columns=['HN1'])

adata = ad.AnnData(M)
adata.obs_names = [f"Cell_{i+1:d}" for i in range(adata.n_obs)]
adata.obs['area'] = C['area'].values
adata.obs['eccentricity'] = C['eccentricity'].values
adata.obs['circularity'] = C['circularity'].values
adata.obs['extent'] = C['extent'].values
adata.obs['num_neighbors'] = N.values[:,0]
adata.obs['total_counts'] = adata.X.sum(1)
adata.write('zombie_bothsamples/count_data.h5ad')

P.to_csv('zombie_bothsamples/perturbation_design.csv',sep=' ')


LFC_doug = pd.read_csv('perturbseq/neweffectsfromDoug.csv',index_col=0)
Q_doug = pd.read_csv('perturbseq/dougsQvalues.csv',index_col=0)
y = LFC_doug.values.flatten()
qd = Q_doug.values.flatten()


covariates,num_permute,rank = 'total_counts',10000,34
_=os.system('python ./run_FR_PerturbnoNorm.py --input-h5ad %s/count_data.h5ad --input-perturbation-matrix %s/perturbation_design.csv --control-perturbation-name Control --out %s/all_data --covariates %s --compute-pval --num-perms %d --rank %d --lambda1 0.0 --lambda2 0.0' % (dir,dir,dir,covariates,num_permute,rank))
res = print_correlations('all_data',dir='zombie_bothsamples',return_values=True)
row_order = hierarchy.leaves_list(hierarchy.linkage(LFC_doug.values,optimal_ordering=True))
col_order = hierarchy.leaves_list(hierarchy.linkage(LFC_doug.values.T,optimal_ordering=True))
_=sns.heatmap(LFC_doug.iloc[row_order,col_order], cmap='vlag',center=0,vmin=-0.4,vmax=0.4)
_=sns.heatmap(res[0].iloc[row_order,col_order], cmap='vlag',center=0,vmin=-0.4,vmax=0.4)


LFC_melt = LFC_doug.stack().reset_index()
LFC_melt.columns = ['Gene','Target','Perturb-Seq effect']
LFC_melt['Perturb-FISH effect'] = res[0].stack().reset_index()[0]
idx_q = (res[3] < 0.1)*(qd < 0.1)
_=sns.scatterplot(data=LFC_melt[idx_q],x='Perturb-Seq effect',y='Perturb-FISH effect')
plt.savefig('%s/LFC_scatter_both_significant2.pdf' % dir)
plt.close()

Z = []
for target in LFC_doug.columns:
	z = LFC_melt[idx_q][LFC_melt[idx_q].Target == target]
	if len(z) > 3:
		c = 1-distance.correlation(z['Perturb-Seq effect'], z['Perturb-FISH effect'])
		print(target,z.shape[0],c)
		Z.append(z)


Z = pd.concat(Z, ignore_index=True)
g = sns.FacetGrid(Z, col="Target",  col_wrap=3)
_=g.map(sns.scatterplot, "Perturb-Seq effect", "Perturb-FISH effect")
plt.savefig('%s/LFC_scatter_both_significant_individual_genes2.pdf' % dir)
plt.close()


#trying splits of data
Corr_1 = defaultdict(list)
Corr_both = defaultdict(list)
for cov in adata.obs.columns:
	for _ in range(4):
		run_random_split_withoutnorm('zombie_bothsamples',covariates=cov)
		res = check_consistency('zombie_bothsamples')
		Corr_1[cov].append(res[0])
		Corr_both[cov].append(res[1])

for cov in ['total_counts,circularity','total_counts,num_neighbors','area,circularity']:
	for _ in range(4):
		run_random_split_withoutnorm('zombie_bothsamples',covariates=cov)
		res = check_consistency('zombie_bothsamples')
		Corr_1[cov].append(res[0])
		Corr_both[cov].append(res[1])


for rank in np.arange(32,36,2):
	for _ in range(4):
		run_random_split_withoutnorm('zombie_bothsamples',rank=rank)
		res = check_consistency('zombie_bothsamples')
		Corr_1[rank].append(res[0])
		Corr_both[rank].append(res[1])

for key in Corr_1.keys():
	print(key,np.average(Corr_1[key]),np.average(Corr_both[key]))



#density on pooled

M = pd.read_csv('zombie_bothsamples/merfishpooleddata.csv') # cells x genes
P = pd.read_csv('zombie_bothsamples/intrinsicdesignpooleddata.csv', index_col=0) # target x cells
C = pd.read_csv('zombie_bothsamples/covariates.csv') # cells x covariates
N = pd.read_csv('zombie_bothsamples/covnumberofneighborspooled.csv', header=None)

# 'HN1' is not in doug data
M = M.drop(columns=['HN1'])

adata = ad.AnnData(M)
adata.obs_names = [f"Cell_{i+1:d}" for i in range(adata.n_obs)]
adata.obs['area'] = C['area'].values
adata.obs['eccentricity'] = C['eccentricity'].values
adata.obs['circularity'] = C['circularity'].values
adata.obs['extent'] = C['extent'].values
adata.obs['num_neighbors'] = N.values[:,0]
adata.obs['total_counts'] = adata.X.sum(1)
adata.write('zombie_bothsamples/count_data.h5ad')

P.to_csv('zombie_bothsamples/perturbation_design.csv',sep=' ')


LFC_doug = pd.read_csv('perturbseq/neweffectsfromDoug.csv',index_col=0)
Q_doug = pd.read_csv('perturbseq/dougsQvalues.csv',index_col=0)
y = LFC_doug.values.flatten()
qd = Q_doug.values.flatten()



n = N.values
n[n>7] = 7
Z1 = M[['TNF']] + 1
Z1['Neighbors'] = N
Z1.columns = ['Expression','Neighbors']
Z1['Gene'] = 'TNF'
Z1=Z1[Z1.Expression<20]
Z2 = M[['IL1A']] + 1
Z2['Neighbors'] = N
Z2.columns = ['Expression','Neighbors']
Z2['Gene'] = 'IL1A'
Z2=Z2[Z2.Expression<20]
Z = pd.concat([Z1,Z2], ignore_index=True)
g = sns.FacetGrid(Z, col='Gene')
_=g.map(sns.violinplot, "Neighbors", "Expression", log_scale=True)
plt.savefig('%s/density_violins_TNF_IL1A.pdf' % dir)
plt.close()



idx1 = np.where(N[0].values < 3)[0]
idx2 = np.where(N[0].values >= 3)[0]
dir = 'zombie_bothsamples'
covariates,num_permute,rank = 'total_counts',10000,34
write_data_subset_nocovs(idx1,'%s/split1_' % dir)
write_data_subset_nocovs(idx2,'%s/split2_' % dir)
_=os.system('python ./run_FR_PerturbnoNorm.py --input-h5ad %s/split1_count_data.h5ad --input-perturbation-matrix %s/split1_perturbation_design.csv --control-perturbation-name Control --out %s/split1 --covariates %s --compute-pval --num-perms %d --rank %d --lambda1 0.0 --lambda2 0.0' % (dir,dir,dir,covariates,num_permute,rank))
_=os.system('python ./run_FR_PerturbnoNorm.py --input-h5ad %s/split2_count_data.h5ad --input-perturbation-matrix %s/split2_perturbation_design.csv --control-perturbation-name Control --out %s/split2 --covariates %s --compute-pval --num-perms %d --rank %d --lambda1 0.0 --lambda2 0.0' % (dir,dir,dir,covariates,num_permute,rank))
res = check_consistency(dir)
LFC_low, Q_low = res[2:4]
LFC_high,Q_high = res[4:6]
Q_sum = (Q_low < 0.2) + (Q_high < 0.2)
idx_target = np.where(Q_sum.values.sum(0) > 0)[0]
idx_exp = np.where(Q_sum.values.sum(1) > 0)[0]
LFC_diff = LFC_low.iloc[idx_exp,idx_target] - LFC_high.iloc[idx_exp,idx_target]

LFC_diff.stack().index[np.argsort(abs(LFC_diff.values.flatten()))[-20:]]

density_violins('EBI3','MAP3K7',savepath='%s/density_violins_EBI3_MAP3K7.pdf' % dir)
density_violins('TNF','JUN',savepath='%s/density_violins_TNF_JUN.pdf' % dir)
density_violins('TNF','NFKB1',savepath='%s/density_violins_TNF_NFKB1.pdf' % dir)
density_violins_eps('IL1A','TRAF6',savepath='%s/density_violins_IL1A_TRAF6.eps' % dir)

_=sns.clustermap(LFC_diff, cmap="vlag", center=0)
plt.savefig('%s/density_violins_LFC_diff.pdf' % dir)
plt.close()
