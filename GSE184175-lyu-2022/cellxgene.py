import pandas as pd
import os
import anndata as ad

os.chdir('GSE184175-lyu-2022')

meta_data = pd.read_csv('results/meta-data.csv', header=0, index_col=0)
meta_data.rename(columns = {'percent.MT': 'MtFrac_RNA'}, inplace = True)
meta_data['Clusters'] = pd.Categorical(meta_data['Clusters'])
header = ['nCount_RNA', 'nFeature_RNA', 'MtFrac_RNA', 'S.Score', 'G2M.Score', 
          'Phase', 'Clusters']
meta_data = meta_data.loc[:, header]
# meta_data = meta_data.iloc[:, [1, 2, 7, 9, 10, 11, 6, 8]]

UMAP = pd.read_csv('results/UMAP.csv', header=0, index_col=0)

# unimputed
expr = pd.read_csv('results/unimputed-expr.csv', header=0, index_col=0).transpose()
expr.index = [ ind.replace(".", "-") for ind in expr.index ]
ad_obj = ad.AnnData(X = expr, obs = meta_data, obsm = {'X_umap' : UMAP.values})
ad_obj.write_h5ad('results/unimputed.h5ad')

# imputed
expr = pd.read_csv('results/imputed-expr.csv', header=0, index_col=0).transpose()
expr.index = [ ind.lstrip("X") for ind in expr.index ]
ad_obj = ad.AnnData(X = expr, obs = meta_data, obsm = {'X_umap' : UMAP.values})
ad_obj.write_h5ad('results/imputed.h5ad')

# md = meta_data[:3000]
# umap = UMAP.loc[md.index]
# expr2 = expr.loc[md.index]
# ad_obj = ad.AnnData(X = expr2, obs = md, obsm = {'X_umap' : umap.values})
# ad_obj.write_h5ad('results/imputed-small.h5ad')