import numpy as np
import pandas as pd
import math
import matplotlib
import matplotlib.pyplot as plt
import palantir as pl
from matplotlib.backends.backend_pdf import PdfPages
import pyreadr as pr
import re
import os
matplotlib.rcParams['font.family'] = ['serif']

os.chdir('/Users/parkt1/0-workspace/CCR7_DC/TC_all_LN')

def plot_gene_expr(expr, vis, dim1, dim2, genes, file, n_cols, s=3, 
                   newGeneName=None):
    cmap = matplotlib.cm.Spectral_r
    fig = pl.plot.FigureGrid(len(genes), n_cols)
    for g, ax in zip(genes, fig):
        c = expr.loc[vis.index, g.upper()]
        ax.scatter(vis.loc[:, dim1], vis.loc[:, dim2], s=s, c=c, cmap=cmap)
        ax.set_axis_off()
        if newGeneName and g in newGeneName:
            figTitle = newGeneName[g]
        else:
            figTitle = g
        ax.set_title(figTitle)
        normalize = matplotlib.colors.Normalize(vmin=np.min(c), vmax=np.max(c))
        cax, _ = matplotlib.colorbar.make_axes(ax)
        cb = matplotlib.colorbar.ColorbarBase(cax, norm=normalize, cmap=cmap)
        #cb.ax.set_title('Expression')
        cb.set_label('Expression')
    n_rows = math.ceil(len(genes) / n_cols)
    fig.figure.set_size_inches(7 * n_cols, 5 * n_rows)
    file.savefig(fig.figure, bbox_inches='tight')
    return 0

md = pd.read_csv('results/meta-data.csv', header=0, index_col=0)
# umap = pd.read_csv('results/UMAP.csv', header=0, index_col=0)
umap_s = pd.read_csv('results/UMAP-subset.csv', header=0, index_col=0)
imp_df = pr.read_r('results/imputed-expr.rds')[None]
imp_df = imp_df.transpose()
unimp_df = pr.read_r('results/unimputed-expr.rds')[None]
unimp_df = unimp_df.transpose()

genes = imp_df.columns
os.mkdir("plots/gene-expr")
############################################################################
geneList = ['ALDH1A2']

for gene in geneList:
    if gene not in genes:
        print(gene)

title = 'RALDH2'
ncols = 1
ci = umap_s.index
cellgroup = ''
    
file = PdfPages('plots/gene-expr/UMAP-{}gene-expr-imputed-{}.pdf'.format(cellgroup, title))
plot_gene_expr(expr=imp_df.loc[ci, :], vis=umap_s.loc[ci, :], dim1='UMAP_1', dim2='UMAP_2', genes=geneList,
                file=file, n_cols=ncols, s=8)
file.close()
file = PdfPages('plots/gene-expr/UMAP-{}gene-expr-unimputed-{}.pdf'.format(cellgroup, title))
plot_gene_expr(expr=unimp_df.loc[ci, :], vis=umap_s.loc[ci, :], dim1='UMAP_1', dim2='UMAP_2', genes=geneList,
                file=file, n_cols=ncols, s=8)
file.close()
