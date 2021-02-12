# This data was downloaded from https://zenodo.org/record/3711134#.YBNCrndKhTY
# As part of the paper “A cell atlas of human thymic development defines T
# cell repertoire formation”
# The best way to access cell types was found in the jupyter notebook located:
# /Users/wellskr/Documents/Analysis/Holger_Russ/mtec_organoid_multi/files/thymus_annotated_matrix_files/Data_navigator.ipynb
# opened by navigating to the folder and running jupyter notebook
# Run this using the scanpy envrionment

import scanpy as sc
import scipy.sparse
import pandas
import anndata2ri
from rpy2.robjects import r
from rpy2.robjects.conversion import localconverter
sc.settings.verbosity = 3

data_dir = "/Users/wellskr/Documents/Analysis/Holger_Russ/mtec_organoid_multi/files/thymus_annotated_matrix_files/"

fig1_data = data_dir + 'HTA08.v01.A05.Science_human_fig1.h5ad'
hsc_data = data_dir + 'HTA08.v01.A05.Science_human_hsc.h5ad'
tcell_data = data_dir + 'HTA08.v01.A06.Science_human_tcells.h5ad'
human_epi_data = data_dir + 'HTA08.v01.A05.Science_human_epi.h5ad'
mouse_stromal_data = data_dir + 'HTA08.v02.A04.Science_mouse_stromal.h5ad'
mouse_total_data = data_dir + 'HTA08.v02.A04.Science_mouse_total.h5ad'

meta_data_file = data_dir + "all_cell_metadata.csv"
counts_file = data_dir + "all_cell_counts.csv"
genes_file = data_dir + "all_cell_genes.csv"
cells_file = data_dir + "all_cell_cells.csv"

fig1 = sc.read(fig1_data)

# They plotted using "Anno_level_fig1" There are many anno levels in their obs, but I will
# grab this one.
# They are all different levels, so I will actually use them all.
#sc.pl.umap(fig1,color = 'Anno_level_fig1'.split(','))
meta_df = pandas.DataFrame(fig1.obs)
meta_df.to_csv(meta_data_file)

count_df = pandas.DataFrame.sparse.from_spmatrix(fig1.X)
count_df.to_csv(counts_file)
fig1.var.to_csv(genes_file)
fig1.obs.to_csv(cells_file)


