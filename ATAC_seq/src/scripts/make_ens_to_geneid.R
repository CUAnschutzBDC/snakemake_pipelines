library(here)
library(biomaRt)

mart <- useEnsembl(biomart = 'genes', 
                   dataset = 'mmusculus_gene_ensembl',
                   version = 102)

G_list <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol"),
                mart = mart)

saveRDS(G_list, here("files/ens_to_geneid.rds"))