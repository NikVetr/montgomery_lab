library(biomaRt)
library(KEGGREST)
library(org.Hs.eg.db)
library(limma)

# kegg_info <- keggGet(c("hsa:2947"))
# str(kegg_info)
# kegg_info[[1]]$PATHWAY
# p = keggLink("hsa", "pathway")
# map = split(names(p), unname(p))  # gene -> pathway map
# map[["hsa:574537"]]

load("~/data/smontgom/node_metadata_list.RData")
ensembl_genes <- node_metadata_list$`8w`$human_ensembl_gene

tab <- getGeneKEGGLinks(species="hsa")
tab$ENSEMBL <- mapIds(org.Hs.eg.db, tab$GeneID, column="ENSEMBL", keytype="ENTREZID")
tab$SYMBOL <- mapIds(org.Hs.eg.db, tab$GeneID, column="SYMBOL", keytype="ENTREZID")

# length(unique(tab$ENSEMBL))
# length(map)

ensembl_to_pathway <- lapply(setNames(ensembl_genes, ensembl_genes), function(eg) tab$PathwayID[which(eg == tab$ENSEMBL)])

km <- c(read.csv("~/data/kegg_pathways.csv", header = F))[[1]]
l1is <- c(grep("\\. ", km), length(km) + 1)
l1l <- lapply(2:length(l1is), function(l1i) {
  l1n <- km[l1is[l1i-1]]
  l1p <- km[(l1is[l1i-1]+1):(l1is[l1i]-1)]
  setNames(list(l1p), l1n)
  })

kegg_map <- do.call(rbind, lapply(1:length(l1l), function(l1li) {
  l2dat <- l1l[[l1li]][[1]]
  l2is <- c(grep("\\.", l2dat), length(l2dat) + 1)
  l2lsub <- lapply(2:length(l2is), function(l2i) {
    l2n <- l2dat[l2is[l2i-1]]
    l2p <- l2dat[(l2is[l2i-1]+1):(l2is[l2i]-1)]
    l2p <- as.data.frame(matrix(l2p, ncol = 2, byrow = T))
    colnames(l2p) <- c("kegg_ID", "pathway")
    l2p$subgroup <- l2n
    l2p$group <- names(l1l[[l1li]])
    l2p
  })
  l2lsub <- do.call(rbind, l2lsub)
}))
str(kegg_map)

kegg_map$kegg_ID <- as.numeric(sapply(kegg_map$kegg_ID, function(kid) strsplit(kid, "\\s+")[[1]][1]))
kegg_map$subgroup <- sapply(kegg_map$subgroup, function(sg) paste0(strsplit(sg, "\\s+")[[1]][-1], collapse = " "))
kegg_map$group <- sapply(kegg_map$group, function(gn) paste0(strsplit(gn, "\\s+")[[1]][-1], collapse = " "))


genepath_mat <- as.data.frame(matrix(0, nrow = length(ensembl_to_pathway), ncol = nrow(kegg_map), dimnames = list(names(ensembl_to_pathway), kegg_map$kegg_ID)))
for(gi in 1:length(ensembl_to_pathway)){
  if(gi %% 100 == 0) cat(paste0(round(gi / length(ensembl_to_pathway) * 100, 2), " "))
  gene <- names(ensembl_to_pathway)[gi]
  paths <- ensembl_to_pathway[gene][[1]]
  paths <- as.numeric(do.call(rbind, strsplit(paths, "hsa"))[,2])
  paths <- unique(paths)
  genepath_mat[gene, as.character(paths)] <- 1
}

table(unlist(genepath_mat))
table(apply(genepath_mat, 2, sum))
length(unique(paste0(names(unlist(ensembl_to_pathway)), "-", unlist(ensembl_to_pathway))))
