#!/usr/bin/env Rscript
library(dplyr)
library(httr)
library(MAST)
library(unixtools)
library(ggplot2)
library(cowplot)
library(future)
#library(ggmin)
library(pheatmap)
library(ggrepel)
library(shiny)
library(shinythemes)
library(plotly)
library(magrittr)
library(rlang)
library(tidyr)
library(tibble)
library(tidyverse)
library(grid)
library(data.table)
library(shinycssloaders)
library(DT)
library(shinydashboard)
library(shinyjs)

library(Matrix)
library(slingshot)
library(tradeSeq)
library(SingleCellExperiment)
library(RColorBrewer)
library(gam)
library(scales)
library(gridExtra)
library(destiny)
library(ComplexHeatmap)
library(phateR)
library(monocle)
library(monocle3)
# library(SeuratWrappers)
library(patchwork)
library(circlize)
library(cowplot)
theme_set(theme_cowplot())
library(Seurat)




# adapt to your path
#setwd("/home/jr345y/tcell_cd18KO_compact/tcell_cd18KO_compact/")

##Reading in the list of precomputed Seurat objects (Resolution 0.15 0.25&0.35, 0.45&0.55 computed separately read in and combined then written to disk and )
tcells_combined_umap_list_res_skinny<-readRDS("tcells_combined_umap_list_res_skinny.rds")

##Reading in the list of precomputed table of cluster markers
tcells_combined_clusters_tables_res<-readRDS("tcells_combined_clusters_tables_res.rds")

##Reading in the list of precomputed table of differential expressed (DE) genes
tcells_combined_de_tables = readRDS("tcells_combined_de_tables.rds")

##Reading in the list of precomputed table of cluster markers
tcells_combined_de_ggplots_table = readRDS("tcells_combined_de_ggplots_table.rds")

##Reading in all the genes that are present in WT1, WT2 and KO raw objects
all_genes_common_in_all_groups = readRDS("all_genes_common_in_all_groups.rds")

##Reading in and processing table of links to PlasmoDB.. 
uniprot_info = fread("Pberghei.Anntation.txt", stringsAsFactors = F)

##Reading in Pseudotime
#UMAP
sds = readRDS("sds.rds")
sce = readRDS("sce.rds")
counts = readRDS("counts.rds")
clusters = readRDS("clusters.rds")

# PHATE
sdsPhate = readRDS("sdsPhate.rds")
scePhate = readRDS("scePhate.rds")
countsPhate = readRDS("countsPhate.rds")
clustersPhate = readRDS("clustersPhate.rds")

cell_colors_clust = readRDS("cell_colors_clust.rds")

rerunPhate = readRDS("rerunPhate.rds")
slingUMAPHeat_ALL = readRDS("slingUMAPHeat_ALL.rds")
slingPHATEHeat_ALL = readRDS("slingPHATEHeat_ALL.rds")

tradeUMAPHeat_ALL = readRDS("tradeUMAPHeat_ALL.rds")
tradePHATEHeat_ALL = readRDS("tradePHATEHeat_ALL.rds")


# MONOCLE 2
cds2 = readRDS("cds2.rds")
new_monocle2_heatmap = readRDS("new_monocle2_heatmap.rds")

# MONOCLE3
cds3 = readRDS("cds3.rds")
new_monocle3_heatmap = readRDS("new_monocle2_heatmap.rds")





##Declaring and assigning variables
dim=24
res1 = 0.51
res2 = 0.51
diff_res = 0.2
#cluster_names = c("Early Rings", "Mid Rings", "Late Rings", "Early Trophs 1", 
#                  "Early Trophs 2", "Mid Trophs", "Late Trophs","Schizonts","Early Gams","Males","Females");
cluster_names = c("Early Rings", "Mid Rings", "Late Rings", "Early Trophs 1", 
                  "Early Trophs 2", "Mid Trophs", "Late Trophs","Schizonts","Early Gams","Males","Females");
#
#Early Trophs 2 Mid Rings      Mid Trophs     Early Rings    Early Trophs 1 Late Rings     Late Trophs    Early Gams    Females        Males          Schizonts 
fav_genes = c("PBANKA-1437500.1","PBANKA-1334300.1","PBANKA-0515000.1")
conditions = c("B", "BM","S")
umap_names = c("Early Rings", "Mid Rings", "Late Rings", "Early Trophs 1", 
               "Early Trophs 2", "Mid Trophs", "Late Trophs","Schizonts","Early Gams","Males","Females")

pairwise <- combn(conditions, 2)
conds = lapply(1:ncol(pairwise), function(x) paste(pairwise[,x], collapse = " VS ")) %>% unlist()
cluster.colours <-c("#A42537","dodgerblue2","#95E949","#FF1E1A")
group.cols <- c("red", "deepskyblue","green","blue")
names(group.cols) <- conditions
choice_gene = "PBANKA-1437500.1"
cond = "Mouse tissues"





if (!(all(exists("tcells_combined_umap_list_res_skinny"), exists("tcells_combined_clusters_tables_res"), exists("tcells_combined_de_tables")))) {
  
  
  
  
  Combined.filt<-readRDS("Pb.combined.noMCA.rds")
  DimPlot(Combined.filt)
  ## transfer the correct variables
  Combined.filt$sample<-Combined.filt$orig.ident; # orig.ident is the default name for a sample in Seurat
  
  ## in this case the group is the organ which is S for spleen BM for Bone marrow and B for blood
  Combined.filt$group<-Combined.filt$Organ
  Combined.filt$groups<-Combined.filt$Organ
  
  ## Potential bug, as some genes are not in all objects!
  all_genes_common_in_all_groups = rownames(Combined.filt@assays$RNA); #Reduce(intersect,list(all_genes_ko,all_genes_wt1,all_genes_wt3))
  saveRDS(all_genes_common_in_all_groups, "all_genes_common_in_all_groups.rds")
  
  Combined.filt <- FindNeighbors(Combined.filt)
  DefaultAssay(object = Combined.filt) <- "integrated"
  
  Combined.filt <-FindClusters(Combined.filt)
  ##Precomputing and saving the list of Seurat objects with different clusters through adjusting of resolution from 0.15, 0.25, 0.35, 0.45 & 0.55
  tcells_combined_umap_list_res <- lapply(seq(res1, res2, by = diff_res), function(x) FindClusters(Combined.filt, resolution = x))
  tcells_combined_umap_list_res[[2]] <- Combined.filt;
  # tcells_combined_umap_list_res <-readRDS("tcells_combined_umap_list_res.rds")
  remove(tcells_combined_umap_list_res[[2]])
  tcells_combined_umap_list_res[[2]]=NULL
  
  ##Precomputing and saving the conserved markers
  tcells_combined_clusters_tables_res = lapply(tcells_combined_umap_list_res, function(x) { 
    DefaultAssay(x) = "RNA"
    lapply(0:(length(unique(x$seurat_clusters))-1), function(y) {
      FindConservedMarkers(x, ident.1 = y, grouping.var = "group")
      
    })
  })
  
  saveRDS(tcells_combined_clusters_tables_res, "tcells_combined_clusters_tables_res.rds")
  
  
  ##Generating pairwise list for all DE group comparisons per cluster
  grps = unique(tcells_combined_umap_list_res[[1]]@meta.data$group)
  pairwise <- combn(grps, 2)
  
  ##Precomputing and saving the list of tables of DE genes per cluster
  tcells_combined_de_tables = lapply(tcells_combined_umap_list_res, function(x) { 
    DefaultAssay(x) = "RNA"
    x$celltype.group <- paste(Idents(x), x$group, sep = "_")
    x$celltype <- Idents(x)
    Idents(x) <- "celltype.group"
    
    lapply(0:(length(unique(x$seurat_clusters))-1), function(y) {
    # lapply(unique(x$Stages), function(y) {
      lapply(1:ncol(pairwise), function(z) {
        #     print(y)
        #      print(z)
        #      print (paste(y, pairwise[1,z], sep = "_"))
        tryCatch(FindMarkers(x, ident.1 = paste(y, pairwise[1,z], sep = "_"), ident.2 = paste(y, pairwise[2,z], sep = "_"), 
                             verbose = T, min.cells.group = 3, assay = "RNA"), error=function(e) NULL)
        
      })
    })
  })
  
  saveRDS(tcells_combined_de_tables, "tcells_combined_de_tables.rds")
  
  ##Precomputing and saving the list of ggplot per cluster for all resolutions
  tcells_combined_de_ggplots_table = lapply(tcells_combined_umap_list_res, function(x) { 
    DefaultAssay(x) = "RNA"
    
    #    lapply(0:(length(unique(x$seurat_clusters))-1), function(y) {
    # the seurat_cluster are just niumbers but we need the cluster names... 
    lapply(0:(length(unique(x$seurat_clusters))-1), function(y) {
    # lapply(unique(x$Stages), function(y) {
      cells_type <- subset(x, idents = y)
      #Idents(cells_type) <- "sample"
      Idents(cells_type) <- "group"
      avg.cells <- log1p(AverageExpression(cells_type, verbose = FALSE)$RNA)
      avg.cells$gene <- rownames(avg.cells)
      # avg.cells <- avg.cells %>% filter(!grepl("^mt-", gene)) %>% dplyr::left_join(x = ., y = uniprot_info, by = c("gene" = "Gene")) %>% dplyr::distinct(., gene, .keep_all = T)%>% select(gene, `Protein name`, PlasmoDB, S, BM, B)
      avg.cells <- avg.cells %>% filter(!grepl("^mt-", gene))
    })
  })
  
  saveRDS(tcells_combined_de_ggplots_table, "tcells_combined_de_ggplots_table.rds")
  
  saveRDS(tcells_combined_umap_list_res,"tcells_combined_umap_list_res.rds")
  
  
  ##PSEUDOTIME####
  # Sys.time()
  ##PHATE####
  # to add PHATE reduction 
  for (i in 1:length(tcells_combined_umap_list_res)) {
    phate_output <- as.matrix(phate(t(tcells_combined_umap_list_res[[i]]@assays$integrated@data)))
    colnames(x= phate_output) <- paste0("PHATE_", 1:ncol(x = phate_output))
    phate.reduction <- CreateDimReducObject(
      embeddings = phate_output,
      key = "PHATE_",
      assay = "integrated"
    )
    tcells_combined_umap_list_res[[i]]@reductions$phate = phate.reduction
  }
  
  ##Slingshot####
  genes_to_test = vector(mode = "list")
  forHeatPhateSce = vector(mode = "list")
  forHeatSce = vector(mode = "list")
  sdsPhate = vector(mode = "list")
  sds = vector(mode = "list")
  sce = vector(mode = "list")
  scePhate = vector(mode = "list")
  slingUMAPHeat_ALL = vector(mode = "list")
  slingPHATEHeat_ALL = vector(mode = "list")
  
  for (i in 1:length(tcells_combined_umap_list_res)) {
    sceLoop <- as.SingleCellExperiment(tcells_combined_umap_list_res[[i]])
    reducedDim(sceLoop) <- reducedDim(sceLoop)[, 1:10]
    sce[[i]] = slingshot(sceLoop, reducedDim = 'UMAP', clusterLabels = 'seurat_clusters')
    scePhate[[i]] = slingshot(sceLoop, reducedDim = 'PHATE', clusterLabels = 'seurat_clusters')
    sds[[i]] = SlingshotDataSet(sce[[i]])
    sdsPhate[[i]] = SlingshotDataSet(scePhate[[i]])
    forHeatSce[[i]] = sce[[i]] #for the heatmaps
    forHeatPhateSce[[i]] = scePhate[[i]] #for the heatmaps 
  }
  remove(sceLoop)
  
  Sys.time()
  #heatmap:
  #UMAP
  for (i in 1:length(tcells_combined_umap_list_res)) {
    genes_to_test <- VariableFeatures(tcells_combined_umap_list_res[[i]])[1:700]
    lineage_cells <- colnames(forHeatSce[[i]])[!is.na(genes_to_test)]
    cnts <- logcounts(forHeatSce[[i]])[genes_to_test, lineage_cells]
    
    ptime = colData(forHeatSce[[i]])[c(paste("slingPseudotime_", 1:(length(sds[[i]]@lineages)), sep = ""))]
    # ptime = ptime[!is.na(ptime)]
    
    gam.pval <- apply(cnts, 1, function(z){
      d <- data.frame(z = z, ptime = ptime)
      tmp <- suppressWarnings(gam(z ~ lo(ptime), data=d))
      p <- summary(tmp)[4][[1]][1, 5]
      p
    })
    res <- tibble(
      id = names(gam.pval),
      pvals = gam.pval,
      qval = p.adjust(gam.pval, method = "fdr")) %>% 
      arrange(qval)
    
    
    to_plot <- as.matrix(forHeatSce[[i]]@assays@data@listData$logcounts[res$id[1:40], lineage_cells])
    ptime_order <- colnames(to_plot)[order(ptime)]
    remove(ptime)
    annotations <- colData(forHeatSce[[i]])[lineage_cells, 
                                            c(c(paste("slingPseudotime_", 1:(length(sds[[i]]@lineages)), sep = "")),
                                              "seurat_clusters")] %>% as.data.frame()
    
    
    col = hue_pal()(length(unique(tcells_combined_umap_list_res[[i]]$seurat_clusters)))[unique(annotations[,"seurat_clusters"])]
    col = as.matrix(col)
    rownames(col) = unique(annotations[,"seurat_clusters"])
    
    
    ha <- HeatmapAnnotation(df = annotations, col = list(seurat_clusters = col[,1]))
    slingUMAPHeat_ALL[[i]]= Heatmap(to_plot,
                                    column_order = ptime_order,
                                    show_column_names = FALSE,
                                    show_row_names = TRUE,
                                    top_annotation = ha,
                                    row_names_gp = grid::gpar(fontsize = 8))
  }
  
  #PHATE
  for (i in 1:length(tcells_combined_umap_list_res)) {
    genes_to_test <- VariableFeatures(tcells_combined_umap_list_res[[i]])[1:700]
    lineage_cells <- colnames(forHeatPhateSce[[i]])[!is.na(genes_to_test)]
    cnts <- logcounts(forHeatPhateSce[[i]])[genes_to_test, lineage_cells]
    
    ptime = colData(forHeatPhateSce[[i]])[c(paste("slingPseudotime_", 1:(length(sdsPhate[[i]]@lineages)), sep = ""))]
    
    gam.pval <- apply(cnts, 1, function(z){
      d <- data.frame(z = z, ptime = ptime)
      tmp <- suppressWarnings(gam(z ~ lo(ptime), data=d))
      p <- summary(tmp)[4][[1]][1, 5]
      p
    })
    
    res <- tibble(
      id = names(gam.pval),
      pvals = gam.pval,
      qval = p.adjust(gam.pval, method = "fdr")) %>% 
      arrange(qval)
    
    
    to_plot <- as.matrix(forHeatPhateSce[[i]]@assays@data@listData$logcounts[res$id[1:40], lineage_cells])
    ptime_order <- colnames(to_plot)[order(ptime)]
    annotations <- colData(forHeatPhateSce[[i]])[lineage_cells, 
                                                 c(c(paste("slingPseudotime_", 1:(length(sdsPhate[[i]]@lineages)), sep = "")),
                                                   "seurat_clusters")] %>% as.data.frame()
    
    
    col = hue_pal()(length(unique(tcells_combined_umap_list_res[[i]]$seurat_clusters)))[unique(annotations[,"seurat_clusters"])]
    col = as.matrix(col)
    rownames(col) = unique(annotations[,"seurat_clusters"])
    
    
    ha <- HeatmapAnnotation(df = annotations, col = list(seurat_clusters = col[,1]))
    slingPHATEHeat_ALL[[i]]= Heatmap(to_plot,
                                     column_order = ptime_order,
                                     show_column_names = FALSE,
                                     show_row_names = TRUE,
                                     top_annotation = ha,
                                     row_names_gp = grid::gpar(fontsize = 8))
  }
  
  saveRDS(slingUMAPHeat_ALL, "slingUMAPHeat_ALL.rds")
  saveRDS(slingPHATEHeat_ALL, "slingPHATEHeat_ALL.rds")
  Sys.time()
  saveRDS(sds, "sds.rds")
  saveRDS(sdsPhate, "sdsPhate.rds")
  
  
  
  ##tradeSeq####
  
  tradeUMAPHeat_ALL = vector(mode = "list")
  tradePHATEHeat_ALL = vector(mode = "list")
  clusters = vector(mode = "list")
  clustersPhate = vector(mode = "list")
  counts = vector(mode = "list")
  
  for (i in 1:length(sds)) {
    print(i)
    
    clusters[[i]] <- forHeatSce[[i]]$seurat_clusters
    clustersPhate[[i]] <- forHeatPhateSce[[i]]$seurat_clusters
    counts[[i]] = as.matrix(tcells_combined_umap_list_res[[i]]@assays$RNA@counts)
    keep <- rowSums(counts[[i]] > 1) >= 120
    counts[[i]] = counts[[i]][keep,]
    
    # to choose the number of knots: 
    # set.seed(5)
    # icMat[[i]] <- evaluateK(counts = counts[[i]], sds = sds[[i]], k = 3:10, nGenes = 200, verbose = T, plot = T, ylim=c(3,9))
    
    set.seed(7)
    sce[[i]] <- fitGAM(counts = counts[[i]], sds = sds[[i]], nknots = 5, verbose = TRUE, sce = TRUE)
    scePhate[[i]] <- fitGAM(counts = counts[[i]], sds = sdsPhate[[i]], nknots = 5, verbose = TRUE, sce = TRUE)
  }
  
  # heatmap:
  #UMAP
  for (i in 1:length(sce)) {
    ATres <- associationTest(sce[[i]])
    ATres = na.omit(ATres)
    ATres= ATres[order(ATres$pvalue), ][1:40,]
    ATres = na.omit(ATres)
    
    heatdata <- as.matrix(assays(sce[[i]])$counts[rownames(ATres), ])
    heatdata = log1p(heatdata)
    heatdata = heatdata[apply(heatdata[,-1], 1, function(x) !all(x==0)),] #removes values of 0
    
    annotations = colData(forHeatSce[[i]])[colnames(heatdata), 
                                           c(c(paste("slingPseudotime_", 1:(length(sds[[i]]@lineages)), sep = "")),
                                             "seurat_clusters")] %>% as.data.frame()
    
    col = hue_pal()(length(unique(tcells_combined_umap_list_res[[i]]$seurat_clusters)))[unique(annotations[,"seurat_clusters"])]
    col = as.matrix(col)
    rownames(col) = unique(annotations[,"seurat_clusters"])
    ha <- HeatmapAnnotation(df = annotations, col = list(seurat_clusters = col[,1]))
    
    
    tradeUMAPHeat_ALL[[i]] = Heatmap(heatdata,
                                     show_column_names = FALSE,
                                     show_row_names = TRUE,
                                     top_annotation = ha,
                                     show_column_dend = FALSE,
                                     row_names_gp = grid::gpar(fontsize = 8))
  }
  
  for (i in 1:length(scePhate)) {
    ATres <- associationTest(scePhate[[i]])
    ATres = na.omit(ATres)
    ATres= ATres[order(ATres$pvalue), ][1:40,]
    ATres = na.omit(ATres)
    
    heatdata <- as.matrix(assays(scePhate[[i]])$counts[rownames(ATres), ])
    heatdata = log1p(heatdata)
    heatdata = heatdata[apply(heatdata[,-1], 1, function(x) !all(x==0)),] #removes values of 0
    
    annotations <- colData(forHeatPhateSce[[i]])[colnames(heatdata), 
                                                 c(c(paste("slingPseudotime_", 1:(length(sdsPhate[[i]]@lineages)), sep = "")),
                                                   "seurat_clusters")] %>% as.data.frame()
    
    col = hue_pal()(length(unique(tcells_combined_umap_list_res[[i]]$seurat_clusters)))[unique(annotations[,"seurat_clusters"])]
    col = as.matrix(col)
    rownames(col) = unique(annotations[,"seurat_clusters"])
    ha <- HeatmapAnnotation(df = annotations, col = list(seurat_clusters = col[,1]))
    
    
    tradePHATEHeat_ALL[[i]] = Heatmap(heatdata,
                                      show_column_names = FALSE,
                                      show_row_names = TRUE,
                                      top_annotation = ha,
                                      show_column_dend = FALSE,
                                      row_names_gp = grid::gpar(fontsize = 8))
  }
  
  
  saveRDS(sce, "sce.rds")
  saveRDS(scePhate, "scePhate.rds")
  saveRDS(tradeUMAPHeat_ALL, "tradeUMAPHeat_ALL.rds")
  saveRDS(tradePHATEHeat_ALL, "tradePHATEHeat_ALL.rds")
  saveRDS(clustersPhate, "clustersPhate.rds")
  saveRDS(clusters, "clusters.rds")
  saveRDS(counts, "counts.rds")
  
  
  #function to add colours to clusters 
  cell_pal <- function(cell_vars, pal_fun,...) {
    if (is.numeric(cell_vars)) {
      pal <- pal_fun(100, ...)
      return(pal[cut(cell_vars, breaks = 100)])
    } else {
      categories <- sort(unique(cell_vars))
      pal <- setNames(pal_fun(length(categories), ...), categories)
      return(pal[cell_vars])
    }
  }
  
  
  #MONOCLE 3 
  cds3  = vector(mode = "list")
  set.seed(1234)
  for (i in 1:length(tcells_combined_umap_list_res)) {
    cds3[[i]] = as.cell_data_set(tcells_combined_umap_list_res[[i]])
    cds3[[i]] = cluster_cells(cds3[[i]])
    
    integrated.sub <- subset(as.Seurat(cds3[[i]]), monocle3_partitions == 1)
    cds3[[i]] <- as.cell_data_set(integrated.sub)
    cds3[[i]] <- learn_graph(cds3[[i]])
    
    #to find the root clusters
    max.clus <- which.max(unlist(FetchData(integrated.sub, "seurat_clusters")))
    max.clus <- colnames(integrated.sub)[max.clus]
    cds3[[i]] <- order_cells(cds3[[i]], root_cells = max.clus)
  }
  saveRDS(cds3, "cds3.rds")
  
  #heatmap
  new_monocle3_heatmap = vector(mode = "list")
  
  for(i in 1:length(cds3)){
    modulated_genes <- graph_test(cds3[[i]], neighbor_graph = "principal_graph", cores = 4)
    genes <- row.names(modulated_genes)
    genes = genes[1:40]
    
    pt.matrix <- exprs(cds3[[i]])[match(genes,rownames(rowData(cds3[[i]]))),order(modulated_genes$q_value)]
    pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
    pt.matrix = pt.matrix[1:40,]
    
    annotations <- colData(cds3[[i]])[colnames(pt.matrix),
                                      "seurat_clusters"] %>% as.data.frame()
    
    
    annotations["pseudotime"] = as.data.frame(pseudotime(cds3[[i]])[colnames(pt.matrix)])
    names(annotations)[1] = "seurat_clusters"
    
    col = hue_pal()(length(unique(tcells_combined_umap_list_res[[i]]$seurat_clusters)))[unique(annotations[,"seurat_clusters"])]
    col = as.matrix(col)
    rownames(col) = unique(annotations[,"seurat_clusters"])
    ha <- HeatmapAnnotation(df = annotations, col = list(seurat_clusters = col[,1]))
    
    new_monocle3_heatmap[[i]] = Heatmap(pt.matrix,
                                        show_column_names = FALSE,
                                        show_row_names = TRUE,
                                        top_annotation = ha,
                                        show_column_dend = FALSE,
                                        row_names_gp = grid::gpar(fontsize = 8))
  }

  saveRDS(new_monocle3_heatmap, "new_monocle3_heatmap.rds")
  
  
  #MONOCLE 2
  
  counts.data = vector(mode = "list")
  pheno.data = vector(mode = "list")
  feature.data = vector(mode = "list")
  cds2 = vector(mode = "list")
  
  for(i in 1:length(tcells_combined_umap_list_res)){
    counts.data <- as(as.matrix(tcells_combined_umap_list_res[[i]]@assays$RNA@data), 'sparseMatrix')
    pheno.data <- new('AnnotatedDataFrame', data = tcells_combined_umap_list_res[[i]]@meta.data)
    feature.data <- data.frame(gene_short_name = row.names(counts.data), row.names = row.names(counts.data))
    feature.data <- new('AnnotatedDataFrame', data = feature.data)
    #Construct a CellDataSet
    cds2[[i]] <- newCellDataSet(counts.data, phenoData = pheno.data, featureData = feature.data, lowerDetectionLimit = 0.5, expressionFamily = negbinomial.size())
    
    cds2[[i]] <- estimateSizeFactors(cds2[[i]])
    
    cds2[[i]] = estimateDispersions(cds2[[i]])
    disp_table <- dispersionTable(cds2[[i]])
    unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
    cds2[[i]] <- setOrderingFilter(cds2[[i]] , unsup_clustering_genes$gene_id)
    plot_ordering_genes(cds2[[i]])

    cds2[[i]] <- reduceDimension(cds2[[i]], reduction_method = "DDRTree", pseudo_expr = 1)

    cds2[[i]] <- orderCells(cds2[[i]], reverse = FALSE)
  }
  
  saveRDS(cds2, "cds2.rds")
  
  
  monocle2_heatmap = vector(mode = "list")
  new_monocle2_heatmap = vector(mode = "list")
  
  
  for(i in 1:length(cds2)){
    my_pseudotime_de <- differentialGeneTest(cds2[[i]],
                                             fullModelFormulaStr = "~sm.ns(Pseudotime)",
                                             cores = 8)
    
    my_pseudotime_de = my_pseudotime_de %>% arrange(qval) %>% head(40)
    my_pseudotime_de %>% arrange(qval) %>% head(40) %>% select(gene_short_name) -> gene_to_cluster
    gene_to_cluster <- gene_to_cluster$gene_short_name
    
    
    x = as.matrix(cds2[[i]]@assayData[["exprs"]][gene_to_cluster,])
    
    
    annotations <- cds2[[i]]@phenoData@data[colnames(x),
                                            c("seurat_clusters", "Pseudotime")] %>% as.data.frame()
    
    col = hue_pal()(length(unique(tcells_combined_umap_list_res[[i]]$seurat_clusters)))[unique(annotations[,"seurat_clusters"])]
    col = as.matrix(col)
    rownames(col) = unique(annotations[,"seurat_clusters"])
    ha <- HeatmapAnnotation(df = annotations, col = list(seurat_clusters = col[,1]))
    
    new_monocle2_heatmap[[i]] = Heatmap(x,
                                      show_column_names = FALSE,
                                      show_row_names = TRUE,
                                      top_annotation = ha,
                                      show_column_dend = FALSE,
                                      row_names_gp = grid::gpar(fontsize = 8))
  }
  
  saveRDS(new_monocle2_heatmap, "new_monocle2_heatmap.rds")
  
  
  
  cell_colors_clust = lapply(tcells_combined_umap_list_res, function(x){
    cell_pal(x$seurat_clusters, hue_pal())
  }) 
  saveRDS(cell_colors_clust, "cell_colors_clust.rds")

  
  #to have different starting clusters for slingshot:
  newsce = vector(mode = "list")    
  diffstartclus = vector(mode = "list")
  for(t in 1:length(tcells_combined_umap_list_res)){
    newsce[[t]] <- as.SingleCellExperiment(tcells_combined_umap_list_res[[t]])
    reducedDim(newsce[[t]]) <- reducedDim(newsce[[t]])[, 1:10]
    for(i in 1:length(unique(tcells_combined_umap_list_res[[t]]@meta.data$seurat_clusters))){
      diffstartclus[[i]] = vector(mode="list")
      diffstartclus[[i]][["UMAP"]] = slingshot(newsce[[t]], reducedDim = 'UMAP', clusterLabels = 'seurat_clusters',start.clus = i)
      diffstartclus[[i]][["UMAP"]] = SlingshotDataSet(diffstartclus[[i]][["UMAP"]]) 
      diffstartclus[[i]][["PHATE"]] = slingshot(newsce[[t]], reducedDim = 'PHATE', clusterLabels = 'seurat_clusters',start.clus = i)
      diffstartclus[[i]][["PHATE"]] = SlingshotDataSet(diffstartclus[[i]][["PHATE"]]) 
      
    }
  }
  saveRDS(diffstartclus, "diffstartclus.rds")
  
  
  tcells_combined_umap_list_res_skinny = tcells_combined_umap_list_res
  tcells_combined_umap_list_res_skinny[[1]]@assays$integrated = NULL
  tcells_combined_umap_list_res_skinny[[1]]@assays$SCT = NULL
  tcells_combined_umap_list_res_skinny[[2]] =NULL
  saveRDS(tcells_combined_umap_list_res_skinny, "tcells_combined_umap_list_res_skinny.rds")
  
}

##functions
##function to make sidebar menu expanded my default
modify_stop_propagation <- function(x) {
  x$children[[1]]$attribs$onclick = "event.stopPropagation()"
  x
}
