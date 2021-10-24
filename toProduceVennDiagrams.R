
#venn diagrams produced for the thesis 
library(Venn)
l = 1
c = list("Slingshot UMAP"=slingUMAPHeat_ALL[[l]]@row_names_param[["labels"]], "Slingshot PHATE" = slingPHATEHeat_ALL[[l]]@row_names_param[["labels"]],
         "Monocle 2" = new_monocle2_heatmap[[l]]@row_names_param$labels, "Monocle 3" = new_monocle3_heatmap[[l]]@row_names_param$labels, 
         "tradeSeq UMAP"=tradeUMAPHeat_ALL[[l]]@row_names_param[["labels"]], "tradeSeq PHATE"=tradePHATEHeat_ALL[[l]]@row_names_param[["labels"]])
venn(c, counts, snames = "", ilabels = FALSE, ellipse = FALSE, zcolor = c("green", "lightgreen", "red", "yellow", "navyblue", "lightblue"), 
     opacity = 0.3, plotsize = 15, ilcs = 1, sncs = 1, borders = TRUE, box = FALSE, par = TRUE, ggplot = FALSE)
