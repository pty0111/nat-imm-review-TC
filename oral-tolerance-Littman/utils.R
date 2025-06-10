'%!in%' <- function(x,y)!('%in%'(x,y))
'%ni%' <- function(x,y)!('%in%'(x,y))

cols <- rev(RColorBrewer::brewer.pal(11,"Spectral"))
transitions <- c(0, 25, 30, 35, 40, 50, 60, 65, 70, 75, 100)
scaled_transitions <- scales::rescale(transitions, from = c(0, 100), to = c(0, 1))
scale.color <- scale_color_gradientn(colors = cols, values = scaled_transitions)

plot.QC.violin <- function(cell.md, ident, feature, yintercept = NULL, br, Log10 = F, pal = NULL, fontsize.xlab = 10){
  if(feature == "nCount_RNA"){ name <- "number of transcripts"
  } else if(feature == "nFeature_RNA"){ name <- "number of detected genes"
  } else if(feature == "MtFrac_RNA"){ name <- "percentage of mitochondrial transcripts"
  } else if(feature == "percent.MT"){ name <- "percentage of mitochondrial transcripts"
  } else if(feature == "S.Score"){ name <- "S-phase score"
  } else if(feature == "G2M.Score") {name <- "G2/M-phase score"
  } else if(feature == "ATAC.nFrags") {name <- "ATAC.nFrags"
  } else if(feature == "ATAC.TSSEnrichment") {name <- "ATAC.TSSEnrichment"
  } else if(feature == "ATAC.DoubletEnrichment") {name <- "ATAC.DoubletEnrichment"
  } else if(feature == "ATAC.DoubletScore") {name <- "ATAC.DoubletScore"
  } else if(feature == "nFrags") {name <- "nFrags"
  } else if(feature == "TSSEnrichment") {name <- "TSSEnrichment"
  } else if(feature == "DoubletEnrichment") {name <- "DoubletEnrichment"
  }else if(feature == "DoubletScore") {name <- "DoubletScore"
  } else name <- feature
  ft <- cell.md[, feature]
  if(Log10){
    ft <- log10(ft)
    if(!is.null(yintercept)) yintercept <- log10(yintercept)
  }
  gp <- ggplot(data.frame(id = cell.md[, ident], ft = ft),
               aes(x = id, y = ft, color = id, fill = id)) +
    geom_violin(show.legend = F) +
    geom_hline(yintercept = yintercept, linetype = 2, color = "violetred4") +
    scale_y_continuous(breaks = seq(0, max(ft), by = br)) +
    theme_classic() +
    theme(axis.text.x = element_text(size=fontsize.xlab)) +
    labs(x = ident, y = "", title = ifelse(Log10, paste0("log10(", name, ")"), name))
  if(!is.null(pal))
    gp <- gp + scale_color_manual(values = pal) + scale_fill_manual(values = pal)
  return(gp)
}

plot.all.QC <- function(obj = NULL, obj.type = "Seurat", ident = 'sample', 
                            thr = NULL, col = NULL, CC = F){
  if (obj.type == "ArchR"){
    cell.md <- data.frame(obj@cellColData)
    cell.md$log10nFrags <- log10(cell.md$nFrags)
    curr.assay <- "ATAC"
    L <- labs(color = ident, x = "number of fragments (log10)", y = "TSSenrichment")
  }else if (obj.type == "Seurat"){
    cell.md <- obj@meta.data
    curr.assay <- "RNA"
    L <- labs(color = ident, x = "number of transcripts", y = "number of detected genes")
  } else if (obj.type == "CellSpace"){
    cell.md <- data.frame(obj@meta.data)
    cell.md$log10nFrags <- log10(cell.md$nFrags)
    curr.assay <- "ATAC"
    L <- labs(color = ident, x = "number of fragments (log10)", y = "TSSenrichment")
  } else{
    print("Invalid input")
    return(1)
  }
  if (curr.assay=='RNA'){
    gp1 <- FeatureScatter(sro, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = ident) + theme_classic() + L
    if(!is.null(col)){
      names(col) <- levels(cell.md[, ident])
      gp1 <- gp1 + scale_color_manual(values = col)
    }
    gp2 <- ggplot(cell.md, aes(x = nCount_RNA, y = nFeature_RNA)) + geom_density_2d_filled() + theme_classic() + L
    
    pl <- list(
      plot.QC.violin(cell.md, ident = ident, feature = "nCount_RNA", Log10 = T, yintercept = c(thr$nc.min, thr$nc.max), br = 0.2, pal = col),
      plot.QC.violin(cell.md, ident = ident, feature = "nFeature_RNA", yintercept = c(thr$nf.min, thr$nf.max), br = 500, pal = col),
      plot.QC.violin(cell.md, ident = ident, feature = "MtFrac_RNA", yintercept = thr$mp.max, br = 2, pal = col)
    )
    if(is.null(thr)){
      pl[[4]] <- gp1
      pl[[5]] <- gp2
    } else {
      hl <- geom_hline(yintercept = thr$nf.min, linetype = 2, color = "violetred4")
      vl <- geom_vline(xintercept = thr$nc.min, linetype = 2, color = "violetred4")
      pl[[4]] <- gp1 + hl + vl
      pl[[5]] <- gp2 + hl + vl
    }
    if (CC) {
      pl[[6]] <- plot.QC.violin(cell.md, ident = ident, feature = "S.Score", br = 0.2, pal = col)
      pl[[7]] <- plot.QC.violin(cell.md, ident = ident, feature = "G2M.Score", br = 0.2, pal = col)
    }
  } else {
    gp1 <- ggplot(cell.md, aes(x = log10nFrags, y = TSSEnrichment, color = .data[[ident]])) + geom_point() + theme_classic() + L
    if(!is.null(col)){
      names(col) <- levels(cell.md[, ident])
      gp1 <- gp1 + scale_color_manual(values = col)
    }
    gp2 <- ggplot(cell.md, aes(x = log10nFrags, y = TSSEnrichment)) + geom_density_2d_filled() + theme_classic() + L
    pl <-list(
      plot.QC.violin(cell.md, ident = ident, feature = "nFrags", Log10 = T, yintercept = c(thr$nfrags.min, thr$nfrags.max), br = 0.2, pal = col),
      plot.QC.violin(cell.md, ident = ident, feature = "TSSEnrichment", yintercept = thr$tss.min, br = 5, pal = col),
      plot.QC.violin(cell.md, ident = ident, feature = "DoubletEnrichment", yintercept = thr$doublet.enrichment.max, br = 5, pal = col)
    )
    if(is.null(thr)){
      pl[[4]] <- gp1
      pl[[5]] <- gp2
    } else {
      thr$log10nfrags.min <- log10(thr$nfrags.min)
      hl <- geom_hline(yintercept = thr$tss.min, linetype = 2, color = "violetred4")
      vl <- geom_vline(xintercept = thr$log10nfrags.min, linetype = 2, color = "violetred4")
      pl[[4]] <- gp1 + hl + vl
      pl[[5]] <- gp2 + hl + vl
    }
  }
  return(pl)
}

# black border
plot.clusters <- function(
    SRO = NULL, idx = NULL,
    vis = SRO@reductions$umap@cell.embeddings, 
    axis.titles = c("UMAP1", "UMAP2"),
    groups = SRO$Clusters, clusters.col = 'Clusters', 
    col = pal[[clusters.col]], new.levels = NULL,
    pref.C = T, labels = T, line = T,
    point.size = 1.5, point.alpha = 0.7,
    label.size = point.size * 2, label.pad = 1, border.size = 0.5, 
    legend.size = 3, show_legend = T
){
  if(class(groups) == "factor"){ 
    clusters <- groups
  } else {
    clusters <- factor(groups, levels = sort(unique(groups)))
  }
  gp <- ggplot() + theme_classic() +
    labs(x = axis.titles[1], y = axis.titles[2], color = clusters.col) +
    theme(axis.ticks = element_blank(), axis.text = element_blank())
  if(!line) gp <- gp + theme(axis.line = element_blank())
  if(!is.null(idx)){
    vis <- vis[idx, ]
    if(is.null(new.levels)){
      new.levels <- as.character(sort(unique(clusters[idx])))
    }
    clusters <- factor(as.character(clusters[idx]), levels = new.levels)
    if(!is.null(col)) col <- col[new.levels]
  }
  gp <- gp + geom_point(
    mapping = aes(x = vis[, 1], y = vis[, 2], color = clusters),
    size = point.size, alpha = point.alpha, show.legend = show_legend
  )
  if(!is.null(col)) gp <- gp + scale_color_manual(values = col)
  if(labels){
    C <- levels(clusters)
    CL <- paste0(ifelse(pref.C, "C", ""), levels(clusters))
    vis.cl <- data.frame(matrix(NA, nrow = length(C), ncol = 2, dimnames = list(C, NULL)))
    for(c in C){
      vis.cl[c, 1] <- median(vis[clusters == c, 1], na.rm = T)
      vis.cl[c, 2] <- median(vis[clusters == c, 2], na.rm = T)
    }
    gp + geom_label_repel(
      aes(x = vis.cl[, 1], vis.cl[, 2], color = C, label = CL),
      color = 'black', seed = 1,
      label.size = border.size,
      size = label.size, label.padding = unit(label.pad, "mm"), show.legend = F
    ) +
      geom_label_repel(
        aes(x = vis.cl[, 1], vis.cl[, 2], color = C, label = CL),
        size = label.size, seed = 1,
        label.size = NA, label.padding = unit(label.pad, "mm"), show.legend = F
      )+ guides(colour = guide_legend(override.aes = list(size = legend.size)))
  } else gp + guides(colour = guide_legend(override.aes = list(size = legend.size)))
}

plot.clusters.highlight.one <- function(
    SRO = NULL, idx1 = NULL, idx2 = NULL,
    vis = SRO@reductions$umap@cell.embeddings, 
    axis.titles = c("UMAP1", "UMAP2"),
    groups = SRO$Clusters, clusters.col = 'Clusters', 
    col = pal[[clusters.col]], new.levels = NULL,
    pref.C = T, labels = T, line = T,
    point.size = 1.5, point.alpha = 0.7,
    label.size = point.size * 2, label.pad = 1
){
  if(class(groups) == "factor"){ 
    clusters <- groups
  } else {
    clusters <- factor(groups, levels = sort(unique(groups)))
  }
  gp <- ggplot() + theme_classic() +
    labs(x = axis.titles[1], y = axis.titles[2], color = clusters.col) +
    theme(axis.ticks = element_blank(), axis.text = element_blank())
  if(!line) gp <- gp + theme(axis.line = element_blank())
  if(!is.null(idx1)){
    vis1 <- vis[idx1, ]
    if(is.null(new.levels)){
      new.levels <- as.character(sort(unique(clusters[idx1])))
    }
    clusters1 <- factor(as.character(clusters[idx1]), levels = new.levels)
    if(!is.null(col)) col <- col[new.levels]
  }
  if(!is.null(idx2)){
    vis2 <- vis[idx2, ]
    if(is.null(new.levels)){
      new.levels <- as.character(sort(unique(clusters[idx2])))
    }
    clusters2 <- factor(as.character(clusters[idx2]), levels = new.levels)
  }
  gp <- gp + geom_point(
    mapping = aes(x = vis1[, 1], y = vis1[, 2]), color = 'lightgray',
    size = point.size, alpha = point.alpha
  ) + geom_point(
    mapping = aes(x = vis2[, 1], y = vis2[, 2], color = clusters2),
    size = point.size, alpha = point.alpha, show.legend = T
  )
  if(!is.null(col)) gp <- gp + scale_color_manual(values = col)
  if(labels){
    C <- levels(clusters2)
    CL <- paste0(ifelse(pref.C, "C", ""), levels(clusters2))
    vis.cl <- data.frame(matrix(NA, nrow = length(C), ncol = 2, dimnames = list(C, NULL)))
    for(c in C){
      vis.cl[c, 1] <- median(vis2[clusters2 == c, 1], na.rm = T)
      vis.cl[c, 2] <- median(vis2[clusters2 == c, 2], na.rm = T)
    }
    gp + geom_label_repel(
      aes(x = vis.cl[, 1], vis.cl[, 2], label = CL), 
      color = 'black',
      seed=1,
      size = label.size, 
      label.size = 0.5,
      label.padding = unit(label.pad, "mm"), show.legend = F
    )
  } else gp
}

plot.clusters.label.dots <- function(
    SRO = NULL, idx = NULL,
    vis = SRO@reductions$umap@cell.embeddings,
    cells = NULL, label.color = 'black',
    clusters.col = "Clusters", col = pal[[clusters.col]],
    pref.C = T, labels = T, line = F,
    point.size = 1.5, point.alpha = 0.7,
    label.size = point.size * 2, label.pad = 1,
    label.shift = 0.3, dist.thr = 1
){
  if(class(SRO$Clusters) == "factor"){ clusters <- SRO@meta.data[, clusters.col]
  } else clusters <- factor(SRO@meta.data[, clusters.col], levels = names(col))
  gp <- ggplot() + theme_classic() + labs(color = clusters.col) +
    theme(axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank())
  if(!line) gp <- gp + theme(axis.line = element_blank())
  if(!is.null(idx)){
    vis <- vis[idx, ]
    new.levels <- as.character(sort(unique(clusters[idx])))
    clusters <- factor(as.character(clusters[idx]), levels = new.levels)
    col <- col[new.levels]
  }
  gp <- gp + geom_point(
    mapping = aes(x = vis[, 1], y = vis[, 2], color = clusters),
    size = point.size, alpha = point.alpha
  )
  if(!is.null(col)) gp <- gp + scale_color_manual(values = col)
  if(labels){
    set.seed(0)
    
    C <- levels(clusters)
    CL <- paste0(ifelse(pref.C, "C", ""), levels(clusters))
    vis.cl <- data.frame(matrix(NA, nrow = length(C), ncol = 2, dimnames = list(C, NULL)))
    for(c in C){
      vis.cl[c, 1] <- median(vis[clusters == c, 1]) + label.shift * sample(c(1, -1), 1)
      vis.cl[c, 2] <- median(vis[clusters == c, 2]) + label.shift * sample(c(1, -1), 1)
    }
    
    while(T){
      label.dist <- reshape2::melt(as.matrix(dist(vis.cl))) %>%
        subset(Var1 != Var2 & value < dist.thr * sqrt(2))
      if(nrow(label.dist) == 0) break
      dummy <- apply(label.dist, 1, function(pair){
        dim <- sample(1:2, size = 1)
        shift.dim <- vis.cl[pair[1], dim] - vis.cl[pair[2], dim]
        if(shift.dim < dist.thr)
          vis.cl[pair[1], dim] <<- vis.cl[pair[1], dim] + shift.dim / 2
      })
    }
    
    gp <- gp + geom_label(
      aes(x = vis.cl[, 1], vis.cl[, 2], color = C, label = CL),
      size = label.size, label.padding = unit(label.pad, "mm"), show.legend = F
    )
  }
  if (!is.null(cells)){
    vis.cl2 <- data.frame(matrix(NA, nrow = length(cells), ncol = 2))
    rownames(vis.cl2) <- cells
    for(c in cells){
      vis.cl2[c, 1] <- vis[rownames(vis) == c, 1] + label.shift * sample(c(1, -1), 1)
      vis.cl2[c, 2] <- vis[rownames(vis) == c, 2] + label.shift * sample(c(1, -1), 1)
    }
    gp + geom_label_repel(
      aes(x = vis.cl2[, 1], vis.cl2[, 2], label = rownames(vis.cl2)), 
      arrow = arrow(length = unit(0.015, "npc")),
      color = label.color,
      size = label.size, label.padding = unit(label.pad, "mm"), show.legend = F, box.padding = 0.5
    )
  } else gp
}

plot.continuous.value <- function(
    sro = NULL, idx = NULL,
    vis = sro@reductions$umap@cell.embeddings,
    val, val.name, axis.blank = T, axis.titles = c("UMAP1", "UMAP2"),
    scale.color=scale_color_distiller(palette = "Spectral"), point.size = 1
){
  gp <- ggplot() +
    geom_point(
      mapping = aes(x = vis[idx, 1], y = vis[idx, 2], color = val[idx]),
      size = point.size, alpha = 0.8
    ) +
    scale.color + theme_classic() + labs(x = axis.titles[1], y = axis.titles[2], color = val.name)
  if (axis.blank){
    gp <- gp + theme(axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank())
  } else{
    gp
  }
  return(gp)
}

break.down.bar.plot <- function(md, group1, group2, col = pal[[group2]], 
                                axis.text.x.size = 12, axis.text.y.size = 12,
                                label.groups = F, threshold = 0, rotate.x.label = F, ...){
  df <- md[, c(group1, group2)] %>% table() %>%
    as.data.frame() %>% subset(Freq > 0)
  colnames(df)[1:2] <- c("group1", "group2")
  sum.count <- table(md[, group1])
  df$Perc <- round(df$Freq / as.numeric(sum.count[df$group1]), 3) * 100
  df <- df[order(df$Perc, decreasing = F), ]
  df$group2.lab <- ifelse(df$Perc > threshold, gsub(" ", "\n", df$group2), NA)
  
  gp <- ggplot(df, aes(x = group1, y = Perc)) +
    geom_bar(aes(fill = group2), stat = "identity", position = "stack") +
    geom_label(aes(label = paste0(Perc, "%"), color = group2), show.legend = F,
               stat = "identity", position = "stack", ...) +
    scale_color_manual(values = col) +
    scale_fill_manual(values = col) +
    scale_y_continuous(breaks = seq(0, 100, by = 10)) +
    labs(x = group1, y = "percentage of cells", fill = group2) + theme_bw() +
    theme(panel.grid.major = element_blank(), 
          axis.text.x = element_text(size = axis.text.x.size), 
          axis.text.y = element_text(size = axis.text.y.size)
    )
  if (label.groups){
    gp <- gp + geom_label(aes(label = group2.lab, color = group2), show.legend = F,
                          stat = "identity", position = "stack", vjust = 1.5)
  }
  if (rotate.x.label){
    gp <- gp + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1))
  }
  return(gp)
}

break.down.tile.plot <- function(md, group1, group2, label.cells = T, 
                                 nonzero.only = T, 
                                 label.nonzero.only = T,
                                 label.size = 3, scale.color = scale_fill_distiller(palette = "Spectral")){
  eb <- element_blank()
  if (nonzero.only){
    df <- md[, c(group1, group2)] %>% table() %>% as.data.frame() %>% subset(Freq > 0)
  } else {
    df <- md[, c(group1, group2)] %>% table() %>% as.data.frame()
  }
  colnames(df)[1:2] <- c("group1", "group2")
  if(class(md[, group1]) == "factor"){
    group1.levels <- levels(md[, group1])
  } else {
    group1.levels <- sort(unique(md[, group1]))
  }
  if(class(md[, group2]) == "factor"){
    group2.levels <- levels(md[, group2])
  } else {
    group2.levels <- sort(unique(md[, group2]))
  }
  df$group1 <- factor(df$group1, levels = group1.levels)
  df$group2 <- factor(df$group2, levels = group2.levels)
  sum.count <- table(md[, group1])
  df$Perc <- round(df$Freq / as.numeric(sum.count[df$group1]), 3) * 100
  df <- df[order(df$Perc, decreasing = F), ]
  df2 <- df[df$Perc > 0,]
  gp <- ggplot(df, aes(x = group1, y = group2)) +
    geom_tile(aes(fill = Perc), color = "black")
  if (label.cells){
    if (label.nonzero.only){
      gp <- gp + geom_label(data = df2, aes(label = round(Perc, 3), color = Perc), size = label.size)
    } else {
      gp <- gp + geom_label(aes(label = round(Perc, 3), color = Perc), size = label.size)
    }
  }
  gp <- gp + scale.color +
    labs(x = group1, y = group2, title = "Composition of cells across groups") + 
    theme_classic() + theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16)) + theme(legend.position = "none")
  return(gp)
}
# 
# break.down.tile.plot <- function(md, group1, group2, nonzero.only = T, 
#                                  label.nonzero.only = T,
#                                  label.size = 3, ...){
#   eb <- element_blank()
#   if (nonzero.only){
#     df <- md[, c(group1, group2)] %>% table() %>% as.data.frame() %>% subset(Freq > 0)
#   } else {
#     df <- md[, c(group1, group2)] %>% table() %>% as.data.frame()
#   }
#   colnames(df)[1:2] <- c("group1", "group2")
#   if(class(md[, group1]) == "factor"){
#     group1.levels <- levels(md[, group1])
#   } else {
#     group1.levels <- sort(unique(md[, group1]))
#   }
#   if(class(md[, group2]) == "factor"){
#     group2.levels <- levels(md[, group2])
#   } else {
#     group2.levels <- sort(unique(md[, group2]))
#   }
#   df$group1 <- factor(df$group1, levels = group1.levels)
#   df$group2 <- factor(df$group2, levels = group2.levels)
#   sum.count <- table(md[, group1])
#   df$Perc <- round(df$Freq / as.numeric(sum.count[df$group1]), 3) * 100
#   df <- df[order(df$Perc, decreasing = F), ]
#   df2 <- df[df$Perc > 0,]
#   gp <- ggplot(df, aes(x = group1, y = group2)) +
#     geom_tile(aes(fill = Perc), color = "black") 
#   if (label.nonzero.only){
#     gp <- gp + geom_label(data = df2, aes(label = round(Perc, 3), color = Perc), size = label.size)
#   } else {
#     gp <- gp + geom_label(aes(label = round(Perc, 3), color = Perc), size = label.size)
#   }
#   gp <- gp + scale_fill_distiller(palette = "Spectral") +
#     labs(x = group1, y = group2, title = "Composition of cells across groups") +
#     theme_classic() + theme(legend.position = "none")
#   return(gp)
# }


plot.histogram <- function(df, x.col, x.title = x.col,  title.lab = "", 
                           nbins = 50, color = 'white',
                           fill = 'lightblue', title.size = 16,
                           axis.title.size = 14, axis.text.size = 12){
  ggplot() + geom_histogram(data = df, aes(x = .data[[x.col]]), 
                            bins = nbins, 
                            color = color, fill = fill) + 
    theme_bw() + labs(title = title.lab, x = x.title) +
    theme(axis.title = element_text(size = axis.title.size),
          axis.text = element_text(size = axis.text.size), 
          plot.title = element_text(size = title.size)
    )
}

run.PhenoGraph <- function(sro, npcs, k){
  sro@misc[["PhenoGraph"]] <- Rphenograph(sro@reductions$pca@cell.embeddings[, 1:npcs], k = k)
  adj <- as_adjacency_matrix(sro@misc$PhenoGraph[[1]], sparse = F)
  rownames(adj) <- colnames(sro); colnames(adj) <- colnames(sro)
  sro@graphs[["PhenoGraph"]] <- as.Graph(adj)
  sro$Clusters <- sro@misc$PhenoGraph[[2]]$membership
  sro$Clusters <- factor(sro$Clusters, levels = min(sro$Clusters):max(sro$Clusters))
  return(sro)
}

run.fdl <- function(sobject,
                    graph = "RNA_snn",
                    weighted = TRUE, 
                    dims = 2) {
  if (!graph %in% names(sobject@graphs)) {
    stop(graph, " graph not found in Seurat object")
  }
  
  graph <-
    igraph::graph_from_adjacency_matrix(
      adjmatrix = sobject[[graph]],
      mode = "undirected",
      weighted = weighted,
      add.colnames = TRUE
    )
  
  fdl <- igraph::layout_with_fr(graph, grid = "nogrid", dim = dims)
  
  rownames(fdl) <- colnames(sobject)
  colnames(fdl) <- paste0("fdr_", 1:dims)
  
  sobject[["fdl"]] <-
    SeuratObject::CreateDimReducObject(
      embeddings = fdl,
      key = "fdl_",
      assay = SeuratObject::DefaultAssay(sobject)
    )
  
  return(sobject)
}

# MAST ####
MAST.DE <- function(fn, sro, output.dir = 'results/markers/', ...){
  m <- suppressWarnings(FindMarkers(object = sro, test.use = "MAST", logfc.threshold = 0, ...))
  m$gene_id <- sro@assays$RNA@meta.features[rownames(m), "ID"]
  m$gene_name <- sro@assays$RNA@meta.features[rownames(m), "symbol.unique"]
  m$gene <- rownames(m)
  write.csv(m, file = paste0(output.dir, fn, ".csv"), row.names = T, quote = F)
}

MAST.DE.multiple <- function(fn, sro, idents = "Clusters", output.dir = 'results/markers/', ...){
  Idents(sro) <- sro@meta.data[, idents]
  m <- FindAllMarkers(object = sro, test.use = "MAST", logfc.threshold = 0, ...)
  m$gene_id <- sro@assays$RNA@meta.features[m$gene, "ID"]
  m$gene_name <- sro@assays$RNA@meta.features[m$gene, "symbol.unique"]
  write.csv(m, file = paste0(output.dir, fn, ".csv"), row.names = F, quote = F)
  return(m)
}

select.markers <- function(fn, output.dir = 'results/markers/', pairwise = F, 
                           fc.thr = 1.5, apv.thr = 0.01, n = Inf, 
                           remove.MT = T, remove.duplicated = T){
  markers <- read.csv(file = paste0(output.dir, fn, ".csv"), header = T, stringsAsFactors = F)
  if (remove.MT){
    markers <- markers %>% subset(!grepl("(^MT-)|(^RPS)|(^RPL)|(^MRPL)|(^MRPS)", toupper(gene)))
  }
  if(pairwise){
    markers <- subset(markers, abs(avg_log2FC) > log2(fc.thr) & p_val_adj < apv.thr)
    markers$cluster <- ifelse(markers$avg_log2FC > 0, "up", "down")
    if(!is.infinite(n)) markers <- markers %>% group_by(cluster) %>% top_n(n, abs(avg_log2FC))
    markers <- markers[order(markers$avg_log2FC), ]
  } else{
    markers <- subset(markers, avg_log2FC > log2(fc.thr) & p_val_adj < apv.thr)
    if (remove.duplicated){
      if(any(duplicated(markers$gene))) markers <- markers %>% group_by(gene) %>% top_n(1, avg_log2FC)
    }
    if(!is.infinite(n)) markers <- markers %>% group_by(cluster) %>% top_n(n, avg_log2FC)
    markers <- markers[order(markers$cluster, -markers$avg_log2FC), ]
  }
  return(markers)
}

select.archr.markers <- function(fn, pairwise = F, lfc.thr = 0, apv.thr = 0.01, 
                                 lfc.colname = 'MeanDiff', apv.colname = 'FDR',
                                 n = Inf){
  markers <- read.csv(file = fn, header = T, stringsAsFactors = F) %>%
    subset(!grepl("(^MT-)|(^RPS)|(^RPL)|(^MRPL)|(^MRPS)", toupper(name)))
  if(pairwise){
    markers <- markers[(abs(markers[lfc.colname]) > lfc.thr) & (markers[apv.colname] < apv.thr), ] # markers[lfc.colname] does not need abs
    if(!is.infinite(n)) markers <- markers %>% group_by(group_name) %>% top_n(n, abs(.data[[lfc.colname]]))
    markers <- markers %>% dplyr::arrange(group_name, desc(.data[[lfc.colname]]))
  } else{
    markers <- markers[(markers[lfc.colname] > lfc.thr) & (markers[apv.colname] < apv.thr), ]
    if(any(duplicated(markers$name))) markers <- markers %>% group_by(name) %>% top_n(1, .data[[lfc.colname]])
    if(!is.infinite(n)) markers <- markers %>% group_by(group) %>% top_n(n, abs(.data[[lfc.colname]]))
    markers <- markers %>% dplyr::arrange(group, desc(.data[[lfc.colname]]))
  }
  return(markers)
}

DE.heatmap <- function(
    sro = NULL, expr = sro@assays$RNA@data, md = sro@meta.data, 
    markers, cells = rownames(md), split = "Clusters", 
    cluster_columns = T, cluster_column_slices = T, 
    column_title_rot = 0, row_title_rot = 90,
    rowsplit = NULL, cluster_rows = F, cluster_row_slices = F,
    legend.name = "scaled\nimputed\nexpression", 
    scale = T,
    ca.col = list(Sample = pal$sample, Clusters = pal$Clusters, celltype=pal$celltype),
    split.columns = names(ca.col),
    rasterize = T,
    ...
){
  annot.df <- data.frame(md[cells, split.columns])
  colnames(annot.df) <- split.columns
  ca <- columnAnnotation(
    df = annot.df,
    col = ca.col
  )
  col.ramp <- colorRamp2(c(-10, -5, seq(-3, 3, 1), 5, 10),
                         c("#173c68", "#173c68", "#1c5785", "#166db0", "#a6d1e3", "#ffffff",
                           "#f09173", "#d46052", "#b32339", "#690822", "#690822"))
  expr <- expr[markers, cells]
  if (scale) {expr <- t(scale(t(expr)))}
  Heatmap(
    matrix = expr, name = legend.name, col = col.ramp,
    show_column_names = F, 
    column_split = md[cells, split], top_annotation = ca,
    cluster_rows = cluster_rows, cluster_row_slices = cluster_row_slices,
    row_names_side = "right", row_split = rowsplit, 
    column_title_rot = column_title_rot, row_title_rot = row_title_rot,
    cluster_columns = cluster_columns, cluster_column_slices = cluster_column_slices,
    row_names_gp = gpar(fontface = "italic", fontsize = 10), use_raster = rasterize
  )
}

DE.heatmap.ra <- function(
    sro = NULL, expr = sro@assays$RNA@data, md = sro@meta.data, marker.info, marker.colname = "Family",
    markers, cells = rownames(md), split = "Clusters", row_labels = markers,
    cluster_columns = T, cluster_column_slices = T, column_title_rot = 0, show_column_names = F,
    rowsplit = NULL, cluster_rows = F, cluster_row_slices = F, show_row_names = T,
    legend.name = "scaled\nimputed\nexpression", scale = T,
    split.columns = c("Sample", "Clusters", 'celltype'),
    ca.col = list(Sample = pal$sample, Clusters = pal$Clusters, celltype=pal$celltype),
    ra.col = NULL, colramp = 1,
    ...
){
  annot.df <- data.frame(md[cells, split.columns])
  colnames(annot.df) <- split.columns
  ca <- columnAnnotation(
    df = annot.df,
    col = ca.col
  )
  rowannot.df <- data.frame(marker.info[c(marker.colname)])
  ra = rowAnnotation(
    df = rowannot.df,
    col = ra.col
  )
  if (colramp == 1){
    col.ramp <- colorRamp2(c(-10, -5, seq(-3, 3, 1), 5, 10),
                           c("#173c68", "#173c68", "#1c5785", "#166db0", "#a6d1e3", "#ffffff",
                             "#f09173", "#d46052", "#b32339", "#690822", "#690822"))
  } else {
    zs.cols <- c(RColorBrewer::brewer.pal(name = "RdBu", n = 7)[7:4],
                 RColorBrewer::brewer.pal(name = "YlOrRd", n = 7)[1:2])
    col.ramp <- colorRamp2(c(-10, seq(-5, 5, length.out = 9), 10),
                           c("#1c1768", rep(zs.cols[1], 2), rep(zs.cols[2], 2), zs.cols[4],
                             rep(zs.cols[6], 2), rep("#ffbb03", 3)))
  }
  expr <- expr[markers, cells]
  if (scale) {expr <- t(scale(t(expr[markers, cells])))} else {expr <- as.matrix(expr[markers, cells])}
  Heatmap(
    matrix = expr, name = legend.name, col = col.ramp, row_labels = row_labels,
    show_column_names = show_column_names, show_row_names = show_row_names,
    column_split = md[cells, split], top_annotation = ca, left_annotation = ra,
    cluster_rows = cluster_rows, cluster_row_slices = cluster_row_slices,
    row_names_side = "right", row_split = rowsplit, column_title_rot = column_title_rot,
    cluster_columns = cluster_columns, cluster_column_slices = cluster_column_slices,
    row_names_gp = gpar(fontface = "italic", fontsize = 10), use_raster = T, ...
  )
}

average.score.mat <- function(expr, cell.md, 
                              features = rownames(expr), 
                              group = "annotations", gene = 'gene'){
  # expr (dataframe): feature by cell score matrix
  expr <- as.data.frame(t(expr[features, ]))
  expr$Run <- rownames(expr)
  melted <- reshape2::melt(expr, id = "Run")
  colnames(melted) <- c("Run", gene, "score")
  group_cols <- c(gene, group)
  melted.2 <- dplyr::inner_join(melted, cell.md, by='Run')
  melted.2 <- melted.2 %>% 
    dplyr::group_by(across(all_of(group_cols))) %>% 
    dplyr::summarise(avg_score = mean(score))
  wide.mtx <- reshape2::dcast(melted.2, paste0(gene, " ~ ", group))
  rownames(wide.mtx) <- wide.mtx$gene
  wide.mtx$gene <- NULL
  wide.mtx <- as.matrix(wide.mtx)
  return(wide.mtx)
}

DE.heatmap.pseudobulk <- function(
    sro = NULL, expr = NULL, assay = "RNA", md = sro@meta.data, 
    markers, marker.info = NULL,
    curr.ident = "Cluster", 
    column_title_rot = 0, row_title_rot = 90,
    cluster_columns = T, 
    scale = T,
    rowsplit = NULL, cluster_rows = F, cluster_row_slices = F,
    legend.name = "scaled\nimputed\nexpression",
    split.columns = c("Clusters"),
    ca.col = list(Clusters = pal$Clusters),
    ra.col = NULL,
    ...
){
  if (assay == "RNA"){
    avg.expr <- AverageExpression(sro, assays = assay, features = markers, group.by = curr.ident)[[1]]
  } else if (assay == "ATAC"){
    avg.expr <- average.score.mat(expr = expr, cell.md = md, features = markers, group = curr.ident)
  }
  md <- data.frame(md %>% group_by(.data[[curr.ident]]) %>% dplyr::filter(row_number()==1))
  annot.df <- data.frame(md[, split.columns])
  colnames(annot.df) <- split.columns
  ca <- columnAnnotation(
    df = annot.df,
    col = ca.col
  )
  if (!is.null(marker.info)){
    rowannot.df <- data.frame(marker.info[names(ra.col)])
    ra = rowAnnotation(
      df = rowannot.df,
      col = ra.col
    )
  }else ra <- NULL
  col.ramp <- colorRamp2(c(-10, -5, seq(-3, 3, 1), 5, 10),
                         c("#173c68", "#173c68", "#1c5785", "#166db0", "#a6d1e3", "#ffffff",
                           "#f09173", "#d46052", "#b32339", "#690822", "#690822"))
  if (scale) {avg.expr <- t(scale(t(avg.expr[markers, ])))} else {avg.expr <- avg.expr[markers, ]}
  Heatmap(
    matrix = avg.expr, name = legend.name, col = col.ramp,
    show_column_names = T, top_annotation = ca, left_annotation = ra,
    cluster_rows = cluster_rows, cluster_row_slices = cluster_row_slices,
    row_names_side = "right", row_split = rowsplit, 
    column_title_rot = column_title_rot, row_title_rot = row_title_rot,
    cluster_columns = cluster_columns, 
    row_names_gp = gpar(fontface = "italic", fontsize = 8), use_raster = T, ...
  )
}

plot.dot <- function(sro, curr.ident, genes, gene.labels = genes, scale = T, cluster.idents = F){
  p <- DotPlot(sro, assay = "RNA", features = genes, group.by = curr.ident, scale = scale,
               cols = c("blue", "red"), col.min = -3, col.max = 3, cluster.idents = cluster.idents) +
    labs(x = "", y = "Cluster") +
    scale_x_discrete(breaks = genes, labels = gene.labels) +
    # scale_y_discrete(breaks = levels(sro@active.ident)) +
    scale_size_continuous(limits = c(0, 100), range = c(0.1, 10)) +
    theme(
      axis.text.x.bottom = element_text(size = 15, angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y.left = element_text(size = 15),
      legend.position = "bottom", legend.justification = "center",
      legend.key.width = unit(0.5, "inch"), legend.box = "vertical",
      legend.title = element_text(vjust = 1, hjust = 1)
    )
  return(p)
}


# CellSpace ####
# extract cell by variable tile accessibility matrix:
# Re-order cell barcodes according to archr.obj@cellColData
get.tile.matrix.sep <- function(archr.obj, genome){
  var.tiles <- archr.obj@reducedDims$IterativeLSI$LSIFeatures[, -3]
  var.tiles.gr <- GRanges(
    seqinfo = seqinfo(genome),
    seqnames = var.tiles$seqnames, strand = "+",
    ranges = IRanges(
      start = var.tiles$start,
      width = archr.obj@reducedDims$IterativeLSI$tileSize,
      names = paste0("tile", var.tiles$idx)
    )
  )
  matrices <- c()
  for (chrom in getSeqnames(archr.obj)){
    print(chrom)
    chr.mtx <- getMatrixFromProject(archr.obj, useMatrix = "TileMatrix", useSeqnames = chrom, binarize = T)
    chr.var.tile.mtx = assays(chr.mtx)$TileMatrix[na.omit(match(var.tiles, chr.mtx@elementMetadata)), ]
    rm(chr.mtx)
    matrices <- c(matrices, chr.var.tile.mtx)
  }
  combined.mtx = do.call(rbind, matrices)
  # sort matrix columns (cell)
  ci <- match(rownames(archr.obj@cellColData), colnames(combined.mtx))
  sorted.mtx <- combined.mtx[, ci]
  
  return(list(var.tiles = var.tiles, 
              var.tiles.gr = var.tiles.gr, 
              genome = genome,
              var.tile.mtx = sorted.mtx))
}
# extract cell by variable tile accessibility matrix:

plot.loss <- function(results.dir, job.epoch, filename = "stdout.txt", 
                      x.ticks.interval = 50, plot.title.size = 34, axis.title.size = 28, axis.text.size = 24){
  stdout <- readLines(paste0(results.dir, filename))
  training.errors <- c()
  
  i <- 0
  for (line in stdout){
    if (startsWith(line, " ---+++ ")){
      words <- str_split(line, pattern=' ')[[1]]
      if (i<10){
        idx <- 26
      } else if (i<100){
        idx <- 25
      } else{
        idx <- 24
      }
      
      i <- i+1
      training.errors <- c(training.errors, as.numeric(words[idx]))
      
    }
  }
  
  training.df <- data.frame(epoch = seq(1,length(training.errors)), training.error = training.errors)
  
  p <- ggplot(data=training.df, aes(x=epoch, y=training.error)) +
    geom_line() +
    geom_point() + 
    scale_x_continuous(breaks = seq(0, job.epoch, by = x.ticks.interval)) +
    theme_classic() +
    theme(plot.title = element_text(size = plot.title.size),
          axis.title = element_text(size = axis.title.size),
          axis.text = element_text(size = axis.text.size),
          axis.text.x = element_text(angle=45, vjust=1, hjust=1)
    )
  return(p)
}



# misc helpers ####
add.custom.umap <- function(sro, umap, slot.name, key, 
                            global = T, assay = "RNA"){
  sro[[slot.name]] <- CreateDimReducObject(embeddings = as.matrix(umap), 
                                           key = key, global = global, 
                                           assay = assay)
  return(sro)
}

RenameGenesSeurat <- function(obj, newnames) { # Replace gene names in different slots of a Seurat object. Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.
  print("Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.")
  RNA <- obj@assays$RNA
  
  if (nrow(RNA) == length(newnames)) {
    if (length(RNA@counts)) RNA@counts@Dimnames[[1]] <- newnames
    if (length(RNA@data)) RNA@data@Dimnames[[1]] <- newnames
    if (length(RNA@scale.data)) rownames(RNA@scale.data) <- newnames
    if (length(RNA@meta.features)) {
      RNA@meta.features$prev.names <- rownames(RNA@meta.features)
      rownames(RNA@meta.features) <- newnames
    }
    if (length(RNA@var.features)){
      RNA@var.features <- rownames(RNA@meta.features[match(RNA@var.features, RNA@meta.features$prev.names), ])
    }
  } else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"}
  obj@assays$RNA <- RNA
  return(obj)
}

add.sro.to.archr <- function(archr.obj, rna.md) { 
  rna.md <- rna.md %>% dplyr::select(-starts_with("ATAC"))
  rownames(rna.md) <- gsub("_", "#", rownames(rna.md))
  archr.md <- DataFrame(data.frame(archr.obj@cellColData) %>% dplyr::select(-starts_with("RNA")))
  archr.obj@cellColData <- cbind(archr.md, RNA = rna.md[rownames(archr.obj@cellColData), ])
  return(archr.obj)
}

add.archr.to.sro <- function(sro, archr.md) { 
  archr.md <- data.frame(archr.md) %>% dplyr::select(-starts_with("RNA"))
  # rownames(archr.md) <- gsub("#", "_", rownames(archr.md))
  rna.md <- sro@meta.data %>% dplyr::select(-starts_with("ATAC"))
  sro@meta.data <- cbind(rna.md, ATAC = archr.md[rownames(sro@meta.data), ])
  return(sro)
}

add.cso.to.sro <- function(sro, cs.md) { 
  cs.md <- data.frame(cs.md) %>% dplyr::select(-starts_with("RNA"))
  rownames(cs.md) <- gsub("#", "_", rownames(cs.md))
  rna.md <- sro@meta.data %>% dplyr::select(-starts_with("CS"))
  sro@meta.data <- cbind(rna.md, CS = cs.md[rownames(sro@meta.data), ])
  return(sro)
}

add.cso.to.archr <- function(archr.obj, cs.md) { 
  archr.md <- DataFrame(data.frame(archr.obj@cellColData) %>% dplyr::select(-starts_with("CS")))
  archr.obj@cellColData <- cbind(archr.md, CS = cs.md[rownames(archr.obj@cellColData), ])
  return(archr.obj)
}

get.named.vector <- function(md, colname){
  val.vec <- md[,colname]; names(val.vec) <- rownames(md)
  return(val.vec)
}

get.named.vector.cso <- function(cso, colname){
  val.vec <- cso@meta.data[, colname]; names(val.vec) <- rownames(cso@meta.data)
  return(val.vec)
}

get.named.vector.sro <- function(sro, colname){
  val.vec <- sro@meta.data[, colname]; names(val.vec) <- rownames(sro@meta.data)
  return(val.vec)
}






# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}


# palantir ####
gene.trend.heatmap <- function(
    expr, md,
    pal.res, trends,
    branch, genes, 
    cluster_rows, pal = NULL, row_split = NULL, 
    trend.name = "scaled\nchromVARscore\ntrend", expr.name = "scaled\nchromVAR scores"
){
  ca1 <- columnAnnotation(
    pseudotime = seq(min(pal.res$pseudotime), max(pal.res$pseudotime), length.out = ncol(trends)),
    col = list(pseudotime = colorRamp2(c(0, 0.5, 1), c("#fcff9c", "#ff4900", "#6e0000"))),
    show_annotation_name = F
  )
  cell.order <- order(
    round(pal.res$pseudotime, digits = 2),
    md[rownames(pal.res), "Clusters.res_0.6.merged"],
    decreasing = F
  )
  branch.cells <- rownames(pal.res)[cell.order]
  ca2 <- columnAnnotation(
    df = cbind(
      md[branch.cells, c("RNA.final.annot", "Clusters.res_0.6.merged")],
      pseudotime = pal.res[branch.cells, "pseudotime"]
    ),
    col = list(
      RNA.final.annot = pal$M_annotations, Clusters.res_0.6.merged = pal$Clusters.res_0.6.merged,
      pseudotime = colorRamp2(c(0, 0.5, 1), c("#fcff9c", "#ff4900", "#6e0000"))
    )
  )
  trends.scaled <- t(scale(t(trends[genes, ])))
  max.val <- ceiling(max(abs(c(trends.scaled))))
  col.ramp1 <- colorRamp2(
    breaks = round(seq(-max.val, max.val, length.out = 11)),
    colors = rev(RColorBrewer::brewer.pal(11, "PuOr"))
  )
  col.ramp2 <- colorRamp2(
    breaks = c(-10, -5, seq(-3, 3, 1), 5, 10),
    colors = c(
      "#173c68", "#173c68", "#1c5785", "#166db0", "#a6d1e3", "#ffffff",
      "#f09173", "#d46052", "#b32339", "#690822", "#690822"
    )
  )
  hm1 <- Heatmap(
    matrix = trends.scaled,
    name = trend.name, col = col.ramp1,
    show_column_names = F, cluster_columns = F, top_annotation = ca1,
    cluster_rows = cluster_rows, row_split = row_split, row_names_side = "right",
    row_names_gp = gpar(fontface = "italic", fontsize = 7),
    use_raster = T, width = unit(6, "inch")
  )
  hm2 <- Heatmap(
    matrix = t(scale(t(expr[genes, branch.cells]))),
    name = expr.name, col = col.ramp2,
    show_column_names = F, cluster_columns = F, top_annotation = ca2,
    cluster_rows = cluster_rows, row_split = row_split, row_names_side = "right",
    row_names_gp = gpar(fontface = "italic", fontsize = 7),
    use_raster = F, width = unit(6, "inch")
  )
  return(list(hm1, hm2))
  # draw(
  #   hm1 + hm2,
  #   ht_gap = unit(c(0.1, 0.1), "inch"),
  #   column_title = paste("branch ending at", branch)
  # )
}

gene.trend.heatmap.custom.row.labels <- function(
    expr, md,
    pal.res, trends,
    branch, genes, 
    cluster_rows, row.labels = genes, label.size = 18,
    pal = NULL, row_split = NULL,
    trend.name = "scaled\nchromVARscore\ntrend", expr.name = "scaled\nchromVAR scores"
){
  ca1 <- columnAnnotation(
    pseudotime = seq(min(pal.res$pseudotime), max(pal.res$pseudotime), length.out = ncol(trends)),
    col = list(pseudotime = colorRamp2(c(0, 0.5, 1), c("#fcff9c", "#ff4900", "#6e0000"))),
    show_annotation_name = F, show_legend=FALSE
  )
  cell.order <- order(
    round(pal.res$pseudotime, digits = 2),
    md[rownames(pal.res), "annotations"],
    decreasing = F
  )
  branch.cells <- rownames(pal.res)[cell.order]
  ca2 <- columnAnnotation(
    df = cbind(
      md[branch.cells, c("Sample", "annotations")],
      pseudotime = pal.res[branch.cells, "pseudotime"]
    ),
    col = c(
      pal[c("Sample", "annotations")],
      pseudotime = colorRamp2(c(0, 0.5, 1), c("#fcff9c", "#ff4900", "#6e0000"))
    ), show_legend=FALSE
  )
  trends.scaled <- t(scale(t(trends[genes, ])))
  max.val <- ceiling(max(abs(c(trends.scaled))))
  col.ramp1 <- colorRamp2(
    breaks = round(seq(-max.val, max.val, length.out = 11)),
    colors = rev(RColorBrewer::brewer.pal(11, "PuOr"))
  )
  col.ramp2 <- colorRamp2(
    breaks = c(-10, -5, seq(-3, 3, 1), 5, 10),
    colors = c(
      "#173c68", "#173c68", "#1c5785", "#166db0", "#a6d1e3", "#ffffff",
      "#f09173", "#d46052", "#b32339", "#690822", "#690822"
    )
  )
  hm1 <- Heatmap(
    matrix = trends.scaled,
    name = trend.name, col = col.ramp1,
    show_column_names = F, cluster_columns = F, top_annotation = ca1,
    cluster_rows = cluster_rows, row_split = row_split, row_names_side = "right",
    row_names_gp = gpar(fontface = "italic", fontsize = 7), 
    use_raster = T, width = unit(6, "inch")
  )
  ha = rowAnnotation(foo = anno_mark(at = match(row.labels, genes), 
                                     labels = row.labels, labels_gp = gpar(fontsize=label.size)))
  
  hm2 <- Heatmap(
    matrix = t(scale(t(expr[genes, branch.cells]))),
    name = expr.name, col = col.ramp2,
    show_column_names = F, cluster_columns = F, top_annotation = ca2,
    cluster_rows = cluster_rows, row_split = row_split, 
    show_row_names = F, row_names_side = "right",
    row_names_gp = gpar(fontface = "italic", fontsize = 7), 
    use_raster = F, width = unit(6, "inch"), right_annotation = ha,
  )
  return(list(hm1, hm2))
  # draw(
  #   hm1 + hm2,
  #   ht_gap = unit(c(0.1, 0.1), "inch"),
  #   column_title = paste("branch ending at", branch),
  #   annotation_legend_side = "bottom", show_heatmap_legend = FALSE
  # )
}


# Peaks ####
# annot <- read.csv(file = "chipseeker/AireWT/annotation.tsv", sep = "\t") # 1-indexed
# annot.gr <- create.GR(annot) # 1-indexed
# annot.summit <- create.GR.summit(annot) # 1-indexed
# dir.create("bed-files/summits")
# export.bed(annot.summit, "bed-files/summits/all-summits.bed") # this function assumes annot object is 1-indexed

# read narrowPeak from idr or macs2
# start is 0-indexed
read.narrowPeak <- function(fn, format = 'idr'){
  if (format == 'idr'){
    col.names = c("seqnames", "start", 'end', 'name', 'score', 'strand', 'signalValue',
                  'pval', 'qval', 'summit', 'localIDR', 'globalIDR',
                  'rep1_start', 'rep1_end', 'rep1_signalValue', 'rep1_summit',
                  'rep2_start', 'rep2_end', 'rep2_signalValue', 'rep2_summit')
  } else if (format == 'narrowPeak') {
    col.names = c("seqnames", "start", 'end', 'name', 'score', 'strand', 'signalValue',
                  'pval', 'qval', 'summit')
  } else if (format == "bed"){
    col.names = c("seqnames", "start", 'end', 'name', 'score', 'strand')
  }
  peaks <- read.table(fn, header = F, col.names = col.names)
  return(peaks)
}

# annot is chipseeker output
create.GR <- function(annot, input.is.0.indexed = F){
  mapper <- setNames(c("-", "+", "*", "*"), c("-", "+", ".", "*"))
  if (input.is.0.indexed){
    annot$start <- annot$start + 1
  }
  gr <- GRanges(seqnames = annot$seqnames, 
                IRanges(annot$start, annot$end),
                strand = mapper[annot$strand],
                annot[,!(colnames(annot) %in% 
                           c("seqnames", "ranges", "strand", "seqlevels", 
                             "seqlengths", "isCircular", "start", "end", 
                             "width", "element"))])
  return(gr)
}

# annot is chipseeker output
create.GR.summit <- function(annot, width = 200, input.is.0.indexed = F){
  mapper <- setNames(c("-", "+", "*", "*"), c("-", "+", ".", "*"))
  width.2 <- floor(width / 2)
  if (input.is.0.indexed){
    annot$start <- annot$start + 1
  }
  annot$summit.pos <- annot$start + annot$summit
  gr <- GRanges(seqnames = annot$seqnames, 
                IRanges(start = annot$summit.pos - width.2, 
                        end = annot$summit.pos + width.2 - 1),
                strand = mapper[annot$strand], 
                annot[,!(colnames(annot) %in% 
                           c("seqnames", "ranges", "strand", "seqlevels", 
                             "seqlengths", "isCircular", "start", "end", 
                             "width", "element"))])
  return(gr)
}

loadGMT <- function(gmtFile) {
  lines <- readLines(gmtFile)  # Read lines from the GMT file
  geneSetList <- lapply(lines, function(line) {
    parts <- strsplit(line, "\t")[[1]]  # Split each line by tab
    df <- data.frame(name = parts[1], description = parts[2], genes = parts[-c(1,2)])
    return(df)
  })
  geneSetDF <- do.call(rbind, geneSetList)
  return(geneSetDF)
}
