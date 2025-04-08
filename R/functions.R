#' @import Seurat
#' @import SingleR
#' @import edgeR
#' @import Matrix
NULL
#> NULL
#> 

.clusterReorder = function(cts) {
  d = Embeddings(cts,"umap")
  d = as.data.frame(d)
  d$cluster = cts$seurat_clusters
  cls.pos = aggregate(.~ cluster, FUN=mean, data = d)
  cls = cls.pos$cluster
  
  cls.dis = lapply(1:nrow(cls.pos), function(i) {
    r = t(cls.pos[,2:3]) - as.numeric(cls.pos[i,2:3])
    colnames(r) = cls.pos$cluster
    r = r^2
    r = colSums(r)
    return(sort(r))
  })
  names(cls.dis) = cls
  
  cls.ord = as.character(cls.pos$cluster[order(cls.pos$umap_1)][1])
  
  for (i in 1:(length(cls)-1)) {
    cls.left = cls.dis[[as.character(cls.ord[i])]]
    cls.left = cls.left[setdiff(cls, cls.ord)]
    cls.left = sort(cls.left)
    cls.ord = c(cls.ord, names(cls.left)[1])
  }
  
  cls.reorder = setNames(1:length(cls),cls.ord)
  cts$cluster = as.character(cls.reorder[as.character(cts$seurat_clusters)])
  cts$cluster = as.numeric(cts$cluster)
  return(cts)
}


.scProc <- function(sc, n = 2000, pc =30, res =0.5, normalize = TRUE, hvg = TRUE, pca = TRUE, umap = TRUE, cluster = TRUE) {
  DefaultAssay(sc) = "RNA"
  if (normalize) {
    sc <- Seurat::NormalizeData(sc, verbose = FALSE)
  } 
  if (hvg) {
    sc <- Seurat::FindVariableFeatures(sc, selection.method = "vst", nfeatures = n, verbose = FALSE)
  }
  if (pca) {
    sc <- Seurat::ScaleData(sc, verbose = FALSE)
    sc <- Seurat::RunPCA(sc, npcs = pc, verbose = FALSE)
  }
  if (umap) {
    sc <- Seurat::RunUMAP(sc, reduction = "pca", dims = 1:pc)
  }
  
  if (cluster) {
    sc <- Seurat::FindNeighbors(sc, dims=1:pc)
    sc <- Seurat::FindClusters(sc, resolution = res)
    sc <- .clusterReorder(cts = sc)
  }
  return(sc)
}

.getPseudobulk <- function(obj, cluster, keep.all.genes = FALSE, min.cell = 5, assay = "RNA") {
  
  require(edgeR)
  obj$cluster = paste0("Cluster_", obj@meta.data[, cluster])
  cls.no = table(obj$cluster)
  cls.no = cls.no[cls.no>=min.cell]
  obj = obj[, obj$cluster %in% names(cls.no)]
  psk = Seurat::AggregateExpression(object = obj, assays = assay, return.seurat = F, group.by = "cluster")
  psk = psk[[1]]
  colnames(psk) = names(cls.no)
  if (keep.all.genes) {
    psk = DGEList(counts = psk)
  } else {
    psk = DGEList(counts = psk[rowSums(psk)>0, ])
  }
  psk = calcNormFactors(psk)
  logtpm = log2(edgeR::cpm(psk)+1)
  
  return(list(dgelist = psk, logtpm = logtpm))
}

#' List available tissue types for reference
#'
#' @export
list.tissues <- function() {
  print(vbls$tissues)
}
#' Annotation for single cell data at cluster level
#'
#' Annotate single cell at cluster level using reference built from Tabula sapiens, human protein atlas and public datasets. The annotation algorithm used is SingleR.
#' @param query.srt query single cell data in Seurat object format
#' @param ref.tissue reference tissue type, default is "all". Other options include "immune" and specific tissue(s). To see list of available tissues, use list.tissues(). Multiple tissues separated by "|"
#' @param cluster user provided cluster assignment to use. it should be a column name in meta.data. if NULL, will use Seurat::FindCluster to find clusters
#' @param resolution resolution used to identify cell clusters, Default is 0.5 if cluster is not set.
#' @return Seurat object with "labels","pruned.labels" and "score" added to metadata
#' @export
cellAnno <- function(query.srt, ref.tissue = "all", cluster = NULL, resolution = 0.5) {
  if (ref.tissue == "all") {
    ref = all.ref
  } else if (ref.tissue == "immune") {
    ref = immgen.ref
  } else {
    ref = all.ref[,grepl(ref.tissue, colnames(all.ref))| colnames(all.ref) %in% c(vbls$stromal.cells, vbls$immune.cells)]
  }
  if (is.null(cluster)) {
    query.srt = .scProc(query.srt, res = resolution)
    cls = paste0("Cluster_", query.srt@meta.data[, "seurat_clusters"])
    psk = .getPseudobulk(obj = query.srt,  cluster = "seurat_clusters", min.cell = 2)
  } else {
    cls = paste0("Cluster_", query.srt@meta.data[, cluster])
    psk = .getPseudobulk(obj = query.srt,  cluster = cluster, min.cell = 2)
  }
  
  pred = SingleR::SingleR(test = psk$logtpm, ref = ref, assay.type.ref = "logcounts", 
                          assay.type.test = "logcounts", labels = colnames(ref), fine.tune = FALSE)
  cluster.map = data.frame(cluster = colnames(psk$logtpm),labels= pred$labels,
                           pruned.labels = pred$pruned.labels, score = apply(pred$scores, 1, max),
                           stringsAsFactors = F)
  cluster.meta = cluster.map[match(cls, cluster.map$cluster),]
  rownames(cluster.meta) = colnames(query.srt)
  query.srt = AddMetaData(query.srt, cluster.meta)
  return(query.srt)
}
#' Cell-of-Origin prediction for cancer bulk RNA-seq data
#'
#' Predict COO using reference built from Tabula sapiens, human protein atlas and public datasets. Right now it only support prediction for solid tumor type
#' @param logtpm.matrix query data matrix in log2-TPM
#' @param ref reference to use, default is "prostate".To see list of available tissues, use list.tissues(). Multiple tissues separated by "|". Set to "all" if not sure (for samples from metastatic sites, it is possible top hit is cell from the biopsy site). 
#' @return output of SingleR. "labels" is the predicted cell type. 
#' @export
cooPredict <- function(logtpm.matrix, ref = "all", exclude.tissues = NULL) {
  if (ref=="all") {
    refdat = all.ref[,1:184]
  } else {
    refdat = all.ref[,grepl(ref, colnames(all.ref))]
  }
  if (!is.null(exclude.tissues)) {
    refdat = refdat[, !grepl(exclude.tissues, colnames(refdat))]
  }
  pred = SingleR::SingleR(test = logtpm.matrix, ref = refdat, assay.type.ref = "logcounts",
                          assay.type.test = "logcounts",
                          labels = colnames(refdat))
  return(pred)
}

