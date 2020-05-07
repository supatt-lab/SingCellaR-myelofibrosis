#Author: Supat Thongjuea, MRC Molecular Haematology Unit, Weatherall Institute of Molecular Medicine, University of Oxford, UK 
#Contact email : supat.thongjuea@imm.ox.ac.uk 
#Maintainer: Supat Thongjuea
#Title: Single-Cell Analyses Reveal Aberrant Pathways for Megakaryocyte-Biased Hematopoiesis in Myelofibrosis and Identify Mutant Clone-Specific Targets
#Journal : Molecular Cell
#Year : 2020
###############################################################################
setGeneric(
		name="get_dir_path_10x_matrix",
		def=function(object){
			standardGeneric("get_dir_path_10x_matrix")
		})
setMethod("get_dir_path_10x_matrix",
		signature(object = "SingCellaR"),
		function (object){
			object@dir_path_10x_matrix
		})
########
########
setGeneric(
		name="get_sample_uniq_id",
		def=function(object){
			standardGeneric("get_sample_uniq_id")
		})
setMethod("get_sample_uniq_id",
		signature(object = "SingCellaR"),
		function (object){
			object@sample_uniq_id
		})
########
########
setGeneric(
		name="get_umi_count",
		def=function(object){
			standardGeneric("get_umi_count")
		})
setMethod("get_umi_count",
		signature(object = "SingCellaR"),
		function (object){
			assay(object,"counts")
		})

########
########
setGeneric(
		name="get_normalized_umi",
		def=function(object){
			standardGeneric("get_normalized_umi")
		})
setMethod("get_normalized_umi",
		signature(object = "SingCellaR"),
		function (object){
		  assay(object,"normalized.umi")
		})
########
########
setGeneric(
  name="get_regressout_log_data",
  def=function(object){
    standardGeneric("get_regressout_log_data")
  })
setMethod("get_regressout_log_data",
          signature(object = "SingCellaR"),
          function (object){
            object@regressout_matrix
          })
########
########
setGeneric(
		name="get_cells_annotation",
		def=function(object){
			standardGeneric("get_cells_annotation")
		})
setGeneric(
		name="get_cells_annotation<-",
		def=function(object,value){
			standardGeneric("get_cells_annotation<-")
		})
setMethod("get_cells_annotation",
		signature(object = "SingCellaR"),
		function (object){
			object@meta.data
		})
setReplaceMethod(
		f="get_cells_annotation",
		signature="SingCellaR",
		definition=function(object,value){
			initialize(object,meta.data=value)
		})
########
########
setGeneric(
		name="get_genes_metadata",
		def=function(object){
			standardGeneric("get_genes_metadata")
		})
setGeneric(
		name="get_genes_metadata<-",
		def=function(object,value){
			standardGeneric("get_genes_metadata<-")
		})
setMethod("get_genes_metadata",
		signature(object = "SingCellaR"),
		function (object){
			object@genes.info
		})
setReplaceMethod(
		f="get_genes_metadata",
		signature="SingCellaR",
		definition=function(object,value){
			initialize(object,genes.info=value)
		})
########
########
setGeneric(
		name="get_pca.result",
		def=function(object){
			standardGeneric("get_pca.result")
		})
setGeneric(
		name="get_pca.result<-",
		def=function(object,value){
			standardGeneric("get_pca.result<-")
		})
setMethod("get_pca.result",
		signature(object = "SingCellaR"),
		function (object){
			object@pca.result
		})
setReplaceMethod(
		f="get_pca.result",
		signature="SingCellaR",
		definition=function(object,value){
			initialize(object,pca.result=value)
		})
########
########
setGeneric(
  name="get_nnmf.result",
  def=function(object){
    standardGeneric("get_nnmf.result")
  })
setGeneric(
  name="get_nnmf.result<-",
  def=function(object,value){
    standardGeneric("get_nnmf.result<-")
  })
setMethod("get_nnmf.result",
          signature(object = "SingCellaR"),
          function (object){
            object@nnmf.result
          })
setReplaceMethod(
  f="get_nnmf.result",
  signature="SingCellaR",
  definition=function(object,value){
    initialize(object,nnmf.result=value)
  })
########
########
setGeneric(
		name="get_tsne.result",
		def=function(object){
			standardGeneric("get_tsne.result")
		})
setGeneric(
		name="get_tsne.result<-",
		def=function(object,value){
			standardGeneric("get_tsne.result<-")
		})
setMethod("get_tsne.result",
		signature(object = "SingCellaR"),
		function (object){
			object@tsne.result
		})
setReplaceMethod(
		f="get_tsne.result",
		signature="SingCellaR",
		definition=function(object,value){
			initialize(object,tsne.result=value)
		})

########
########
setGeneric(
		name="get_dfm.result",
		def=function(object){
			standardGeneric("get_dfm.result")
		})
setGeneric(
		name="get_dfm.result<-",
		def=function(object,value){
			standardGeneric("get_dfm.result<-")
		})
setMethod("get_dfm.result",
		signature(object = "SingCellaR"),
		function (object){
			object@diffusionmap.result
		})
setReplaceMethod(
		f="get_dfm.result",
		signature="SingCellaR",
		definition=function(object,value){
			initialize(object,diffusionmap.result=value)
		})

########
########
setGeneric(
		name="get_umap.result",
		def=function(object){
			standardGeneric("get_umap.result")
		})
setGeneric(
		name="get_umap.result<-",
		def=function(object,value){
			standardGeneric("get_umap.result<-")
		})
setMethod("get_umap.result",
		signature(object = "SingCellaR"),
		function (object){
			object@umap.result
		})
setReplaceMethod(
		f="get_umap.result",
		signature="SingCellaR",
		definition=function(object,value){
			initialize(object,umap.result=value)
		})
########
########
setGeneric(
		name="get_tsne.result.high.contour",
		def=function(object){
			standardGeneric("get_tsne.result.high.contour")
		})
setGeneric(
		name="get_tsne.result.high.contour<-",
		def=function(object,value){
			standardGeneric("get_tsne.result.high.contour<-")
		})
setMethod("get_tsne.result.high.contour",
		signature(object = "SingCellaR"),
		function (object){
			object@tsne.result.high.contour
		})
setReplaceMethod(
		f="get_tsne.result.high.contour",
		signature="SingCellaR",
		definition=function(object,value){
			initialize(object,tsne.result.high.contour=value)
		})
########
########
setGeneric(
		name="get_knn_graph.graph",
		def=function(object){
			standardGeneric("get_knn_graph.graph")
		})
setGeneric(
		name="get_knn_graph.graph<-",
		def=function(object,value){
			standardGeneric("get_knn_graph.graph<-")
		})
setMethod("get_knn_graph.graph",
		signature(object = "SingCellaR"),
		function (object){
			object@knn_graph.graph
		})
setReplaceMethod(
		f="get_knn_graph.graph",
		signature="SingCellaR",
		definition=function(object,value){
			initialize(object,knn_graph.graph=value)
		})
########
########
setGeneric(
		name="get_knn_graph.layout",
		def=function(object){
			standardGeneric("get_knn_graph.layout")
		})
setGeneric(
		name="get_knn_graph.layout<-",
		def=function(object,value){
			standardGeneric("get_knn_graph.layout<-")
		})
setMethod("get_knn_graph.layout",
		signature(object = "SingCellaR"),
		function (object){
			object@knn_graph.layout
		})
setReplaceMethod(
		f="get_knn_graph.layout",
		signature="SingCellaR",
		definition=function(object,value){
			initialize(object,knn_graph.layout=value)
		})
########
########
setGeneric(
		name="get_knn_graph.kmeans.cluster",
		def=function(object){
			standardGeneric("get_knn_graph.kmeans.cluster")
		})
setGeneric(
		name="get_knn_graph.kmeans.cluster<-",
		def=function(object,value){
			standardGeneric("get_knn_graph.kmeans.cluster<-")
		})
setMethod("get_knn_graph.kmeans.cluster",
		signature(object = "SingCellaR"),
		function (object){
			object@knn_graph.kmeans.cluster
		})
setReplaceMethod(
		f="get_knn_graph.kmeans.cluster",
		signature="SingCellaR",
		definition=function(object,value){
			initialize(object,knn_graph.kmeans.cluster=value)
		})
########
########
setGeneric(
		name="get_clusters",
		def=function(object){
			standardGeneric("get_clusters")
		})
setGeneric(
		name="get_clusters<-",
		def=function(object,value){
			standardGeneric("get_clusters<-")
		})
setMethod("get_clusters",
		signature(object = "SingCellaR"),
		function (object){
			object@sc.clusters
		})
setReplaceMethod(
		f="get_clusters",
		signature="SingCellaR",
		definition=function(object,value){
			initialize(object,sc.clusters=value)
		})

########
setGeneric(
		name="get_marker_genes",
		def=function(object){
			standardGeneric("get_marker_genes")
		})
setGeneric(
		name="get_marker_genes<-",
		def=function(object,value){
			standardGeneric("get_marker_genes<-")
		})
setMethod("get_marker_genes",
		signature(object = "SingCellaR"),
		function (object){
			object@marker.genes
		})
setReplaceMethod(
		f="get_marker_genes",
		signature="SingCellaR",
		definition=function(object,value){
			initialize(object,marker.genes=value)
		})

########
setGeneric(
		name="get_igraph.graph",
		def=function(object){
			standardGeneric("get_igraph.graph")
		})
setGeneric(
		name="get_igraph.graph<-",
		def=function(object,value){
			standardGeneric("get_igraph.graph<-")
		})
setMethod("get_igraph.graph",
		signature(object = "SingCellaR"),
		function (object){
			object@igraph.graph
		})
setReplaceMethod(
		f="get_igraph.graph",
		signature="SingCellaR",
		definition=function(object,value){
			initialize(object,igraph.graph=value)
		})
########
########
setGeneric(
		name="get_fa2_graph.layout",
		def=function(object){
			standardGeneric("get_fa2_graph.layout")
		})
setGeneric(
		name="get_fa2_graph.layout<-",
		def=function(object,value){
			standardGeneric("get_fa2_graph.layout<-")
		})
setMethod("get_fa2_graph.layout",
		signature(object = "SingCellaR"),
		function (object){
			object@fa2_graph.layout
		})
setReplaceMethod(
		f="get_fa2_graph.layout",
		signature="SingCellaR",
		definition=function(object,value){
			initialize(object,fa2_graph.layout=value)
		})
########
