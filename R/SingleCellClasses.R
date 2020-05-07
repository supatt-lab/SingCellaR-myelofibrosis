#Author: Supat Thongjuea, MRC Molecular Haematology Unit, Weatherall Institute of Molecular Medicine, University of Oxford, UK 
#Contact email : supat.thongjuea@imm.ox.ac.uk 
#Maintainer: Supat Thongjuea
#Title: Single-Cell Analyses Reveal Aberrant Pathways for Megakaryocyte-Biased Hematopoiesis in Myelofibrosis and Identify Mutant Clone-Specific Targets
#Journal : Molecular Cell
#Year : 2020
###############################################################################
setClass("SingCellaR",
		contains="SingleCellExperiment",
		slots = c(
				dir_path_10x_matrix = "character",
				sample_uniq_id="character",
				genes.info = "data.frame",
				meta.data  = "data.frame",
				regressout_matrix = "ANY",
				pca.result = "ANY",
				nnmf.result= "ANY",
				tsne.result = "ANY",
				umap.result ="ANY",
				knn_graph.graph = "ANY",
				knn_graph.layout = "ANY",
				knn_graph.kmeans.cluster = "ANY",
				igraph.graph="ANY",
				fa2_graph.layout = "ANY",
				sc.clusters="ANY",
				marker.genes="ANY"
		)
)

setClass("SingCellaR_Int",
		contain="SingCellaR",
		slots=c(
				dir_path_SingCellR_object_files="character",
				SingCellR_object_files="character",
				Variable.genes="ANY",
				SingCellaR.individual.clusters="ANY",
				GSEA.aligned.clusters="ANY",
				Scanorama.integrative.matrix="ANY",
				Harmony.embeddings="ANY",
				SupervisedHarmony.embeddings="ANY",
				Liger.embeddings="ANY",
				Seurat.embeddings="ANY",
				Combat.embeddings="ANY"
		)
)


# show method for scRNAseq
#' @param object A scRNAseq object
#' @name show
#' @aliases show,SingCellaR-method
#' @docType methods
#' @rdname show-methods
#' 
setMethod(
		f ="show",
		signature ="SingCellaR",
		definition = function(object) {
			cat(
					"An object of class",
					class(object),
					"with a matrix of :",
					dim(object)[1],
					"genes across",
					dim(object)[2],
					"samples.\n"
			)
			invisible(x = NULL)
		}
)
