#Author: Supat Thongjuea, MRC Molecular Haematology Unit, Weatherall Institute of Molecular Medicine, University of Oxford, UK 
#Contact email : supat.thongjuea@imm.ox.ac.uk 
#Maintainer: Supat Thongjuea
#Title: Single-Cell Analyses Reveal Aberrant Pathways for Megakaryocyte-Biased Hematopoiesis in Myelofibrosis and Identify Mutant Clone-Specific Targets
#Journal : Molecular Cell
#Year : 2020
###############################################################################
library(SingleCellExperiment)
library(Matrix)
library(matrixStats)
library(ggplot2)
library(gridExtra)
library(igraph)
library(ComplexHeatmap)
library(circlize)
library(DT)


source("R/SingleCellClasses.R")
source("R/SingleCellGenerics.R")
source("R/SingleCellPlots.R")

shinyUI(fluidPage(
  h1(span("SingCellaR : Visualisation tool for scRNA-seq data", style = "font-weight: 50"), 
     style = "font-family: 'Arial';
     color: yellow; text-align: center;font-size: 400%;
     background-image: url('bg.png');
     padding: 60px"),
     br(),
				sidebarLayout(
						sidebarPanel(width = 3,
								htmlOutput('list_of_projects'),
								br(),
								htmlOutput('list_of_genes')
						),
						mainPanel(
								tabsetPanel(
										# Plot
										tabPanel("Plot By Gene",
												plotOutput("plot_by_gene",height = "650px")
										),
										# PCA analysis 
										tabPanel("Plot By Cluster",
												#scatterplotThreeOutput("plot_by_cluster",height = "650px")
												plotOutput("plot_by_cluster",height = "650px")

										),
										# show marker genes
										#tabPanel("Heatmap of marker genes per cluster",
										#         #scatterplotThreeOutput("plot_by_cluster",height = "650px")
										#         numericInput('top_genes', 'Select top genes', 10,
										#                      min = 5, max = 50),
										#         actionButton("submit_top_genes", "Plot"),
										#         plotOutput("heatmap_marker_genes",height = "800px")
										#         
										#),
										tabPanel("Table of marker genes per cluster",
										         htmlOutput('list_of_clusters'),
										         dataTableOutput('MarkerGenesTable'),HTML('<hr>')
										         
										),
										id="tabs1"
								)
						
						)

				)
		)
)