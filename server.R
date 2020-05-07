#Author: Supat Thongjuea, MRC Molecular Haematology Unit, Weatherall Institute of Molecular Medicine, University of Oxford, UK 
#Contact email : supat.thongjuea@imm.ox.ac.uk 
#Maintainer: Supat Thongjuea
#Title: Single-Cell Analyses Reveal Aberrant Pathways for Megakaryocyte-Biased Hematopoiesis in Myelofibrosis and Identify Mutant Clone-Specific Targets
#Journal : Molecular Cell
#Year : 2020
###############################################################################
shinyServer(function(input, output){
			##***Interface for loading projects*****
			output$list_of_projects <- renderUI({
						project.files<-list.files("projects")
						
						if(length(project.files) > 0){
							fluidRow(
									column(6,selectInput("My_project", h5("Choose SingCellaR's object:"),choices=as.character(project.files)),
							        actionButton("selected_project","Load")
									)
							)
						} 
					})
			gene_list <- eventReactive(input$selected_project, {
						withProgress(message = 'Loading genes ...', value = 0.8, {
						Sys.sleep(0.25)
						genes_info<-get_genes_metadata(project.obj())
						selected.genes<-subset(genes_info,IsExpress==TRUE)
						return(rownames(selected.genes))	
						Sys.sleep(0.25)
					})
			})
			output$list_of_genes <- renderUI({
						
						genes<-gene_list()
						if(length(genes) > 0){
							fluidPage(
									fluidRow(
											radioButtons("embedding.method", "Select Method:",
											                    c("UMAP" = "umap","Force-directed Graph"="fa2")),
											br()),
									    column(8,selectInput("My_gene",h5("Choose gene:"),genes,selectize=TRUE),
											actionButton("selected_gene","Visualize")
										
								   )		
							)
						}
					})
			project.obj <- eventReactive(input$selected_project, {
						object<-local(get(load(paste("projects/",input$My_project,sep=""))))
						return(object)	
					})
			####################Color Selection########################
			#select_color <- eventReactive(input$my.g.col, {
			#  return(input$my.g.col)	
			#})
			##################plot by gene#############################
			f_plot_by_gene <- eventReactive(input$selected_gene, {
			      
			      withProgress(message = 'Plotting...', value = 0.8, {
			        my_obj<-project.obj()
							Sys.sleep(0.25)
							if(input$embedding.method=="umap"){
			          g<-plot_umap_label_by_genes(my_obj,gene_list=c(input$My_gene))
							  Sys.sleep(0.25)
							}else if(input$embedding.method=="fa2"){
							  g<-plot_forceDirectedGraph_label_by_genes(my_obj,gene_list=c(input$My_gene))
							  Sys.sleep(0.25)
							}
						})
			  return(g)
			})

			output$plot_by_gene<-renderPlot({
			 f_plot_by_gene()
			})
			###########Plot by clusters############################
			f_plot_by_clusters <- eventReactive(input$embedding.method, {
			  
			  withProgress(message = 'Plotting...', value = 0.8, {
			    my_obj<-project.obj()
			    Sys.sleep(0.25)
			    if(input$embedding.method=="umap"){
			      g<-plot_umap_label_by_clusters(my_obj,show_method=c("louvain"))
			      Sys.sleep(0.25)
				  
			    }else if(input$embedding.method=="fa2"){
			      g<-plot_forceDirectedGraph_label_by_clusters(my_obj,show_method=c("louvain"))
			      Sys.sleep(0.25)
			    }
			  })
			  return(g)
			})
			
			output$plot_by_cluster<-renderPlot({
			  f_plot_by_clusters()
			})
			##############Heatmap#####################################
			heatmap_plot <- eventReactive(input$submit_top_genes,{
			  withProgress(message = 'Plotting...', value = 0.8, {
			    my_obj<-project.obj()
			    Sys.sleep(0.25)
			    g<-plot_heatmap_for_marker_genes(my_obj,cluster.type=c("louvain"),n.TopGenes=input$top_genes)
			  })
			  return(g)
			})
			
			output$heatmap_marker_genes<-renderPlot({
			  heatmap_plot()
			})
			##############Table#########################################
			get_cluster_ids <- eventReactive(input$selected_project, {
			  withProgress(message = 'Loading genes ...', value = 0.8, {
			    my_obj<-project.obj()
			    Sys.sleep(0.25)
			    clusters.ids<-names(my_obj@marker.genes$louvain)
			    return(clusters.ids)	
			    Sys.sleep(0.25)
			  })
			})
			output$list_of_clusters <- renderUI({
			  
			  cluster_ids<-get_cluster_ids()
			  if(length(cluster_ids) > 0){
			    fluidPage(
			      fluidRow(
			        br()),
			      column(8,selectInput("My_cluster",h5("Choose cluster:"),cluster_ids,selectize=TRUE),
			             actionButton("selected_cluster","Show Table")
			             
			      )		
			    )
			  }
			})
			get_cluster_info <- eventReactive(input$selected_cluster, {
			  withProgress(message = 'Loading genes ...', value = 0.8, {
			    Sys.sleep(0.25)
			    clusters.info<-getMarkerGenesInfo(project.obj(),cluster.type=c("louvain"),cluster_id=input$My_cluster)
			    return(clusters.info)
			    #return(isolate(input$My_cluster))
			    Sys.sleep(0.25)
			  })
			})
			output$MarkerGenesTable <- renderDataTable({
			  get_cluster_info()
			}, options = list(orderClasses = TRUE))
})
