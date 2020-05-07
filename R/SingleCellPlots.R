#Author: Supat Thongjuea, MRC Molecular Haematology Unit, Weatherall Institute of Molecular Medicine, University of Oxford, UK 
#Contact email : supat.thongjuea@imm.ox.ac.uk 
#Maintainer: Supat Thongjuea
#Title: Single-Cell Analyses Reveal Aberrant Pathways for Megakaryocyte-Biased Hematopoiesis in Myelofibrosis and Identify Mutant Clone-Specific Targets
#Journal : Molecular Cell
#Year : 2020
###############################################################################
plot_umap_label_by_genes<-function(object,gene_list=c(),point.size=0.3,point.color1="blue",point.color2="red"){
	
	if(!is(object,"SingCellaR")){
		stop("Need to initialize the scRNAseq object")
	}
	if(length(gene_list)==0){
		stop("Required list of genes!")
	}
	#####
	res.umap<-get_umap.result(object)
	#####
	umi<-get_umi_count(object)
	umi.used<-umi[,as.character(res.umap$Cell)]
	#####
	my.dat<-res.umap[,c("Cell","UMAP1","UMAP2")]
	for(i in 1:length(gene_list)){
		genes.x<-gene_list[i]
		g.ind<-which(rownames(umi.used)==genes.x)
		my.sub.umi<-as.data.frame(log1p(umi.used[g.ind,]))
		colnames(my.sub.umi)<-genes.x
		#################################
		my.dat<-cbind(my.dat,my.sub.umi)
	}
	X.colnames<-colnames(my.dat)
	X.colnames<-gsub("-","_",X.colnames)
	X.colnames<-gsub("\\+","_",X.colnames)
	colnames(my.dat)<-X.colnames
	
	my.fg.names<-colnames(my.dat)
	my.fg.names<-my.fg.names[4:length(my.fg.names)]
	
	p.x <- list()
	for(j in 1:length(my.fg.names)){
		genes.x<-my.fg.names[j]
		genes.x<-gsub("-","_",genes.x)
		genes.x<-gsub("\\+","_",genes.x)
		p.x[[j]] <- ggplot(my.dat,aes(UMAP1,UMAP2)) + geom_point(aes_string(colour = genes.x),size=point.size) + scale_colour_gradientn(colours = c("gray87",point.color1,point.color2),values=c(0,0.01,1))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
	}
	do.call(grid.arrange,p.x)
}

plot_umap_label_by_clusters<-function(object,show_method=c("walktrap","louvain","kmeans","merged_walktrap","merged_louvain","merged_kmeans"),point.size=2,
		mark.clusters=TRUE,mark.font.size=10,mark.font.color="yellow",
		add.cluster.in.cells=FALSE,cluster.font.size=1,cluster.font.color="yellow"){
	
	if(!is(object,"SingCellaR")){
		stop("Need to initialize the SingCellaR object")
	}
	
	show_method <- match.arg(show_method)
	###
	res.umap<-get_umap.result(object)
	clusters.info<-get_clusters(object)
	
	if(show_method=="walktrap"){
		clusters.info<-clusters.info[,c("Cell","walktrap_cluster","walktrap_cluster_color")]
		colnames(clusters.info)<-c("Cell","cluster","cluster_color")
	}else if(show_method=="louvain"){
		clusters.info<-clusters.info[,c("Cell","louvain_cluster","louvain_cluster_color")]
		colnames(clusters.info)<-c("Cell","cluster","cluster_color")
	}else if(show_method=="kmeans"){
		clusters.info<-get_knn_graph.kmeans.cluster(object)
		colnames(clusters.info)<-c("Cell","cluster","cluster_color")
	}else if(show_method=="merged_walktrap"){
		clusters.info<-clusters.info[,c("Cell","merged_walktrap","merged_walktrap_color")]
		colnames(clusters.info)<-c("Cell","cluster","cluster_color")
	}else if(show_method=="merged_louvain"){
		clusters.info<-clusters.info[,c("Cell","merged_louvain","merged_louvain_color")]
		colnames(clusters.info)<-c("Cell","cluster","cluster_color")
	}else if(show_method=="merged_kmeans"){
		clusters.info<-get_knn_graph.kmeans.cluster(object)[,c("Cell","merged_kmeans","merged_kmeans_color")]
		colnames(clusters.info)<-c("Cell","cluster","cluster_color")
	}else{
		stop("Need community detection or clustering-method names!")
	}
	res.umap<-merge(res.umap,clusters.info)
	cluster.id<-gsub("cl","",res.umap$cluster)
	
	if(add.cluster.in.cells==TRUE){
		p<- qplot(UMAP1,UMAP2, data=res.umap,colour=cluster)+
				geom_point(aes(fill=cluster),colour="black",pch=21,size=point.size)+ 
				geom_text(aes(label = add.cluster.in.cells),size = cluster.font.size,colour=cluster.font.color)+
				theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
				theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
		
		if(mark.clusters==TRUE){
			edata <- res.umap[,c("UMAP1","UMAP2")]
			edata$cluster.id<-cluster.id
			colnames(edata) <- c('x', "y", "z")
			center <- aggregate(cbind(x,y) ~ z, data = edata, median)
			p<-p+geom_text(data=center, aes_string(x = "x", y = "y", label = "z"), size = mark.font.size, colour = mark.font.color)
		}
		
	}else{
		p<- qplot(UMAP1,UMAP2, data=res.umap,colour=cluster)+
				geom_point(aes(fill=cluster),colour="black",pch=21,size=point.size)+
				theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
				theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
		
		if(mark.clusters==TRUE){
			edata <- res.umap[,c("UMAP1","UMAP2")]
			edata$cluster.id<-cluster.id
			colnames(edata) <- c('x', "y", "z")
			center <- aggregate(cbind(x,y) ~ z, data = edata, median)
			p<-p+geom_text(data=center, aes_string(x = "x", y = "y", label = "z"), size = mark.font.size, colour = mark.font.color)
		}
	}
	p<-p+ggtitle(show_method)
	return(p)
}

plot_heatmap_for_marker_genes<-function(object,cluster.type=c("walktrap","louvain","kmeans","merged_louvain","merged_walktrap","merged_kmeans"),
		n.TopGenes=5,min.log2FC=0.5,min.expFraction=0.3,isClusterByRow=F,col.scaled.min = -2.5,col.scaled.max = 2.5,
		col.low="blue",col.mid = "black",col.high = "red",rowFont.size=6,split.line.col="white", split.line.type=1,split.line.lwd=1){
	
	if(!is(object,"SingCellaR")){
		stop("Need to initialize the SingCellaR object")
	}
	
	clusters.info<-get_clusters(object)
	cluster.type <- match.arg(cluster.type)
	####################################
	if(cluster.type=="walktrap"){
		clusters.info<-clusters.info[,c("Cell","walktrap_cluster","walktrap_cluster_color")]
	}else if(cluster.type=="louvain"){
		clusters.info<-clusters.info[,c("Cell","louvain_cluster","louvain_cluster_color")]
	}else if(cluster.type=="kmeans"){
		clusters.info<-get_knn_graph.kmeans.cluster(object)
	}else if(cluster.type=="merged_walktrap"){
		clusters.info<-clusters.info[,c("Cell","merged_walktrap","merged_walktrap_color")]
	}else if(cluster.type=="merged_louvain"){
		clusters.info<-clusters.info[,c("Cell","merged_louvain","merged_louvain_color")]
	}else if(cluster.type=="merged_kmeans"){
		clusters.info<-get_knn_graph.kmeans.cluster(object)[,c("Cell","merged_kmeans","merged_kmeans_color")]
	}else{
		stop("Need a clustering-method name!")
	}
	clusters.info$id<-0
	clusters.info$id<-as.numeric(gsub("cl","",clusters.info[,2]))
	
	clusters.ids<-table(clusters.info[,4])
	TopGenes.info<-data.frame()
	for(k in 1:length(clusters.ids)){
		my.cluster<-paste("cl",names(clusters.ids)[k],sep="")
		##################
		sig.genes<-getMarkerGenesInfo(object,cluster.type = cluster.type,cluster_id = my.cluster)
		if(nrow(sig.genes) >0){
			sig.genes<-subset(sig.genes,fishers.pval< 0.05 
							& wilcoxon.pval < 0.05 
							& ExpFractionA > ExpFractionB 
							& ExpFractionA > min.expFraction 
							& log2FC >= min.log2FC)
		}
		if(nrow(sig.genes) > 0){
			sig.genes<-sig.genes[order(sig.genes$FoldChange,decreasing = T),]	
			top.genes<-head(sig.genes,n=n.TopGenes)
			top.genes$cluster<-my.cluster
			TopGenes.info<-rbind(TopGenes.info,top.genes)
		}else{
			message<-paste("There is no gene passed the cut-off for:",my.cluster,"!.",sept="")
			print(message)
		}
	}
	topGenes<-unique(as.character(TopGenes.info$Gene))
	umi.dat<-get_normalized_umi(object)
	###############
	clusters.info.ranked<-clusters.info[order(as.numeric(clusters.info$id)),]
	umi.dat<-umi.dat[as.character(topGenes),as.character(clusters.info.ranked$Cell)]
	umi.dat<-log1p(umi.dat)
	###############
	mat_scaled = t(apply(umi.dat, 1, scale))
	df<-data.frame(cluster = as.character(clusters.info.ranked[,2]))
	###get color###
	clusters.info.part<-clusters.info.ranked[,c(2:3)]
	clusters.info.part<-unique(clusters.info.part)
	my.cols<-clusters.info.part[,2]
	names(my.cols)<-clusters.info.part[,1]
	cluster.cols<-list(cluster=my.cols)
	###############
	cluster.annotation = HeatmapAnnotation(df = df,col=cluster.cols)
	###############
	ht1<-Heatmap(mat_scaled, name ="expression",col = colorRamp2(c(col.scaled.min, 0, col.scaled.max), c(col.low,col.mid,col.high)), 
			cluster_rows = isClusterByRow, cluster_columns = F,show_row_dend = F,clustering_distance_rows = "euclidean",show_row_names = T,show_column_names = F,
			top_annotation = cluster.annotation,row_names_gp = gpar(fontsize = rowFont.size))
	
	draw(ht1)
	loc<-cumsum(as.numeric(clusters.ids))
	decorate_heatmap_body("expression", {
				for(l in 1:length(loc)){
					i = loc[l]
					x = i/ncol(mat_scaled)
					grid.lines(c(x, x), c(0, 1), gp = gpar(col = split.line.col, lty = split.line.type,lwd=split.line.lwd))
				}
			})
}

plot_forceDirectedGraph_label_by_genes<-function(object,gene_list=c(),vertex.size=0.5,vertex.color1="blue",
		vertex.color2="red",edge.size=0.2,edge.color="gray",vertex.base.color="gray35",background.color="white"){
	
	if(!is(object,"SingCellaR")){
		stop("Need to initialize the SingCellaR object")
	}
	if(length(gene_list)==0){
		stop("Required list of genes!")
	}
	#####
	fa2.layout<-as.data.frame(get_fa2_graph.layout(object))
	colnames(fa2.layout) <- c("X1","X2")
	#####
	umi<-get_umi_count(object)
	umi.used<-umi[,as.character(rownames(fa2.layout))]
	#####
	rownames(fa2.layout) <- 1:nrow(fa2.layout)
	#####
	my.dat<-fa2.layout[,(1:2)]
	
	for(i in 1:length(gene_list)){
		genes.x<-gene_list[i]
		g.ind<-which(rownames(umi.used)==genes.x)
		my.sub.umi<-as.data.frame(log1p(umi.used[g.ind,]))
		colnames(my.sub.umi)<-genes.x
		#################################
		my.dat<-cbind(my.dat,my.sub.umi)
	}
	
	X.colnames<-colnames(my.dat)
	X.colnames<-gsub("-","_",X.colnames)
	X.colnames<-gsub("\\+","_",X.colnames)
	colnames(my.dat)<-X.colnames
	
	my.fg.names<-colnames(my.dat)
	my.fg.names<-my.fg.names[3:length(my.fg.names)]
	#########
	edgelist <- get.edgelist(get_igraph.graph(object))
	edges <- data.frame(fa2.layout[edgelist[,1],], fa2.layout[edgelist[,2],])
	colnames(edges) <- c("X1","Y1","X2","Y2")
	
	p.x <- list()
	for(j in 1:length(my.fg.names)){
		genes.x<-my.fg.names[j]
		genes.x<-gsub("-","_",genes.x)
		genes.x<-gsub("\\+","_",genes.x)
		
		p.x[[j]] <- ggplot(my.dat,aes(X1,X2)) + 
				geom_segment(aes(x=X1, y=Y1, xend = X2, yend = Y2), data=edges, size = edge.size, colour=edge.color)+
				geom_point(aes_string(colour = genes.x),size=vertex.size) + 
				scale_colour_gradientn(colours = c(vertex.base.color,vertex.color1,vertex.color2),values=c(0,0.01,1))+
				theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
				theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y=element_blank(),axis.title.x=element_blank())+
				theme(plot.background = element_rect(fill = background.color))
		
	}
	do.call(grid.arrange,p.x)
}

plot_forceDirectedGraph_label_by_clusters<-function(object,show_method=c("walktrap","louvain","kmeans","merged_walktrap","merged_louvain","merged_kmeans"),vertex.size=2,
		mark.clusters=TRUE,mark.font.size=10,mark.font.color="yellow",
		add.cluster.in.cells=FALSE,cluster.font.size=1,cluster.font.color="yellow",
		edge.size=0.2,edge.color="gray",background.color="white"){
	
	if(!is(object,"SingCellaR")){
		stop("Need to initialize the SingCellaR object")
	}
	show_method <- match.arg(show_method)
	###
	fa2.layout<-as.data.frame(get_fa2_graph.layout(object))
	clusters.info<-get_clusters(object)
	
	if(show_method=="walktrap"){
		clusters.info<-clusters.info[,c("Cell","walktrap_cluster","walktrap_cluster_color")]
		colnames(clusters.info)<-c("Cell","cluster","cluster_color")
	}else if(show_method=="louvain"){
		clusters.info<-clusters.info[,c("Cell","louvain_cluster","louvain_cluster_color")]
		colnames(clusters.info)<-c("Cell","cluster","cluster_color")
	}else if(show_method=="kmeans"){
		clusters.info<-get_knn_graph.kmeans.cluster(object)
		colnames(clusters.info)<-c("Cell","cluster","cluster_color")
	}else if(show_method=="merged_walktrap"){
		clusters.info<-clusters.info[,c("Cell","merged_walktrap","merged_walktrap_color")]
		colnames(clusters.info)<-c("Cell","cluster","cluster_color")
	}else if(show_method=="merged_louvain"){
		clusters.info<-clusters.info[,c("Cell","merged_louvain","merged_louvain_color")]
		colnames(clusters.info)<-c("Cell","cluster","cluster_color")
	}else if(show_method=="merged_kmeans"){
		clusters.info<-get_knn_graph.kmeans.cluster(object)[,c("Cell","merged_kmeans","merged_kmeans_color")]
		colnames(clusters.info)<-c("Cell","cluster","cluster_color")
	}else{
		stop("Need community detection or clustering-method names!")
	}
	colnames(fa2.layout)<-c("X1","X2")
	fa2.layout$Cell<-rownames(fa2.layout)
	
	res.fa2<-merge(fa2.layout,clusters.info)
	cluster.id<-gsub("cl","",res.fa2$cluster)
	####################################
	rownames(fa2.layout) <- 1:nrow(fa2.layout)
	fa2.layout<-fa2.layout[,c(1:2)]
	####################################
	edgelist <- get.edgelist(get_igraph.graph(object))
	edges <- data.frame(fa2.layout[edgelist[,1],], fa2.layout[edgelist[,2],])
	colnames(edges) <- c("X1","Y1","X2","Y2")
	####################################
	if(add.cluster.in.cells==TRUE){
		p<- qplot(X1,X2, data=res.fa2,colour=cluster)+
				geom_segment(aes(x=X1, y=Y1, xend = X2, yend = Y2), data=edges, size = edge.size, colour=edge.color)+
				geom_point(aes(fill=cluster),colour="black",pch=21,size=vertex.size)+ 
				geom_text(aes(label = add.cluster.in.cells),size = cluster.font.size,colour=cluster.font.color)+
				theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
				theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank())+
				theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y=element_blank())+
				theme(plot.background = element_rect(fill = background.color))
		
		if(mark.clusters==TRUE){
			edata <- res.fa2[,c("X1","X2")]
			edata$cluster.id<-cluster.id
			colnames(edata) <- c('x', "y", "z")
			center <- aggregate(cbind(x,y) ~ z, data = edata, median)
			p<-p+geom_text(data=center, aes_string(x = "x", y = "y", label = "z"), size = mark.font.size, colour = mark.font.color)
		}
		
	}else{
		p<- qplot(X1,X2, data=res.fa2,colour=cluster)+
				geom_segment(aes(x=X1, y=Y1, xend = X2, yend = Y2), data=edges, size = edge.size, colour=edge.color)+
				geom_point(aes(fill=cluster),colour="black",pch=21,size=vertex.size)+
				theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
				theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank())+
				theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y=element_blank())+
				theme(plot.background = element_rect(fill = background.color))
		
		if(mark.clusters==TRUE){
			edata <- res.fa2[,c("X1","X2")]
			edata$cluster.id<-cluster.id
			colnames(edata) <- c('x', "y", "z")
			center <- aggregate(cbind(x,y) ~ z, data = edata, median)
			p<-p+geom_text(data=center, aes_string(x = "x", y = "y", label = "z"), size = mark.font.size, colour = mark.font.color)
		}
	}
	p<-p+ggtitle(show_method)
	return(p)
}

getMarkerGenesInfo <- function(object,cluster.type=c("louvain","walktrap","kmeans","merged_louvain","merged_walktrap","merged_kmeans"),cluster_id="cl1"){
	objName <- deparse(substitute(object))
	if(!is(object,"SingCellaR")){
		stop("Need to initialize the SingCellaR object")
	}
	cluster.type <- match.arg(cluster.type)
	
	if(cluster.type=="walktrap"){
		markers_genes<-get_marker_genes(object)
		markers_genes<-markers_genes$walktrap
		cluster.genes<-as.data.frame(markers_genes[[cluster_id]])	
		return(cluster.genes)
	}else if(cluster.type=="louvain"){
		markers_genes<-get_marker_genes(object)
		markers_genes<-markers_genes$louvain
		cluster.genes<-as.data.frame(markers_genes[[cluster_id]])
		return(cluster.genes)
	}else if(cluster.type=="kmeans"){
		markers_genes<-get_marker_genes(object)
		markers_genes<-markers_genes$kmeans
		cluster.genes<-as.data.frame(markers_genes[[cluster_id]])
		return(cluster.genes)
	}else if(cluster.type=="merged_louvain"){
		markers_genes<-get_marker_genes(object)
		markers_genes<-markers_genes$merged_louvain
		cluster.genes<-as.data.frame(markers_genes[[cluster_id]])
		return(cluster.genes)
	}else if(cluster.type=="merged_walktrap"){
		markers_genes<-get_marker_genes(object)
		markers_genes<-markers_genes$merged_walktrap
		cluster.genes<-as.data.frame(markers_genes[[cluster_id]])
		return(cluster.genes)
	}else if(cluster.type=="merged_kmeans"){
		markers_genes<-get_marker_genes(object)
		markers_genes<-markers_genes$merged_kmeans
		cluster.genes<-as.data.frame(markers_genes[[cluster_id]])
		return(cluster.genes)
	}else{
		stop("Need a clustering-method name!")
	}	
}