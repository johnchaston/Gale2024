# bc_plot <- pcoa_2var(
# 	folder_name = "Efemale/bray_curtis", 
# 	var1 = "time_point", 
# 	var2 = "genotype", 
# 	source_var1 = c(1,2,3), 
# 	var1_colors = c("red","blue","black"), 
# 	var2_shape = c(21,22,23,24), 
# 	title_name_add = "Bray Curtis", 
# 	legend_position = "none")

# 0.0.2 - allows use of variable mapper name
#rm(folder_name, var1, var2, source_var1, var1)
# 
# pco1 = 3
# pco2 = 1
# 

# 
# folder_name = paste0(file_path,"/unweighted_unifrac")
# shape_var = "experiment" # which variable is the classifier
# color_var = "experiment"
# fill_var = "experiment"
# circle_var = "experiment"
# plot_shapes <- c(15,15,15,15,15,15,15)
# plot_fills <- c("green","blue","black","red","cyan","yellow","magenta")
# plot_colors <- c("green","blue","black","red","cyan","yellow","magenta")
# plot_ellipses <- c(1,1)
# linetype_var <- "experiment"
# circle_color <- c("green","blue","black","red","cyan","yellow","magenta")
# i = "bray_curtis"
# title_name_add = title_name_add
# legend_position = legend_position
# pco1=1
# pco2=2
# legend_name = NULL
# legend_label = NULL

pcoa_manyvar <- function(folder_name = folder_name, mapper_file = mapper_file, title_name_add = title_name_add, legend_position = legend_position,pco1=1,pco2=2,shape_var = shape_var,color_var = color_var,fill_var = fill_var,circle_var = circle_var, linetype_var = linetype_var, plot_shapes = plot_shapes, plot_fills = plot_fills, plot_colors = plot_colors, plot_ellipses = plot_ellipses, circle_color = "black", legend_name = NULL, legend_label = NULL) {
  
  varnames = c(shape_var, color_var, fill_var, circle_var,linetype_var)
  varnames <- varnames[!duplicated(varnames)]
  	wfpc <- head(read.table(paste('core-metrics-results-',folder_name,'_pcoa_results/ordination.txt',sep=""), sep="\t", fill=T, skip = 9,blank.lines.skip = T, header=F),-2)
	pc2 <- head(read.table(paste('core-metrics-results-',folder_name,'_pcoa_results/ordination.txt', sep = ""), sep="\t", fill=T, skip = 3, header=F, colClasses = "character"),2)
	pc_values <- as.numeric(as.character(pc2[2,]))*100

		wfpc <- wfpc[,c(1,pco1+1, pco2+1)] 
	
	fut2 <- read.table(paste('core-metrics-results-',folder_name,'_distance_matrix/distance-matrix.tsv',sep=""), header=T, sep="\t", stringsAsFactors = T) %>% 
	  dplyr::select(X) %>% 
	  mutate(X=as.character(X)) %>% 
	  left_join(read.table(paste0(mapper_file),comment.char = "", header=T, sep="\t") %>% dplyr::select(all_of(varnames),X.SampleID), by=c("X"="X.SampleID"))
	
	uwmpc_all <- wfpc %>% 
	  inner_join(fut2, by=c("V1"="X")) %>%
	  mutate(shape_var = as.factor(as.character(get(shape_var))),
	         color_var = as.factor(as.character(get(color_var))),
	         circle_var = as.factor(as.character(get(circle_var))),
	         fill_var = as.factor(as.character(get(fill_var))),
	         linetype_var = as.factor(as.character(get(linetype_var)))
	  ) %>% droplevels()
	
	shape_table <- table(list(uwmpc_all[,shape_var]))
	fill_table <- table(list(uwmpc_all[,fill_var]))
	color_table <- table(list(uwmpc_all[,color_var]))
	circle_table <- table(list(uwmpc_all[,circle_var]))
	linetype_table <- table(list(uwmpc_all[,linetype_var]))
	
	cat("Shapes:",length(shape_table),"values")
	print(shape_table)
	cat("Fill:",length(fill_table),"values")
	print(fill_table)
	cat("Color:",length(color_table),"values")
	print(color_table)
	cat("Ellipse:",length(circle_table),"values")
	print(circle_table)
	cat("Linetype:",length(circle_table),"values")
	print(linetype_table)
	# 
	# plot_shapes <- c("square","circle")
	# plot_fills <- c("green","green")
	# plot_colors <- c("red","black")
	# plot_ellipses <- c(1:38)
	# 
	# legend_position = "bottom"
	#str(uwmpc_all)
	
	uwmpc_all$axis1 <- round(uwmpc_all[,paste0("V",pco1+1)],10)
	uwmpc_all$axis2 <- round(uwmpc_all[,paste0("V",pco2+1)],10)
	
	if(is.null(legend_name)){
	  legend_name = c(color_var, fill_var, shape_var, linetype_var)
	}
	
	if(is.null(legend_label)){
	  legend_label = list(colors = names(color_table)[color_table!=0], fills = names(fill_table)[fill_table!=0], shapes = names(shape_table)[shape_table!=0], linetypes = names(linetype_table)[linetype_table!=0])
	}
	
	ggplot(uwmpc_all, aes(x = axis1, y = axis2, colour=color_var, shape=shape_var, fill=fill_var, lty = linetype_var)) + 
		geom_point(alpha = 1, size=2.5) + #, shape=plot_shapes, size=3, col=plot_colors) + 
		scale_color_manual(name=legend_name[1], labels=legend_label[[1]], values=plot_colors) + 
		scale_fill_manual(name=legend_name[2], labels=legend_label[[2]], values=plot_fills) + 
		scale_shape_manual(name=legend_name[3],labels=legend_label[[3]], values=plot_shapes) + 
#	 scale_linetype_manual(name=legend_name[4],labels=legend_label[[4]], values=plot_ellipses) +
	  #		scale_fill_manual(name="Legend", values=c("red","black","white")) +
		theme(panel.background = element_blank(), 
					axis.line = element_line(), 
					axis.ticks=element_blank(), 
					axis.title=element_text(size=14),	
					#           title=element_text(size=16),
					legend.position=legend_position,
					plot.title = element_text(size=16, hjust=0),
					axis.text = element_blank())+
		labs(x = paste("PCo",pco1," ( ",round(as.numeric(pc_values[pco1]),1),"% )",sep=""), 
				 y = paste("PCo",pco2," ( ",round(as.numeric(pc_values[pco2]),1),"% )",sep=""), 
				 title = title_name_add) + 
		stat_ellipse(aes(x = axis1, y = axis2, group= circle_var), show.legend = T, type = "t", geom = "polygon", level = 0.95, alpha = 0, inherit.aes=T) 
	  

}



pcoa_manyvar_dashtype <- function(folder_name = folder_name, mapper_file = mapper_file, title_name_add = title_name_add, legend_position = legend_position,pco1=1,pco2=2,shape_var = shape_var,color_var = color_var,fill_var = fill_var,circle_var = circle_var, linetype_var = linetype_var, plot_shapes = plot_shapes, plot_fills = plot_fills, plot_colors = plot_colors, plot_ellipses = plot_ellipses, circle_color = "black") {
  
  varnames = c(shape_var, color_var, fill_var, circle_var,linetype_var,"dashtype")
  varnames <- varnames[!duplicated(varnames)]
  wfpc <- head(read.table(paste('core-metrics-results-',folder_name,'_pcoa_results/ordination.txt',sep=""), sep="\t", fill=T, skip = 9,blank.lines.skip = T, header=F),-2)
  pc2 <- head(read.table(paste('core-metrics-results-',folder_name,'_pcoa_results/ordination.txt', sep = ""), sep="\t", fill=T, skip = 3, header=F, colClasses = "character"),2)
  pc_values <- as.numeric(as.character(pc2[2,]))*100
  
  wfpc <- wfpc[,c(1,pco1+1, pco2+1)] 
  
  fut2 <- read.table(paste('core-metrics-results-',folder_name,'_distance_matrix/distance-matrix.tsv',sep=""), header=T, sep="\t", stringsAsFactors = T) %>% 
    dplyr::select(X) %>% 
    mutate(X=as.character(X)) %>% 
    left_join(read.table(paste0(mapper_file),comment.char = "", header=T, sep="\t") %>% dplyr::select(all_of(varnames),X.SampleID), by=c("X"="X.SampleID"))
  
  uwmpc_all <- wfpc %>% 
    inner_join(fut2, by=c("V1"="X")) %>%
    mutate(shape_var = as.factor(as.character(get(shape_var))),
           color_var = as.factor(as.character(get(color_var))),
           circle_var = as.factor(as.character(get(circle_var))),
           fill_var = as.factor(as.character(get(fill_var))),
           linetype_var = as.factor(as.character(get(linetype_var)))
    ) %>% droplevels()
  
  shape_table <- table(list(uwmpc_all[,shape_var]))
  fill_table <- table(list(uwmpc_all[,fill_var]))
  color_table <- table(list(uwmpc_all[,color_var]))
  circle_table <- table(list(uwmpc_all[,circle_var]))
  linetype_table <- table(list(uwmpc_all[,linetype_var]))
  
  cat("Shapes:",length(shape_table),"values")
  print(shape_table)
  cat("Fill:",length(fill_table),"values")
  print(fill_table)
  cat("Color:",length(color_table),"values")
  print(color_table)
  cat("Ellipse:",length(circle_table),"values")
  print(circle_table)
  cat("Linetype:",length(circle_table),"values")
  print(linetype_table)
  # 
  # plot_shapes <- c("square","circle")
  # plot_fills <- c("green","green")
  # plot_colors <- c("red","black")
  # plot_ellipses <- c(1:38)
  # 
  # legend_position = "bottom"
  #str(uwmpc_all)
  
  ggplot(uwmpc_all, aes(x = get(paste0("V",pco1+1)), y = get(paste0("V",pco2+1)), colour=color_var, shape=shape_var, fill=fill_var, lty = linetype_var)) + 
    geom_point(alpha = 1, size=2.5) + #, shape=plot_shapes, size=3, col=plot_colors) + 
    scale_color_manual(name=color_var, labels=names(color_table)[color_table!=0], values=plot_colors) + 
    scale_fill_manual(name=fill_var, labels=names(fill_table)[fill_table!=0], values=plot_fills) + 
    scale_shape_manual(name=shape_var,labels=names(shape_table)[shape_table!=0], values=plot_shapes) + 
    scale_linetype_manual(name=linetype_var,labels=names(linetype_table)[linetype_table!=0], values=plot_ellipses) +
    #		scale_fill_manual(name="Legend", values=c("red","black","white")) +
    theme(axis.text=element_text(size=12), 
          panel.background = element_blank(), 
          axis.line = element_line(), 
          axis.ticks=element_line(), 
          axis.title=element_text(size=14),	
          #           title=element_text(size=16),
          legend.position=legend_position,
          plot.title = element_text(size=16, hjust=0))+
    labs(x = paste("PCo1 ( ",round(as.numeric(pc_values[1]),1),"% )",sep=""), 
         y = paste("PCo2 ( ",round(as.numeric(pc_values[2]),1),"% )",sep=""), 
         title = title_name_add) + 
    stat_ellipse(aes(x=V2, y=V3, group= circle_var), show.legend = T, type = "t", geom = "polygon", level = 0.95, alpha = 0, inherit.aes=T) 
  
  
}
# 
# linetype_var_ellipse = "st2"
# color_var_ellipse = "st2"
# fill_var_ellipse = "st2"
# linetype_vars_ellipse = c(1,2,3,4)
# color_vars_ellipse = c("grey","magenta","black","yellow")
# fill_vars_ellipse = c("black","grey","black","grey")
# shape_var = "site" # which variable is the classifier
# color_var = "site"
# fill_var = "site"
# circle_var = "st2"
# plot_shapes <- c(15,16,17,18,19,20,21) # starved sqaure ---- unstarved circle ---- fruit triangle ----- soil diamond
# plot_fills <- c("orange","violet","midnight blue","yellow","red","green1","blue","magenta")
# plot_colors <- c("orange","violet","midnight blue","yellow","red","green1","blue","magenta")
# plot_ellipses <- c(rep(1,4))
# linetype_var <- "st2" 
# circle_color <- c("orange","violet","midnight blue","yellow","red","green1","blue","magenta")
# i = "bray_curtis"
# 
# 
# pcoa_manyvar <- function(folder_name = folder_name, mapper_file = mapper_file, title_name_add = title_name_add, legend_position = legend_position,pco1=1,pco2=2,shape_var = shape_var,color_var = color_var,fill_var = fill_var,circle_var = circle_var, linetype_var = linetype_var, plot_shapes = plot_shapes, plot_fills = plot_fills, plot_colors = plot_colors, plot_ellipses = plot_ellipses, circle_color = "black", linetype_var_ellipse = linetype_var_ellipse, color_var_ellipse = color_var_ellipse, fill_var_ellipse = fill_var_ellipse, linetype_vars_ellipse = linetype_vars_ellipse, color_vars_ellipse = color_vars_ellipse, fill_vars_ellipse = fill_vars_ellipse) {
#   
#   varnames = c(shape_var, color_var, fill_var, circle_var,linetype_var)
#   varnames <- varnames[!duplicated(varnames)]
#   wfpc <- head(read.table(paste('core-metrics-results-',folder_name,'_pcoa_results/ordination.txt',sep=""), sep="\t", fill=T, skip = 9,blank.lines.skip = T, header=F),-2)
#   pc2 <- head(read.table(paste('core-metrics-results-',folder_name,'_pcoa_results/ordination.txt', sep = ""), sep="\t", fill=T, skip = 3, header=F, colClasses = "character"),2)
#   pc_values <- as.numeric(as.character(pc2[2,]))*100
#   
#   wfpc <- wfpc[,c(1,pco1+1, pco2+1)] 
#   
#   fut2 <- read.table(paste('core-metrics-results-',folder_name,'_distance_matrix/distance-matrix.tsv',sep=""), header=T, sep="\t", stringsAsFactors = T) %>% 
#     dplyr::select(X) %>% 
#     mutate(X=as.character(X)) %>% 
#     left_join(read.table(paste0(mapper_file),comment.char = "", header=T, sep="\t") %>% dplyr::select(all_of(varnames),X.SampleID), by=c("X"="X.SampleID"))
#   
#   uwmpc_all <- wfpc %>% 
#     inner_join(fut2, by=c("V1"="X")) %>%
#     mutate(shape_var = as.factor(as.character(get(shape_var))),
#            color_var = as.factor(as.character(get(color_var))),
#            circle_var = as.factor(as.character(get(circle_var))),
#            fill_var = as.factor(as.character(get(fill_var))),
#            linetype_var = as.factor(as.character(get(linetype_var))),
#            linetype_var_ellipse = as.factor(as.character(get(linetype_var_ellipse))),
#            color_var_ellipse = as.factor(as.character(get(color_var_ellipse))),
#            fill_var_ellipse = as.factor(as.character(get(fill_var_ellipse)))
#     ) %>% droplevels()
#   
#   shape_table <- table(list(uwmpc_all[,shape_var]))
#   fill_table <- table(list(uwmpc_all[,fill_var]))
#   color_table <- table(list(uwmpc_all[,color_var]))
#   circle_table <- table(list(uwmpc_all[,circle_var]))
#   linetype_table <- table(list(uwmpc_all[,linetype_var]))
#   
#   cat("Shapes:",length(shape_table),"values")
#   print(shape_table)
#   cat("Fill:",length(fill_table),"values")
#   print(fill_table)
#   cat("Color:",length(color_table),"values")
#   print(color_table)
#   cat("Ellipse:",length(circle_table),"values")
#   print(circle_table)
#   cat("Linetype:",length(circle_table),"values")
#   print(linetype_table)
#   # 
#   # plot_shapes <- c("square","circle")
#   # plot_fills <- c("green","green")
#   # plot_colors <- c("red","black")
#   # plot_ellipses <- c(1:38)
#   # 
#   # legend_position = "bottom"
#   #str(uwmpc_all)
#   
#   uwmpc_all$axis1 <- round(uwmpc_all[,paste0("V",pco1+1)],10)
#   uwmpc_all$axis2 <- round(uwmpc_all[,paste0("V",pco2+1)],10)
#   
#   ggplot(uwmpc_all, aes(x = axis1, y = axis2, colour=color_var, shape=shape_var, fill=fill_var)) + 
#     geom_point(alpha = 1, size=2.5) + #, shape=plot_shapes, size=3, col=plot_colors) + 
#     scale_color_manual(name=color_var, labels=names(color_table)[color_table!=0], values=plot_colors) + 
#     scale_fill_manual(name=fill_var, labels=names(fill_table)[fill_table!=0], values=plot_fills) + 
#     scale_shape_manual(name=shape_var,labels=names(shape_table)[shape_table!=0], values=plot_shapes) + 
#     #	  scale_linetype_manual(name=linetype_var,labels=names(linetype_table)[linetype_table!=0], values=plot_ellipses) +
#     #		scale_fill_manual(name="Legend", values=c("red","black","white")) +
#     theme(panel.background = element_blank(), 
#           axis.line = element_line(), 
#           axis.ticks=element_blank(), 
#           axis.title=element_text(size=14),	
#           #           title=element_text(size=16),
#           legend.position=legend_position,
#           plot.title = element_text(size=16, hjust=0),
#           axis.text = element_blank())+
#     labs(x = paste("PCo",pco1," ( ",round(as.numeric(pc_values[pco1]),1),"% )",sep=""), 
#          y = paste("PCo",pco2," ( ",round(as.numeric(pc_values[pco2]),1),"% )",sep=""), 
#          title = title_name_add) + 
#     stat_ellipse(aes(x = axis1, y = axis2, group= circle_var, color = color_var_ellipse, lty = linetype_var_ellipse, fill = fill_var_ellipse), show.legend = T, type = "t", geom = "polygon", level = 0.95, alpha = 0, inherit.aes=F) + 
#     scale_color_manual(name=color_var_ellipse, values=color_vars_ellipse) + 
#     scale_fill_manual(name=fill_var_ellipse, values=fill_vars_ellipse) + 
#     scale_linetype_manual(name=linetype_var_ellipse, values=linetype_vars_ellipse)
#   


pcoa_manyvar_solidlineplot <- function(folder_name = folder_name, mapper_file = mapper_file, title_name_add = title_name_add, legend_position = legend_position,pco1=1,pco2=2,shape_var = shape_var,color_var = color_var,fill_var = fill_var,circle_var = circle_var, linetype_var = linetype_var, plot_shapes = plot_shapes, plot_fills = plot_fills, plot_colors = plot_colors, plot_ellipses = plot_ellipses, circle_color = "black", legend_name = NULL, legend_label = NULL) {
  
  varnames = c(shape_var, color_var, fill_var, circle_var,linetype_var)
  varnames <- varnames[!duplicated(varnames)]
  	wfpc <- head(read.table(paste('core-metrics-results-',folder_name,'_pcoa_results/ordination.txt',sep=""), sep="\t", fill=T, skip = 9,blank.lines.skip = T, header=F),-2)
	pc2 <- head(read.table(paste('core-metrics-results-',folder_name,'_pcoa_results/ordination.txt', sep = ""), sep="\t", fill=T, skip = 3, header=F, colClasses = "character"),2)
	pc_values <- as.numeric(as.character(pc2[2,]))*100

		wfpc <- wfpc[,c(1,pco1+1, pco2+1)] 
	
	fut2 <- read.table(paste('core-metrics-results-',folder_name,'_distance_matrix/distance-matrix.tsv',sep=""), header=T, sep="\t", stringsAsFactors = T) %>% 
	  dplyr::select(X) %>% 
	  mutate(X=as.character(X)) %>% 
	  left_join(read.table(paste0(mapper_file),comment.char = "", header=T, sep="\t") %>% dplyr::select(all_of(varnames),X.SampleID), by=c("X"="X.SampleID"))
	
	uwmpc_all <- wfpc %>% 
	  inner_join(fut2, by=c("V1"="X")) %>%
	  mutate(shape_var = as.factor(as.character(get(shape_var))),
	         color_var = as.factor(as.character(get(color_var))),
	         circle_var = as.factor(as.character(get(circle_var))),
	         fill_var = as.factor(as.character(get(fill_var))),
	         linetype_var = as.factor(as.character(get(linetype_var)))
	  ) %>% droplevels()
	
	shape_table <- table(list(uwmpc_all[,shape_var]))
	fill_table <- table(list(uwmpc_all[,fill_var]))
	color_table <- table(list(uwmpc_all[,color_var]))
	circle_table <- table(list(uwmpc_all[,circle_var]))
	linetype_table <- table(list(uwmpc_all[,linetype_var]))
	
	cat("Shapes:",length(shape_table),"values")
	print(shape_table)
	cat("Fill:",length(fill_table),"values")
	print(fill_table)
	cat("Color:",length(color_table),"values")
	print(color_table)
	cat("Ellipse:",length(circle_table),"values")
	print(circle_table)
	cat("Linetype:",length(circle_table),"values")
	print(linetype_table)
	# 
	# plot_shapes <- c("square","circle")
	# plot_fills <- c("green","green")
	# plot_colors <- c("red","black")
	# plot_ellipses <- c(1:38)
	# 
	# legend_position = "bottom"
	#str(uwmpc_all)
	
	uwmpc_all$axis1 <- round(uwmpc_all[,paste0("V",pco1+1)],10)
	uwmpc_all$axis2 <- round(uwmpc_all[,paste0("V",pco2+1)],10)
	
	if(is.null(legend_name)){
	  legend_name = c(color_var, fill_var, shape_var, linetype_var)
	}
	
	if(is.null(legend_label)){
	  legend_label = list(colors = names(color_table)[color_table!=0], fills = names(fill_table)[fill_table!=0], shapes = names(shape_table)[shape_table!=0], linetypes = names(linetype_table)[linetype_table!=0])
	}
	
	ggplot(uwmpc_all, aes(x = axis1, y = axis2, colour=color_var, shape=shape_var, fill=fill_var, lty = linetype_var)) + 
		geom_point(alpha = 1, size=2.5) + #, shape=plot_shapes, size=3, col=plot_colors) + 
		scale_color_manual(name=legend_name[1], labels=legend_label[[1]], values=plot_colors) + 
		scale_fill_manual(name=legend_name[2], labels=legend_label[[2]], values=plot_fills) + 
		scale_shape_manual(name=legend_name[3],labels=legend_label[[3]], values=plot_shapes) + 
#	 scale_linetype_manual(name=legend_name[4],labels=legend_label[[4]], values=plot_ellipses) +
	  #		scale_fill_manual(name="Legend", values=c("red","black","white")) +
		theme(panel.background = element_blank(), 
					axis.line = element_line(), 
					axis.ticks=element_blank(), 
					axis.title=element_text(size=14),	
					#           title=element_text(size=16),
					legend.position=legend_position,
					plot.title = element_text(size=16, hjust=0),
					axis.text = element_blank())+
		labs(x = paste("PCo",pco1," ( ",round(as.numeric(pc_values[pco1]),1),"% )",sep=""), 
				 y = paste("PCo",pco2," ( ",round(as.numeric(pc_values[pco2]),1),"% )",sep=""), 
				 title = title_name_add) + 
		stat_ellipse(aes(x = axis1, y = axis2, group= circle_var, linetype = circle_var), show.legend = T, type = "t", geom = "polygon", level = 0.95, alpha = 0, inherit.aes=T, lty = 1) 
	  

}
