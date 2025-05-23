
make_taxonomy_forR <- function(taxon_file_path) {
  taxa_table <- read.table(paste0(taxon_file_path,"/taxonomy.tsv"), sep = "\t", header = T) %>% 
    tidyr::separate(col = Taxon, into = c("kingdom","phylum","class","order","family","genus","species"), remove = T, sep = ";") %>%
    dplyr::select(-Confidence)
  write.csv(taxa_table, file = paste0(taxon_file_path,"/taxonomy_forR.csv"), quote = F, sep = ",", row.names = F)
}

# folder_path = "seasonE473_noC/bray_curtis"
# var1 <- "year"
# comparisons <- c("2014_2014","2017_2017")
# the_xid <- "XID"


# 0.0.2 - with mapper file path
# seasonE473_noC/bray_curtis","year",comparisons = c("2014_2014","2014_2017
# folder_path = "Efemale/bray_curtis",var1 = "time_point",comparisons = c("1_1","2_2","3_3")


# 0.0.3 - updated 2022-06-17

# folder_path = "bray_curtis"
# var1 = "plotorder"
# title_name_add = "hi"
# the_xid = "X"
# comparisons = c('fruit_NA_fruit_NA','fruit_NA_soil_NA','flies_unstarved_fruit_NA','flies_starved_fruit_NA')


# var1="experiment"
# comparisons = "east coast_east coast"
make_distance_plot_1var <- function (file_path = file_path, folder_path = "weighted_unifrac", mapper_file = mapper_file, var1,comparisons, title_name_add="", the_xid="X", plot_type = "dot") {
   
	flies_unifrac_table <- read.table(paste('core-metrics-results-',file_path,'/',folder_path,'_distance_matrix/distance-matrix.tsv',sep=""), header=T, sep="\t") %>% 
	  dplyr::select(X) %>% 
	  left_join(read.table(mapper_file,comment.char = "", header=T, sep = "\t"), by=c("X"="X.SampleID")) %>% 
	  mutate(the_xid=gsub(x = X,pattern = "-",replacement = ".")) %>%
	  dplyr::select_(.dots=list(the_xid, var1)) 
	
	flies_unifrac_table
	
	flies_dm <- read.table(paste('core-metrics-results-',file_path,'/',folder_path,'_distance_matrix/distance-matrix.tsv',sep=""), header=T, sep="\t")
	fs2 <- as.dist(flies_dm[,2:dim(flies_dm)[2]])
	
	fs3 <- as.matrix(fs2); fs3[lower.tri(fs3)] <- NA
	dim(fs3)
	print(flies_unifrac_table[1,])
	
	colnames(fs3)
	## melt the table and get out just the right comparisons
	melt_table <- data.frame(melt(fs3)) %>% 
		filter(is.na(value)==F, value!=0) %>% 
	  mutate(Var1 = gsub("\\.","",Var1), Va2 = gsub("\\.","",Var2)) %>%
	  inner_join(flies_unifrac_table, by=c("Var1"=the_xid)) %>% 
		inner_join(flies_unifrac_table, by=c("Var2"=the_xid)) %>%  
		mutate(gs.x = as.character(get(paste0(var1,".x"))), gs.y=as.character(get(paste0(var1,".y")))) %>% 
		mutate(gs=ifelse(gs.x < gs.y, paste(gs.x,gs.y, sep="_"), paste(gs.y, gs.x, sep="_"))) %>% 
		mutate(gs2 = as.factor(as.character(gs))) 
	
	print(table(list(melt_table$gs2)))
	
	melt_table <- melt_table %>% 
		filter(gs2%in%comparisons) %>% droplevels()
	
	## run statistics
	def <- kruskal.test(value ~ gs2, melt_table);   #print(def$p.value)
	efg <- dunn.test(melt_table$value,melt_table$gs2, method = "bh", table = F, kw = F)
	ghi <- cldList(comparison = efg$comparisons, p.value = efg$P.adjusted, threshold = 0.05);   print(ghi)
	
	## make the summary statistics
	melt_table_plot <- melt_table %>% 
		group_by(gs2) %>% 
		dplyr::summarize(mean=mean(value), sem=sd(value)/sqrt(length(value))) %>%
	  ungroup()
	
	  
	## make the plot
	if (plot_type == "bar") {  
	  ggplot(melt_table_plot, aes(x=gs2, y=mean)) + 
				 	geom_bar(position=position_dodge(), stat="identity") +
				 	geom_errorbar(aes(ymin=mean+sem, ymax=mean-sem),
				 								width=.2,                    # Width of the error bars
				 								position=position_dodge(.9)) +
				 	coord_cartesian(ylim=c((min(melt_table_plot$mean-melt_table_plot$sem)*.9),max(melt_table_plot$mean+melt_table_plot$sem)*1.1)) + #scale_y_continuous(limits=c(.25,.4)) + 
				 	#	scale_x_discrete(labels=c("Base","Time1","Time2","Time3")) +
				 	theme(axis.text=element_text(size=14), 
				 				panel.background = element_blank(),
				 				axis.line = element_line(), 
				 				axis.ticks=element_line(), 
				 				axis.title=element_text(size=16), plot.title=element_text(hjust=0)) +
				 	labs(x="Season",
				 			 y="Index distance",
				 			 title=title_name_add) +
				 	geom_text(aes(label=ghi$Letter, y=mean+sem, vjust=-1.5), size=6)
	} else {
	  ggplot(melt_table, aes(x= gs2, y = value)) +
	    geom_jitter(width = 0.2) + 
	    geom_crossbar(data=melt_table_plot, aes(x = gs2, y = mean, ymin = mean, ymax = mean), size=.5, col="blue", width = .35)
	} 	  
	  
}

make_distance_plot_1var2023 <- function(file_path = file_path, folder_path = "weighted_unifrac", mapper_file = mapper_file, var1,comparisons, title_name_add="", the_xid="X", plot_type = "dot", alpha = 1) {
  
  flies_unifrac_table <- read.table(paste('core-metrics-results-',file_path,'/',folder_path,'_distance_matrix/distance-matrix.tsv',sep=""), header=T, sep="\t") %>% 
    dplyr::select(X) %>% 
    left_join(read.table(mapper_file,comment.char = "", header=T, sep = "\t"), by=c("X"="X.SampleID")) %>% 
    mutate(the_xid=gsub(x = X,pattern = "-",replacement = ".")) %>%
    dplyr::select("the_xid", var1) 
  
  flies_dm <- read.table(paste('core-metrics-results-',file_path,'/',folder_path,'_distance_matrix/distance-matrix.tsv',sep=""), header=T, sep="\t")
  fs2 <- as.dist(flies_dm[,2:dim(flies_dm)[2]])
  
  fs3 <- as.matrix(fs2); fs3[lower.tri(fs3)] <- NA

    ## melt the table and get out just the right comparisons
  melt_table <- data.frame(melt(fs3)) %>% 
    filter(is.na(value)==F, value!=0) %>% 
    #	  mutate(Var1 = gsub(".","",Var1), Va2 = gsub(".","",Var2)) %>%
    inner_join(flies_unifrac_table, by=c("Var1"="the_xid")) %>% 
    inner_join(flies_unifrac_table, by=c("Var2"="the_xid")) %>%  
    mutate(gs.x = as.character(get(paste0(var1,".x"))), gs.y=as.character(get(paste0(var1,".y")))) %>% 
    mutate(gs=ifelse(gs.x < gs.y, paste(gs.x,gs.y, sep="_"), paste(gs.y, gs.x, sep="_"))) %>% 
    mutate(gs2 = as.factor(as.character(gs))) 
  
  print(table(list(melt_table$gs2)))
  
  melt_table <- melt_table %>% 
    filter(gs2%in%comparisons) %>% droplevels()
  
  ## run statistics
  def <- kruskal.test(value ~ gs2, melt_table);   #print(def$p.value)
  efg <- dunn.test(melt_table$value,melt_table$gs2, method = "bh", table = F, kw = F)
  ghi <- cldList(comparison = efg$comparisons, p.value = efg$P.adjusted, threshold = 0.05);   print(ghi)
  
  ## make the summary statistics
  melt_table_plot <- melt_table %>% 
    group_by(gs2) %>% 
    dplyr::summarize(mean=mean(value), sem=sd(value)/sqrt(length(value))) %>%
    ungroup()
  
  
  ## make the plot
  if (plot_type == "bar") {  
    ggplot(melt_table_plot, aes(x=gs2, y=mean)) + 
      geom_bar(position=position_dodge(), stat="identity") +
      geom_errorbar(aes(ymin=mean+sem, ymax=mean-sem),
                    width=.2,                    # Width of the error bars
                    position=position_dodge(.9)) +
      coord_cartesian(ylim=c((min(melt_table_plot$mean-melt_table_plot$sem)*.9),max(melt_table_plot$mean+melt_table_plot$sem)*1.1)) + #scale_y_continuous(limits=c(.25,.4)) + 
      #	scale_x_discrete(labels=c("Base","Time1","Time2","Time3")) +
      theme(axis.text=element_text(size=14), 
            panel.background = element_blank(),
            axis.line = element_line(), 
            axis.ticks=element_line(), 
            axis.title=element_text(size=16), plot.title=element_text(hjust=0)) +
      labs(x="Season",
           y="Index distance",
           title=title_name_add) +
      geom_text(aes(label=ghi$Letter, y=mean+sem, vjust=-1.5), size=6)
  } else {
    ggplot(melt_table, aes(x= gs2, y = value)) +
      geom_jitter(width = 0.2, alpha = alpha) + 
      geom_crossbar(data=melt_table_plot, aes(x = gs2, y = mean, ymin = mean, ymax = mean), size=.5, col="yellow", width = .35)
  } 	  
  
  
}

# var2="treatment"
# folder_path <- "season4_noC/weighted_unifrac"
#folder_path <- "bray_curtis"
make_distance_plot_2var <- function (folder_path, mapper_file = mapper_file, var1, var2, comparisons, title_name_add="", the_xid="XID",anglenum = anglenum) {
  flies_unifrac_table <- read.table(paste('core-metrics-results-',folder_path,'_distance_matrix/distance-matrix.tsv',sep=""), header=T, sep="\t") %>% dplyr::select(X) %>% left_join(read.table(mapper_file,comment.char = "", sep="\t", header=T), by=c("X"="X.SampleID")) %>% mutate(XID=gsub(x = X,pattern = "-",replacement = ".")) %>% dplyr::select_(.dots=list(the_xid, var1, var2))
  
  flies_dm <- read.table(paste('core-metrics-results-',folder_path,'_distance_matrix/distance-matrix.tsv',sep=""), header=T, sep="\t")
  fs2 <- as.dist(flies_dm[,2:dim(flies_dm)[2]])
  
  fs3 <- as.matrix(fs2); fs3[lower.tri(fs3)] <- NA
  
  print(flies_unifrac_table[1,])
  ## melt the table and get out just the right comparisons
  fs3[1,]
  melt_table <- data.frame(melt(fs3)) %>% 
    filter(is.na(value)==F, value!=0) %>% 
    inner_join(flies_unifrac_table, by=c("Var1"="XID")) %>% 
    inner_join(flies_unifrac_table, by=c("Var2"="XID")) %>%  
    mutate(gs.x = paste0(as.character(get(paste0(var1,".x"))),as.character(get(paste0(var2,".x")))), gs.y=paste0(as.character(get(paste0(var1,".y"))),as.character(get(paste0(var2,".y"))))) %>% 
    mutate(gs=ifelse(gs.x < gs.y, paste(gs.x,gs.y, sep="_"), paste(gs.y, gs.x, sep="_"))) %>% 
    mutate(gs2 = as.factor(as.character(gs))) 
  
  str(melt_table)
  print(table(list(melt_table$gs2)))
  
  melt_table <- melt_table %>% 
    filter(gs2%in%comparisons) %>% droplevels()
  
  melt_table$gs2 <- factor(melt_table$gs2,levels=comparisons)
  
  table(list(melt_table$gs2))
  ## run statistics
  def <- kruskal.test(value ~ gs2, melt_table);   #print(def$p.value)
  print(def)
  efg <- dunn.test(melt_table$value,melt_table$gs2, method = "bh", table = F, kw = F)
  ghi <- cldList(comparison = efg$comparisons, p.value = efg$P.adjusted, threshold = 0.05);   print(ghi)
  
  #ghi$Group2 <- factor(ghi$Group, levels=reorder(comparisons))
  ghi <- ghi %>% slice(match(comparisons %>% gsub(pattern = " ",replacement = "",x = comparisons), Group))
  ghi
  
  ## make the summary statistics
  melt_table_plot <- melt_table %>% 
    group_by(gs2) %>% 
    summarize(mean=mean(value), sem=sd(value)/sqrt(length(value)))
  
  ## make the plot
  ggplot(melt_table_plot, aes(x=gs2, y=mean)) + 
    geom_bar(position=position_dodge(), stat="identity") +
    geom_errorbar(aes(ymin=mean+sem, ymax=mean-sem),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9)) +
    coord_cartesian(ylim=c((min(melt_table_plot$mean-melt_table_plot$sem)*.9),max(melt_table_plot$mean+melt_table_plot$sem)*1.1)) + #scale_y_continuous(limits=c(.25,.4)) + 
    #	scale_x_discrete(labels=c("Base","Time1","Time2","Time3")) +
    theme(axis.text=element_text(size=14),
          axis.text.x=element_text(angle=anglenum),
          panel.background = element_blank(),
          axis.line = element_line(), 
          axis.ticks=element_line(), 
          axis.title=element_text(size=16), plot.title=element_text(hjust=0)) +
    labs(x="Season",
         y="Index distance",
         title=title_name_add) +
    geom_text(aes(label=ghi$Letter, y=mean+sem, vjust=-1.5), size=6)
}
make_distance_plot_2var_with0 <- function (folder_path, mapper_file = mapper_file, var1, var2, comparisons, comparisons2, title_name_add="", the_xid="XID",anglenum = anglenum) {
	flies_unifrac_table <- read.table(paste('core-metrics-results-',folder_path,'_distance_matrix/distance-matrix.tsv',sep=""), header=T, sep="\t") %>% dplyr::select(X) %>% left_join(read.table(mapper_file,comment.char = "", sep="\t", header=T), by=c("X"="X.SampleID")) %>% mutate(XID=gsub(x = X,pattern = "-",replacement = ".")) %>% dplyr::select_(.dots=list(the_xid, var1, var2))
	
	flies_dm <- read.table(paste('core-metrics-results-',folder_path,'_distance_matrix/distance-matrix.tsv',sep=""), header=T, sep="\t")
	fs2 <- as.dist(flies_dm[,2:dim(flies_dm)[2]])
	
	fs3 <- as.matrix(fs2); fs3[lower.tri(fs3)] <- NA
	
	print(flies_unifrac_table[1,])
	## melt the table and get out just the right comparisons
	fs3[1,]
	melt_table <- data.frame(melt(fs3)) %>% 
		filter(is.na(value)==F, value!=0) %>% 
		inner_join(flies_unifrac_table, by=c("Var1"="XID")) %>% 
		inner_join(flies_unifrac_table, by=c("Var2"="XID")) %>%  
		mutate(gs.x = paste0(as.character(get(paste0(var1,".x"))),as.character(get(paste0(var2,".x")))), gs.y=paste0(as.character(get(paste0(var1,".y"))),as.character(get(paste0(var2,".y"))))) %>% 
		mutate(gs=ifelse(gs.x < gs.y, paste(gs.x,gs.y, sep="_"), paste(gs.y, gs.x, sep="_"))) %>% 
		mutate(gs2 = as.factor(as.character(gs))) 
	
	str(melt_table)
	print(table(list(melt_table$gs2)))
	
	melt_table <- melt_table %>% 
		filter(gs2%in%comparisons) %>% droplevels()
	
	melt_table$gs2 <- factor(melt_table$gs2,levels=comparisons)
	
	table(list(melt_table$gs2))
	## run statistics
	def <- kruskal.test(value ~ gs2, melt_table);   #print(def$p.value)
	print(def)
	efg <- dunn.test(melt_table$value,melt_table$gs2, method = "bh", table = F, kw = F)
	ghi <- cldList(comparison = efg$comparisons, p.value = efg$P.adjusted, threshold = 0.05);   print(ghi)
	#ghi$Group2 <- factor(ghi$Group, levels=reorder(comparisons))
	ghi <- ghi %>% slice(match(gsub(pattern = " ",replacement = "",x = comparisons2), Group))
	
	## make the summary statistics
	melt_table_plot <- melt_table %>% 
		group_by(gs2) %>% 
		summarize(mean=mean(value), sem=sd(value)/sqrt(length(value)))
	
	## make the plot
	ggplot(melt_table_plot, aes(x=gs2, y=mean)) + 
		geom_bar(position=position_dodge(), stat="identity") +
		geom_errorbar(aes(ymin=mean+sem, ymax=mean-sem),
									width=.2,                    # Width of the error bars
									position=position_dodge(.9)) +
		coord_cartesian(ylim=c((min(melt_table_plot$mean-melt_table_plot$sem)*.9),max(melt_table_plot$mean+melt_table_plot$sem)*1.1)) + #scale_y_continuous(limits=c(.25,.4)) + 
		#	scale_x_discrete(labels=c("Base","Time1","Time2","Time3")) +
		theme(axis.text=element_text(size=14),
					axis.text.x=element_text(angle=anglenum),
					panel.background = element_blank(),
					axis.line = element_line(), 
					axis.ticks=element_line(), 
					axis.title=element_text(size=16), plot.title=element_text(hjust=0)) +
		labs(x="Season",
				 y="Index distance",
				 title=title_name_add) +
		geom_text(aes(label=ghi$Letter, y=mean+sem, vjust=-1.5), size=6)
}

make_distance_plot_2var_with0_1comp <- function (folder_path, mapper_file = mapper_file, var1, var2, comparisons, comparisons2, title_name_add="", the_xid="XID",anglenum = anglenum) {
  flies_unifrac_table <- read.table(paste('core-metrics-results-',folder_path,'_distance_matrix/distance-matrix.tsv',sep=""), header=T, sep="\t") %>% dplyr::select(X) %>% left_join(read.table(mapper_file,comment.char = "", sep="\t", header=T), by=c("X"="X.SampleID")) %>% mutate(XID=gsub(x = X,pattern = "-",replacement = ".")) %>% dplyr::select_(.dots=list(the_xid, var1, var2))
  
  flies_dm <- read.table(paste('core-metrics-results-',folder_path,'_distance_matrix/distance-matrix.tsv',sep=""), header=T, sep="\t")
  fs2 <- as.dist(flies_dm[,2:dim(flies_dm)[2]])
  
  fs3 <- as.matrix(fs2); fs3[lower.tri(fs3)] <- NA
  
  print(flies_unifrac_table[1,])
  ## melt the table and get out just the right comparisons
  fs3[1,]
  melt_table <- data.frame(melt(fs3)) %>% 
    filter(is.na(value)==F, value!=0) %>% 
    inner_join(flies_unifrac_table, by=c("Var1"="XID")) %>% 
    inner_join(flies_unifrac_table, by=c("Var2"="XID")) %>%  
    mutate(gs.x = paste0(as.character(get(paste0(var1,".x"))),as.character(get(paste0(var2,".x")))), gs.y=paste0(as.character(get(paste0(var1,".y"))),as.character(get(paste0(var2,".y"))))) %>% 
    mutate(gs=ifelse(gs.x < gs.y, paste(gs.x,gs.y, sep="_"), paste(gs.y, gs.x, sep="_"))) %>% 
    mutate(gs2 = as.factor(as.character(gs))) 
  
  str(melt_table)
  print(table(list(melt_table$gs2)))
  
  melt_table <- melt_table %>% 
    filter(gs2%in%comparisons) %>% droplevels()
  
  melt_table$gs2 <- factor(melt_table$gs2,levels=comparisons)
  
  table(list(melt_table$gs2))

  ## make the summary statistics
  melt_table_plot <- melt_table %>% 
    group_by(gs2) %>% 
    summarize(mean=mean(value), sem=sd(value)/sqrt(length(value)))
  
  ## make the plot
  ggplot(melt_table_plot, aes(x=gs2, y=mean)) + 
    geom_bar(position=position_dodge(), stat="identity") +
    geom_errorbar(aes(ymin=mean+sem, ymax=mean-sem),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9)) +
    coord_cartesian(ylim=c((min(melt_table_plot$mean-melt_table_plot$sem)*.9),max(melt_table_plot$mean+melt_table_plot$sem)*1.1)) + #scale_y_continuous(limits=c(.25,.4)) + 
    #	scale_x_discrete(labels=c("Base","Time1","Time2","Time3")) +
    theme(axis.text=element_text(size=14),
          axis.text.x=element_text(angle=anglenum),
          panel.background = element_blank(),
          axis.line = element_line(), 
          axis.ticks=element_line(), 
          axis.title=element_text(size=16), plot.title=element_text(hjust=0)) +
    labs(x="Season",
         y="Index distance",
         title=title_name_add)
}

get_distance_comparisons <- function (file_path, folder_path, mapper_file = mapper_file, var1, var2, comparisons, title_name_add="", the_xid="XID",anglenum = anglenum) {
  flies_unifrac_table <- read.table(paste('core-metrics-results-',file_path,'/',folder_path,'_distance_matrix/distance-matrix.tsv',sep=""), header=T, sep="\t") %>% dplyr::select(X) %>% left_join(read.table(mapper_file,comment.char = "", sep="\t", header=T), by=c("X"="X.SampleID")) %>% mutate(XID=gsub(x = X,pattern = "-",replacement = ".")) %>% dplyr::select_(.dots=list(the_xid, var1, var2))
  
  flies_dm <- read.table(paste('core-metrics-results-',file_path,'/',folder_path,'_distance_matrix/distance-matrix.tsv',sep=""), header=T, sep="\t")
  fs2 <- as.dist(flies_dm[,2:dim(flies_dm)[2]])
  
  fs3 <- as.matrix(fs2); fs3[lower.tri(fs3)] <- NA
  
  print(flies_unifrac_table[1,])
  ## melt the table and get out just the right comparisons
  fs3[1,]
  melt_table <- data.frame(melt(fs3)) %>% 
    filter(is.na(value)==F, value!=0) %>% 
    inner_join(flies_unifrac_table, by=c("Var1"="XID")) %>% 
    inner_join(flies_unifrac_table, by=c("Var2"="XID")) %>%  
    mutate(gs.x = paste0(as.character(get(paste0(var1,".x"))),as.character(get(paste0(var2,".x")))), gs.y=paste0(as.character(get(paste0(var1,".y"))),as.character(get(paste0(var2,".y"))))) %>% 
    mutate(gs=ifelse(gs.x < gs.y, paste(gs.x,gs.y, sep="_"), paste(gs.y, gs.x, sep="_"))) %>% 
    mutate(gs2 = as.factor(as.character(gs))) 
  
  str(melt_table)
  print(table(list(melt_table$gs2)))
  return(table(list(melt_table$gs2)))
}


# make_distance_plot_1var <- function (file_path, var1, the_xid="XID",comparisons) {
# 	flies_unifrac_table <- read.table(paste('core-metrics-results-',file_path,'/unweighted_unifrac_distance_matrix/distance-matrix.tsv',sep=""), header=T, sep="\t") %>% dplyr::select(X) %>% left_join(read.table("metadata.tsv",comment.char = "", header=T), by=c("X"="X.SampleID")) %>% mutate(XID=gsub(x = X,pattern = "-",replacement = ".")) %>% dplyr::select_(.dots=list(the_xid, var1))
# 	
# 	flies_unifrac_table[1,]
# 	## melt the table and get out just the right comparisons
# 	melt_table <- data.frame(melt(fs3)) %>% 
# 		filter(is.na(value)==F, value!=0) %>% 
# 		inner_join(flies_unifrac_table, by=c("X1"="XID")) %>% 
# 		inner_join(flies_unifrac_table, by=c("X2"="XID")) %>%  
# 		mutate(gs.x = as.character(get(paste0(var1,".x"))), gs.y=as.character(get(paste0(var1,".y")))) %>% 
# 		mutate(gs=ifelse(gs.x < gs.y, paste(gs.x,gs.y, sep="_"), paste(gs.y, gs.x, sep="_"))) %>% 
# 		mutate(gs2 = as.factor(as.character(gs))) %>% 
# 		filter(gs2%in%comparisons) %>% droplevels()
# 	
# 	## run statistics
# 	def <- kruskal.test(value ~ gs2, melt_table);   #print(def$p.value)
# 	efg <- dunn.test(melt_table$value,melt_table$gs2, method = "bh", table = F, kw = F)
# 	ghi <- cldList(comparison = efg$comparisons, p.value = efg$P.adjusted, threshold = 0.05);   print(ghi)
# 	
# 	## make the summary statistics
# 	melt_table_plot <- melt_table %>% 
# 		group_by(gs2) %>% 
# 		summarize(mean=mean(value), sem=sd(value)/sqrt(length(value)))
# 	
# 	melt_table
# 	## make the plot
# 	return(ggplot(melt_table_plot, aes(x=gs2, y=mean)) + 
# 				 	geom_bar(position=position_dodge(), stat="identity") +
# 				 	geom_errorbar(aes(ymin=mean+sem, ymax=mean-sem),
# 				 								width=.2,                    # Width of the error bars
# 				 								position=position_dodge(.9)) +
# 				 	coord_cartesian(ylim=c((min(melt_table_plot$mean-melt_table_plot$sem)*.9),max(melt_table_plot$mean+melt_table_plot$sem)*1.1)) + #scale_y_continuous(limits=c(.25,.4)) + 
# 				 	#	scale_x_discrete(labels=c("Base","Time1","Time2","Time3")) +
# 				 	theme(axis.text=element_text(size=14), 
# 				 				panel.background = element_blank(),
# 				 				axis.line = element_line(), 
# 				 				axis.ticks=element_line(), 
# 				 				axis.title=element_text(size=16)) +
# 				 	labs(x="Season",
# 				 			 y="weighted Unifrac distance") +
# 				 	geom_text(aes(label=ghi$Letter, y=mean+sem, vjust=-1.5), size=6))
# }



make_distance_plot_1var_nostat <- function (folder_path, var1,comparisons, title_name_add="", the_xid="XID") {
	flies_unifrac_table <- read.table(paste('core-metrics-results-',folder_path,'_distance_matrix/distance-matrix.tsv',sep=""), header=T, sep="\t") %>% dplyr::select(X) %>% left_join(read.table("metadata.tsv",comment.char = "", header=T), by=c("X"="X.SampleID")) %>% mutate(XID=gsub(x = X,pattern = "-",replacement = ".")) %>% dplyr::select_(.dots=list(the_xid, var1))
	
	
	flies_dm <- read.table(paste('core-metrics-results-',folder_path,'_distance_matrix/distance-matrix.tsv',sep=""), header=T, sep="\t")
	fs2 <- as.dist(flies_dm[,2:dim(flies_dm)[2]])
	
	fs3 <- as.matrix(fs2); fs3[lower.tri(fs3)] <- NA
	
	print(flies_unifrac_table[1,])
	## melt the table and get out just the right comparisons
	melt_table <- data.frame(melt(fs3)) %>% 
		filter(is.na(value)==F, value!=0) %>% 
		inner_join(flies_unifrac_table, by=c("X1"="XID")) %>% 
		inner_join(flies_unifrac_table, by=c("X2"="XID")) %>%  
		mutate(gs.x = as.character(get(paste0(var1,".x"))), gs.y=as.character(get(paste0(var1,".y")))) %>% 
		mutate(gs=ifelse(gs.x < gs.y, paste(gs.x,gs.y, sep="_"), paste(gs.y, gs.x, sep="_"))) %>% 
		mutate(gs2 = as.factor(as.character(gs))) 
	
	print(table(list(melt_table$gs2)))
	
	melt_table <- melt_table %>% 
		filter(gs2%in%comparisons) %>% droplevels()
	
	## make the summary statistics
	melt_table_plot <- melt_table %>% 
		group_by(gs2) %>% 
		summarize(mean=mean(value), sem=sd(value)/sqrt(length(value)))
	
	## make the plot
	ggplot(melt_table_plot, aes(x=gs2, y=mean)) + 
		geom_bar(position=position_dodge(), stat="identity") +
		geom_errorbar(aes(ymin=mean+sem, ymax=mean-sem),
									width=.2,                    # Width of the error bars
									position=position_dodge(.9)) +
		coord_cartesian(ylim=c((min(melt_table_plot$mean-melt_table_plot$sem)*.9),max(melt_table_plot$mean+melt_table_plot$sem)*1.1)) + #scale_y_continuous(limits=c(.25,.4)) + 
		#	scale_x_discrete(labels=c("Base","Time1","Time2","Time3")) +
		theme(axis.text=element_text(size=14), 
					panel.background = element_blank(),
					axis.line = element_line(), 
					axis.ticks=element_line(), 
					axis.title=element_text(size=16), plot.title=element_text(hjust=0)) +
		labs(x="Season",
				 y="Index distance",
				 title=title_name_add) 
}

