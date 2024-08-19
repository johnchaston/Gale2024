library(dplyr)
#taxonomic_level = "family"
make_ancom2_plots <- function(file_path, map_path="", mapper_file=mapper_file, taxonomic_level=taxonomic_level,var1=var1,var2=var2,var3=var3, var4 = var4,var5 = NULL, newcol="newcol", the_id = "X.SampleID", main.var=main.var, adj.formula=adj.formula, repeat.var = repeat.var, multcorr=multcorr, sig=sig, prev.cut=prev.cut, Group = Group, random.formula=NULL, return_plots = FALSE, return_tables = FALSE) {
  
  rm(otu_table, ref_names, otu2, map2, phyl2, phyl3, col1drop, col2drop, phyl4, ancom.OTU, detected_taxa, otu3, otu4, otu5, phyl5, plist, row_list)

  otu_table <- read.table(paste('core-metrics-results-',file_path,'/rarefied_table.txt',sep=""), comment.char="", header=T, sep="\t", fill=T, skip=1) %>% 
    left_join(read.csv(paste0('taxonomy',map_path,'/taxonomy_forR.csv')), by=c("X.OTU.ID"="Feature.ID"))
  
  ref_names <- subset(colnames(otu_table), (colnames(otu_table)%in%c("X.OTU.ID","kingdom","phylum","class","order","family","genus","species")==F))
  
#  otu_table$X.OTU.ID <- as.factor(otu_table$X.OTU.ID)
  
  #  taxonomic_level <- "X.OTU.ID"
    ## Cluster rows by taxonomy hierarchically 
  taxon_vector <- c("kingdom","phylum","class","order","family","genus","species","X.OTU.ID")
  taxon_number <- which(taxon_vector==taxonomic_level)
  taxon_vector2 <- c()
  for (i in 1:taxon_number) {
    taxon_vector2 <- c(taxon_vector2,taxon_vector[i])
  }

  if(taxonomic_level!="X.OTU.ID") {
    otu2 <- otu_table[,c("X.OTU.ID",taxon_vector2)] %>% 
      tidyr::unite(clustered_taxonomy, (taxon_vector[1:length(taxon_vector2)])) %>%
      mutate(clustered_taxonomy = as.character(clustered_taxonomy)) %>%
      mutate(clustered_taxonomy = gsub(pattern = " ",replacement = "",x = clustered_taxonomy)) %>%
      mutate(clustered_taxonomy = gsub(pattern = "^__",replacement = "__unassigned",x = clustered_taxonomy)) %>%
      mutate(clustered_taxonomy = gsub(pattern = "___",replacement = "__unassigned_",x = clustered_taxonomy)) %>%
      mutate(clustered_taxonomy = gsub(pattern = "__$",replacement = "__unassigned",x = clustered_taxonomy))
  } else if (taxonomic_level == "X.OTU.ID") {
    otu2 <- otu_table[,c("X.OTU.ID",taxon_vector2)] %>% 
      mutate(OTUID2 = X.OTU.ID) %>%
      tidyr::unite(clustered_taxonomy, c(taxon_vector[1:7],"OTUID2")) %>%
      mutate(clustered_taxonomy = as.character(clustered_taxonomy)) %>%
      mutate(clustered_taxonomy = gsub(pattern = " ",replacement = "",x = clustered_taxonomy)) %>%
      mutate(clustered_taxonomy = gsub(pattern = "^__",replacement = "__unassigned",x = clustered_taxonomy)) %>%
      mutate(clustered_taxonomy = gsub(pattern = "___",replacement = "__unassigned_",x = clustered_taxonomy)) %>%
      mutate(clustered_taxonomy = gsub(pattern = "__$",replacement = "__unassigned",x = clustered_taxonomy)) %>%
      dplyr::select(-X.OTU.ID.1)
        
    } 

  if(is.null(var5)) {
    map2 <- read.table(mapper_file,comment.char = "", header=T, fill=T,sep="\t") %>% 
      dplyr::select_(.dots = list(the_id,var1,var2,var3,var4, Group))
  } else {
    map2 <- read.table(mapper_file,comment.char = "", header=T, fill=T,sep="\t") %>% 
      dplyr::select_(.dots = list(the_id,var1,var2,var3,var4, var5,Group))
  }
  #print(head(map2))
    
  phyl2 <- otu_table %>% 
    inner_join(otu2, by = "X.OTU.ID") %>%
    group_by(clustered_taxonomy) %>% 
    summarize_at(ref_names,sum, na.rm=T)
  rownames(phyl2) <- as.character(unlist(phyl2$clustered_taxonomy))
  
  phyl3 <- phyl2 %>% dplyr::select(-clustered_taxonomy) %>% t() %>% data.frame() 
  colnames(phyl3) <- rownames(phyl2)
  phyl3 <- mutate(phyl3, X.SampleID=rownames(phyl3))
  phyl3$Sample.ID <- gsub("\\.","", phyl3$X.SampleID)
  phyl3$Sample.ID <- gsub("_","", phyl3$Sample.ID)
  map2$Sample.ID <- gsub("_","",map2$X.SampleID)
  map2$Sample.ID <- gsub("-","",map2$Sample.ID)
  map2$Sample.ID <- gsub("\\.","",map2$Sample.ID) # added 2021-05-05
  
  phyl5 <- phyl3 %>% dplyr::select(Sample.ID, everything()) %>% dplyr::select(-X.SampleID)
  colnames(phyl5)[1] <- "Sample.ID"
  
  map3 <- map2 %>% dplyr::select(Sample.ID, everything()) %>% dplyr::select(-X.SampleID) %>% filter(Sample.ID %in% phyl5$Sample.ID)
  colnames(map3)[1] <- "Sample.ID"
  
  phyl4 <- phyl3 %>% inner_join(map2, by=c("Sample.ID")) %>% mutate(Group = get(Group)) %>% droplevels()#%>% dplyr::select_(.dots = list(col1drop, col2drop, paste0("-",var1))) %>% dplyr::select(-sample2)
  rm(comparison_test)
  
  if (taxonomic_level == "X.OTU.ID") {
    colnames(phyl5) <- paste0("X",colnames(phyl5))
    colnames(phyl5)[1] <- "Sample.ID"
  }
  
  write.csv(phyl5,"phyl5b.csv")
  comparison_test = ANCOM(otu_data = phyl5, 
              meta_data = map3, 
              main_var = main.var,  
              zero_cut = prev.cut, 
              p_adjust_method = multcorr, 
              alpha = sig, 
              adj_formula = adj.formula, 
              rand_formula = random.formula)
  
  #print(comparison_test)
  write.csv(comparison_test, paste("ancomexcel_",main.var,"_",file_path,"_",taxonomic_level,".csv",sep=""))
  dt1 <- data.frame(comparison_test) %>% filter(detected_0.9==T) %>% dplyr::select(otu_id)
  detected_taxa <- data.frame(OTU=unlist(sapply(X = dt1$otu_id, function(x) ifelse(substring(x,1,1)=="X",substring(x,2),substring(x,1)))))	%>% unlist() %>% unname() %>% as.character()
  
  otu3 <- otu_table[,c("X.OTU.ID","kingdom","phylum","class","order","family","genus","species")]
  otu3$OTU <- as.character(otu3$X.OTU.ID)
  
  if(taxonomic_level!="OTU") {
   otu3[,paste(taxonomic_level)] = gsub(" ","",otu3[,paste(taxonomic_level)])
  }

  otu5 <- otu3 %>% inner_join(otu2, by = "X.OTU.ID") %>% filter(clustered_taxonomy%in%detected_taxa) %>% distinct(get(taxonomic_level), .keep_all = T)

  extra_cols2 <- names(table(list(phyl4$Group)))
  extra_cols <- unname(c(extra_cols2,sapply(extra_cols2, function(x) paste0(x,"_sem"))))
  otu5[extra_cols] <- "NA"
  
  phyl3$sample2 <- gsub("\\.","", phyl3$X.SampleID)
  phyl3$sample2 <- gsub("_","", phyl3$sample2)
  
  i
  for (i in 1:length(detected_taxa)) {
    
    rm(phyl6, column_name)
    column_name <- detected_taxa[i]
    try(phyl6 <-phyl4 %>% mutate(rabun = get(column_name)/sum(otu_table[,2])) %>% dplyr::select(rabun, Group) %>% group_by(Group) %>% dplyr::summarize(mean=mean(rabun), sem=sd(rabun)/sqrt(length(rabun))),T)
    if(exists("phyl6")==F) {
      try(phyl6 <- phyl4 %>% mutate(rabun = get(paste("X",column_name, sep=""))/sum(otu_table[,2])) %>% dplyr::select(rabun, Group) %>% group_by(Group) %>% dplyr::summarize(mean=mean(rabun), sem=sd(rabun)/sqrt(length(rabun))))
    }
    
    phyl6
    rm(t)
    for (t in 1:length(extra_cols2)) {
      otu5[i,paste0(extra_cols2[t])] <- phyl6 %>% filter(Group==extra_cols2[t]) %>% dplyr::select(mean) %>% unlist()
      otu5[i,paste0(extra_cols2[t],"_sem")] <- phyl6 %>% filter(Group==extra_cols2[t]) %>% dplyr::select(sem) %>% unlist()
    }
    
    assign(paste("p",i,sep="_"),value = ggplot(phyl6, aes(x=Group, y=mean)) + 
             geom_bar(position=position_dodge(), stat="identity") +
             geom_errorbar(aes(ymin=mean+sem, ymax=mean-sem),
                           width=.2,                    # Width of the error bars
                           position=position_dodge(.9)) +
             coord_cartesian(ylim=c((min(phyl6$mean-phyl6$sem)*.9),max(phyl6$mean+phyl6$sem)*1.3)) + #scale_y_continuous(limits=c(.25,.4)) + 
             theme(axis.text=element_text(size=14), 
                   panel.background = element_blank(),
                   axis.line = element_line(), 
                   axis.ticks=element_line(), 
                   axis.title=element_text(size=16),
                   title=element_text(size=13)) +
             labs(y="relative abundance",x=var1,title=detected_taxa[i]))
  }
  
  ## make a list of the plots
  if(length(detected_taxa)>1) {
    plist <- list(p_1)
    for(q in 2:length(detected_taxa)) {
      plist[[q]] <- get(paste0("p_",q))
    }
  } else if (length(detected_taxa)>0) {
    plist <- list(p_1)
  }
  
  write.csv(otu5, paste("ancom_",var1,"_",file_path,"_",taxonomic_level,".csv",sep=""))
  
  jpeg(h=800*1.25, w=1600*1.25, paste("ancom_",var1,"_",file_path,"_",taxonomic_level,".jpg",sep=""), units = "px", quality = 0.9)
  do.call("grid.arrange", c(plist))
  dev.off()
  
  if(return_plots == T) {
  return(plist)
  }
  
  if(return_tables == T) {
    for (i in 1:length(detected_taxa)) {
      
      rm(phyl6, phyl7, column_name)
      column_name <- detected_taxa[i]
      try(phyl6 <-phyl4 %>% mutate(rabun = get(column_name)/sum(otu_table[,2])) %>% dplyr::select(rabun, Group) %>% group_by(Group) %>% dplyr::summarize(mean=mean(rabun), sem=sd(rabun)/sqrt(length(rabun))),T)
      try(phyl7 <-phyl4 %>% mutate(rabun = get(column_name)/sum(otu_table[,2])) %>% dplyr::select(rabun, Group),T)
      if(exists("phyl6")==F) {
        try(phyl6 <- phyl4 %>% mutate(rabun = get(paste("X",column_name, sep=""))/sum(otu_table[,2])) %>% dplyr::select(rabun, Group) %>% group_by(Group) %>% dplyr::summarize(mean=mean(rabun), sem=sd(rabun)/sqrt(length(rabun))))
      }
      if(exists("phyl7")==F) {
        try(phyl7 <- phyl4 %>% mutate(rabun = get(paste("X",column_name, sep=""))/sum(otu_table[,2])) %>% dplyr::select(rabun, Group),T)
      }
      phyl6
      phyl7
      rm(t)
      for (t in 1:length(extra_cols2)) {
        otu5[i,paste0(extra_cols2[t])] <- phyl6 %>% filter(Group==extra_cols2[t]) %>% dplyr::select(mean) %>% unlist()
        otu5[i,paste0(extra_cols2[t],"_sem")] <- phyl6 %>% filter(Group==extra_cols2[t]) %>% dplyr::select(sem) %>% unlist()
      }
      
      assign(paste("table",i,sep="_"),value = phyl6)
      assign(paste("data",i,sep="_"),value = phyl7)
    }
    
    ## make a list of the plots
    if(length(detected_taxa)>1) {
      tablelist <- list(table_1)
      datalist <- list(data_1)
      for(q in 2:length(detected_taxa)) {
        tablelist[[q]] <- get(paste0("table_",q))
        datalist[[q]] <- get(paste0("data_",q))
      }
    } else if (length(detected_taxa)>0) {
      tablelist <- list(table_1)
      datalist <- list(data_1)
    }
    return(list(tablelist, datalist))
  }
}


make_abundance_plot <- function(file_path, map_path="", mapper_file=mapper_file, taxonomic_level=taxonomic_level,var1=var1,var2=var2,var3=var3, var4 = var4,var5 = NULL, newcol="newcol", the_id = "X.SampleID", main.var=main.var, adj.formula=adj.formula, repeat.var = repeat.var, multcorr=multcorr, sig=sig, prev.cut=prev.cut, Group = Group, random.formula=NULL, return_plots = FALSE, return_tables = FALSE, specified_group = specified_group, return_table = F) {
  
  rm(otu_table, ref_names, otu2, map2, phyl2, phyl3, col1drop, col2drop, phyl4, ancom.OTU, detected_taxa, otu3, otu4, otu5, phyl5, plist, row_list)
  
  otu_table <- read.table(paste('core-metrics-results-',file_path,'/rarefied_table.txt',sep=""), comment.char="", header=T, sep="\t", fill=T, skip=1) %>% 
    left_join(read.csv(paste0('taxonomy',map_path,'/taxonomy_forR.csv')), by=c("X.OTU.ID"="Feature.ID"))
  
  ref_names <- subset(colnames(otu_table), (colnames(otu_table)%in%c("X.OTU.ID","kingdom","phylum","class","order","family","genus","species")==F))
  
  #  otu_table$X.OTU.ID <- as.factor(otu_table$X.OTU.ID)
  
  #  taxonomic_level <- "X.OTU.ID"
  ## Cluster rows by taxonomy hierarchically 
  taxon_vector <- c("kingdom","phylum","class","order","family","genus","species","X.OTU.ID")
  taxon_number <- which(taxon_vector==taxonomic_level)
  taxon_vector2 <- c()
  for (i in 1:taxon_number) {
    taxon_vector2 <- c(taxon_vector2,taxon_vector[i])
  }
  
  if(taxonomic_level!="X.OTU.ID") {
    otu2 <- otu_table[,c("X.OTU.ID",taxon_vector2)] %>% 
      tidyr::unite(clustered_taxonomy, (taxon_vector[1:length(taxon_vector2)])) %>%
      mutate(clustered_taxonomy = as.character(clustered_taxonomy)) %>%
      mutate(clustered_taxonomy = gsub(pattern = " ",replacement = "",x = clustered_taxonomy)) %>%
      mutate(clustered_taxonomy = gsub(pattern = "^__",replacement = "__unassigned",x = clustered_taxonomy)) %>%
      mutate(clustered_taxonomy = gsub(pattern = "___",replacement = "__unassigned_",x = clustered_taxonomy)) %>%
      mutate(clustered_taxonomy = gsub(pattern = "__$",replacement = "__unassigned",x = clustered_taxonomy))
  } else if (taxonomic_level == "X.OTU.ID") {
    otu2 <- otu_table[,c("X.OTU.ID",taxon_vector2)] %>% 
      mutate(OTUID2 = X.OTU.ID) %>%
      tidyr::unite(clustered_taxonomy, c(taxon_vector[1:7],"OTUID2")) %>%
      mutate(clustered_taxonomy = as.character(clustered_taxonomy)) %>%
      mutate(clustered_taxonomy = gsub(pattern = " ",replacement = "",x = clustered_taxonomy)) %>%
      mutate(clustered_taxonomy = gsub(pattern = "^__",replacement = "__unassigned",x = clustered_taxonomy)) %>%
      mutate(clustered_taxonomy = gsub(pattern = "___",replacement = "__unassigned_",x = clustered_taxonomy)) %>%
      mutate(clustered_taxonomy = gsub(pattern = "__$",replacement = "__unassigned",x = clustered_taxonomy)) %>%
      dplyr::select(-X.OTU.ID.1)
    
  } 
  
  if(is.null(var5)) {
    map2 <- read.table(mapper_file,comment.char = "", header=T, fill=T,sep="\t") %>% 
      dplyr::select_(.dots = list(the_id,var1,var2,var3,var4, Group))
  } else {
    map2 <- read.table(mapper_file,comment.char = "", header=T, fill=T,sep="\t") %>% 
      dplyr::select_(.dots = list(the_id,var1,var2,var3,var4, var5,Group))
  }
  #print(head(map2))
  
  phyl2 <- otu_table %>% 
    inner_join(otu2, by = "X.OTU.ID") %>%
    group_by(clustered_taxonomy) %>% 
    summarize_at(ref_names,sum, na.rm=T)
  rownames(phyl2) <- as.character(unlist(phyl2$clustered_taxonomy))
  
  phyl3 <- phyl2 %>% dplyr::select(-clustered_taxonomy) %>% t() %>% data.frame() 
  colnames(phyl3) <- rownames(phyl2)
  phyl3 <- mutate(phyl3, X.SampleID=rownames(phyl3))
  phyl3$Sample.ID <- gsub("\\.","", phyl3$X.SampleID)
  phyl3$Sample.ID <- gsub("_","", phyl3$Sample.ID)
  map2$Sample.ID <- gsub("_","",map2$X.SampleID)
  map2$Sample.ID <- gsub("-","",map2$Sample.ID)
  map2$Sample.ID <- gsub("\\.","",map2$Sample.ID) # added 2021-05-05
  
  phyl5 <- phyl3 %>% dplyr::select(Sample.ID, everything()) %>% dplyr::select(-X.SampleID)
  colnames(phyl5)[1] <- "Sample.ID"
  
  map3 <- map2 %>% dplyr::select(Sample.ID, everything()) %>% dplyr::select(-X.SampleID) %>% filter(Sample.ID %in% phyl5$Sample.ID)
  colnames(map3)[1] <- "Sample.ID"
  
  phyl4 <- phyl3 %>% inner_join(map2, by=c("Sample.ID")) %>% mutate(Group = get(Group)) %>% droplevels()#%>% dplyr::select_(.dots = list(col1drop, col2drop, paste0("-",var1))) %>% dplyr::select(-sample2)
  rm(comparison_test)
  
  if (taxonomic_level == "X.OTU.ID") {
    colnames(phyl5) <- paste0("X",colnames(phyl5))
    colnames(phyl5)[1] <- "Sample.ID"
  }
  
  detected_taxa <- specified_group
  
  otu3 <- otu_table[,c("X.OTU.ID","kingdom","phylum","class","order","family","genus","species")]
  otu3$OTU <- as.character(otu3$X.OTU.ID)
  
  if(taxonomic_level!="OTU") {
    otu3[,paste(taxonomic_level)] = gsub(" ","",otu3[,paste(taxonomic_level)])
  }
  
# if(need_titles == T) {
#     otu5 <- otu3 %>% inner_join(otu2, by = "X.OTU.ID")%>% distinct(get(taxonomic_level), .keep_all=T)
#     return(otu5$clustered_taxonomy)
#   } else {
  
    otu5 <- otu3 %>% inner_join(otu2, by = "X.OTU.ID") %>% filter(clustered_taxonomy%in%detected_taxa) %>% distinct(get(taxonomic_level), .keep_all = T)
    
    extra_cols2 <- names(table(list(phyl4$Group)))
    extra_cols <- unname(c(extra_cols2,sapply(extra_cols2, function(x) paste0(x,"_sem"))))
    otu5[extra_cols] <- "NA"
    
    phyl3$sample2 <- gsub("\\.","", phyl3$X.SampleID)
    phyl3$sample2 <- gsub("_","", phyl3$sample2)
    
       
      rm(phyl6, column_name)
      column_name <- detected_taxa
      try(phyl6 <-phyl4 %>% mutate(rabun = get(column_name)/sum(otu_table[,2])) %>% dplyr::select(rabun, Group) %>% group_by(Group) %>% dplyr::summarize(mean=mean(rabun), sem=sd(rabun)/sqrt(length(rabun))),T)
      if(exists("phyl6")==F) {
        try(phyl6 <- phyl4 %>% mutate(rabun = get(paste("X",column_name, sep=""))/sum(otu_table[,2])) %>% dplyr::select(rabun, Group) %>% group_by(Group) %>% dplyr::summarize(mean=mean(rabun), sem=sd(rabun)/sqrt(length(rabun))))
      }
      
      rm(t)
      for (t in 1:length(extra_cols2)) {
        otu5[i,paste0(extra_cols2[t])] <- phyl6 %>% filter(Group==extra_cols2[t]) %>% dplyr::select(mean) %>% unlist()
        otu5[i,paste0(extra_cols2[t],"_sem")] <- phyl6 %>% filter(Group==extra_cols2[t]) %>% dplyr::select(sem) %>% unlist()
      }
      
      assign("outplot",value = ggplot(phyl6, aes(x=Group, y=mean)) + 
               geom_bar(position=position_dodge(), stat="identity") +
               geom_errorbar(aes(ymin=mean+sem, ymax=mean-sem),
                             width=.2,                    # Width of the error bars
                             position=position_dodge(.9)) +
               coord_cartesian(ylim=c((min(phyl6$mean-phyl6$sem)*.9),max(phyl6$mean+phyl6$sem)*1.3)) + #scale_y_continuous(limits=c(.25,.4)) + 
               theme(axis.text=element_text(size=14), 
                     panel.background = element_blank(),
                     axis.line = element_line(), 
                     axis.ticks=element_line(), 
                     axis.title=element_text(size=16),
                     title=element_text(size=13)) +
               labs(y="relative abundance",x=var1,title=detected_taxa))
    
      if(return_table == T) {
        return(phyl6)
      } else {
        return(outplot)
        
      }
    
  }




make_ancom2_plots_randomize <- function (file_path, map_path="", mapper_file=mapper_file, taxonomic_level=taxonomic_level,var1=var1,var2=var2,var3=var3 ,correction_level=3,newcol="newcol", the_id = "X.SampleID", main.var=main.var, adj.formula=adj.formula, repeat.var = repeat.var, multcorr=multcorr, sig=sig, prev.cut=prev.cut, Group = Group, random.formula=NULL, counter = counter) {
  
  rm(otu_table, ref_names, otu2, map2, phyl2, phyl3, col1drop, col2drop, phyl4, ancom.OTU, detected_taxa, otu3, otu4, otu5, phyl5, plist, row_list)
  
  otu_table <- read.table(paste('core-metrics-results-',file_path,'/rarefied_table.txt',sep=""), comment.char="", header=T, sep="\t", fill=T, skip=1) %>% 
    left_join(read.csv(paste0('taxonomy',map_path,'/taxonomy_forR.csv')), by=c("X.OTU.ID"="Feature.ID"))
  
  ref_names <- subset(colnames(otu_table), (colnames(otu_table)%in%c("X.OTU.ID","kingdom","phylum","class","order","family","genus","species")==F))
  
  #  taxonomic_level <- "X.OTU.ID"
  ## Cluster rows by taxonomy hierarchically 
  taxon_vector <- c("kingdom","phylum","class","order","family","genus","species","X.OTU.ID")
  taxon_number <- which(taxon_vector==taxonomic_level)
  taxon_vector2 <- c()
  for (i in 1:taxon_number) {
    taxon_vector2 <- c(taxon_vector2,taxon_vector[i])
  }
  
  if(taxonomic_level!="X.OTU.ID") {
    otu2 <- otu_table[,c("X.OTU.ID",taxon_vector2)] %>% 
      tidyr::unite(clustered_taxonomy, (taxon_vector[1:length(taxon_vector2)])) %>%
      mutate(clustered_taxonomy = as.character(clustered_taxonomy)) %>%
      mutate(clustered_taxonomy = gsub(pattern = " ",replacement = "",x = clustered_taxonomy)) %>%
      mutate(clustered_taxonomy = gsub(pattern = "^__",replacement = "__unassigned",x = clustered_taxonomy)) %>%
      mutate(clustered_taxonomy = gsub(pattern = "___",replacement = "__unassigned_",x = clustered_taxonomy)) %>%
      mutate(clustered_taxonomy = gsub(pattern = "__$",replacement = "__unassigned",x = clustered_taxonomy))
  } else if (taxonomic_level == "X.OTU.ID") {
    otu2 <- otu_table[,c("X.OTU.ID",taxon_vector2)] %>% 
      mutate(OTUID2 = X.OTU.ID) %>%
      tidyr::unite(clustered_taxonomy, c(taxon_vector[1:7],"OTUID2")) %>%
      mutate(clustered_taxonomy = as.character(clustered_taxonomy)) %>%
      mutate(clustered_taxonomy = gsub(pattern = " ",replacement = "",x = clustered_taxonomy)) %>%
      mutate(clustered_taxonomy = gsub(pattern = "^__",replacement = "__unassigned",x = clustered_taxonomy)) %>%
      mutate(clustered_taxonomy = gsub(pattern = "___",replacement = "__unassigned_",x = clustered_taxonomy)) %>%
      mutate(clustered_taxonomy = gsub(pattern = "__$",replacement = "__unassigned",x = clustered_taxonomy)) %>%
      dplyr::select(-X.OTU.ID.1)
    
  } 
  
  map2 <- read.table(mapper_file,comment.char = "", header=T, fill=T,sep="\t") %>% 
    dplyr::select_(.dots = list(the_id,var1,var2,var3,Group))
  
  phyl2 <- otu_table %>% 
    inner_join(otu2, by = "X.OTU.ID") %>%
    group_by(clustered_taxonomy) %>% 
    summarize_at(ref_names,sum, na.rm=T)
  rownames(phyl2) <- as.character(unlist(phyl2$clustered_taxonomy))
  
  phyl3 <- phyl2 %>% dplyr::select(-clustered_taxonomy) %>% t() %>% data.frame() 
  colnames(phyl3) <- rownames(phyl2)
  phyl3 <- mutate(phyl3, X.SampleID=rownames(phyl3))
  phyl3$Sample.ID <- gsub("\\.","", phyl3$X.SampleID)
  phyl3$Sample.ID <- gsub("_","", phyl3$Sample.ID)
  map2$Sample.ID <- gsub("_","",map2$X.SampleID)
  map2$Sample.ID <- gsub("-","",map2$Sample.ID)
  map2$Sample.ID <- gsub("\\.","",map2$Sample.ID) # added 2021-05-05
  
  phyl5 <- phyl3 %>% dplyr::select(Sample.ID, everything()) %>% dplyr::select(-X.SampleID)
  colnames(phyl5)[1] <- "Sample.ID"
  
  map3 <- map2 %>% dplyr::select(Sample.ID, everything()) %>% dplyr::select(-X.SampleID) %>% filter(Sample.ID %in% phyl5$Sample.ID)
  colnames(map3)[1] <- "Sample.ID"

  ## new for random  
  # random_lookup <- phyl4 %>%
  #   group_by(Sample) %>% 
  #   dplyr::summarize(count = dplyr::n()) %>%
  #   data.frame() %>%
  #   rowwise() %>%
  #   mutate(rand = runif(1,0,1)) %>%
  #   arrange(desc(rand)) %>% 
  #   ungroup() %>%
  #   slice_head(n = 5) %>%
  #   pull(Sample) %>%
  #   as.character()
  # 
  ## new for random()

  phyl4 <- phyl3 %>% 
    inner_join(map2, by=c("Sample.ID")) %>% 
    mutate(Group = get(Group)) %>% droplevels()

  random_lookup <- phyl4 %>%
    group_by(Sample) %>% 
    dplyr::summarize(count = dplyr::n()) %>%
    data.frame() %>%
    sample_n(size = 5) %>% 
    pull(Sample) %>%
    as.character()
  #print(random_lookup)
    for(i in 1:length(map3[,var1])) {
    if(map3[i,"Sample"] %in% random_lookup) {
      map3[i,var1] <- "E"
    } else {map3[i,var1] <- "C"  }}
 #random_lookup <- c("46J", "4D",  "28J", "7H",  "33R") 
 #random_lookup <- c("31N", "40A",  "42C", "11L",  "20A")
   if (taxonomic_level == "X.OTU.ID") {
    colnames(phyl5) <- paste0("X",colnames(phyl5))
    colnames(phyl5)[1] <- "Sample.ID"
  }
  
  write.csv(phyl5,"phyl5b.csv")
  comparison_test = ANCOM(otu_data = phyl5, 
                          meta_data = map3, 
                          main_var = main.var,  
                          zero_cut = prev.cut, 
                          p_adjust_method = multcorr, 
                          alpha = sig, 
                          adj_formula = adj.formula, 
                          rand_formula = random.formula)
  names(table(map3 %>% filter(ASD2=="E") %>% select(Sample) %>% droplevels()))
# print(comparison_test)
#  write.csv(comparison_test, paste("ancomexcel_",main.var,"_",file_path,"_",taxonomic_level,".csv",sep=""))
  dt1 <- data.frame(comparison_test) %>% filter(detected_0.9==T) %>% dplyr::select(otu_id)
  detected_taxa <- data.frame(OTU=unlist(sapply(X = dt1$otu_id, function(x) ifelse(substring(x,1,1)=="X",substring(x,2),substring(x,1)))))	%>% unlist() %>% unname() %>% as.character()
  #print(detected_taxa)
  otu3 <- otu_table[,c("X.OTU.ID","kingdom","phylum","class","order","family","genus","species")]
  otu3$OTU <- as.character(otu3$X.OTU.ID)
  
  if(taxonomic_level!="OTU") {
    otu3[,paste(taxonomic_level)] = gsub(" ","",otu3[,paste(taxonomic_level)])
  }
  
  otu5 <- otu3 %>% inner_join(otu2, by = "X.OTU.ID") %>% filter(clustered_taxonomy%in%detected_taxa) %>% distinct(get(taxonomic_level), .keep_all = T) %>% mutate(randrun = counter)

  extra_cols2 <- names(table(list(phyl4$Group)))
  extra_cols <- unname(c(extra_cols2,sapply(extra_cols2, function(x) paste0(x,"_sem"))))
  otu5[extra_cols] <- "NA"

  phyl3$sample2 <- gsub("\\.","", phyl3$X.SampleID)
  phyl3$sample2 <- gsub("_","", phyl3$sample2)

  for (i in 1:length(detected_taxa)) {

    rm(phyl6, column_name)
    column_name <- detected_taxa[i]
    try(phyl6 <-phyl4 %>% mutate(rabun = get(column_name)/sum(otu_table[,2])) %>% dplyr::select(rabun, Group) %>% group_by(Group) %>% summarize(mean=mean(rabun), sem=sd(rabun)/sqrt(length(rabun))),T)
    if(exists("phyl6")==F) {
      try(phyl6 <- phyl4 %>% mutate(rabun = get(paste("X",column_name, sep=""))/sum(otu_table[,2])) %>% dplyr::select(rabun, Group) %>% group_by(Group) %>% summarize(mean=mean(rabun), sem=sd(rabun)/sqrt(length(rabun))))
    }

    phyl6
    rm(t)
    for (t in 1:length(extra_cols2)) {
      otu5[i,paste0(extra_cols2[t])] <- phyl6 %>% filter(Group==extra_cols2[t]) %>% dplyr::select(mean) %>% unlist()
      otu5[i,paste0(extra_cols2[t],"_sem")] <- phyl6 %>% filter(Group==extra_cols2[t]) %>% dplyr::select(sem) %>% unlist()
    }

    # assign(paste("p",i,sep="_"),value = ggplot(phyl6, aes(x=Group, y=mean)) +
    #          geom_bar(position=position_dodge(), stat="identity") +
    #          geom_errorbar(aes(ymin=mean+sem, ymax=mean-sem),
    #                        width=.2,                    # Width of the error bars
    #                        position=position_dodge(.9)) +
    #          coord_cartesian(ylim=c((min(phyl6$mean-phyl6$sem)*.9),max(phyl6$mean+phyl6$sem)*1.3)) + #scale_y_continuous(limits=c(.25,.4)) +
    #          theme(axis.text=element_text(size=14),
    #                panel.background = element_blank(),
    #                axis.line = element_line(),
    #                axis.ticks=element_line(),
    #                axis.title=element_text(size=16),
    #                title=element_text(size=13)) +
    #          labs(y="relative abundance",x=var1,title=detected_taxa[i]))
  }

  ## make a list of the plots
  # if(length(detected_taxa)>1) {
  #   plist <- list(p_1)
  #   for(q in 2:length(detected_taxa)) {
  #     plist[[q]] <- get(paste0("p_",q))
  #   }
  # } else if (length(detected_taxa)>0) {
  #   plist <- list(p_1)
  # }

  # otu6 <- read.csv(paste("randancom_",var1,"_",file_path,"_",taxonomic_level,".csv",sep="")) %>% dplyr::select(-X)
  # otu7 <- rbind(otu6,otu5)
  # print(otu5)
  # write.csv(otu7, paste("randancom_",var1,"_",file_path,"_",taxonomic_level,".csv",sep=""))
  write.table(counter, paste("randancom_",var1,"_",file_path,"_",taxonomic_level,".csv",sep=""), sep = "\t",append = T)
  write.table(random_lookup, paste("randancom_",var1,"_",file_path,"_",taxonomic_level,".csv",sep=""), sep = "\t",append = T)
  write.table(otu5, paste("randancom_",var1,"_",file_path,"_",taxonomic_level,".csv",sep=""), sep = "\t",append = T)
  
  # jpeg(h=800*1.25, w=1600*1.25, paste("ancom_",var1,"_",file_path,"_",taxonomic_level,".jpg",sep=""), units = "px", quality = 0.9)
  # do.call("grid.arrange", c(plist))
  # dev.off()
}

make_ancom2_plots_jafari <- function (file_path, map_path="", mapper_file=mapper_file, taxonomic_level=taxonomic_level,var1=var1,var2=var2,var3=var3, var4 = var4, correction_level=3,newcol="newcol", the_id = "X.SampleID", main.var=main.var, adj.formula=adj.formula, repeat.var = repeat.var, multcorr=multcorr, sig=sig, prev.cut=prev.cut, Group = Group, random.formula=NULL) {
  
  rm(otu_table, ref_names, otu2, map2, phyl2, phyl3, col1drop, col2drop, phyl4, ancom.OTU, detected_taxa, otu3, otu4, otu5, phyl5, plist, row_list)
  
  otu_table <- read.table(paste('core-metrics-results-',file_path,'/rarefied_table.txt',sep=""), comment.char="", header=T, sep="\t", fill=T, skip=1) %>% 
    left_join(read.csv(paste0('taxonomy',map_path,'/taxonomy_forR.csv')), by=c("X.OTU.ID"="Feature.ID"))
  
  ref_names <- subset(colnames(otu_table), (colnames(otu_table)%in%c("X.OTU.ID","kingdom","phylum","class","order","family","genus","species")==F))
  
  #  taxonomic_level <- "X.OTU.ID"
  ## Cluster rows by taxonomy hierarchically 
  taxon_vector <- c("kingdom","phylum","class","order","family","genus","species","X.OTU.ID")
  taxon_number <- which(taxon_vector==taxonomic_level)
  taxon_vector2 <- c()
  for (i in 1:taxon_number) {
    taxon_vector2 <- c(taxon_vector2,taxon_vector[i])
  }
  
  if(taxonomic_level!="X.OTU.ID") {
    otu2 <- otu_table[,c("X.OTU.ID",taxon_vector2)] %>% 
      tidyr::unite(clustered_taxonomy, (taxon_vector[1:length(taxon_vector2)])) %>%
      mutate(clustered_taxonomy = as.character(clustered_taxonomy)) %>%
      mutate(clustered_taxonomy = gsub(pattern = " ",replacement = "",x = clustered_taxonomy)) %>%
      mutate(clustered_taxonomy = gsub(pattern = "^__",replacement = "__unassigned",x = clustered_taxonomy)) %>%
      mutate(clustered_taxonomy = gsub(pattern = "___",replacement = "__unassigned_",x = clustered_taxonomy)) %>%
      mutate(clustered_taxonomy = gsub(pattern = "__$",replacement = "__unassigned",x = clustered_taxonomy))
  } else if (taxonomic_level == "X.OTU.ID") {
    otu2 <- otu_table[,c("X.OTU.ID",taxon_vector2)] %>% 
      mutate(OTUID2 = X.OTU.ID) %>%
      tidyr::unite(clustered_taxonomy, c(taxon_vector[1:7],"OTUID2")) %>%
      mutate(clustered_taxonomy = as.character(clustered_taxonomy)) %>%
      mutate(clustered_taxonomy = gsub(pattern = " ",replacement = "",x = clustered_taxonomy)) %>%
      mutate(clustered_taxonomy = gsub(pattern = "^__",replacement = "__unassigned",x = clustered_taxonomy)) %>%
      mutate(clustered_taxonomy = gsub(pattern = "___",replacement = "__unassigned_",x = clustered_taxonomy)) %>%
      mutate(clustered_taxonomy = gsub(pattern = "__$",replacement = "__unassigned",x = clustered_taxonomy)) %>%
      dplyr::select(-X.OTU.ID.1)
    
  } 
  
  map2 <- read.table(mapper_file,comment.char = "", header=T, fill=T,sep="\t") %>% 
    dplyr::select_(.dots = list(the_id,var1,var2,var3,var4, Group))
  
  phyl2 <- otu_table %>% 
    inner_join(otu2, by = "X.OTU.ID") %>%
    group_by(clustered_taxonomy) %>% 
    summarize_at(ref_names,sum, na.rm=T)
  rownames(phyl2) <- as.character(unlist(phyl2$clustered_taxonomy))
  
  phyl3 <- phyl2 %>% dplyr::select(-clustered_taxonomy) %>% t() %>% data.frame() 
  colnames(phyl3) <- rownames(phyl2)
  phyl3 <- mutate(phyl3, X.SampleID=rownames(phyl3))
  phyl3$Sample.ID <- gsub("\\.","_", phyl3$X.SampleID)
  phyl3$Sample.ID <- gsub("X","s", phyl3$Sample.ID)
  map2$Sample.ID <- map2$X.SampleID

  phyl5 <- phyl3 %>% dplyr::select(Sample.ID, everything()) %>% dplyr::select(-X.SampleID) %>% mutate()
  colnames(phyl5)[1] <- "Sample.ID"
  
  map3 <- map2 %>% dplyr::select(Sample.ID, everything()) %>% dplyr::select(-X.SampleID) %>% filter(Sample.ID %in% phyl5$Sample.ID)
  colnames(map3)[1] <- "Sample.ID"
  
  phyl4 <- phyl3 %>% inner_join(map2, by=c("Sample.ID")) %>% mutate(Group = get(Group)) %>% droplevels()#%>% dplyr::select_(.dots = list(col1drop, col2drop, paste0("-",var1))) %>% dplyr::select(-sample2)
  rm(comparison_test)
  
  if (taxonomic_level == "X.OTU.ID") {
    colnames(phyl5) <- paste0("X",colnames(phyl5))
    colnames(phyl5)[1] <- "Sample.ID"
  }
  
  phyl5[1,]
  map3[1,]
  print(dim(phyl5))
  write.csv(phyl5,"phyl5b.csv")
  comparison_test = ANCOM(otu_data = phyl5, 
                          meta_data = map3, 
                          main_var = main.var,  
                          zero_cut = prev.cut, 
                          p_adjust_method = multcorr, 
                          alpha = sig, 
                          adj_formula = adj.formula, 
                          rand_formula = random.formula)
  
  #print(comparison_test)
  write.csv(comparison_test, paste("ancomexcel_",main.var,"_",file_path,"_",taxonomic_level,".csv",sep=""))
  dt1 <- data.frame(comparison_test) %>% filter(detected_0.9==T) %>% dplyr::select(otu_id)
  detected_taxa <- data.frame(OTU=unlist(sapply(X = dt1$otu_id, function(x) ifelse(substring(x,1,1)=="X",substring(x,2),substring(x,1)))))	%>% unlist() %>% unname() %>% as.character()
  
  otu3 <- otu_table[,c("X.OTU.ID","kingdom","phylum","class","order","family","genus","species")]
  otu3$OTU <- as.character(otu3$X.OTU.ID)
  
  if(taxonomic_level!="OTU") {
    otu3[,paste(taxonomic_level)] = gsub(" ","",otu3[,paste(taxonomic_level)])
  }
  
  otu5 <- otu3 %>% inner_join(otu2, by = "X.OTU.ID") %>% filter(clustered_taxonomy%in%detected_taxa) %>% distinct(get(taxonomic_level), .keep_all = T)
  
  extra_cols2 <- names(table(list(phyl4$Group)))
  extra_cols <- unname(c(extra_cols2,sapply(extra_cols2, function(x) paste0(x,"_sem"))))
  otu5[extra_cols] <- "NA"
  
  phyl3$sample2 <- gsub("\\.","", phyl3$X.SampleID)
  phyl3$sample2 <- gsub("_","", phyl3$sample2)
  
  for (i in 1:length(detected_taxa)) {
    
    rm(phyl6, column_name)
    column_name <- detected_taxa[i]
    try(phyl6 <-phyl4 %>% mutate(rabun = get(column_name)/sum(otu_table[,2])) %>% dplyr::select(rabun, Group) %>% group_by(Group) %>% summarize(mean=mean(rabun), sem=sd(rabun)/sqrt(length(rabun))),T)
    if(exists("phyl6")==F) {
      try(phyl6 <- phyl4 %>% mutate(rabun = get(paste("X",column_name, sep=""))/sum(otu_table[,2])) %>% dplyr::select(rabun, Group) %>% group_by(Group) %>% summarize(mean=mean(rabun), sem=sd(rabun)/sqrt(length(rabun))))
    }
    
    phyl6
    rm(t)
    for (t in 1:length(extra_cols2)) {
      otu5[i,paste0(extra_cols2[t])] <- phyl6 %>% filter(Group==extra_cols2[t]) %>% dplyr::select(mean) %>% unlist()
      otu5[i,paste0(extra_cols2[t],"_sem")] <- phyl6 %>% filter(Group==extra_cols2[t]) %>% dplyr::select(sem) %>% unlist()
    }
    
    assign(paste("p",i,sep="_"),value = ggplot(phyl6, aes(x=Group, y=mean)) + 
             geom_bar(position=position_dodge(), stat="identity") +
             geom_errorbar(aes(ymin=mean+sem, ymax=mean-sem),
                           width=.2,                    # Width of the error bars
                           position=position_dodge(.9)) +
             coord_cartesian(ylim=c((min(phyl6$mean-phyl6$sem)*.9),max(phyl6$mean+phyl6$sem)*1.3)) + #scale_y_continuous(limits=c(.25,.4)) + 
             theme(axis.text=element_text(size=14), 
                   panel.background = element_blank(),
                   axis.line = element_line(), 
                   axis.ticks=element_line(), 
                   axis.title=element_text(size=16),
                   title=element_text(size=13)) +
             labs(y="relative abundance",x=var1,title=detected_taxa[i]))
  }
  
  ## make a list of the plots
  if(length(detected_taxa)>1) {
    plist <- list(p_1)
    for(q in 2:length(detected_taxa)) {
      plist[[q]] <- get(paste0("p_",q))
    }
  } else if (length(detected_taxa)>0) {
    plist <- list(p_1)
  }
  
  write.csv(otu5, paste("ancom_",var1,"_",file_path,"_",taxonomic_level,".csv",sep=""))
  
  jpeg(h=800*1.25, w=1600*1.25, paste("ancom_",var1,"_",file_path,"_",taxonomic_level,".jpg",sep=""), units = "px", quality = 0.9)
  do.call("grid.arrange", c(plist))
  dev.off()
}

make_ancom2_plots_seasonal <- function(file_path, map_path="", mapper_file=mapper_file, taxonomic_level=taxonomic_level,var1=var1,var2=var2,var3=var3, var4 = var4,newcol="newcol", the_id = "X.SampleID", main.var=main.var, adj.formula=adj.formula, repeat.var = repeat.var, multcorr=multcorr, sig=sig, prev.cut=prev.cut, Group = Group, random.formula=NULL) {
  
  rm(otu_table, ref_names, otu2, map2, phyl2, phyl3, col1drop, col2drop, phyl4, ancom.OTU, detected_taxa, otu3, otu4, otu5, phyl5, plist, row_list)
  
  otu_table <- read.table(paste('core-metrics-results-',file_path,'/rarefied_table.txt',sep=""), comment.char="", header=T, sep="\t", fill=T, skip=1) %>% 
    left_join(read.csv(paste0('taxonomy',map_path,'/taxonomy_forR.csv')), by=c("X.OTU.ID"="Feature.ID"))
  
  ref_names <- subset(colnames(otu_table), (colnames(otu_table)%in%c("X.OTU.ID","kingdom","phylum","class","order","family","genus","species")==F))
  
  #  otu_table$X.OTU.ID <- as.factor(otu_table$X.OTU.ID)
  
  #  taxonomic_level <- "X.OTU.ID"
  ## Cluster rows by taxonomy hierarchically 
  taxon_vector <- c("kingdom","phylum","class","order","family","genus","species","X.OTU.ID")
  taxon_number <- which(taxon_vector==taxonomic_level)
  taxon_vector2 <- c()
  for (i in 1:taxon_number) {
    taxon_vector2 <- c(taxon_vector2,taxon_vector[i])
  }
  
  if(taxonomic_level!="X.OTU.ID") {
    otu2 <- otu_table[,c("X.OTU.ID",taxon_vector2)] %>% 
      tidyr::unite(clustered_taxonomy, (taxon_vector[1:length(taxon_vector2)])) %>%
      mutate(clustered_taxonomy = as.character(clustered_taxonomy)) %>%
      mutate(clustered_taxonomy = gsub(pattern = " ",replacement = "",x = clustered_taxonomy)) %>%
      mutate(clustered_taxonomy = gsub(pattern = "^__",replacement = "__unassigned",x = clustered_taxonomy)) %>%
      mutate(clustered_taxonomy = gsub(pattern = "___",replacement = "__unassigned_",x = clustered_taxonomy)) %>%
      mutate(clustered_taxonomy = gsub(pattern = "__$",replacement = "__unassigned",x = clustered_taxonomy))
  } else if (taxonomic_level == "X.OTU.ID") {
    otu2 <- otu_table[,c("X.OTU.ID",taxon_vector2)] %>% 
      mutate(OTUID2 = X.OTU.ID) %>%
      tidyr::unite(clustered_taxonomy, c(taxon_vector[1:7],"OTUID2")) %>%
      mutate(clustered_taxonomy = as.character(clustered_taxonomy)) %>%
      mutate(clustered_taxonomy = gsub(pattern = " ",replacement = "",x = clustered_taxonomy)) %>%
      mutate(clustered_taxonomy = gsub(pattern = "^__",replacement = "__unassigned",x = clustered_taxonomy)) %>%
      mutate(clustered_taxonomy = gsub(pattern = "___",replacement = "__unassigned_",x = clustered_taxonomy)) %>%
      mutate(clustered_taxonomy = gsub(pattern = "__$",replacement = "__unassigned",x = clustered_taxonomy)) %>%
      dplyr::select(-X.OTU.ID.1)
    
  } 
  
  map2 <- read.table(mapper_file,comment.char = "", header=T, fill=T,sep="\t") %>% 
    dplyr::select_(.dots = list(the_id,var1,var2,var3,var4, Group)) %>%
    dplyr::mutate(X.SampleID = paste0("X",X.SampleID))
  
  phyl2 <- otu_table %>% 
    inner_join(otu2, by = "X.OTU.ID") %>%
    group_by(clustered_taxonomy) %>% 
    summarize_at(ref_names,sum, na.rm=T)
  rownames(phyl2) <- as.character(unlist(phyl2$clustered_taxonomy))
  
  phyl3 <- phyl2 %>% dplyr::select(-clustered_taxonomy) %>% t() %>% data.frame() 
  colnames(phyl3) <- rownames(phyl2)
  phyl3 <- mutate(phyl3, X.SampleID=rownames(phyl3))
  phyl3$Sample.ID <- gsub("\\.","", phyl3$X.SampleID)
  phyl3$Sample.ID <- gsub("_","", phyl3$Sample.ID)
  map2$Sample.ID <- gsub("_","",map2$X.SampleID)
  map2$Sample.ID <- gsub("-","",map2$Sample.ID)
  map2$Sample.ID <- gsub("\\.","",map2$Sample.ID) # added 2021-05-05
  
  phyl5 <- phyl3 %>% dplyr::select(Sample.ID, everything()) %>% dplyr::select(-X.SampleID)
  colnames(phyl5)[1] <- "Sample.ID"
  
  map3 <- map2 %>% dplyr::select(Sample.ID, everything()) %>% dplyr::select(-X.SampleID) %>% filter(Sample.ID %in% phyl5$Sample.ID)
  map3$date3 <- as.Date(as.character(map3$date2), "%m/%d/%y")
  
  colnames(map3)[1] <- "Sample.ID"
  
  phyl4 <- phyl3 %>% inner_join(map2, by=c("Sample.ID")) %>% mutate(Group = get(Group)) %>% droplevels()#%>% dplyr::select_(.dots = list(col1drop, col2drop, paste0("-",var1))) %>% dplyr::select(-sample2)
  rm(comparison_test)
  
  if (taxonomic_level == "X.OTU.ID") {
    colnames(phyl5) <- paste0("X",colnames(phyl5))
    colnames(phyl5)[1] <- "Sample.ID"
  }
  
  write.csv(phyl5,"phyl5b.csv")
  
  comparison_test = ANCOM(otu_data = phyl5, 
                          meta_data = map3, 
                          main_var = main.var,  
                          zero_cut = prev.cut, 
                          p_adjust_method = multcorr,
                          alpha = sig, 
                          adj_formula = adj.formula, 
                          rand_formula = random.formula)
  
  
  # print(comparison_test)
  write.csv(comparison_test, paste("ancomexcel_",main.var,"_",file_path,"_",taxonomic_level,".csv",sep=""))
  dt1 <- data.frame(comparison_test) %>% filter(detected_0.9==T) %>% dplyr::select(otu_id)
  detected_taxa <- data.frame(OTU=unlist(sapply(X = dt1$otu_id, function(x) ifelse(substring(x,1,1)=="X",substring(x,2),substring(x,1)))))	%>% unlist() %>% unname() %>% as.character()
  
  otu3 <- otu_table[,c("X.OTU.ID","kingdom","phylum","class","order","family","genus","species")]
  otu3$OTU <- as.character(otu3$X.OTU.ID)
  
  if(taxonomic_level!="OTU") {
    otu3[,paste(taxonomic_level)] = gsub(" ","",otu3[,paste(taxonomic_level)])
  }
  
  otu2$clustered_taxonomy
  otu5 <- otu3 %>% inner_join(otu2, by = "X.OTU.ID") %>% filter(clustered_taxonomy%in%detected_taxa) %>% distinct(get(taxonomic_level), .keep_all = T)
  
  extra_cols2 <- names(table(list(phyl4$Group)))
  extra_cols <- unname(c(extra_cols2,sapply(extra_cols2, function(x) paste0(x,"_sem"))))
  otu5[extra_cols] <- "NA"
  
  phyl3$sample2 <- gsub("\\.","", phyl3$X.SampleID)
  phyl3$sample2 <- gsub("_","", phyl3$sample2)
  
  
  for (i in 1:length(detected_taxa)) {
    
    rm(phyl6, column_name)
    column_name <- detected_taxa[i]
    try(phyl6 <-phyl4 %>% mutate(rabun = get(column_name)/sum(otu_table[,2])) %>% dplyr::select(rabun, Group) %>% group_by(Group) %>% dplyr::summarize(mean=mean(rabun), sem=sd(rabun)/sqrt(length(rabun))),T)
    if(exists("phyl6")==F) {
      try(phyl6 <- phyl4 %>% mutate(rabun = get(paste("X",column_name, sep=""))/sum(otu_table[,2])) %>% dplyr::select(rabun, Group) %>% group_by(Group) %>% dplyr::summarize(mean=mean(rabun), sem=sd(rabun)/sqrt(length(rabun))))
    }
    
    phyl6
    rm(t)
    for (t in 1:length(extra_cols2)) {
      otu5[i,paste0(extra_cols2[t])] <- phyl6 %>% filter(Group==extra_cols2[t]) %>% dplyr::select(mean) %>% unlist()
      otu5[i,paste0(extra_cols2[t],"_sem")] <- phyl6 %>% filter(Group==extra_cols2[t]) %>% dplyr::select(sem) %>% unlist()
    }
    
    assign(paste("p",i,sep="_"),value = ggplot(phyl6, aes(x=Group, y=mean)) + 
             geom_bar(position=position_dodge(), stat="identity") +
             geom_errorbar(aes(ymin=mean+sem, ymax=mean-sem),
                           width=.2,                    # Width of the error bars
                           position=position_dodge(.9)) +
             coord_cartesian(ylim=c((min(phyl6$mean-phyl6$sem)*.9),max(phyl6$mean+phyl6$sem)*1.3)) + #scale_y_continuous(limits=c(.25,.4)) + 
             theme(axis.text=element_text(size=14), 
                   panel.background = element_blank(),
                   axis.line = element_line(), 
                   axis.ticks=element_line(), 
                   axis.title=element_text(size=16),
                   title=element_text(size=13)) +
             labs(y="relative abundance",x=var1,title=detected_taxa[i]))
  }
  
  ## make a list of the plots
  if(length(detected_taxa)>1) {
    plist <- list(p_1)
    for(q in 2:length(detected_taxa)) {
      plist[[q]] <- get(paste0("p_",q))
    }
  } else if (length(detected_taxa)>0) {
    plist <- list(p_1)
  }
  
  write.csv(otu5, paste("ancom_",var1,"_",file_path,"_",taxonomic_level,".csv",sep=""))
  
  jpeg(h=800*1.25, w=1600*1.25, paste("ancom_",var1,"_",file_path,"_",taxonomic_level,".jpg",sep=""), units = "px", quality = 0.9)
  do.call("grid.arrange", c(plist))
  dev.off()
}

# 0.0.2 Aug 2019 make compatible with variable taxonomy and mapper names

make_ancom_plots_0.0.2 <- function (file_path, map_path="", mapper_file=mapper_file, taxonomic_level="X.OTU.ID", name_levels=taxonomic_level,var1,var2, correction_level=3,newcol="newcol", the_id = "X.SampleID") {
  
  rm(otu_table, ref_names, otu2, map2, phyl2, phyl3, col1drop, col2drop, phyl4, ancom.OTU, detected_taxa, otu3, otu4, otu5, phyl5, plist, row_list)
  otu_table <- read.table(paste('core-metrics-results-',file_path,'/rarefied_table.txt',sep=""), comment.char="", header=T, sep="\t", fill=T, skip=1) %>% left_join(read.csv(paste0('taxonomy',map_path,'/taxonomy_forR.csv')), by=c("X.OTU.ID"="Feature.ID"))
  #otu_table$OTU <- as.character(otu_table$X.OTU.ID)
  ref_names <- subset(colnames(otu_table), (colnames(otu_table)%in%c("X.OTU.ID","kingdom","phylum","class","order","family","genus","species")==F))
  otu2 <- otu_table[,c("X.OTU.ID","kingdom","phylum","class","order","family","genus","species")]
  otu2$OTU <- as.character(otu2$X.OTU.ID)
  
  map2 <- read.table(mapper_file,comment.char = "", header=T, fill=T) %>% dplyr::select_(.dots = list(the_id,var1))
  
  ## hardcoded
  phyl2 <- otu_table %>% group_by_(.dots=taxonomic_level) %>% summarize_at(ref_names,sum, na.rm=T)
  row_list <- as.character(unlist(phyl2[,paste(taxonomic_level)])) 
  rownames(phyl2) <- gsub(" ","", row_list)
  phyl3 <- phyl2 %>% dplyr::select_(.dots=paste0("-",taxonomic_level)) %>% t() %>% data.frame() 
  phyl3 <- mutate(phyl3, X.SampleID=rownames(phyl3))
  phyl3$sample2 <- gsub("\\.","", phyl3$X.SampleID)
  phyl3$sample2 <- gsub("_","", phyl3$sample2)
  map2$sample2 <- gsub("_","",map2$X.SampleID)
  map2$sample2 <- gsub("-","",map2$sample2)
  
  ## OTU - by inocula
  col1drop <- paste0("-",the_id,".x"); col1drop
  col2drop <- paste0("-",the_id,".y"); col2drop
  
  phyl4 <- phyl3 %>% inner_join(map2, by=c("sample2")) %>% mutate(Group=as.factor(as.character(get(var1)))) %>% dplyr::select_(.dots = list(col1drop, col2drop, paste0("-",var1))) %>% dplyr::select(-sample2)
  
  ancom.OTU <- ANCOM( phyl4, sig = 0.05,multcorr = correction_level, repeated=FALSE )
  
  detected_taxa <- data.frame(OTU=unlist(sapply(X = ancom.OTU$detected, function(x) ifelse(substring(x,1,1)=="X",substring(x,2),x))))%>% filter(OTU!="V1")
  detected_taxa$OTU <- as.character(detected_taxa$OTU)
  
  if(taxonomic_level!="OTU") {
    otu2[,paste(taxonomic_level)] = gsub(" ","",otu2[,paste(taxonomic_level)])
  }
  
  ## filter to OTU table by the data frame
  otu3 <- otu2 %>% mutate(newcol = get(name_levels[1]))#%>% dplyr::select_(.dots=c(taxonomic_level)) %>% inner_join(detected_taxa)
  
  otu5 <- otu2 %>% filter(get(taxonomic_level)%in%detected_taxa$OTU) %>% distinct(get(taxonomic_level), .keep_all = T)
  
  extra_cols2 <- names(table(list(phyl4$Group)))
  extra_cols <- unname(c(extra_cols2,sapply(extra_cols2, function(x) paste0(x,"_sem"))))
  otu5[extra_cols] <- "NA"
  
  otu5
  if(length(name_levels)>1) {
    for(i in 2:length(name_levels)) {
      otu3 <- otu3 %>% mutate(newcol = paste(newcol, get(name_levels[i]), sep="_"))
    }
  }
  
  phyl3$sample2 <- gsub("\\.","", phyl3$X.SampleID)
  phyl3$sample2 <- gsub("_","", phyl3$sample2)
  
  otu3[,paste0(taxonomic_level)]
  otu4 <- otu3 %>% mutate(!!taxonomic_level:=gsub("\\[",".", get(taxonomic_level))) %>% mutate(!!taxonomic_level:=gsub("\\]",".", get(taxonomic_level))) %>% dplyr::select_(.dots=c(taxonomic_level, newcol)) %>% inner_join(detected_taxa, by = structure(names=taxonomic_level, "OTU")) %>% distinct_(.dots=taxonomic_level, .keep_all=T)
  otu3 <- otu4
  
  rm(i)
  ## run the filtered taxa to make plots
  for (i in 1:dim(otu3)[1]) {
    
    rm(phyl5, column_name)
    column_name <- otu3[,taxonomic_level][i]
    try(phyl5 <-phyl4 %>% mutate(rabun = get(column_name)/sum(otu_table[,2])) %>% dplyr::select(rabun, Group) %>% group_by(Group) %>% summarize(mean=mean(rabun), sem=sd(rabun)/sqrt(length(rabun))),T)
    if(exists("phyl5")==F) {
      try(phyl5 <- phyl4 %>% mutate(rabun = get(paste("X",column_name, sep=""))/sum(otu_table[,2])) %>% dplyr::select(rabun, Group) %>% group_by(Group) %>% summarize(mean=mean(rabun), sem=sd(rabun)/sqrt(length(rabun))))
    }
    
    rm(t)
    for (t in 1:length(extra_cols2)) {
      otu5[i,paste0(extra_cols2[t])] <- phyl5 %>% filter(Group==extra_cols2[t]) %>% dplyr::select(mean) %>% unlist()
      otu5[i,paste0(extra_cols2[t],"_sem")] <- phyl5 %>% filter(Group==extra_cols2[t]) %>% dplyr::select(sem) %>% unlist()
    }
    
    assign(paste("p",i,sep="_"),value = ggplot(phyl5, aes(x=Group, y=mean)) + 
             geom_bar(position=position_dodge(), stat="identity") +
             geom_errorbar(aes(ymin=mean+sem, ymax=mean-sem),
                           width=.2,                    # Width of the error bars
                           position=position_dodge(.9)) +
             coord_cartesian(ylim=c((min(phyl5$mean-phyl5$sem)*.9),max(phyl5$mean+phyl5$sem)*1.3)) + #scale_y_continuous(limits=c(.25,.4)) + 
             theme(axis.text=element_text(size=14), 
                   panel.background = element_blank(),
                   axis.line = element_line(), 
                   axis.ticks=element_line(), 
                   axis.title=element_text(size=16),
                   title=element_text(size=13)) +
             labs(y="relative abundance",x=var1,title=otu3$newcol[i])) 
  }
  
  ## make a list of the plots
  if(dim(detected_taxa)[1]>1) {
    plist <- list(p_1)
    for(q in 2:dim(detected_taxa)[1]) {
      plist[[q]] <- get(paste0("p_",q))
    }
  } else if (dim(detected_taxa)[1]>0) {
    plist <- list(p_1)
  }
  write.csv(otu5, paste("ancom_",var1,"_",file_path,"_",taxonomic_level,".csv",sep=""))
  
  
  ##plot them out
  # tiff(h=800*1.25, w=1600*1.25, paste("ancom_",var1,"_",file_path,"_",taxonomic_level,".tiff",sep=""), units = "px", pointsize = 12, compression = "none", res = 100)
  # do.call(grid.arrange, c(plist))
  # dev.off()
  
  jpeg(h=800*1.25, w=1600*1.25, paste("ancom_",var1,"_",file_path,"_",taxonomic_level,".jpg",sep=""), units = "px", quality = 0.9)
  do.call(grid.arrange, c(plist))
  dev.off()
}

make_ancom_covar_plots_0.0.2 <- function (file_path, map_path="", mapper_file=mapper_file, taxonomic_level=taxonomic_level,var1=var1,var2=var2,var3=var3 ,correction_level=3,newcol="newcol", the_id = "X.SampleID", adjusted = adjusted, repeated = repeated, main.var=main.var, adj.formula=adj.formula, repeat.var = repeat.var, longitudinal = latitudinal, multcorr=multcorr, sig=sig, prev.cut=prev.cut, Group = Group, name_levels = name_levels, random.formula=NULL) {
  
  rm(otu_table, ref_names, otu2, map2, phyl2, phyl3, col1drop, col2drop, phyl4, ancom.OTU, detected_taxa, otu3, otu4, otu5, phyl5, plist, row_list)
  otu_table <- read.table(paste('core-metrics-results-',file_path,'/rarefied_table.txt',sep=""), comment.char="", header=T, sep="\t", fill=T, skip=1) %>% left_join(read.csv(paste0('taxonomy',map_path,'/taxonomy_forR.csv')), by=c("X.OTU.ID"="Feature.ID"))
  #otu_table$OTU <- as.character(otu_table$X.OTU.ID)
  ref_names <- subset(colnames(otu_table), (colnames(otu_table)%in%c("X.OTU.ID","kingdom","phylum","class","order","family","genus","species")==F))
  otu2 <- otu_table[,c("X.OTU.ID","kingdom","phylum","class","order","family","genus","species")]
  otu2$OTU <- as.character(otu2$X.OTU.ID)
  
  map2 <- read.table(mapper_file,comment.char = "", header=T, fill=T,sep="\t") %>% dplyr::select_(.dots = list(the_id,var1,var2,var3,Group))
  map2[1,]
  ## hardcoded
  phyl2 <- otu_table %>% group_by_(.dots=taxonomic_level) %>% summarize_at(ref_names,sum, na.rm=T)
  row_list <- as.character(unlist(phyl2[,paste(taxonomic_level)])) 
  rownames(phyl2) <- gsub(" ","", row_list)
  for (i in 1:length(rownames(phyl2))) {
    if (rownames(phyl2)[i] == "" ) {
      print(i)
      rownames(phyl2)[i] <- paste0("o__unassigned",i)
    }
  }
  phyl3 <- phyl2 %>% dplyr::select_(.dots=paste0("-",taxonomic_level)) %>% t() %>% data.frame() 
  colnames(phyl3) <- rownames(phyl2)
  phyl3 <- mutate(phyl3, X.SampleID=rownames(phyl3))
  phyl3$Sample.ID <- gsub("\\.","", phyl3$X.SampleID)
  phyl3$Sample.ID <- gsub("_","", phyl3$Sample.ID)
  map2$Sample.ID <- gsub("_","",map2$X.SampleID)
  map2$Sample.ID <- gsub("-","",map2$Sample.ID)
  
  ## OTU - by inocula
  # col1drop <- paste0("-",the_id,".x"); col1drop
  # col2drop <- paste0("-",the_id,".y"); col2drop
  
  phyl5 <- phyl3 %>% dplyr::select(Sample.ID, everything()) %>% dplyr::select(-X.SampleID)
  colnames(phyl5)[1] <- "Sample.ID"
  
  map3 <- map2 %>% dplyr::select(Sample.ID, everything()) %>% dplyr::select(-X.SampleID) %>% filter(Sample.ID %in% phyl5$Sample.ID)
  colnames(map3)[1] <- "Sample.ID"
  
  phyl4 <- phyl3 %>% inner_join(map2, by=c("Sample.ID")) %>% mutate(Group = get(Group)) %>% droplevels()#%>% dplyr::select_(.dots = list(col1drop, col2drop, paste0("-",var1))) %>% dplyr::select(-sample2)
  rm(comparison_test)
  
  if (taxonomic_level == "X.OTU.ID") {
    colnames(phyl5) <- paste0("X",colnames(phyl5))
    colnames(phyl5)[1] <- "Sample.ID"
  }
  
  phyl5[1,]
  map3[1,]
  # comparison_test=ANCOM.main(OTUdat=phyl5,
  #                            Vardat=map3,
  #                            adjusted=T,
  #                            repeated=F,
  #                            main.var="genotype2",
  #                            adj.formula="time",
  #                            repeat.var=NULL,
  #                            longitudinal = FALSE,
  #                            multcorr=3,
  #                            random.formula="~1|genotype2/cage",
  #                            sig=0.05,
  #                            prev.cut=.5)
  # 
  comparison_test=ANCOM.main(OTUdat=phyl5,
                             Vardat=map3,
                             adjusted=adjusted,
                             repeated=repeated,
                             main.var=main.var,
                             adj.formula=adj.formula,
                             repeat.var=repeat.var,
                             longitudinal = longitudinal,
                             multcorr=multcorr,
                             random.formula=random.formula,
                             sig=sig,
                             prev.cut=prev.cut)
  
  write.csv(comparison_test$W.taxa, paste("ancomexcel_",var1,"_",file_path,"_",taxonomic_level,".csv",sep=""))
  dt1 <- data.frame(comparison_test$W.taxa) %>% filter(detected_0.9==T) %>% dplyr::select(otu.names)
  detected_taxa <- data.frame(OTU=unlist(sapply(X = dt1$otu.names, function(x) ifelse(substring(x,1,1)=="X",substring(x,2),substring(x,1)))))	%>% unlist() %>% unname() %>% as.character()
  
  if(taxonomic_level!="OTU") {
    otu2[,paste(taxonomic_level)] = gsub(" ","",otu2[,paste(taxonomic_level)])
  }
  
  ## filter to OTU table by the data frame
  otu3 <- otu2 %>% mutate(newcol = get(name_levels[1]))#%>% dplyr::select_(.dots=c(taxonomic_level)) %>% inner_join(detected_taxa)
  
  otu5 <- otu2 %>% filter(get(taxonomic_level)%in%detected_taxa) %>% distinct(get(taxonomic_level), .keep_all = T)
  
  extra_cols2 <- names(table(list(phyl4$Group)))
  extra_cols <- unname(c(extra_cols2,sapply(extra_cols2, function(x) paste0(x,"_sem"))))
  otu5[extra_cols] <- "NA"
  
  if(length(name_levels)>1) {
    for(i in 2:length(name_levels)) {
      otu3 <- otu3 %>% mutate(newcol = paste(newcol, get(name_levels[i]), sep="_"))
    }
  }
  
  phyl3$sample2 <- gsub("\\.","", phyl3$X.SampleID)
  phyl3$sample2 <- gsub("_","", phyl3$sample2)
  
  otu4 <- otu3 %>% mutate(!!taxonomic_level:=gsub("\\[",".", get(taxonomic_level))) %>% mutate(!!taxonomic_level:=gsub("\\]",".", get(taxonomic_level))) %>% dplyr::select_(.dots=c(taxonomic_level, newcol)) %>% filter(get(taxonomic_level) %in% detected_taxa) %>% distinct_(.dots=taxonomic_level, .keep_all=T)
  
  ## run the filtered taxa to make plots
  # i=1
  # for (i in 1:dim(otu3)[1]) {
  # 	
  # 	rm(phyl6, column_name)
  # 	column_name <- otu3[,taxonomic_level][i]
  # 	try(phyl6 <-phyl4 %>% mutate(rabun = get(column_name)/sum(otu_table[,2])) %>% dplyr::select(rabun, Group) %>% group_by(Group) %>% summarize(mean=mean(rabun), sem=sd(rabun)/sqrt(length(rabun))),T)
  # 	if(exists("phyl6")==F) {
  # 		try(phyl6 <- phyl4 %>% mutate(rabun = get(paste("X",column_name, sep=""))/sum(otu_table[,2])) %>% dplyr::select(rabun, Group) %>% group_by(Group) %>% summarize(mean=mean(rabun), sem=sd(rabun)/sqrt(length(rabun))))
  # 	}
  # 	
  # 	phyl6
  # 	rm(t)
  # 	for (t in 1:length(extra_cols2)) {
  # 		otu5[i,paste0(extra_cols2[t])] <- phyl6 %>% filter(Group==extra_cols2[t]) %>% dplyr::select(mean) %>% unlist()
  # 		otu5[i,paste0(extra_cols2[t],"_sem")] <- phyl6 %>% filter(Group==extra_cols2[t]) %>% dplyr::select(sem) %>% unlist()
  # 	}
  # 	
  # 	assign(paste("p",i,sep="_"),value = ggplot(phyl6, aes(x=Group, y=mean)) + 
  # 				 	geom_bar(position=position_dodge(), stat="identity") +
  # 				 	geom_errorbar(aes(ymin=mean+sem, ymax=mean-sem),
  # 				 								width=.2,                    # Width of the error bars
  # 				 								position=position_dodge(.9)) +
  # 				 	coord_cartesian(ylim=c((min(phyl6$mean-phyl6$sem)*.9),max(phyl6$mean+phyl6$sem)*1.3)) + #scale_y_continuous(limits=c(.25,.4)) + 
  # 				 	theme(axis.text=element_text(size=14), 
  # 				 				panel.background = element_blank(),
  # 				 				axis.line = element_line(), 
  # 				 				axis.ticks=element_line(), 
  # 				 				axis.title=element_text(size=16),
  # 				 				title=element_text(size=13)) +
  # 				 	labs(y="relative abundance",x=var1,title=otu3$newcol[i])) 
  # }
  # 
  # 
  
  
  otu3[1,]
  
  
  for (i in 1:length(detected_taxa)) {
    
    rm(phyl6, column_name)
    column_name <- detected_taxa[i]
    try(phyl6 <-phyl4 %>% mutate(rabun = get(column_name)/sum(otu_table[,2])) %>% dplyr::select(rabun, Group) %>% group_by(Group) %>% summarize(mean=mean(rabun), sem=sd(rabun)/sqrt(length(rabun))),T)
    if(exists("phyl6")==F) {
      try(phyl6 <- phyl4 %>% mutate(rabun = get(paste("X",column_name, sep=""))/sum(otu_table[,2])) %>% dplyr::select(rabun, Group) %>% group_by(Group) %>% summarize(mean=mean(rabun), sem=sd(rabun)/sqrt(length(rabun))))
    }
    
    phyl6
    rm(t)
    for (t in 1:length(extra_cols2)) {
      otu5[i,paste0(extra_cols2[t])] <- phyl6 %>% filter(Group==extra_cols2[t]) %>% dplyr::select(mean) %>% unlist()
      otu5[i,paste0(extra_cols2[t],"_sem")] <- phyl6 %>% filter(Group==extra_cols2[t]) %>% dplyr::select(sem) %>% unlist()
    }
    
    assign(paste("p",i,sep="_"),value = ggplot(phyl6, aes(x=Group, y=mean)) + 
             geom_bar(position=position_dodge(), stat="identity") +
             geom_errorbar(aes(ymin=mean+sem, ymax=mean-sem),
                           width=.2,                    # Width of the error bars
                           position=position_dodge(.9)) +
             coord_cartesian(ylim=c((min(phyl6$mean-phyl6$sem)*.9),max(phyl6$mean+phyl6$sem)*1.3)) + #scale_y_continuous(limits=c(.25,.4)) + 
             theme(axis.text=element_text(size=14), 
                   panel.background = element_blank(),
                   axis.line = element_line(), 
                   axis.ticks=element_line(), 
                   axis.title=element_text(size=16),
                   title=element_text(size=13)) +
             labs(y="relative abundance",x=var1,title=	otu3 %>% filter(get(taxonomic_level)== paste0(column_name)) %>% dplyr::select(newcol) %>% unlist() %>% unname()))
  }
  
  p_1
  ## make a list of the plots
  if(length(detected_taxa)>1) {
    plist <- list(p_1)
    for(q in 2:length(detected_taxa)) {
      plist[[q]] <- get(paste0("p_",q))
    }
  } else if (length(detected_taxa)>0) {
    plist <- list(p_1)
  }
  
  write.csv(otu5, paste("ancom_",var1,"_",file_path,"_",taxonomic_level,".csv",sep=""))
  
  
  ##plot them out
  # tiff(h=800*1.25, w=1600*1.25, paste("ancom_",var1,"_",file_path,"_",taxonomic_level,".tiff",sep=""), units = "px", pointsize = 12, compression = "none", res = 100)
  # do.call(grid.arrange, c(plist))
  # dev.off()
  
  jpeg(h=800*1.25, w=1600*1.25, paste("ancom_",var1,"_",file_path,"_",taxonomic_level,".jpg",sep=""), units = "px", quality = 0.9)
  do.call(grid.arrange, c(plist))
  dev.off()
}

make_ancom_covar_plots_0.0.2_date <- function (file_path, map_path="", mapper_file=mapper_file, taxonomic_level=taxonomic_level,var1=var1,var2=var2,var3=var3 ,correction_level=3,newcol="newcol", the_id = "X.SampleID", adjusted = adjusted, repeated = repeated, main.var=main.var, adj.formula=adj.formula, repeat.var = repeat.var, longitudinal = latitudinal, multcorr=multcorr, sig=sig, prev.cut=prev.cut, Group = Group, name_levels = name_levels, random.formula=NULL) {
  
  rm(otu_table, ref_names, otu2, map2, phyl2, phyl3, col1drop, col2drop, phyl4, ancom.OTU, detected_taxa, otu3, otu4, otu5, phyl5, plist, row_list)
  otu_table <- read.table(paste('core-metrics-results-',file_path,'/rarefied_table.txt',sep=""), comment.char="", header=T, sep="\t", fill=T, skip=1) %>% left_join(read.csv(paste0('taxonomy',map_path,'/taxonomy_forR.csv')), by=c("X.OTU.ID"="Feature.ID"))
  #otu_table$OTU <- as.character(otu_table$X.OTU.ID)
  ref_names <- subset(colnames(otu_table), (colnames(otu_table)%in%c("X.OTU.ID","kingdom","phylum","class","order","family","genus","species")==F))
  otu2 <- otu_table[,c("X.OTU.ID","kingdom","phylum","class","order","family","genus","species")]
  otu2$OTU <- as.character(otu2$X.OTU.ID)
  
  map2 <- read.table(mapper_file,comment.char = "", header=T, fill=T,sep="\t") %>% dplyr::select_(.dots = list(the_id,var1,var2,var3,Group))
  map2[1,]
  ## hardcoded
  phyl2 <- otu_table %>% group_by_(.dots=taxonomic_level) %>% summarize_at(ref_names,sum, na.rm=T)
  row_list <- as.character(unlist(phyl2[,paste(taxonomic_level)])) 
  rownames(phyl2) <- gsub(" ","", row_list)
  for (i in 1:length(rownames(phyl2))) {
    if (rownames(phyl2)[i] == "" ) {
      print(i)
      rownames(phyl2)[i] <- paste0("o__unassigned",i)
    }
  }
  phyl3 <- phyl2 %>% dplyr::select_(.dots=paste0("-",taxonomic_level)) %>% t() %>% data.frame() 
  colnames(phyl3) <- rownames(phyl2)
  phyl3 <- mutate(phyl3, X.SampleID=rownames(phyl3))
  phyl3$Sample.ID <- gsub("\\.","", phyl3$X.SampleID)
  phyl3$Sample.ID <- gsub("_","", phyl3$Sample.ID)
  map2$Sample.ID <- gsub("_","",map2$X.SampleID)
  map2$Sample.ID <- gsub("-","",map2$Sample.ID)
  
  ## OTU - by inocula
  # col1drop <- paste0("-",the_id,".x"); col1drop
  # col2drop <- paste0("-",the_id,".y"); col2drop
  
  phyl5 <- phyl3 %>% dplyr::select(Sample.ID, everything()) %>% dplyr::select(-X.SampleID)
  colnames(phyl5)[1] <- "Sample.ID"
  
  map3 <- map2 %>% dplyr::select(Sample.ID, everything()) %>% dplyr::select(-X.SampleID) %>% filter(Sample.ID %in% phyl5$Sample.ID)
  colnames(map3)[1] <- "Sample.ID"
  
  map3[1,]
  map3 <- map3 %>% mutate(time = ifelse(time=="base",0,ifelse(time=="T1",1,ifelse(time=="T2",2,ifelse(time=="T3",3,4))))) 
  try(map3 <- map3 %>% mutate(year2 = gsub(pattern = "y", replacement = "", x = year)),T)
  try(map3$year <- as.Date(as.character(map3$year2), "%Y"),T)
  map3[1,]
  
  phyl4 <- phyl3 %>% inner_join(map2, by=c("Sample.ID")) %>% mutate(Group = get(Group)) %>% droplevels()#%>% dplyr::select_(.dots = list(col1drop, col2drop, paste0("-",var1))) %>% dplyr::select(-sample2)
  rm(comparison_test)
  
  if (taxonomic_level == "X.OTU.ID") {
    colnames(phyl5) <- paste0("X",colnames(phyl5))
    colnames(phyl5)[1] <- "Sample.ID"
  }
  
  phyl5[1,]
  map3[1,]
  # comparison_test=ANCOM.main(OTUdat=phyl5,
  #                            Vardat=map3,
  #                            adjusted=T,
  #                            repeated=F,
  #                            main.var="genotype2",
  #                            adj.formula="time",
  #                            repeat.var=NULL,
  #                            longitudinal = FALSE,
  #                            multcorr=3,
  #                            random.formula="~1|genotype2/cage",
  #                            sig=0.05,
  #                            prev.cut=.5)
  # 
  comparison_test=ANCOM.main(OTUdat=phyl5,
                             Vardat=map3,
                             adjusted=adjusted,
                             repeated=repeated,
                             main.var=main.var,
                             adj.formula=adj.formula,
                             repeat.var=repeat.var,
                             longitudinal = longitudinal,
                             multcorr=multcorr,
                             random.formula=random.formula,
                             sig=sig,
                             prev.cut=prev.cut)
  
  write.csv(comparison_test$W.taxa, paste("ancomexcel_",var1,"_",file_path,"_",taxonomic_level,".csv",sep=""))
  dt1 <- data.frame(comparison_test$W.taxa) %>% filter(detected_0.9==T) %>% dplyr::select(otu.names)
  detected_taxa <- data.frame(OTU=unlist(sapply(X = dt1$otu.names, function(x) ifelse(substring(x,1,1)=="X",substring(x,2),substring(x,1)))))	%>% unlist() %>% unname() %>% as.character()
  
  if(taxonomic_level!="OTU") {
    otu2[,paste(taxonomic_level)] = gsub(" ","",otu2[,paste(taxonomic_level)])
  }
  
  ## filter to OTU table by the data frame
  otu3 <- otu2 %>% mutate(newcol = get(name_levels[1]))#%>% dplyr::select_(.dots=c(taxonomic_level)) %>% inner_join(detected_taxa)
  
  otu5 <- otu2 %>% filter(get(taxonomic_level)%in%detected_taxa) %>% distinct(get(taxonomic_level), .keep_all = T)
  
  extra_cols2 <- names(table(list(phyl4$Group)))
  extra_cols <- unname(c(extra_cols2,sapply(extra_cols2, function(x) paste0(x,"_sem"))))
  otu5[extra_cols] <- "NA"
  
  if(length(name_levels)>1) {
    for(i in 2:length(name_levels)) {
      otu3 <- otu3 %>% mutate(newcol = paste(newcol, get(name_levels[i]), sep="_"))
    }
  }
  
  phyl3$sample2 <- gsub("\\.","", phyl3$X.SampleID)
  phyl3$sample2 <- gsub("_","", phyl3$sample2)
  
  otu4 <- otu3 %>% mutate(!!taxonomic_level:=gsub("\\[",".", get(taxonomic_level))) %>% mutate(!!taxonomic_level:=gsub("\\]",".", get(taxonomic_level))) %>% dplyr::select_(.dots=c(taxonomic_level, newcol)) %>% filter(get(taxonomic_level) %in% detected_taxa) %>% distinct_(.dots=taxonomic_level, .keep_all=T)
  
  ## run the filtered taxa to make plots
  # i=1
  # for (i in 1:dim(otu3)[1]) {
  # 	
  # 	rm(phyl6, column_name)
  # 	column_name <- otu3[,taxonomic_level][i]
  # 	try(phyl6 <-phyl4 %>% mutate(rabun = get(column_name)/sum(otu_table[,2])) %>% dplyr::select(rabun, Group) %>% group_by(Group) %>% summarize(mean=mean(rabun), sem=sd(rabun)/sqrt(length(rabun))),T)
  # 	if(exists("phyl6")==F) {
  # 		try(phyl6 <- phyl4 %>% mutate(rabun = get(paste("X",column_name, sep=""))/sum(otu_table[,2])) %>% dplyr::select(rabun, Group) %>% group_by(Group) %>% summarize(mean=mean(rabun), sem=sd(rabun)/sqrt(length(rabun))))
  # 	}
  # 	
  # 	phyl6
  # 	rm(t)
  # 	for (t in 1:length(extra_cols2)) {
  # 		otu5[i,paste0(extra_cols2[t])] <- phyl6 %>% filter(Group==extra_cols2[t]) %>% dplyr::select(mean) %>% unlist()
  # 		otu5[i,paste0(extra_cols2[t],"_sem")] <- phyl6 %>% filter(Group==extra_cols2[t]) %>% dplyr::select(sem) %>% unlist()
  # 	}
  # 	
  # 	assign(paste("p",i,sep="_"),value = ggplot(phyl6, aes(x=Group, y=mean)) + 
  # 				 	geom_bar(position=position_dodge(), stat="identity") +
  # 				 	geom_errorbar(aes(ymin=mean+sem, ymax=mean-sem),
  # 				 								width=.2,                    # Width of the error bars
  # 				 								position=position_dodge(.9)) +
  # 				 	coord_cartesian(ylim=c((min(phyl6$mean-phyl6$sem)*.9),max(phyl6$mean+phyl6$sem)*1.3)) + #scale_y_continuous(limits=c(.25,.4)) + 
  # 				 	theme(axis.text=element_text(size=14), 
  # 				 				panel.background = element_blank(),
  # 				 				axis.line = element_line(), 
  # 				 				axis.ticks=element_line(), 
  # 				 				axis.title=element_text(size=16),
  # 				 				title=element_text(size=13)) +
  # 				 	labs(y="relative abundance",x=var1,title=otu3$newcol[i])) 
  # }
  # 
  # 
  
  
  otu3[1,]
  
  
  for (i in 1:length(detected_taxa)) {
    
    rm(phyl6, column_name)
    column_name <- detected_taxa[i]
    try(phyl6 <-phyl4 %>% mutate(rabun = get(column_name)/sum(otu_table[,2])) %>% dplyr::select(rabun, Group) %>% group_by(Group) %>% summarize(mean=mean(rabun), sem=sd(rabun)/sqrt(length(rabun))),T)
    if(exists("phyl6")==F) {
      try(phyl6 <- phyl4 %>% mutate(rabun = get(paste("X",column_name, sep=""))/sum(otu_table[,2])) %>% dplyr::select(rabun, Group) %>% group_by(Group) %>% summarize(mean=mean(rabun), sem=sd(rabun)/sqrt(length(rabun))))
    }
    
    phyl6
    rm(t)
    for (t in 1:length(extra_cols2)) {
      otu5[i,paste0(extra_cols2[t])] <- phyl6 %>% filter(Group==extra_cols2[t]) %>% dplyr::select(mean) %>% unlist()
      otu5[i,paste0(extra_cols2[t],"_sem")] <- phyl6 %>% filter(Group==extra_cols2[t]) %>% dplyr::select(sem) %>% unlist()
    }
    
    assign(paste("p",i,sep="_"),value = ggplot(phyl6, aes(x=Group, y=mean)) + 
             geom_bar(position=position_dodge(), stat="identity") +
             geom_errorbar(aes(ymin=mean+sem, ymax=mean-sem),
                           width=.2,                    # Width of the error bars
                           position=position_dodge(.9)) +
             coord_cartesian(ylim=c((min(phyl6$mean-phyl6$sem)*.9),max(phyl6$mean+phyl6$sem)*1.3)) + #scale_y_continuous(limits=c(.25,.4)) + 
             theme(axis.text=element_text(size=14), 
                   panel.background = element_blank(),
                   axis.line = element_line(), 
                   axis.ticks=element_line(), 
                   axis.title=element_text(size=16),
                   title=element_text(size=13)) +
             labs(y="relative abundance",x=var1,title=	otu3 %>% filter(get(taxonomic_level)== paste0(column_name)) %>% dplyr::select(newcol) %>% unlist() %>% unname()))
  }
  
  p_1
  ## make a list of the plots
  if(length(detected_taxa)>1) {
    plist <- list(p_1)
    for(q in 2:length(detected_taxa)) {
      plist[[q]] <- get(paste0("p_",q))
    }
  } else if (length(detected_taxa)>0) {
    plist <- list(p_1)
  }
  
  write.csv(otu5, paste("ancom_",var1,"_",file_path,"_",taxonomic_level,".csv",sep=""))
  
  
  ##plot them out
  # tiff(h=800*1.25, w=1600*1.25, paste("ancom_",var1,"_",file_path,"_",taxonomic_level,".tiff",sep=""), units = "px", pointsize = 12, compression = "none", res = 100)
  # do.call(grid.arrange, c(plist))
  # dev.off()
  
  jpeg(h=800*1.25, w=1600*1.25, paste("ancom_",var1,"_",file_path,"_",taxonomic_level,".jpg",sep=""), units = "px", quality = 0.9)
  do.call(grid.arrange, c(plist))
  dev.off()
}

make_ancom_covar_plots_0.0.2_dateCOVAR <- function (file_path, map_path="", mapper_file=mapper_file, taxonomic_level=taxonomic_level,var1=var1,var2=var2,var3=var3 ,correction_level=3,newcol="newcol", the_id = "X.SampleID", adjusted = adjusted, repeated = repeated, main.var=main.var, adj.formula=adj.formula, repeat.var = repeat.var, longitudinal = latitudinal, multcorr=multcorr, sig=sig, prev.cut=prev.cut, Group = Group, name_levels = name_levels, random.formula=NULL) {
  
  rm(otu_table, ref_names, otu2, map2, phyl2, phyl3, col1drop, col2drop, phyl4, ancom.OTU, detected_taxa, otu3, otu4, otu5, phyl5, plist, row_list)
  otu_table <- read.table(paste('core-metrics-results-',file_path,'/rarefied_table.txt',sep=""), comment.char="", header=T, sep="\t", fill=T, skip=1) %>% left_join(read.csv(paste0('taxonomy',map_path,'/taxonomy_forR.csv')), by=c("X.OTU.ID"="Feature.ID"))
  #otu_table$OTU <- as.character(otu_table$X.OTU.ID)
  ref_names <- subset(colnames(otu_table), (colnames(otu_table)%in%c("X.OTU.ID","kingdom","phylum","class","order","family","genus","species")==F))
  otu2 <- otu_table[,c("X.OTU.ID","kingdom","phylum","class","order","family","genus","species")]
  otu2$OTU <- as.character(otu2$X.OTU.ID)
  
  map2 <- read.table(mapper_file,comment.char = "", header=T, fill=T,sep="\t") %>% dplyr::select_(.dots = list(the_id,var1,var2,var3,Group))
  map2[1,]
  ## hardcoded
  phyl2 <- otu_table %>% group_by_(.dots=taxonomic_level) %>% summarize_at(ref_names,sum, na.rm=T)
  row_list <- as.character(unlist(phyl2[,paste(taxonomic_level)])) 
  rownames(phyl2) <- gsub(" ","", row_list)
  for (i in 1:length(rownames(phyl2))) {
    if (rownames(phyl2)[i] == "" ) {
      print(i)
      rownames(phyl2)[i] <- paste0("o__unassigned",i)
    }
  }
  phyl3 <- phyl2 %>% dplyr::select_(.dots=paste0("-",taxonomic_level)) %>% t() %>% data.frame() 
  colnames(phyl3) <- rownames(phyl2)
  phyl3 <- mutate(phyl3, X.SampleID=rownames(phyl3))
  phyl3$Sample.ID <- gsub("\\.","", phyl3$X.SampleID)
  phyl3$Sample.ID <- gsub("_","", phyl3$Sample.ID)
  map2$Sample.ID <- gsub("_","",map2$X.SampleID)
  map2$Sample.ID <- gsub("-","",map2$Sample.ID)
  
  ## OTU - by inocula
  # col1drop <- paste0("-",the_id,".x"); col1drop
  # col2drop <- paste0("-",the_id,".y"); col2drop
  
  phyl5 <- phyl3 %>% dplyr::select(Sample.ID, everything()) %>% dplyr::select(-X.SampleID)
  colnames(phyl5)[1] <- "Sample.ID"
  
  map3 <- map2 %>% dplyr::select(Sample.ID, everything()) %>% dplyr::select(-X.SampleID) %>% filter(Sample.ID %in% phyl5$Sample.ID)
  colnames(map3)[1] <- "Sample.ID"
  
  map3[1,]
  map3 <- map3 %>% mutate(time = ifelse(time=="base",0,ifelse(time=="t1",1,ifelse(time=="t2",2,ifelse(time=="t3",3,4))))) 
  
  phyl4 <- phyl3 %>% inner_join(map2, by=c("Sample.ID")) %>% mutate(Group = get(Group)) %>% droplevels()#%>% dplyr::select_(.dots = list(col1drop, col2drop, paste0("-",var1))) %>% dplyr::select(-sample2)
  rm(comparison_test)
  
  if (taxonomic_level == "X.OTU.ID") {
    colnames(phyl5) <- paste0("X",colnames(phyl5))
    colnames(phyl5)[1] <- "Sample.ID"
  }
  
  phyl5[1,]
  map3[1,]
  
  comparison_test=ANCOM.main(OTUdat=phyl5,
                             Vardat=map3,
                             adjusted=adjusted,
                             repeated=repeated,
                             main.var=main.var,
                             adj.formula=adj.formula,
                             repeat.var=repeat.var,
                             longitudinal = longitudinal,
                             multcorr=multcorr,
                             random.formula=random.formula,
                             sig=sig,
                             prev.cut=prev.cut)
  
  write.csv(comparison_test$W.taxa, paste("ancomexcel_",var1,"_",file_path,"_",taxonomic_level,".csv",sep=""))
  dt1 <- data.frame(comparison_test$W.taxa) %>% filter(detected_0.9==T) %>% dplyr::select(otu.names)
  detected_taxa <- data.frame(OTU=unlist(sapply(X = dt1$otu.names, function(x) ifelse(substring(x,1,1)=="X",substring(x,2),substring(x,1)))))	%>% unlist() %>% unname() %>% as.character()
  
  if(taxonomic_level!="OTU") {
    otu2[,paste(taxonomic_level)] = gsub(" ","",otu2[,paste(taxonomic_level)])
  }
  
  ## filter to OTU table by the data frame
  otu3 <- otu2 %>% mutate(newcol = get(name_levels[1]))#%>% dplyr::select_(.dots=c(taxonomic_level)) %>% inner_join(detected_taxa)
  
  otu5 <- otu2 %>% filter(get(taxonomic_level)%in%detected_taxa) %>% distinct(get(taxonomic_level), .keep_all = T)
  
  extra_cols2 <- names(table(list(phyl4$Group)))
  extra_cols <- unname(c(extra_cols2,sapply(extra_cols2, function(x) paste0(x,"_sem"))))
  otu5[extra_cols] <- "NA"
  
  if(length(name_levels)>1) {
    for(i in 2:length(name_levels)) {
      otu3 <- otu3 %>% mutate(newcol = paste(newcol, get(name_levels[i]), sep="_"))
    }
  }
  
  phyl3$sample2 <- gsub("\\.","", phyl3$X.SampleID)
  phyl3$sample2 <- gsub("_","", phyl3$sample2)
  
  otu4 <- otu3 %>% mutate(!!taxonomic_level:=gsub("\\[",".", get(taxonomic_level))) %>% mutate(!!taxonomic_level:=gsub("\\]",".", get(taxonomic_level))) %>% dplyr::select_(.dots=c(taxonomic_level, newcol)) %>% filter(get(taxonomic_level) %in% detected_taxa) %>% distinct_(.dots=taxonomic_level, .keep_all=T)
  
  
  
  
  otu3[1,]
  
  
  for (i in 1:length(detected_taxa)) {
    
    rm(phyl6, column_name)
    column_name <- detected_taxa[i]
    try(phyl6 <-phyl4 %>% mutate(rabun = get(column_name)/sum(otu_table[,2])) %>% dplyr::select(rabun, Group) %>% group_by(Group) %>% summarize(mean=mean(rabun), sem=sd(rabun)/sqrt(length(rabun))),T)
    if(exists("phyl6")==F) {
      try(phyl6 <- phyl4 %>% mutate(rabun = get(paste("X",column_name, sep=""))/sum(otu_table[,2])) %>% dplyr::select(rabun, Group) %>% group_by(Group) %>% summarize(mean=mean(rabun), sem=sd(rabun)/sqrt(length(rabun))))
    }
    
    phyl6
    rm(t)
    for (t in 1:length(extra_cols2)) {
      otu5[i,paste0(extra_cols2[t])] <- phyl6 %>% filter(Group==extra_cols2[t]) %>% dplyr::select(mean) %>% unlist()
      otu5[i,paste0(extra_cols2[t],"_sem")] <- phyl6 %>% filter(Group==extra_cols2[t]) %>% dplyr::select(sem) %>% unlist()
    }
    
    assign(paste("p",i,sep="_"),value = ggplot(phyl6, aes(x=Group, y=mean)) + 
             geom_bar(position=position_dodge(), stat="identity") +
             geom_errorbar(aes(ymin=mean+sem, ymax=mean-sem),
                           width=.2,                    # Width of the error bars
                           position=position_dodge(.9)) +
             coord_cartesian(ylim=c((min(phyl6$mean-phyl6$sem)*.9),max(phyl6$mean+phyl6$sem)*1.3)) + #scale_y_continuous(limits=c(.25,.4)) + 
             theme(axis.text=element_text(size=14), 
                   panel.background = element_blank(),
                   axis.line = element_line(), 
                   axis.ticks=element_line(), 
                   axis.title=element_text(size=16),
                   title=element_text(size=13)) +
             labs(y="relative abundance",x=var1,title=	otu3 %>% filter(get(taxonomic_level)== paste0(column_name)) %>% dplyr::select(newcol) %>% unlist() %>% unname()))
  }
  
  p_1
  ## make a list of the plots
  if(length(detected_taxa)>1) {
    plist <- list(p_1)
    for(q in 2:length(detected_taxa)) {
      plist[[q]] <- get(paste0("p_",q))
    }
  } else if (length(detected_taxa)>0) {
    plist <- list(p_1)
  }
  
  write.csv(otu5, paste("ancom_",var1,"_",file_path,"_",taxonomic_level,".csv",sep=""))
  
  
  ##plot them out
  # tiff(h=800*1.25, w=1600*1.25, paste("ancom_",var1,"_",file_path,"_",taxonomic_level,".tiff",sep=""), units = "px", pointsize = 12, compression = "none", res = 100)
  # do.call(grid.arrange, c(plist))
  # dev.off()
  
  jpeg(h=800*1.25, w=1600*1.25, paste("ancom_",var1,"_",file_path,"_",taxonomic_level,".jpg",sep=""), units = "px", quality = 0.9)
  do.call(grid.arrange, c(plist))
  dev.off()
}

make_ancom_covar_plots_4var <- function (file_path, map_path="", mapper_file=mapper_file, taxonomic_level=taxonomic_level,var1=var1,var2=var2,var3=var3, var4=var4, correction_level=3,newcol="newcol", the_id = "X.SampleID", adjusted = adjusted, repeated = repeated, main.var=main.var, adj.formula=adj.formula, repeat.var = repeat.var, longitudinal = latitudinal, multcorr=multcorr, sig=sig, prev.cut=prev.cut, Group = Group, name_levels = name_levels, random.formula=NULL) {
  
  rm(otu_table, ref_names, otu2, map2, phyl2, phyl3, col1drop, col2drop, phyl4, ancom.OTU, detected_taxa, otu3, otu4, otu5, phyl5, plist, row_list)
  otu_table <- read.table(paste('core-metrics-results-',file_path,'/rarefied_table.txt',sep=""), comment.char="", header=T, sep="\t", fill=T, skip=1) %>% left_join(read.csv(paste0('taxonomy',map_path,'/taxonomy_forR.csv')), by=c("X.OTU.ID"="Feature.ID"))
  #otu_table$OTU <- as.character(otu_table$X.OTU.ID)
  ref_names <- subset(colnames(otu_table), (colnames(otu_table)%in%c("X.OTU.ID","kingdom","phylum","class","order","family","genus","species")==F))
  otu2 <- otu_table[,c("X.OTU.ID","kingdom","phylum","class","order","family","genus","species")]
  otu2$OTU <- as.character(otu2$X.OTU.ID)
  
  map2 <- read.table(mapper_file,comment.char = "", header=T, fill=T,sep="\t") %>% dplyr::select_(.dots = list(the_id,var1,var2,var3,var4, Group))
  map2[1,]
  ## hardcoded
  phyl2 <- otu_table %>% group_by_(.dots=taxonomic_level) %>% summarize_at(ref_names,sum, na.rm=T)
  row_list <- as.character(unlist(phyl2[,paste(taxonomic_level)])) 
  rownames(phyl2) <- gsub(" ","", row_list)
  phyl3 <- phyl2 %>% dplyr::select_(.dots=paste0("-",taxonomic_level)) %>% t() %>% data.frame() 
  phyl3 <- mutate(phyl3, X.SampleID=rownames(phyl3))
  phyl3$sample2 <- gsub("\\.","", phyl3$X.SampleID)
  phyl3$sample2 <- gsub("_","", phyl3$sample2)
  map2$sample2 <- gsub("_","",map2$X.SampleID)
  map2$sample2 <- gsub("-","",map2$sample2)
  
  ## OTU - by inocula
  # col1drop <- paste0("-",the_id,".x"); col1drop
  # col2drop <- paste0("-",the_id,".y"); col2drop
  
  phyl5 <- phyl3 %>% dplyr::select(sample2, everything()) %>% dplyr::select(-X.SampleID)
  colnames(phyl5)[1] <- "Sample.ID"
  
  map3 <- map2 %>% dplyr::select(sample2, everything()) %>% dplyr::select(-X.SampleID)
  colnames(map3)[1] <- "Sample.ID"
  
  phyl4 <- phyl3 %>% inner_join(map2, by=c("sample2")) %>% mutate(Group = get(Group)) %>% droplevels()#%>% dplyr::select_(.dots = list(col1drop, col2drop, paste0("-",var1))) %>% dplyr::select(-sample2)
  rm(comparison_test)
  comparison_test=ANCOM.main(OTUdat=phyl5,
                             Vardat=map3,
                             adjusted=adjusted,
                             repeated=repeated,
                             main.var=main.var,
                             adj.formula=adj.formula,
                             repeat.var=repeat.var,
                             longitudinal = longitudinal,
                             multcorr=multcorr,
                             random.formula=random.formula,
                             sig=sig,
                             prev.cut=prev.cut)
  
  write.csv(comparison_test$W.taxa, paste("ancomexcel_",var1,"_",file_path,"_",taxonomic_level,".csv",sep=""))
  dt1 <- data.frame(comparison_test$W.taxa) %>% filter(detected_0.9==T) %>% dplyr::select(otu.names)
  detected_taxa <- data.frame(OTU=unlist(sapply(X = dt1$otu.names, function(x) ifelse(substring(x,1,1)=="X",substring(x,2),substring(x,1)))))	%>% unlist() %>% unname() %>% as.character()
  
  if(taxonomic_level!="OTU") {
    otu2[,paste(taxonomic_level)] = gsub(" ","",otu2[,paste(taxonomic_level)])
  }
  
  ## filter to OTU table by the data frame
  otu3 <- otu2 %>% mutate(newcol = get(name_levels[1]))#%>% dplyr::select_(.dots=c(taxonomic_level)) %>% inner_join(detected_taxa)
  
  otu5 <- otu2 %>% filter(get(taxonomic_level)%in%detected_taxa) %>% distinct(get(taxonomic_level), .keep_all = T)
  
  extra_cols2 <- names(table(list(phyl4$Group)))
  extra_cols <- unname(c(extra_cols2,sapply(extra_cols2, function(x) paste0(x,"_sem"))))
  otu5[extra_cols] <- "NA"
  
  if(length(name_levels)>1) {
    for(i in 2:length(name_levels)) {
      otu3 <- otu3 %>% mutate(newcol = paste(newcol, get(name_levels[i]), sep="_"))
    }
  }
  
  phyl3$sample2 <- gsub("\\.","", phyl3$X.SampleID)
  phyl3$sample2 <- gsub("_","", phyl3$sample2)
  
  otu4 <- otu3 %>% mutate(!!taxonomic_level:=gsub("\\[",".", get(taxonomic_level))) %>% mutate(!!taxonomic_level:=gsub("\\]",".", get(taxonomic_level))) %>% dplyr::select_(.dots=c(taxonomic_level, newcol)) %>% filter(get(taxonomic_level) %in% detected_taxa) %>% distinct_(.dots=taxonomic_level, .keep_all=T)
  
  for (i in 1:length(detected_taxa)) {
    
    rm(phyl6, column_name)
    column_name <- detected_taxa[i]
    try(phyl6 <-phyl4 %>% mutate(rabun = get(column_name)/sum(otu_table[,2])) %>% dplyr::select(rabun, Group) %>% group_by(Group) %>% summarize(mean=mean(rabun), sem=sd(rabun)/sqrt(length(rabun))),T)
    if(exists("phyl6")==F) {
      try(phyl6 <- phyl4 %>% mutate(rabun = get(paste("X",column_name, sep=""))/sum(otu_table[,2])) %>% dplyr::select(rabun, Group) %>% group_by(Group) %>% summarize(mean=mean(rabun), sem=sd(rabun)/sqrt(length(rabun))))
    }
    
    phyl6
    rm(t)
    for (t in 1:length(extra_cols2)) {
      otu5[i,paste0(extra_cols2[t])] <- phyl6 %>% filter(Group==extra_cols2[t]) %>% dplyr::select(mean) %>% unlist()
      otu5[i,paste0(extra_cols2[t],"_sem")] <- phyl6 %>% filter(Group==extra_cols2[t]) %>% dplyr::select(sem) %>% unlist()
    }
    
    assign(paste("p",i,sep="_"),value = ggplot(phyl6, aes(x=Group, y=mean)) + 
             geom_bar(position=position_dodge(), stat="identity") +
             geom_errorbar(aes(ymin=mean+sem, ymax=mean-sem),
                           width=.2,                    # Width of the error bars
                           position=position_dodge(.9)) +
             coord_cartesian(ylim=c((min(phyl6$mean-phyl6$sem)*.9),max(phyl6$mean+phyl6$sem)*1.3)) + #scale_y_continuous(limits=c(.25,.4)) + 
             theme(axis.text=element_text(size=14), 
                   panel.background = element_blank(),
                   axis.line = element_line(), 
                   axis.ticks=element_line(), 
                   axis.title=element_text(size=16),
                   title=element_text(size=13)) +
             labs(y="relative abundance",x=var1,title=	otu3 %>% filter(get(taxonomic_level)== paste0(column_name)) %>% dplyr::select(newcol) %>% unlist() %>% unname()))
  }
  
  p_1
  ## make a list of the plots
  if(length(detected_taxa)>1) {
    plist <- list(p_1)
    for(q in 2:length(detected_taxa)) {
      plist[[q]] <- get(paste0("p_",q))
    }
  } else if (length(detected_taxa)>0) {
    plist <- list(p_1)
  }
  
  write.csv(otu5, paste("ancom_",var1,"_",file_path,"_",taxonomic_level,".csv",sep=""))
  
  
  ##plot them out
  # tiff(h=800*1.25, w=1600*1.25, paste("ancom_",var1,"_",file_path,"_",taxonomic_level,".tiff",sep=""), units = "px", pointsize = 12, compression = "none", res = 100)
  # do.call(grid.arrange, c(plist))
  # dev.off()
  
  jpeg(h=800*1.25, w=1600*1.25, paste("ancom_",var1,"_",file_path,"_",taxonomic_level,".jpg",sep=""), units = "px", quality = 0.9)
  do.call(grid.arrange, c(plist))
  dev.off()
}

make_ancom_covar_plots_4var_periodname <- function (file_path, map_path="", mapper_file=mapper_file, taxonomic_level=taxonomic_level,var1=var1,var2=var2,var3=var3, var4=var4, correction_level=3,newcol="newcol", the_id = "X.SampleID", adjusted = adjusted, repeated = repeated, main.var=main.var, adj.formula=adj.formula, repeat.var = repeat.var, longitudinal = latitudinal, multcorr=multcorr, sig=sig, prev.cut=prev.cut, Group = Group, name_levels = name_levels, random.formula=NULL) {
  
  rm(otu_table, ref_names, otu2, map2, phyl2, phyl3, col1drop, col2drop, phyl4, ancom.OTU, detected_taxa, otu3, otu4, otu5, phyl5, plist, row_list)
  otu_table <- read.table(paste('core-metrics-results-',file_path,'/rarefied_table.txt',sep=""), comment.char="", header=T, sep="\t", fill=T, skip=1) %>% left_join(read.csv(paste0('taxonomy',map_path,'/taxonomy_forR.csv')), by=c("X.OTU.ID"="Feature.ID"))
  #otu_table$OTU <- as.character(otu_table$X.OTU.ID)
  ref_names <- subset(colnames(otu_table), (colnames(otu_table)%in%c("X.OTU.ID","kingdom","phylum","class","order","family","genus","species")==F))
  otu2 <- otu_table[,c("X.OTU.ID","kingdom","phylum","class","order","family","genus","species")]
  otu2$OTU <- as.character(otu2$X.OTU.ID)
  
  map2 <- read.table(mapper_file,comment.char = "", header=T, fill=T,sep="\t") %>% dplyr::select_(.dots = list(the_id,var1,var2,var3,var4, Group))
  map2[1,]
  ## hardcoded
  phyl2 <- otu_table %>% group_by_(.dots=taxonomic_level) %>% summarize_at(ref_names,sum, na.rm=T)
  row_list <- as.character(unlist(phyl2[,paste(taxonomic_level)])) 
  rownames(phyl2) <- gsub(" ","", row_list)
  rownames(phyl2)[1] <- "o__other"
  phyl3 <- phyl2 %>% dplyr::select_(.dots=paste0("-",taxonomic_level)) %>% t() %>% data.frame() 
  colnames(phyl3) <- rownames(phyl2)
  phyl3 <- phyl3 %>% mutate(X.SampleID=rownames(phyl3))
  phyl3
  phyl3$sample2 <- gsub("X","", phyl3$X.SampleID)
  map2$sample2 <- gsub("_","",map2$X.SampleID)
  map2$sample2 <- gsub("-","",map2$sample2)
  
  phyl5 <- phyl3 %>% dplyr::select(sample2, everything()) %>% dplyr::select(-X.SampleID)
  colnames(phyl5)[1] <- "Sample.ID"
  
  map3 <- map2 %>% dplyr::select(sample2, everything()) %>% dplyr::select(-X.SampleID)
  colnames(map3)[1] <- "Sample.ID"
  phyl5[1]
  map3[1,]
  map3$date2 <- as.Date(as.character(map3$date2), "%m/%d/%y")
  
  phyl4 <- phyl5 %>% inner_join(map2, by=c("Sample.ID")) %>% mutate(Group = get(Group)) %>% droplevels()#%>% dplyr::select_(.dots = list(col1drop, col2drop, paste0("-",var1))) %>% dplyr::select(-sample2)
  rm(comparison_test)
  comparison_test=ANCOM.main(OTUdat=phyl5,
                             Vardat=map3,
                             adjusted=adjusted,
                             repeated=repeated,
                             main.var=main.var,
                             adj.formula=adj.formula,
                             repeat.var=repeat.var,
                             longitudinal = longitudinal,
                             multcorr=multcorr,
                             random.formula=random.formula,
                             sig=sig,
                             prev.cut=prev.cut)
  
  write.csv(comparison_test$W.taxa, paste("ancomexcel_",var1,"_",file_path,"_",taxonomic_level,".csv",sep=""))
  dt1 <- data.frame(comparison_test$W.taxa) %>% filter(detected_0.9==T) %>% dplyr::select(otu.names)
  detected_taxa <- data.frame(OTU=unlist(sapply(X = dt1$otu.names, function(x) ifelse(substring(x,1,1)=="X",substring(x,2),substring(x,1)))))	%>% unlist() %>% unname() %>% as.character()
  
  if(taxonomic_level!="OTU") {
    otu2[,paste(taxonomic_level)] = gsub(" ","",otu2[,paste(taxonomic_level)])
  }
  
  ## filter to OTU table by the data frame
  otu3 <- otu2 %>% mutate(newcol = get(name_levels[1]))#%>% dplyr::select_(.dots=c(taxonomic_level)) %>% inner_join(detected_taxa)
  
  otu5 <- otu2 %>% filter(get(taxonomic_level)%in%detected_taxa) %>% distinct(get(taxonomic_level), .keep_all = T)
  
  extra_cols2 <- names(table(list(phyl4$Group)))
  extra_cols <- unname(c(extra_cols2,sapply(extra_cols2, function(x) paste0(x,"_sem"))))
  otu5[extra_cols] <- "NA"
  
  if(length(name_levels)>1) {
    for(i in 2:length(name_levels)) {
      otu3 <- otu3 %>% mutate(newcol = paste(newcol, get(name_levels[i]), sep="_"))
    }
  }
  
  phyl3$sample2 <- gsub("\\.","", phyl3$X.SampleID)
  phyl3$sample2 <- gsub("_","", phyl3$sample2)
  
  otu4 <- otu3 %>% mutate(!!taxonomic_level:=gsub("\\[",".", get(taxonomic_level))) %>% mutate(!!taxonomic_level:=gsub("\\]",".", get(taxonomic_level))) %>% dplyr::select_(.dots=c(taxonomic_level, newcol)) %>% filter(get(taxonomic_level) %in% detected_taxa) %>% distinct_(.dots=taxonomic_level, .keep_all=T)
  
  for (i in 1:length(detected_taxa)) {
    
    rm(phyl6, column_name)
    column_name <- detected_taxa[i]
    try(phyl6 <-phyl4 %>% mutate(rabun = get(column_name)/sum(otu_table[,2])) %>% dplyr::select(rabun, Group) %>% group_by(Group) %>% summarize(mean=mean(rabun), sem=sd(rabun)/sqrt(length(rabun))),T)
    if(exists("phyl6")==F) {
      try(phyl6 <- phyl4 %>% mutate(rabun = get(paste("X",column_name, sep=""))/sum(otu_table[,2])) %>% dplyr::select(rabun, Group) %>% group_by(Group) %>% summarize(mean=mean(rabun), sem=sd(rabun)/sqrt(length(rabun))))
    }
    
    phyl6
    rm(t)
    for (t in 1:length(extra_cols2)) {
      otu5[i,paste0(extra_cols2[t])] <- phyl6 %>% filter(Group==extra_cols2[t]) %>% dplyr::select(mean) %>% unlist()
      otu5[i,paste0(extra_cols2[t],"_sem")] <- phyl6 %>% filter(Group==extra_cols2[t]) %>% dplyr::select(sem) %>% unlist()
    }
    
    assign(paste("p",i,sep="_"),value = ggplot(phyl6, aes(x=Group, y=mean)) + 
             geom_bar(position=position_dodge(), stat="identity") +
             geom_errorbar(aes(ymin=mean+sem, ymax=mean-sem),
                           width=.2,                    # Width of the error bars
                           position=position_dodge(.9)) +
             coord_cartesian(ylim=c((min(phyl6$mean-phyl6$sem)*.9),max(phyl6$mean+phyl6$sem)*1.3)) + #scale_y_continuous(limits=c(.25,.4)) + 
             theme(axis.text=element_text(size=14), 
                   panel.background = element_blank(),
                   axis.line = element_line(), 
                   axis.ticks=element_line(), 
                   axis.title=element_text(size=16),
                   title=element_text(size=13)) +
             labs(y="relative abundance",x=var1,title=	otu3 %>% filter(get(taxonomic_level)== paste0(column_name)) %>% dplyr::select(newcol) %>% unlist() %>% unname()))
  }
  
  p_1
  ## make a list of the plots
  if(length(detected_taxa)>1) {
    plist <- list(p_1)
    for(q in 2:length(detected_taxa)) {
      plist[[q]] <- get(paste0("p_",q))
    }
  } else if (length(detected_taxa)>0) {
    plist <- list(p_1)
  }
  
  write.csv(otu5, paste("ancom_",var1,"_",file_path,"_",taxonomic_level,".csv",sep=""))
  
  
  ##plot them out
  # tiff(h=800*1.25, w=1600*1.25, paste("ancom_",var1,"_",file_path,"_",taxonomic_level,".tiff",sep=""), units = "px", pointsize = 12, compression = "none", res = 100)
  # do.call(grid.arrange, c(plist))
  # dev.off()
  
  jpeg(h=800*1.25, w=1600*1.25, paste("ancom_",var1,"_",file_path,"_",taxonomic_level,".jpg",sep=""), units = "px", quality = 0.9)
  do.call(grid.arrange, c(plist))
  dev.off()
}


make_ancom_covar_plots_periodname_0.0.2 <- function (file_path, map_path="", mapper_file=mapper_file, taxonomic_level=taxonomic_level,var1=var1,var2=var2,var3=var3 ,correction_level=3,newcol="newcol", the_id = "X.SampleID", adjusted = adjusted, repeated = repeated, main.var=main.var, adj.formula=adj.formula, repeat.var = repeat.var, longitudinal = latitudinal, multcorr=multcorr, sig=sig, prev.cut=prev.cut, Group = Group, name_levels = name_levels, random.formula=NULL) {
  
  rm(otu_table, ref_names, otu2, map2, phyl2, phyl3, col1drop, col2drop, phyl4, ancom.OTU, detected_taxa, otu3, otu4, otu5, phyl5, plist, row_list)
  otu_table <- read.table(paste('core-metrics-results-',file_path,'/rarefied_table.txt',sep=""), comment.char="", header=T, sep="\t", fill=T, skip=1) %>% left_join(read.csv(paste0('taxonomy',map_path,'/taxonomy_forR.csv')), by=c("X.OTU.ID"="Feature.ID"))
  #otu_table$OTU <- as.character(otu_table$X.OTU.ID)
  ref_names <- subset(colnames(otu_table), (colnames(otu_table)%in%c("X.OTU.ID","kingdom","phylum","class","order","family","genus","species")==F))
  otu2 <- otu_table[,c("X.OTU.ID","kingdom","phylum","class","order","family","genus","species")]
  otu2$OTU <- as.character(otu2$X.OTU.ID)
  
  map2 <- read.table(mapper_file,comment.char = "", header=T, fill=T,sep="\t") %>% dplyr::select_(.dots = list(the_id,var1,var2,var3,Group))
  map2[1,]
  ## hardcoded
  phyl2 <- otu_table %>% group_by_(.dots=taxonomic_level) %>% summarize_at(ref_names,sum, na.rm=T)
  row_list <- as.character(unlist(phyl2[,paste(taxonomic_level)])) 
  rownames(phyl2) <- gsub(" ","", row_list)
  phyl3 <- phyl2 %>% dplyr::select_(.dots=paste0("-",taxonomic_level)) %>% t() %>% data.frame() 
  phyl3 <- mutate(phyl3, X.SampleID=rownames(phyl3))
  phyl3$sample2 <- gsub("X","", phyl3$X.SampleID)
  map2$sample2 <- gsub("_","",map2$X.SampleID)
  map2$sample2 <- gsub("-","",map2$sample2)
  
  phyl5 <- phyl3 %>% dplyr::select(sample2, everything()) %>% dplyr::select(-X.SampleID)
  colnames(phyl5)[1] <- "Sample.ID"
  
  map3 <- map2 %>% dplyr::select(sample2, everything()) %>% dplyr::select(-X.SampleID)
  colnames(map3)[1] <- "Sample.ID"
  
  phyl4 <- phyl3 %>% inner_join(map2, by=c("sample2")) %>% mutate(Group = get(Group)) %>% droplevels()#%>% dplyr::select_(.dots = list(col1drop, col2drop, paste0("-",var1))) %>% dplyr::select(-sample2)
  rm(comparison_test)
  comparison_test=ANCOM.main(OTUdat=phyl5,
                             Vardat=map3,
                             adjusted=adjusted,
                             repeated=repeated,
                             main.var=main.var,
                             adj.formula=adj.formula,
                             repeat.var=repeat.var,
                             longitudinal = longitudinal,
                             multcorr=multcorr,
                             random.formula=random.formula,
                             sig=sig,
                             prev.cut=prev.cut)
  
  write.csv(comparison_test$W.taxa, paste("ancomexcel_",var1,"_",file_path,"_",taxonomic_level,".csv",sep=""))
  dt1 <- data.frame(comparison_test$W.taxa) %>% filter(detected_0.9==T) %>% dplyr::select(otu.names)
  detected_taxa <- data.frame(OTU=unlist(sapply(X = dt1$otu.names, function(x) ifelse(substring(x,1,1)=="X",substring(x,2),substring(x,1)))))	%>% unlist() %>% unname() %>% as.character()
  
  if(taxonomic_level!="OTU") {
    otu2[,paste(taxonomic_level)] = gsub(" ","",otu2[,paste(taxonomic_level)])
  }
  
  ## filter to OTU table by the data frame
  otu3 <- otu2 %>% mutate(newcol = get(name_levels[1]))#%>% dplyr::select_(.dots=c(taxonomic_level)) %>% inner_join(detected_taxa)
  
  otu5 <- otu2 %>% filter(get(taxonomic_level)%in%detected_taxa) %>% distinct(get(taxonomic_level), .keep_all = T)
  
  extra_cols2 <- names(table(list(phyl4$Group)))
  extra_cols <- unname(c(extra_cols2,sapply(extra_cols2, function(x) paste0(x,"_sem"))))
  otu5[extra_cols] <- "NA"
  
  if(length(name_levels)>1) {
    for(i in 2:length(name_levels)) {
      otu3 <- otu3 %>% mutate(newcol = paste(newcol, get(name_levels[i]), sep="_"))
    }
  }
  
  phyl3$sample2 <- gsub("\\.","", phyl3$X.SampleID)
  phyl3$sample2 <- gsub("_","", phyl3$sample2)
  
  otu4 <- otu3 %>% mutate(!!taxonomic_level:=gsub("\\[",".", get(taxonomic_level))) %>% mutate(!!taxonomic_level:=gsub("\\]",".", get(taxonomic_level))) %>% dplyr::select_(.dots=c(taxonomic_level, newcol)) %>% filter(get(taxonomic_level) %in% detected_taxa) %>% distinct_(.dots=taxonomic_level, .keep_all=T)
  
  ## run the filtered taxa to make plots
  # i=1
  # for (i in 1:dim(otu3)[1]) {
  # 	
  # 	rm(phyl6, column_name)
  # 	column_name <- otu3[,taxonomic_level][i]
  # 	try(phyl6 <-phyl4 %>% mutate(rabun = get(column_name)/sum(otu_table[,2])) %>% dplyr::select(rabun, Group) %>% group_by(Group) %>% summarize(mean=mean(rabun), sem=sd(rabun)/sqrt(length(rabun))),T)
  # 	if(exists("phyl6")==F) {
  # 		try(phyl6 <- phyl4 %>% mutate(rabun = get(paste("X",column_name, sep=""))/sum(otu_table[,2])) %>% dplyr::select(rabun, Group) %>% group_by(Group) %>% summarize(mean=mean(rabun), sem=sd(rabun)/sqrt(length(rabun))))
  # 	}
  # 	
  # 	phyl6
  # 	rm(t)
  # 	for (t in 1:length(extra_cols2)) {
  # 		otu5[i,paste0(extra_cols2[t])] <- phyl6 %>% filter(Group==extra_cols2[t]) %>% dplyr::select(mean) %>% unlist()
  # 		otu5[i,paste0(extra_cols2[t],"_sem")] <- phyl6 %>% filter(Group==extra_cols2[t]) %>% dplyr::select(sem) %>% unlist()
  # 	}
  # 	
  # 	assign(paste("p",i,sep="_"),value = ggplot(phyl6, aes(x=Group, y=mean)) + 
  # 				 	geom_bar(position=position_dodge(), stat="identity") +
  # 				 	geom_errorbar(aes(ymin=mean+sem, ymax=mean-sem),
  # 				 								width=.2,                    # Width of the error bars
  # 				 								position=position_dodge(.9)) +
  # 				 	coord_cartesian(ylim=c((min(phyl6$mean-phyl6$sem)*.9),max(phyl6$mean+phyl6$sem)*1.3)) + #scale_y_continuous(limits=c(.25,.4)) + 
  # 				 	theme(axis.text=element_text(size=14), 
  # 				 				panel.background = element_blank(),
  # 				 				axis.line = element_line(), 
  # 				 				axis.ticks=element_line(), 
  # 				 				axis.title=element_text(size=16),
  # 				 				title=element_text(size=13)) +
  # 				 	labs(y="relative abundance",x=var1,title=otu3$newcol[i])) 
  # }
  # 
  # 
  
  
  otu3[1,]
  
  
  for (i in 1:length(detected_taxa)) {
    
    rm(phyl6, column_name)
    column_name <- detected_taxa[i]
    try(phyl6 <-phyl4 %>% mutate(rabun = get(column_name)/sum(otu_table[,2])) %>% dplyr::select(rabun, Group) %>% group_by(Group) %>% summarize(mean=mean(rabun), sem=sd(rabun)/sqrt(length(rabun))),T)
    if(exists("phyl6")==F) {
      try(phyl6 <- phyl4 %>% mutate(rabun = get(paste("X",column_name, sep=""))/sum(otu_table[,2])) %>% dplyr::select(rabun, Group) %>% group_by(Group) %>% summarize(mean=mean(rabun), sem=sd(rabun)/sqrt(length(rabun))))
    }
    
    phyl6
    rm(t)
    for (t in 1:length(extra_cols2)) {
      otu5[i,paste0(extra_cols2[t])] <- phyl6 %>% filter(Group==extra_cols2[t]) %>% dplyr::select(mean) %>% unlist()
      otu5[i,paste0(extra_cols2[t],"_sem")] <- phyl6 %>% filter(Group==extra_cols2[t]) %>% dplyr::select(sem) %>% unlist()
    }
    
    assign(paste("p",i,sep="_"),value = ggplot(phyl6, aes(x=Group, y=mean)) + 
             geom_bar(position=position_dodge(), stat="identity") +
             geom_errorbar(aes(ymin=mean+sem, ymax=mean-sem),
                           width=.2,                    # Width of the error bars
                           position=position_dodge(.9)) +
             coord_cartesian(ylim=c((min(phyl6$mean-phyl6$sem)*.9),max(phyl6$mean+phyl6$sem)*1.3)) + #scale_y_continuous(limits=c(.25,.4)) + 
             theme(axis.text=element_text(size=14), 
                   panel.background = element_blank(),
                   axis.line = element_line(), 
                   axis.ticks=element_line(), 
                   axis.title=element_text(size=16),
                   title=element_text(size=13)) +
             labs(y="relative abundance",x=var1,title=	otu3 %>% filter(get(taxonomic_level)== paste0(column_name)) %>% dplyr::select(newcol) %>% unlist() %>% unname()))
  }
  
  p_1
  ## make a list of the plots
  if(length(detected_taxa)>1) {
    plist <- list(p_1)
    for(q in 2:length(detected_taxa)) {
      plist[[q]] <- get(paste0("p_",q))
    }
  } else if (length(detected_taxa)>0) {
    plist <- list(p_1)
  }
  
  write.csv(otu5, paste("ancom_",var1,"_",file_path,"_",taxonomic_level,".csv",sep=""))
  
  
  ##plot them out
  # tiff(h=800*1.25, w=1600*1.25, paste("ancom_",var1,"_",file_path,"_",taxonomic_level,".tiff",sep=""), units = "px", pointsize = 12, compression = "none", res = 100)
  # do.call(grid.arrange, c(plist))
  # dev.off()
  
  jpeg(h=800*1.25, w=1600*1.25, paste("ancom_",var1,"_",file_path,"_",taxonomic_level,".jpg",sep=""), units = "px", quality = 0.9)
  do.call(grid.arrange, c(plist))
  dev.off()
}



make_ancom_plots_0.0.1_period_name <- function (file_path, taxonomic_level="X.OTU.ID", name_levels=taxonomic_level,var1, correction_level=3,newcol="newcol", the_id = "X.SampleID") {
  
  rm(otu_table, ref_names, otu2, map2, phyl2, phyl3, col1drop, col2drop, phyl4, ancom.OTU, detected_taxa, otu3, otu4, otu5, phyl5, plist, row_list)
  otu_table <- read.table(paste('core-metrics-results-',file_path,'/rarefied_table.txt',sep=""), comment.char="", header=T, sep="\t", fill=T, skip=1) %>% left_join(read.csv('taxonomy/taxonomy_forR.csv'), by=c("X.OTU.ID"="Feature.ID"))
  #otu_table$OTU <- as.character(otu_table$X.OTU.ID)
  ref_names <- subset(colnames(otu_table), (colnames(otu_table)%in%c("X.OTU.ID","kingdom","phylum","class","order","family","genus","species")==F))
  otu2 <- otu_table[,c("X.OTU.ID","kingdom","phylum","class","order","family","genus","species")]
  otu2$OTU <- as.character(otu2$X.OTU.ID)
  
  map2 <- read.table("metadata.tsv",comment.char = "", header=T) %>% dplyr::select_(.dots = list(the_id,var1))
  
  ## hardcoded
  phyl2 <- otu_table %>% group_by_(.dots=taxonomic_level) %>% summarize_at(ref_names,sum, na.rm=T)
  row_list <- as.character(unlist(phyl2[,paste(taxonomic_level)])) 
  rownames(phyl2) <- gsub(" ","", row_list)
  phyl3 <- phyl2 %>% dplyr::select_(.dots=paste0("-",taxonomic_level)) %>% t() %>% data.frame() 
  phyl3 <- mutate(phyl3, X.SampleID=rownames(phyl3))
  #	phyl3$sample2 <- gsub("\\.","", phyl3$X.SampleID)
  phyl3$sample2 <- gsub("X","", phyl3$X.SampleID)
  map2$sample2 <- gsub("_","",map2$X.SampleID)
  map2$sample2 <- gsub("-","",map2$sample2)
  
  ## OTU - by inocula
  col1drop <- paste0("-",the_id,".x"); col1drop
  col2drop <- paste0("-",the_id,".y"); col2drop
  
  phyl4 <- phyl3 %>% inner_join(map2, by=c("sample2")) %>% mutate(Group=as.factor(as.character(get(var1)))) %>% dplyr::select_(.dots = list(col1drop, col2drop, paste0("-",var1))) %>% dplyr::select(-sample2)
  
  ancom.OTU <- ANCOM( phyl4, sig = 0.05,multcorr = correction_level, repeated=FALSE )
  
  detected_taxa <- data.frame(OTU=unlist(sapply(X = ancom.OTU$detected, function(x) ifelse(substring(x,1,1)=="X",substring(x,2),x))))
  detected_taxa$OTU <- as.character(detected_taxa$OTU)
  print(detected_taxa$OTU)
  if(taxonomic_level!="OTU") {
    otu2[,paste(taxonomic_level)] = gsub(" ","",otu2[,paste(taxonomic_level)])
  }
  
  ## filter to OTU table by the data frame
  otu3 <- otu2 %>% mutate(newcol = get(name_levels[1]))#%>% dplyr::select_(.dots=c(taxonomic_level)) %>% inner_join(detected_taxa)
  
  ## melt OTU table for non-OTU stuff
  # melted_table <- data.frame(OTU_ID=character(),Count=integer(), Sample=character(), stringsAsFactors = F)
  # for (i in 2:(dim(otu_table)[2]-7-1)) {
  # 	rm(new_table)
  # 	new_table <- data.frame(OTU=as.character(otu_table$OTU),Count=as.integer(otu_table[,i]),Sample=as.character(colnames(otu_table)[i]), stringsAsFactors = F)
  # 	melted_table <- rbind(melted_table, new_table)
  # }
  # 
  # ## start merging
  # melted_table_no0 <- melted_table %>% filter(Count>-1)
  # mt0 <- melted_table_no0 %>% inner_join(otu2, by=c("OTU"="X.OTU.ID"))
  # mt0[1,]
  # mt0$sample2 <- gsub("\\.","", mt0$Sample)
  # mt0$sample2 <- gsub("_","", mt0$sample2)
  # map2$sample2 <- gsub("_","",map2$X.SampleID)
  # map2$sample2 <- gsub("-","",map2$sample2)
  # mt0[1,]
  # mt1 <- mt0 %>% inner_join(map2)
  # mt2 <- mt1 %>% arrange(phylum, class,order,family,genus,species,OTU)
  # mt2$sort_order <- c(1:dim(mt2)[1])
  # mt2$class <- paste(mt2$phylum, mt2$class,sep="_")
  # mt2$order <- paste(mt2$class, mt2$order,sep="_")
  # mt2$family <- paste(mt2$order, mt2$family,sep="_")
  # mt2$genus <- paste(mt2$family, mt2$genus,sep="_")
  # mt2$species <- paste(mt2$genus, mt2$species,sep="_")
  # mt2$OTU <- paste(mt2$species, mt2$OTU,sep="_")
  # 
  otu5 <- otu2 %>% filter(get(taxonomic_level)%in%detected_taxa$OTU) %>% distinct(get(taxonomic_level), .keep_all = T)
  
  extra_cols2 <- names(table(list(phyl4$Group)))
  extra_cols <- unname(c(extra_cols2,sapply(extra_cols2, function(x) paste0(x,"_sem"))))
  otu5[extra_cols] <- "NA"
  
  otu5
  if(length(name_levels)>1) {
    for(i in 2:length(name_levels)) {
      otu3 <- otu3 %>% mutate(newcol = paste(newcol, get(name_levels[i]), sep="_"))
    }
  }
  
  phyl3$sample2 <- gsub("\\.","", phyl3$X.SampleID)
  phyl3$sample2 <- gsub("_","", phyl3$sample2)
  
  otu3[,paste0(taxonomic_level)]
  otu4 <- otu3 %>% mutate(!!taxonomic_level:=gsub("\\[",".", get(taxonomic_level))) %>% mutate(!!taxonomic_level:=gsub("\\]",".", get(taxonomic_level))) %>% dplyr::select_(.dots=c(taxonomic_level, newcol)) %>% inner_join(detected_taxa, by = structure(names=taxonomic_level, "OTU")) %>% distinct_(.dots=taxonomic_level, .keep_all=T)
  otu3 <- otu4
  
  rm(i)
  ## run the filtered taxa to make plots
  for (i in 1:dim(otu3)[1]) {
    
    rm(phyl5, column_name)
    column_name <- otu3[,taxonomic_level][i]
    try(phyl5 <-phyl4 %>% mutate(rabun = get(column_name)/sum(otu_table[,2])) %>% dplyr::select(rabun, Group) %>% group_by(Group) %>% summarize(mean=mean(rabun), sem=sd(rabun)/sqrt(length(rabun))),T)
    if(exists("phyl5")==F) {
      try(phyl5 <- phyl4 %>% mutate(rabun = get(paste("X",column_name, sep=""))/sum(otu_table[,2])) %>% dplyr::select(rabun, Group) %>% group_by(Group) %>% summarize(mean=mean(rabun), sem=sd(rabun)/sqrt(length(rabun))))
    }
    
    rm(t)
    for (t in 1:length(extra_cols2)) {
      otu5[i,paste0(extra_cols2[t])] <- phyl5 %>% filter(Group==extra_cols2[t]) %>% dplyr::select(mean) %>% unlist()
      otu5[i,paste0(extra_cols2[t],"_sem")] <- phyl5 %>% filter(Group==extra_cols2[t]) %>% dplyr::select(sem) %>% unlist()
    }
    
    assign(paste("p",i,sep="_"),value = ggplot(phyl5, aes(x=Group, y=mean)) + 
             geom_bar(position=position_dodge(), stat="identity") +
             geom_errorbar(aes(ymin=mean+sem, ymax=mean-sem),
                           width=.2,                    # Width of the error bars
                           position=position_dodge(.9)) +
             coord_cartesian(ylim=c((min(phyl5$mean-phyl5$sem)*.9),max(phyl5$mean+phyl5$sem)*1.3)) + #scale_y_continuous(limits=c(.25,.4)) + 
             theme(axis.text=element_text(size=14), 
                   panel.background = element_blank(),
                   axis.line = element_line(), 
                   axis.ticks=element_line(), 
                   axis.title=element_text(size=16),
                   title=element_text(size=13)) +
             labs(y="relative abundance",x=var1,title=otu3$newcol[i])) 
  }
  
  ## make a list of the plots
  if(dim(detected_taxa)[1]>1) {
    plist <- list(p_1)
    for(q in 2:dim(detected_taxa)[1]) {
      plist[[q]] <- get(paste0("p_",q))
    }
  } else if (dim(detected_taxa)[1]>0) {
    plist <- list(p_1)
  }
  write.csv(otu5, paste("ancom_",var1,"_",file_path,"_",taxonomic_level,".csv",sep=""))
  
  
  ##plot them out
  # tiff(h=800*1.25, w=1600*1.25, paste("ancom_",var1,"_",file_path,"_",taxonomic_level,".tiff",sep=""), units = "px", pointsize = 12, compression = "none", res = 100)
  # do.call(grid.arrange, c(plist))
  # dev.off()
  
  jpeg(h=800*1.25, w=1600*1.25, paste("ancom_",var1,"_",file_path,"_",taxonomic_level,".jpg",sep=""), units = "px", quality = 0.9)
  do.call(grid.arrange, c(plist))
  dev.off()
}
