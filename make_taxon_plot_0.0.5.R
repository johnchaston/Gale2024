# 0.0.2	update Aug 2019 ot make compatible with varying file headers and names in the mapper file and taxonomy folders

prep_taxon_plot_core <- function(
  file_path, 
  map_path="", 
  plotvarnames = plotvarnames, 
  mapper_file = mapper_file, 
  taxonomic_level="ASV", 
  rpa_in_chart = 0.05, 
  read_depth = read_depth, 
  legend_position=legend_position,
  sort_x_axis = sort_x_axis,
  x_cluster = "X.SampleID"
  
  ) {
  
  varnames = c(plotvarnames, x_cluster, sort_x_axis)
  varnames <- varnames[!duplicated(varnames)]
  
  hierarchy <- list(kindgom = "kingdom",
                    phylum = c("kingdom","phylum"),
                    class = c("kingdom","phylum","class"),
                    order = c("kingdom","phylum","class","order"),
                    family = c("kingdom","phylum","class","order","family"),
                    genus = c("kingdom","phylum","class","order","family","genus"),
                    species = c("kingdom","phylum","class","order","family","genus","species"),
                    ASV = c("kingdom","phylum","class","order","family","genus","species","Feature.ID")
  )
  
  tax_level <- unlist(unname(hierarchy[taxonomic_level]))
  
  otu_table <- read.table(paste('core-metrics-results-',file_path,'/rarefied_table.txt',sep=""), comment.char="", header=T, sep="\t", fill=T, skip=1) %>% 
    left_join(read.csv(paste0('taxonomy',map_path,'/taxonomy_forR.csv')) %>%
                tidyr::unite(cluster,all_of(tax_level),sep = "_", remove = F) %>% 
                dplyr::select(Feature.ID, cluster), 
              by=c("X.OTU.ID"="Feature.ID")
              )
otu_table[1,]
otu_table %>% filter(X.OTU.ID %in% "d3ae5e3b180c63f5fca1a99523f2bf93")
  map2 <- read.table(mapper_file,comment.char = "", header=T, fill=T, sep="\t") %>%
    mutate(sample2 = gsub("_","", X.SampleID)) %>%
    mutate(sample2 = gsub("-","", sample2))
  melted_table[1,]
  
  melted_table %>% filter(OTU %in% "d3ae5e3b180c63f5fca1a99523f2bf93")
  ## melt
  melted_table <- otu_table %>%
    reshape2::melt(id.vars = c("X.OTU.ID","cluster")) %>% 
    filter(value>0) %>%
    dplyr::select(OTU = X.OTU.ID, cluster, Count = value, Sample = variable) %>%
    mutate(Sample = as.character(Sample), OTU = as.character(OTU)) %>%
    filter(Count>-1) %>%
    mutate(sample2 = gsub("\\.","", Sample)) %>%
    mutate(sample2 = gsub("_","", sample2)) %>%
    inner_join(map2, by = "sample2") %>%
    arrange(cluster) %>%
    dplyr::select(X.SampleID, OTU = 1, cluster = 2, Count = 3,Sample = 4,  all_of(varnames)) %>%
    tidyr::unite(col = twovar,  all_of(plotvarnames), sep = "_", remove = F) %>% 
    filter()
    
    ## find the rare taxa
  rare_taxa <- data.frame(table(melted_table$cluster)) %>%
    mutate(perc = Freq/sum(Freq)) %>%
    filter(perc < rpa_in_chart)  %>%
    dplyr::select(Var1) %>%
    unlist() %>% unname() %>% as.character()

  rare_taxa  
melted_table %>%  filter(X.SampleID %in% "d3ae5e3b180c63f5fca1a99523f2bf93" )
  
  # find the abundant taxa
  # abun_taxa <- melted_table %>% 
  #   group_by(across(all_of(varnames)),twovar,OTU,X.SampleID) %>% 
  #   filter(cluster%in%rare_taxa==F) %>% 
  #   mutate(cluster = as.factor(cluster)) %>%
  #   data.frame()
  
    ## specify the x-axis order by two values
  axis_order <- melted_table %>% 
    dplyr::select(all_of(x_cluster),all_of(sort_x_axis)) %>% 
    arrange(get(sort_x_axis)) %>% 
    distinct(get(x_cluster)) %>% 
    unlist() %>% unname() %>% droplevels()
  
  abun_taxa <- data.frame(table(melted_table$cluster)) %>%
    mutate(perc = Freq/sum(Freq)) %>%
    filter(perc >= rpa_in_chart) %>%
    mutate(Var1 = as.character(Var1),
           shortname = as.character(""))
  
  print(abun_taxa$Var1)

  for(i in 1:length(abun_taxa$Var1)) {
    try(space_split_vector <- strsplit(abun_taxa$Var1[i], split = " "))
    try(while (paste(tail(strsplit(tail(space_split_vector[[1]], 1),split = "")[[1]], 2), collapse="") == "__") {
      space_split_vector[[1]] <- head(space_split_vector[[1]], -1); space_split_vector[[1]]
    })
    try(abun_taxa$shortname[i] <- gsub(strsplit(tail(space_split_vector[[1]],1), "__")[[1]][2], pattern = "_",replacement = ""))
  }

    ## test if any shortnames occur more than onces
  countvar <- table(abun_taxa$shortname)[table(abun_taxa$shortname)>1]

    ## add a counter to any names that do occur more than once
  if(sum(countvar)>0) {
    for(i in names(countvar)) {
      abun_taxa <- rbind(abun_taxa %>%
        filter(shortname != i) %>% 
          droplevels(),
        abun_taxa %>%
          filter(shortname == i) %>%
          droplevels() %>%
          tibble::rowid_to_column("index") %>%
          mutate(shortname = paste0(shortname,index)) %>%
          dplyr::select(-index)
      )
    }
  }

  ptp <- list(melted_table, rare_taxa, abun_taxa, axis_order)
  print(abun_taxa$shortname)
  return(list(melted_table, rare_taxa, abun_taxa, axis_order))
}


prep_taxon_plot_numbersfirst <- function(
    file_path, 
    map_path="", 
    plotvarnames = plotvarnames, 
    mapper_file = mapper_file, 
    taxonomic_level="ASV", 
    rpa_in_chart = 0.05, 
    read_depth = read_depth, 
    legend_position=legend_position,
    sort_x_axis = sort_x_axis,
    x_cluster = "X.SampleID",
    predefined_taxa = "", 
    predefined_shortname = predefined_shortname
    
) {
  
  varnames = c(plotvarnames, x_cluster, sort_x_axis)
  varnames <- varnames[!duplicated(varnames)]
  
  hierarchy <- list(kindgom = "kingdom",
                    phylum = c("kingdom","phylum"),
                    class = c("kingdom","phylum","class"),
                    order = c("kingdom","phylum","class","order"),
                    family = c("kingdom","phylum","class","order","family"),
                    genus = c("kingdom","phylum","class","order","family","genus"),
                    species = c("kingdom","phylum","class","order","family","genus","species"),
                    ASV = c("kingdom","phylum","class","order","family","genus","species","Feature.ID")
  )
  
  tax_level <- unlist(unname(hierarchy[taxonomic_level]))
  
  otu_table <- read.table(paste('core-metrics-results-',file_path,'/rarefied_table.txt',sep=""), comment.char="", header=T, sep="\t", fill=T, skip=1) %>% 
    left_join(read.csv(paste0('taxonomy',map_path,'/taxonomy_forR.csv')) %>%
                tidyr::unite(cluster,all_of(tax_level),sep = "_", remove = F) %>% 
                dplyr::select(Feature.ID, cluster), 
              by=c("X.OTU.ID"="Feature.ID")
    )
  
  map2 <- read.table(mapper_file,comment.char = "", header=T, fill=T, sep="\t") %>%
    mutate(sample2 = gsub("_","", X.SampleID)) %>%
    mutate(sample2 = gsub("-","", sample2)) %>%
    mutate(sample2 = gsub("\\.","",sample2)) %>% 
    mutate(sample2 = paste0("X",sample2))
  
  
  
  ## melt
  melted_table <- otu_table %>%
    reshape2::melt(id.vars = c("X.OTU.ID","cluster")) %>%
    filter(value>0) %>%
    dplyr::select(OTU = X.OTU.ID, cluster, Count = value, Sample = variable) %>%
    mutate(Sample = as.character(Sample), OTU = as.character(OTU)) %>%
    filter(Count>-1) %>%
    mutate(sample2 = gsub("\\.","", Sample)) %>%
    mutate(sample2 = gsub("_","", sample2)) %>%
    inner_join(map2, by = "sample2") %>%
    arrange(cluster) %>%
    dplyr::select(X.SampleID, OTU = 1, cluster = 2, Count = 3,Sample = 4,  all_of(varnames)) %>%
    tidyr::unite(col = twovar,  all_of(plotvarnames), sep = "_", remove = F) %>%
    droplevels()
  
  rare_taxa <- melted_table %>%
    group_by(cluster) %>%
    dplyr::summarize(total = sum(Count)/length(table(melted_table$X.SampleID))/read_depth) %>%
    filter(total < rpa_in_chart) %>%
    filter(!cluster %in% predefined_taxa) %>%
    dplyr::select(cluster) %>%
    unlist() %>% unname() %>% as.character()
  
  ## specify the x-axis order by two values
  axis_order <- melted_table %>% 
    dplyr::select(all_of(x_cluster),all_of(sort_x_axis)) %>% 
    arrange(get(sort_x_axis)) %>% 
    distinct(get(x_cluster)) %>% 
    unlist() %>% unname() #%>% droplevels()
  
  abun_taxa <- melted_table %>%
    group_by(cluster) %>%
    dplyr::summarize(perc = sum(Count)/length(table(melted_table$X.SampleID))/read_depth) %>%
    #mutate(Freq = perc * length(table(melted_table[,x_cluster])) * read_depth) %>%
    mutate(Freq = perc * length(table(melted_table$X.SampleID)) * read_depth) %>%
    filter(perc >= rpa_in_chart) %>%
    dplyr::select(Var1 = cluster,Freq,perc) %>%
    mutate(shortname = as.character(""))
  
  for(i in 1:length(abun_taxa$Var1)) {
    try(space_split_vector <- strsplit(abun_taxa$Var1[i], split = " "))
    try(while (paste(tail(strsplit(tail(space_split_vector[[1]], 1),split = "")[[1]], 2), collapse="") == "__") {
      space_split_vector[[1]] <- head(space_split_vector[[1]], -1); space_split_vector[[1]]
    })
    try(abun_taxa$shortname[i] <- gsub(strsplit(tail(space_split_vector[[1]],1), "__")[[1]][2], pattern = "_",replacement = ""))
  }
  
  ## test if any shortnames occur more than onces
  countvar <- table(abun_taxa$shortname)[table(abun_taxa$shortname)>1]
  
  ## add a counter to any names that do occur more than once
  if(sum(countvar)>0) {
    for(i in names(countvar)) {
      abun_taxa <- rbind(abun_taxa %>%
                           filter(shortname != i) %>% 
                           droplevels(),
                         abun_taxa %>%
                           filter(shortname == i) %>%
                           droplevels() %>%
                           tibble::rowid_to_column("index") %>%
                           mutate(shortname = paste0(shortname,index)) %>%
                           dplyr::select(-index)
      )
    }
  }
  
  ptp <- list(melted_table, rare_taxa, abun_taxa, axis_order)
  print(abun_taxa$shortname)
  if(stringr::str_c(predefined_taxa, collapse = "") != "") {
    print(abun_taxa$shortname[!abun_taxa$shortname %in% predefined_shortname])
    
  }
  print(abun_taxa$Var1)
  return(list(melted_table, rare_taxa, abun_taxa, axis_order))
}

prep_taxon_plot <- function(
    file_path, 
    map_path="", 
    plotvarnames = plotvarnames, 
    mapper_file = mapper_file, 
    taxonomic_level="ASV", 
    rpa_in_chart = 0.05, 
    read_depth = read_depth, 
    legend_position=legend_position,
    sort_x_axis = sort_x_axis,
    x_cluster = "X.SampleID",
    predefined_taxa = "", 
    predefined_shortname = predefined_shortname
) {
  
  varnames = c(plotvarnames, x_cluster, sort_x_axis)
  varnames <- varnames[!duplicated(varnames)]

  hierarchy <- list(kindgom = "kingdom",
                    phylum = c("kingdom","phylum"),
                    class = c("kingdom","phylum","class"),
                    order = c("kingdom","phylum","class","order"),
                    family = c("kingdom","phylum","class","order","family"),
                    genus = c("kingdom","phylum","class","order","family","genus"),
                    species = c("kingdom","phylum","class","order","family","genus","species"),
                    ASV = c("kingdom","phylum","class","order","family","genus","species","Feature.ID")
  )

  tax_level <- unlist(unname(hierarchy[taxonomic_level]))

  otu_table <- read.table(paste('core-metrics-results-',file_path,'/rarefied_table.txt',sep=""), comment.char="", header=T, sep="\t", fill=T, skip=1) %>%
    left_join(read.csv(paste0('taxonomy',map_path,'/taxonomy_forR.csv')) %>%
                tidyr::unite(cluster,all_of(tax_level),sep = "_", remove = F) %>%
                dplyr::select(Feature.ID, cluster),
              by=c("X.OTU.ID"="Feature.ID")
    )

  map2 <- read.table(mapper_file,comment.char = "", header=T, fill=T, sep="\t") %>%
    mutate(sample2 = gsub("_","", X.SampleID)) %>%
    mutate(sample2 = gsub("-","", sample2)) %>%
    mutate(sample2 = gsub("\\.","",sample2))


  ## melt
  melted_table <- otu_table %>%
    reshape2::melt(id.vars = c("X.OTU.ID","cluster")) %>%
    filter(value>0) %>%
    dplyr::select(OTU = X.OTU.ID, cluster, Count = value, Sample = variable) %>%
    mutate(Sample = as.character(Sample), OTU = as.character(OTU)) %>%
    filter(Count>-1) %>%
    mutate(sample2 = gsub("\\.","", Sample)) %>%
    mutate(sample2 = gsub("_","", sample2)) %>%
    inner_join(map2, by = "sample2") %>%
    arrange(cluster) %>%
    dplyr::select(X.SampleID, OTU = 1, cluster = 2, Count = 3,Sample = 4,  all_of(varnames)) %>%
    tidyr::unite(col = twovar,  all_of(plotvarnames), sep = "_", remove = F) %>%
    droplevels()

  rare_taxa <- melted_table %>%
    group_by(cluster) %>%
    dplyr::summarize(total = sum(Count)/length(table(melted_table$X.SampleID))/read_depth) %>%
    filter(total < rpa_in_chart) %>%
    filter(!cluster %in% predefined_taxa) %>%
    dplyr::select(cluster) %>%
    unlist() %>% unname() %>% as.character()

   ## specify the x-axis order by two values
  axis_order <- melted_table %>%
    dplyr::select(all_of(x_cluster),all_of(sort_x_axis)) %>%
    arrange(get(sort_x_axis)) %>%
    distinct(get(x_cluster)) %>%
    unlist() %>% unname() #%>% droplevels()

  abun_taxa <- melted_table %>%
    group_by(cluster) %>%
    dplyr::summarize(perc = sum(Count)/length(table(melted_table$X.SampleID))/read_depth) %>%
    #mutate(Freq = perc * length(table(melted_table[,x_cluster])) * read_depth) %>%
    mutate(Freq = perc * length(table(melted_table$X.SampleID)) * read_depth) %>%
    filter(perc >= rpa_in_chart | cluster %in% predefined_taxa) %>%
    dplyr::select(Var1 = cluster,Freq,perc) %>%
    mutate(shortname = as.character(""))

  for(i in 1:length(abun_taxa$Var1)) {
    try(space_split_vector <- strsplit(abun_taxa$Var1[i], split = " "))
    try(while (paste(tail(strsplit(tail(space_split_vector[[1]], 1),split = "")[[1]], 2), collapse="") == "__") {
      space_split_vector[[1]] <- head(space_split_vector[[1]], -1); space_split_vector[[1]]
    })
    try(abun_taxa$shortname[i] <- gsub(strsplit(tail(space_split_vector[[1]],1), "__")[[1]][2], pattern = "_",replacement = ""))
  }

  ## test if any shortnames occur more than onces
  countvar <- table(abun_taxa$shortname)[table(abun_taxa$shortname)>1]

  ## add a counter to any names that do occur more than once
  if(sum(countvar)>0) {
    for(i in names(countvar)) {
      abun_taxa <- rbind(abun_taxa %>%
                           filter(shortname != i) %>%
                           droplevels(),
                         abun_taxa %>%
                           filter(shortname == i) %>%
                           droplevels() %>%
                           tibble::rowid_to_column("index") %>%
                           mutate(shortname = paste0(shortname,index)) %>%
                           dplyr::select(-index)
      )
    }
  }

  ptp <- list(melted_table, rare_taxa, abun_taxa, axis_order)
  print(abun_taxa$shortname)
  if(stringr::str_c(predefined_taxa, collapse = "") != "") {
    print(abun_taxa$shortname[!abun_taxa$shortname %in% predefined_shortname])

  }
  print(abun_taxa$Var1)
  return(list(melted_table, rare_taxa, abun_taxa, axis_order))
}



make_taxon_plot <- function (
  ptp = ptp,
  plot_color = plot_color,
  facet_formula = facet_formula,
  sort_x_axis = sort_x_axis, 
  x_cluster = x_cluster,
  plotvarnames = plotvarnames,
  legend_position = legend_position,
  read_depth = read_depth, 
  plot_order = plot_order) {
  
  varnames = c(plotvarnames, x_cluster, sort_x_axis)
  varnames <- varnames[!duplicated(varnames)]
  
  rta <- bind_rows(ptp[[1]] %>% 
    filter((cluster%in%ptp[[2]])==F) %>% 
    group_by(across(all_of(varnames)),twovar, cluster,X.SampleID) %>% 
    dplyr::summarise(phylum.abun = sum(Count)/read_depth)%>%
    inner_join(ptp[[3]], by = c("cluster"="Var1")) %>%
    dplyr::select(-Freq,-perc,-cluster),
  ptp[[1]] %>% 
    filter((cluster%in%ptp[[2]])==T) %>% 
    group_by(across(all_of(varnames)),twovar,X.SampleID) %>% 
    dplyr::summarise(phylum.abun = sum(Count)/read_depth) %>% 
    ungroup() %>%
    mutate(shortname = "other")
  ) %>% ungroup()

    ## set the order of the bars and colors
  if(length(plot_order) == 1) {
    plot_order <- ptp[[3]]$shortname %>% as.character()
  }
  rta$genus2 <- rta$shortname
  rta$genus2 <- factor(rta$genus2,levels=c("other",plot_order))
  
    ## order replicates within a cluster - variable x_cluster
  rta$xorder <- factor(unlist(unname(rta[,x_cluster])),levels=c(as.character(ptp[[4]])))
  
  if(length(plot_color) == 1) {
    plot_color <- c(rep("red", length(table(list(rta$genus2)))))
  }
  
  # p<- ggplot(rta, aes(x = xorder, y = rpa, fill = genus2)) + 
  #   facet_grid(as.formula(facet_formula) , drop = TRUE, space = "free", scales = "free") + 
p<- ggplot(rta, aes(x = xorder, y = phylum.abun, fill = genus2)) + 
    facet_grid(as.formula(facet_formula) , drop = TRUE, space = "free", scales = "free") + 
    geom_bar(stat = "identity", width = 1) +
    scale_fill_manual(values=plot_color) +
    theme_cowplot() + 
    theme(legend.position = legend_position,
          axis.text.x = element_blank()) +
    ylab("fractional abundance")

  
  # plot(p)
  # jpeg(h=800, w=1600, paste("taxon_plot_",file_path,"_",taxonomic_level,".jpg",sep=""), units = "px", quality = 0.9)
  # plot(p)
  # dev.off()
  
  return(p)
}


make_taxon_plot_condensed <- function (
    ptp = ptp,
    plot_color = plot_color,
    facet_formula = facet_formula,
    sort_x_axis = sort_x_axis, 
    x_cluster = x_cluster,
    plotvarnames = plotvarnames,
    legend_position = legend_position,
    repcol = "X.SampleID", plot_order = plot_order, read_depth=1, relevel_col = NULL, relevel_vec = NULL, relevel_row = NULL, relevelrow_vec = NULL, xorder = NULL, column_bar_width = 1, dropTheLevels = F, default_colors = F, predefined_taxa = predefined_taxa, predefined_colors = predefined_colors, predefined_shortname = predefined_shortname, default_order = T, order_numbers = F) {
  
  varnames = c(plotvarnames, x_cluster, sort_x_axis)
  varnames <- varnames[!duplicated(varnames)]
  
  var3 <- c("X.SampleID","sample_num")
  varnames2 <- varnames[!varnames%in%c(var3)]
  
  ptp[[1]]$groupid <- ptp[[1]][x_cluster]

  growing_df <- data.frame(body_site = character(), mucosa_or_lumen = character(), treatment = character(),phylum.abun = numeric(), shortname = character())

  for (i in names(table(ptp[[1]]$groupid))) {
  working_df <- ptp[[1]] %>% 
    filter(groupid == i) %>% droplevels()
  head(working_df)
  num_groups <- dim(working_df %>% group_by(get(repcol)) %>% dplyr::summarize(count = dplyr::n()))[1]
  abun_taxa_df <- working_df %>%
    filter((cluster%in%ptp[[3]]$Var1)) %>% 
    group_by(across(all_of(varnames2)),cluster) %>% 
    dplyr::summarise(phylum.abun = (sum(Count))/(read_depth)/num_groups) %>%
    ungroup() %>% 
    inner_join(ptp[[3]], by = c("cluster"="Var1")) %>%
    dplyr::select(-Freq,-perc,-cluster)
  rare_taxa_df <- working_df %>% 
    filter(!cluster%in%(ptp[[3]])$Var1) %>% 
    group_by(across(all_of(varnames2))) %>% 
    dplyr::summarise(phylum.abun = (sum(Count))/(read_depth)/num_groups) %>%
    ungroup() %>% 
    mutate(shortname = "other")    
  growing_df <- bind_rows(growing_df,abun_taxa_df,rare_taxa_df)  
}
growing_df  

  ## set the order of the bars and colors
if(default_colors == T) {
  remaining_shortname <- ptp[[3]]$shortname[!ptp[[3]]$shortname %in% predefined_shortname]
  if(default_order == T) {
    plot_order <- c(remaining_shortname, predefined_shortname)
    newplot_color <- c(plot_color, predefined_colors)
    plot_color <- newplot_color
  } else {
    plot_order <- ptp[[3]]$shortname[order_numbers] 
    plot_color <- c(plot_color)
  }
} 

rta <- growing_df
rta$genus2 <- factor(rta$shortname)

if(length(plot_order) == 1) {
  plot_order <- ptp[[3]]$shortname %>% as.character()
  plot_color <- c(rep("red", length(table(list(rta$genus2)))))
}

if(sum(is.na(plot_order))> 0) {
  plot_order[is.na(plot_order)] <- "other"
  rta$genus2[is.na(rta$genus2)] <- "other"
  rta$genus2 <- factor(rta$genus2,levels=c("other",plot_order[-which(plot_order == "other")]))
} else {
  rta$genus2 <- factor(rta$genus2,levels=c("other",plot_order))
}

if("" %in% plot_order) {
  plot_order[plot_order == ""] <- "other"
  rta$genus2[rta$genus2 == ""] <- "other"
#  levels(factor(rta$genus2))
  rta$genus2 <- factor(rta$genus2,levels=c("other",plot_order[-which(plot_order == "other")]))
} else {
  rta$genus2 <- factor(rta$genus2,levels=c("other",plot_order))
}

if (length(which(!plot_order %in% droplevels(rta$genus2))) > 0) {
  po2 <- plot_order[!(plot_order == "")]
  po3 <- po2[which(!po2 %in% droplevels(rta$genus2))]
  plot_color <- plot_color[-(which(plot_order %in% po3)+1)]
  plot_order <- plot_order[!plot_order %in% po3]
}

rta$genus2 <- droplevels(rta$genus2)

  ## order replicates within a cluster - variable x_cluster
 if(is.null(xorder)) {
   rta$xorder <- factor(unlist(unname(rta[,x_cluster])),levels=c(as.character(ptp[[4]])))
 } else {
   rta$xorder = rta[[xorder]]
 }

  ## reorder x-axis
  if (!is.null(relevel_col)) {
    rta[,relevel_col] <- factor(unlist(unname(rta[,relevel_col])), levels = relevel_vec)
    rta$xorder <- rta[,relevel_col]
  }

  ## reorder y-axis
  if (!is.null(relevel_row)) {
    rta[,relevel_row] <- factor(unlist(unname(rta[,relevel_row])), levels = relevelrow_vec)
  }

print(levels(rta$genus2))
print(plot_color)
p<- ggplot(rta %>% droplevels(), aes(x = xorder, y = phylum.abun, fill = genus2)) + 
    facet_grid(as.formula(facet_formula) , drop = T, space = "free", scales = "free") + 
    geom_bar(stat = "identity", width = column_bar_width) +
    scale_fill_manual(values=plot_color) +
    theme_cowplot() + 
    theme(legend.position = legend_position,
          axis.text.x = element_blank()) +
    ylab("fractional abundance") 
  

  p
  # plot(p)
  # jpeg(h=800, w=1600, paste("taxon_plot_",file_path,"_",taxonomic_level,".jpg",sep=""), units = "px", quality = 0.9)
  # plot(p)
  # dev.off()
  
  return(p)
}


make_taxon_plot_condensed_nofacet <- function (
    ptp = ptp,
    plot_color = plot_color,
#    facet_formula = facet_formula,
    sort_x_axis = sort_x_axis, 
    x_cluster = x_cluster,
    plotvarnames = plotvarnames,
    legend_position = legend_position,
    repcol = "X.SampleID", plot_order = plot_order, read_depth=1, relevel_col = NULL, relevel_vec = NULL, relevel_row = NULL, relevelrow_vec = NULL, xorder = NULL, column_bar_width = 1) {
  
  varnames = c(plotvarnames, x_cluster, sort_x_axis)
  varnames <- varnames[!duplicated(varnames)]
  
  var3 <- c("X.SampleID","sample_num")
  varnames2 <- varnames[!varnames%in%c(var3)]
  
  
  #   cde <- ptp[[1]] %>% group_by(groupid, sample_num) %>% summarize(count = dplyr::n()) %>% ungroup() %>% group_by(groupid) %>% summarize(count_samplenum = dplyr::n())
  # 
  # i <- names(table(ptp[[1]]$groupid)[2])
  # i
  # 
  # ptp[[1]]$groupid
  
  ptp[[1]]$groupid <- ptp[[1]][x_cluster]
  
  growing_df <- data.frame(body_site = character(), mucosa_or_lumen = character(), treatment = character(),phylum.abun = numeric(), shortname = character())
  for (i in names(table(ptp[[1]]$groupid))) {
    working_df <- ptp[[1]] %>% 
      filter(groupid == i) %>% droplevels()
    head(working_df)
    num_groups <- dim(working_df %>% group_by(get(repcol)) %>% dplyr::summarize(count = dplyr::n()))[1]
    abun_taxa_df <- working_df %>%
      filter((cluster%in%ptp[[3]]$Var1)) %>% 
      group_by(across(all_of(varnames2)),cluster) %>% 
      dplyr::summarise(phylum.abun = (sum(Count))/(read_depth)/num_groups) %>%
      ungroup() %>% 
      inner_join(ptp[[3]], by = c("cluster"="Var1")) %>%
      dplyr::select(-Freq,-perc,-cluster)
    rare_taxa_df <- working_df %>% 
      filter(!cluster%in%(ptp[[3]])$Var1) %>% 
      group_by(across(all_of(varnames2))) %>% 
      dplyr::summarise(phylum.abun = (sum(Count))/(read_depth)/num_groups) %>%
      ungroup() %>% 
      mutate(shortname = "other")    
    growing_df <- bind_rows(growing_df,abun_taxa_df,rare_taxa_df)  
  }
  growing_df  
  
  # rta <- bind_rows(ptp[[1]] %>% 
  #                    filter((cluster%in%ptp[[3]]$Var1)) %>% 
  #                    group_by(across(all_of(varnames2)),cluster, sample_num) %>% 
  #                    dplyr::summarise(mean1 = (sum(Count))/(read_depth)) %>%
  #                    ungroup() %>% 
  #                    group_by(cluster, across(all_of(varnames2))) %>%
  #                    summarize(phylum.abun = mean(mean1), num_count = dplyr::n()) %>%
  #                    ungroup() %>%
  #                    inner_join(ptp[[3]], by = c("cluster"="Var1")) %>%
  #                    dplyr::select(-Freq,-perc,-cluster),
  #                  ptp[[1]] %>% 
  #                    filter(!cluster%in%(ptp[[3]])$Var1) %>% 
  #                    group_by(across(all_of(varnames2)), sample_num) %>% 
  #                    dplyr::summarise(mean1 = (sum(Count))/(read_depth)) %>%
  #                    ungroup() %>% 
  #                    group_by(across(all_of(varnames2))) %>%
  #                    summarize(phylum.abun = mean(mean1)) %>%
  #                    ungroup() %>%
  #                    mutate(shortname = "other", num_count = 30)
  #                  # ptp[[1]] %>% 
  #                  #   filter((cluster%in%ptp[[2]])==T) %>% 
  #                  #   group_by(across(all_of(varnames2))) %>% 
  #                  #   dplyr::summarise(phylum.abun = (sum(Count))/(read_depth)) %>% 
  #                  #   ungroup() %>%
  #                  #   mutate(shortname = "other")
  # ) %>% ungroup()
  
  rta <- growing_df# %>%
  #  mutate(x_order = factor(across(all_of(sort_x_axis)))
  # abc <-ptp[[1]] %>% filter(groupid == "c l 125ppm")
  # abc %>% filter(cluster %in% as.character(ptp[[3]]$Var1)) %>% droplevels() %>% group_by(cluster,groupid,sample_num) %>% summarize(meancount = sum(Count)/read_depth) %>% ungroup() %>% group_by(cluster, groupid) %>% summarize(mean2 = mean(meancount), numcount = dplyr::n())
  # abc %>% filter(cluster %in% as.character(ptp[[3]]$Var1)) %>% droplevels() %>% group_by(sample_num) %>% summarize(num = dplyr::n())
  
  
  
  # abc <-ptp[[1]] %>% filter(groupid == "c l 125ppm")
  # def <- ptp[[1]] %>% 
  #   filter((cluster%in%as.character(ptp[[3]]$Var1))==T)
  # 
  # 
  # def
  # ptp[[2]]
  # 
  # table(abc$sample_num)
  #   rta <- bind_rows(ptp[[1]] %>% 
  #                      filter((cluster%in%as.character(ptp[[3]]$Var1))==T) %>% 
  #                      group_by(groupid,cluster) %>% 
  #                      dplyr::summarise(phylum.abun = mean(Count)) %>%
  #                      inner_join(ptp[[3]], by = c("cluster"="Var1")) %>%
  #                      dplyr::select(-Freq,-perc,-cluster),
  #                    ptp[[1]] %>% 
  #                      filter((cluster%in%ptp[[2]])==T) %>% 
  #                      group_by(groupid) %>% 
  #                      dplyr::summarise(phylum.abun = mean(Count)) %>% 
  #                      ungroup() %>%
  #                      mutate(shortname = "other")
  #   ) %>% ungroup()
  #   
  #   rta  
  ## set the order of the bars and colors
  if(length(plot_order) == 1) {
    plot_order <- ptp[[3]]$shortname %>% as.character()
  }
  rta$genus2 <- rta$shortname
  rta$genus2 <- factor(rta$genus2,levels=c("other",plot_order))
  
  
  ## order replicates within a cluster - variable x_cluster
  if(is.null(xorder)) {
    rta$xorder <- factor(unlist(unname(rta[,x_cluster])),levels=c(as.character(ptp[[4]])))
  } else {
    rta$xorder = rta[[xorder]]
  }
  
  ## reorder x-axis
  #  relevel_col = "site"
  #  relevel_vec = c("Charlottesville, VA", "Churchville, MD", "Media, PA", "Middlefield, CT", "Harvard, MA", "Durham, NH", "Bowdoin, ME", "Etna, ME")
  if (!is.null(relevel_col)) {
    rta[,relevel_col] <- factor(unlist(unname(rta[,relevel_col])), levels = relevel_vec)
  }
  
  ## reorder y-axis
  if (!is.null(relevel_row)) {
    rta[,relevel_row] <- factor(unlist(unname(rta[,relevel_row])), levels = relevelrow_vec)
  }
  
  if(length(plot_color) == 1) {
    plot_color <- c(rep("red", length(table(list(rta$genus2)))))
  }
  
  # p<- ggplot(rta, aes(x = xorder, y = rpa, fill = genus2)) + 
  #   facet_grid(as.formula(facet_formula) , drop = TRUE, space = "free", scales = "free") + 
  p<- ggplot(rta, aes(x = xorder, y = phylum.abun, fill = genus2)) + 
 #   facet_grid(as.formula(facet_formula) , drop = TRUE, space = "free", scales = "free") + 
    geom_bar(stat = "identity", width = column_bar_width) +
    scale_fill_manual(values=plot_color) +
    theme_cowplot() + 
    theme(legend.position = legend_position,
          axis.text.x = element_blank()) +
    ylab("fractional abundance") 
  
  p
  # plot(p)
  # jpeg(h=800, w=1600, paste("taxon_plot_",file_path,"_",taxonomic_level,".jpg",sep=""), units = "px", quality = 0.9)
  # plot(p)
  # dev.off()
  
  return(p)
}


make_taxon_plot_condensed_numericX <- function (
    ptp = ptp,
    plot_color = plot_color,
    facet_formula = facet_formula,
    sort_x_axis = sort_x_axis, 
    x_cluster = x_cluster,
    plotvarnames = plotvarnames,
    legend_position = legend_position,
    repcol = "X.SampleID", plot_order = plot_order, read_depth=1, relevel_col = NULL, relevel_vec = NULL, relevel_row = NULL, relevelrow_vec = NULL) {
  
  varnames = c(plotvarnames, x_cluster, sort_x_axis)
  varnames <- varnames[!duplicated(varnames)]
  
  var3 <- c("X.SampleID","sample_num")
  varnames2 <- varnames[!varnames%in%c(var3)]
  
  ptp[[1]]$groupid <- ptp[[1]][x_cluster]
  
  growing_df <- data.frame(body_site = character(), mucosa_or_lumen = character(), treatment = character(),phylum.abun = numeric(), shortname = character())
  for (i in names(table(ptp[[1]]$groupid))) {
    working_df <- ptp[[1]] %>% 
      filter(groupid == i) %>% droplevels()
    head(working_df)
    num_groups <- dim(working_df %>% group_by(get(repcol)) %>% dplyr::summarize(count = dplyr::n()))[1]
    abun_taxa_df <- working_df %>%
      filter((cluster%in%ptp[[3]]$Var1)) %>% 
      group_by(across(all_of(varnames2)),cluster) %>% 
      dplyr::summarise(phylum.abun = (sum(Count))/(read_depth)/num_groups) %>%
      ungroup() %>% 
      inner_join(ptp[[3]], by = c("cluster"="Var1")) %>%
      dplyr::select(-Freq,-perc,-cluster)
    rare_taxa_df <- working_df %>% 
      filter(!cluster%in%(ptp[[3]])$Var1) %>% 
      group_by(across(all_of(varnames2))) %>% 
      dplyr::summarise(phylum.abun = (sum(Count))/(read_depth)/num_groups) %>%
      ungroup() %>% 
      mutate(shortname = "other")    
    growing_df <- bind_rows(growing_df,abun_taxa_df,rare_taxa_df)  
  }

  rta <- growing_df
  
  if(length(plot_order) == 1) {
    plot_order <- ptp[[3]]$shortname %>% as.character()
  }
  rta$genus2 <- rta$shortname
  rta$genus2 <- factor(rta$genus2,levels=c("other",plot_order))

 if(length(plot_color) == 1) {
    plot_color <- c(rep("red", length(table(list(rta$genus2)))))
  }
 
  p <- ggplot(rta, aes(x = as.numeric(as.character(get(x_cluster))), y = phylum.abun, fill = genus2)) + 
  #  facet_grid(as.formula(facet_formula) , drop = TRUE, space = "free", scales = "free") + 
    geom_bar(stat = "identity", width = 1) +
    scale_fill_manual(values=plot_color) +
    theme_cowplot() + 
    theme(legend.position = legend_position,
          axis.text.x = element_blank()) +
    ylab("fractional abundance")
  return(p)
}


prep_taxon_plot_byOrder <- function(
    file_path, 
    map_path = "",
    plotvarnames = plotvarnames, 
    mapper_file = mapper_file, 
    taxonomic_level="ASV", 
    rpa_in_chart = 0.05, 
    read_depth = read_depth, 
    legend_position=legend_position,
    sort_x_axis = sort_x_axis,
    x_cluster = "X.SampleID",
    taxa_filter = NULL
    
) {
  
  varnames = c(plotvarnames, x_cluster, sort_x_axis)
  varnames <- varnames[!duplicated(varnames)]
  
  hierarchy <- list(kindgom = "kingdom",
                    phylum = c("kingdom","phylum"),
                    class = c("kingdom","phylum","class"),
                    order = c("kingdom","phylum","class","order"),
                    family = c("kingdom","phylum","class","order","family"),
                    genus = c("kingdom","phylum","class","order","family","genus"),
                    species = c("kingdom","phylum","class","order","family","genus","species"),
                    ASV = c("kingdom","phylum","class","order","family","genus","species","Feature.ID")
  )
  
  tax_level <- unlist(unname(hierarchy[taxonomic_level]))

  
  abc <- read.table(paste('core-metrics-results-',file_path,'/rarefied_table.txt',sep=""), comment.char="", header=T, sep="\t", fill=T, skip=1) %>%
    left_join(read.csv(paste0('taxonomy',map_path,'/taxonomy_forR.csv')) %>%
                tidyr::unite(cluster,all_of(tax_level),sep = "_", remove = F) %>% 
                dplyr::select(Feature.ID, cluster), 
              by=c("X.OTU.ID"="Feature.ID")
    )
  
  
    
  otu_table <- read.table(paste('core-metrics-results-',file_path,'/rarefied_table.txt',sep=""), comment.char="", header=T, sep="\t", fill=T, skip=1) %>% 
    full_join(read.csv(paste0('taxonomy',map_path,'/taxonomy_forR.csv')) %>%
                tidyr::unite(cluster,all_of(tax_level),sep = "_", remove = F) %>% 
                dplyr::select(Feature.ID, cluster), 
              by=c("X.OTU.ID"="Feature.ID")
    )
  
  map2 <- read.table(mapper_file,comment.char = "", header=T, fill=T, sep="\t") %>%
    mutate(sample2 = gsub("_","", X.SampleID)) %>%
    mutate(sample2 = gsub("-","", sample2)) %>%
    mutate(sample2 = gsub("\\.","",sample2))
  
  
  ## melt
  melted_table <- otu_table %>%
    reshape2::melt(id.vars = c("X.OTU.ID","cluster")) %>% 
    filter(value>0) %>%
    dplyr::select(OTU = X.OTU.ID, cluster, Count = value, Sample = variable) %>%
    mutate(Sample = as.character(Sample), OTU = as.character(OTU)) %>%
    filter(Count>-1) %>%
    mutate(sample2 = gsub("\\.","", Sample)) %>%
    mutate(sample2 = gsub("_","", sample2)) %>%
    inner_join(map2, by = "sample2") %>%
    arrange(cluster) %>%
    dplyr::select(X.SampleID, OTU = 1, cluster = 2, Count = 3,Sample = 4,  all_of(varnames)) %>%
    tidyr::unite(col = twovar,  all_of(plotvarnames), sep = "_", remove = F) %>% 
    droplevels()
  
  if(is.null(taxa_filter)) {
    rare_taxa <- melted_table %>%
      group_by(cluster) %>%
      dplyr::summarize(total = sum(Count)/length(table(melted_table$X.SampleID))/read_depth) %>%
      filter(total < rpa_in_chart) %>%
      dplyr::select(cluster) %>%
      unlist() %>% unname() %>% as.character()
  } else {
    rare_taxa <- melted_table %>%
      group_by(cluster) %>%
      dplyr::summarize(total = sum(Count)/length(table(melted_table$X.SampleID))/read_depth) %>%
      filter(total < rpa_in_chart) %>%
      filter(grepl(x = cluster, taxa_filter)) %>%
      dplyr::select(cluster) %>%
      unlist() %>% unname() %>% as.character()
  }  
  
  ## specify the x-axis order by two values
  axis_order <- melted_table %>% 
    dplyr::select(all_of(x_cluster),all_of(sort_x_axis)) %>% 
    arrange(get(sort_x_axis)) %>% 
    distinct(get(x_cluster)) %>% 
    unlist() %>% unname() #%>% droplevels()
  
  
  if(is.null(taxa_filter)) {
    abun_taxa <- melted_table %>%
      group_by(cluster) %>%
      dplyr::summarize(perc = sum(Count)/length(table(melted_table$X.SampleID))/read_depth) %>%
      #mutate(Freq = perc * length(table(melted_table[,x_cluster])) * read_depth) %>%
      mutate(Freq = perc * length(table(melted_table$X.SampleID)) * read_depth) %>%
      filter(perc >= rpa_in_chart) %>%
      dplyr::select(Var1 = cluster,Freq,perc) %>%
      mutate(shortname = as.character(""))
  } else {
    abun_taxa <- melted_table %>%
      group_by(cluster) %>%
      dplyr::summarize(perc = sum(Count)/length(table(melted_table$X.SampleID))/read_depth) %>%
      #mutate(Freq = perc * length(table(melted_table[,x_cluster])) * read_depth) %>%
      mutate(Freq = perc * length(table(melted_table$X.SampleID)) * read_depth) %>%
      filter(grepl(x = cluster, taxa_filter)) %>%
      filter(perc >= rpa_in_chart) %>%
      dplyr::select(Var1 = cluster,Freq,perc) %>%
      mutate(shortname = as.character(""))
  } 
  
  for(i in 1:length(abun_taxa$Var1)) {
    try(space_split_vector <- strsplit(abun_taxa$Var1[i], split = " "))
    try(while (paste(tail(strsplit(tail(space_split_vector[[1]], 1),split = "")[[1]], 2), collapse="") == "__") {
      space_split_vector[[1]] <- head(space_split_vector[[1]], -1); space_split_vector[[1]]
    })
    try(abun_taxa$shortname[i] <- gsub(strsplit(tail(space_split_vector[[1]],1), "__")[[1]][2], pattern = "_",replacement = ""))
  }
  
  ## test if any shortnames occur more than onces
  countvar <- table(abun_taxa$shortname)[table(abun_taxa$shortname)>1]
  
  ## add a counter to any names that do occur more than once
  if(sum(countvar)>0) {
    for(i in names(countvar)) {
      abun_taxa <- rbind(abun_taxa %>%
                           filter(shortname != i) %>% 
                           droplevels(),
                         abun_taxa %>%
                           filter(shortname == i) %>%
                           droplevels() %>%
                           tibble::rowid_to_column("index") %>%
                           mutate(shortname = paste0(shortname,index)) %>%
                           dplyr::select(-index)
      )
    }
  }
  
  ptp <- list(melted_table, rare_taxa, abun_taxa, axis_order)
  print(abun_taxa$shortname)
  print(abun_taxa$Var1)
  return(list(melted_table, rare_taxa, abun_taxa, axis_order))
}



make_taxon_plot_condensed_byOrder <- function (
    ptp = ptp,
    plot_color = plot_color,
    facet_formula = facet_formula,
    sort_x_axis = sort_x_axis, 
    x_cluster = x_cluster,
    plotvarnames = plotvarnames,
    legend_position = legend_position,
    repcol = "X.SampleID", plot_order = plot_order, read_depth=1, relevel_col = NULL, relevel_vec = NULL, relevel_row = NULL, relevelrow_vec = NULL, xorder = NULL, column_bar_width = 1, taxa_filter = NULL) {
  
  varnames = c(plotvarnames, x_cluster, sort_x_axis)
  varnames <- varnames[!duplicated(varnames)]
  
  var3 <- c("X.SampleID","sample_num")
  varnames2 <- varnames[!varnames%in%c(var3)]
  
  ptp[[1]]$groupid <- ptp[[1]][x_cluster]
  
  growing_df <- data.frame(body_site = character(), mucosa_or_lumen = character(), treatment = character(),phylum.abun = numeric(), shortname = character())
  for (i in names(table(ptp[[1]]$groupid))) {
    working_df <- ptp[[1]] %>% 
      filter(groupid == i) %>%
      #    filter(grepl(x = cluster, taxa_filter)) %>% 
      droplevels()
    num_groups <- dim(working_df %>% group_by(get(repcol)) %>% dplyr::summarize(count = dplyr::n()))[1]
    abun_taxa_df <- working_df %>%
      filter((cluster%in%ptp[[3]]$Var1)) %>% 
      group_by(across(all_of(varnames2)),cluster) %>% 
      dplyr::summarise(phylum.abun = (sum(Count))/(read_depth)/num_groups) %>%
      ungroup() %>% 
      inner_join(ptp[[3]], by = c("cluster"="Var1")) %>%
      dplyr::select(-Freq,-perc,-cluster)
    rare_taxa_df <- working_df %>% 
      filter(cluster%in%(ptp[[2]])) %>% 
      group_by(across(all_of(varnames2))) %>% 
      dplyr::summarise(phylum.abun = (sum(Count))/(read_depth)/num_groups) %>%
      ungroup() %>% 
      mutate(shortname = "other")    
    
    to_add_df <- bind_rows(abun_taxa_df,rare_taxa_df) 
    
    to_add_df2 <- to_add_df %>%
      mutate(phylum.abun = phylum.abun/sum(phylum.abun))
    
    growing_df <- bind_rows(growing_df,to_add_df2)  
  }
  growing_df  
  
  
  
  
  rta <- growing_df# %>%
  
  ## set the order of the bars and colors
  if(length(plot_order) == 1) {
    plot_order <- ptp[[3]]$shortname %>% as.character()
  }
  rta$genus2 <- rta$shortname
  rta$genus2 <- factor(rta$genus2,levels=c("other",plot_order))
  
  
  ## order replicates within a cluster - variable x_cluster
  if(is.null(xorder)) {
    rta$xorder <- factor(unlist(unname(rta[,x_cluster])),levels=c(as.character(ptp[[4]])))
  } else {
    rta$xorder = rta[[xorder]]
  }
  
  ## reorder x-axis
  #  relevel_col = "site"
  #  relevel_vec = c("Charlottesville, VA", "Churchville, MD", "Media, PA", "Middlefield, CT", "Harvard, MA", "Durham, NH", "Bowdoin, ME", "Etna, ME")
  if (!is.null(relevel_col)) {
    rta[,relevel_col] <- factor(unlist(unname(rta[,relevel_col])), levels = relevel_vec)
  }
  
  ## reorder y-axis
  if (!is.null(relevel_row)) {
    rta[,relevel_row] <- factor(unlist(unname(rta[,relevel_row])), levels = relevelrow_vec)
  }
  
  if(length(plot_color) == 1) {
    plot_color <- c(rep("red", length(table(list(rta$genus2)))))
  }
  
  # p<- ggplot(rta, aes(x = xorder, y = rpa, fill = genus2)) + 
  #   facet_grid(as.formula(facet_formula) , drop = TRUE, space = "free", scales = "free") + 
  p<- ggplot(rta, aes(x = xorder, y = phylum.abun, fill = genus2)) + 
    facet_grid(as.formula(facet_formula) , drop = TRUE, space = "free", scales = "free") + 
    geom_bar(stat = "identity", width = column_bar_width) +
    scale_fill_manual(values=plot_color) +
    theme_cowplot() + 
    theme(legend.position = legend_position,
          axis.text.x = element_blank()) +
    ylab("fractional abundance") 
  
  p
  # plot(p)
  # jpeg(h=800, w=1600, paste("taxon_plot_",file_path,"_",taxonomic_level,".jpg",sep=""), units = "px", quality = 0.9)
  # plot(p)
  # dev.off()
  
  return(p)
}

