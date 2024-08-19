make_plot <- function(file_path, taxonomy_extension, mapper_file, sig_level, make_plots = F) {

  hierarchy <- list(phylum = c("kingdom","phylum"),
                    class = c("kingdom","phylum","class"),
                    order = c("kingdom","phylum","class","order"),
                    family = c("kingdom","phylum","class","order","family"),
                    genus = c("kingdom","phylum","class","order","family","genus"),
                    species = c("kingdom","phylum","class","order","family","genus","species"),
                    ASV = c("kingdom","phylum","class","order","family","genus","species","Feature.ID")
  )

  sig_object <- read.csv(paste0('core-metrics-results-',file_path,"/all_sig_",sig_level,".csv"))
  if (sig_object$taxonomy[1] == "Feature.ID") {
    sig_object$taxonomy <- "ASV"
  }

  j <- hierarchy[[as.character(sig_object$taxonomy[1])]]

  taxonomy <- read.table(paste0('taxonomy-',taxonomy_extension,'/taxonomy_forR.csv'), fill = T, header = T, sep = ",") %>%
    unite(cluster,all_of(j),sep = "_", remove = F) %>%
    dplyr::select(X.OTU.ID = Feature.ID, cluster)

  ## read in OTU table, prep columns and row names
  otutable <- read.table(paste0('core-metrics-results-',file_path,"/rarefied_table.txt"), comment.char = "",skip = 1, fill = T, header = T, sep = "\t") %>%
    inner_join(taxonomy) %>%
    dplyr::select(-X.OTU.ID) %>%
    group_by(cluster) %>%
    summarize_all(sum)

  ## prep row names
  clusterIDs <- otutable$cluster
  ## prep column names
  sampleID <- otutable %>% dplyr::select(-cluster) %>% colnames()

  ## transpose table, including assigning sample names
  otutable2 <- otutable %>%
    dplyr::select(-cluster) %>%
    t() %>%
    data.frame() %>%
    mutate(X.SampleID = sampleID)

  ## assign column names
  colnames(otutable2)[1:(dim(otutable2)[2]-1)] <- clusterIDs

  ## merge with metadata
  cortable <- read.table(mapper_file, comment.char = "", header = T, fill = T, sep = "\t") %>%
    inner_join(otutable2)

  sig <- data.frame(cluster = character(), pval = numeric(), rho = numeric(), sum = numeric(), fdr = numeric(), test = character(), taxonomy = character())

  m = 1
  #cortable2 <- cortable
  #cortable <- cortable2 %>% filter(Season == "fall")
  m=5
  for(m in 1:dim(sig_object)[1]) {

    taxon <- "k__Bacteria_ p__Firmicutes_ c__Bacilli_ o__Lactobacillales"
    taxon <- "k__Bacteria_ p__Proteobacteria_ c__Alphaproteobacteria_ o__Rhodospirillales"
    taxon <- as.character(sig_object[m,"cluster"])
    test <- as.character(sig_object[m,"test"])

    plot_df <- cortable[,c(taxon, test)]
    colnames(plot_df) <- c("taxon","test")
    plot_df$taxon <- (plot_df$taxon/4975)*100

    ctest <- cor.test(plot_df$taxon, plot_df$test, method = "spearman", exact = F)

    ggplot(plot_df, aes(x = test, y = taxon)) +
      geom_jitter(width = 1)+
      theme_cowplot() +
      stat_smooth(method = "lm") +
      geom_label(aes(x = max(test)*.7, y = max(taxon)*.9), label = paste0("R2 = ",round(ctest$estimate^2,2),"\np = ",round(ctest$p.value,4)), label.size = NA, size = 6) +
      labs(y = taxon, x = test)

    ggsave(paste0('core-metrics-results-',file_path,"/file_plot",test,"_",m,"_",tail(j,1),".jpg"), bg="white")
  }
}

# sig_level = "species"
# taxonomy_extension = "fv2"

make_plot2 <- function(file_path, taxonomy_extension, mapper_file, sig_level, make_plots = F, jitter_width = 1, default_stats_position = T, xposition, yposition, hjustposition, vjustposition) {
  
  #cluster = sig_level
  
  hierarchy <- list(phylum = c("kingdom","phylum"),
                    class = c("kingdom","phylum","class"),
                    order = c("kingdom","phylum","class","order"),
                    family = c("kingdom","phylum","class","order","family"),
                    genus = c("kingdom","phylum","class","order","family","genus"),
                    species = c("kingdom","phylum","class","order","family","genus","species"),
                    ASV = c("kingdom","phylum","class","order","family","genus","species","Feature.ID")
  )
  
  sig_object <- read.csv(paste0('core-metrics-results-',file_path,"/all_sig_",sig_level,".csv"))
  if (sig_object$taxonomy[1] == "Feature.ID") {
    sig_object$taxonomy <- "ASV"
  }
  
  j <- hierarchy[[as.character(sig_object$taxonomy[1])]]
  
  taxonomy <- read.table(paste0('taxonomy-',taxonomy_extension,'/taxonomy_forR.csv'), fill = T, header = T, sep = ",") %>%
    unite(cluster,all_of(j),sep = "_", remove = F) %>%
    dplyr::select(X.OTU.ID = Feature.ID, cluster)
  
  ## read in OTU table, prep columns and row names
  otutable <- read.table(paste0('core-metrics-results-',file_path,"/rarefied_table.txt"), comment.char = "",skip = 1, fill = T, header = T, sep = "\t") %>%
    inner_join(taxonomy) %>%
    dplyr::select(-X.OTU.ID) %>%
    group_by(cluster) %>%
    summarize_all(sum)
  
  ## prep row names
  clusterIDs <- otutable$cluster
  ## prep column names
  sampleID <- otutable %>% dplyr::select(-cluster) %>% colnames()
  
  ## transpose table, including assigning sample names
  otutable2 <- otutable %>%
    dplyr::select(-cluster) %>%
    t() %>%
    data.frame() %>%
    mutate(X.SampleID = sampleID)
  
  depth = colSums(otutable[,2:3])
  
  ## assign column names
  colnames(otutable2)[1:(dim(otutable2)[2]-1)] <- clusterIDs
  
  ## merge with metadata
  cortable <- read.table(mapper_file, comment.char = "", header = T, fill = T, sep = "\t") %>%
    inner_join(otutable2)
  
  sig <- data.frame(cluster = character(), pval = numeric(), rho = numeric(), sum = numeric(), fdr = numeric(), test = character(), taxonomy = character())
  
  m = 1
  #cortable2 <- cortable
  #cortable <- cortable2 %>% filter(Season == "fall")
  m=5
  m=5
  for(m in 1:dim(sig_object)[1]) {
    
    taxon <- "k__Bacteria_ p__Firmicutes_ c__Bacilli_ o__Lactobacillales"
    taxon <- "k__Bacteria_ p__Proteobacteria_ c__Alphaproteobacteria_ o__Rhodospirillales"
    taxon <- as.character(sig_object[m,"cluster"])
    test <- as.character(sig_object[m,"test"])
    
    plot_df <- cortable[,c(taxon, test)]
    colnames(plot_df) <- c("taxon","test")
    plot_df$taxon <- (plot_df$taxon/depth)*100
    
    ctest <- cor.test(plot_df$taxon, plot_df$test, method = "spearman", exact = F)
  
    if(default_stats_position == T) {
      ggplot(plot_df, aes(x = test, y = taxon)) +
        geom_jitter(width = jitter_width)+
        theme_cowplot() +
        geom_label(aes(x = max(test)*.7, y = max(taxon)*.9), label = paste0("R2 = ",round(ctest$estimate^2,2),"\np = ",round(ctest$p.value,4)), label.size = NA, size = 6) +
        labs(y = taxon, x = test) + 
        stat_smooth(method = "lm")
    } else {
      ggplot(plot_df, aes(x = test, y = taxon)) +
        geom_jitter(width = jitter_width)+
        theme_cowplot() +
        geom_label(aes(x = xposition, y = yposition, hjust = hjustposition, vjust = vjustposition), label = paste0("R2 = ",round(ctest$estimate^2,2),"\np = ",round(ctest$p.value,4)), label.size = NA, size = 6) +
        labs(y = taxon, x = test) + 
        stat_smooth(method = "lm")
    }  
    ggsave(paste0('core-metrics-results-',file_path,"/file_plot",test,"_",m,"_",tail(j,1),".jpg"), bg="white")
  }
}


# file_path = "fruit"
#taxonomy_extension <- "DEST"
# frac_samples <- .25
#averaged = F
#
# make_plot("fruit","eastcoast",mapper_file,"genus")
# make_plots = T
# calculate_correlations("fruit","eastcoast",mapper_file, .5, averaged = c("orchard"))
taxonomy_extension = "fv2"
calculate_correlations <- function(file_path, taxonomy_extension, mapper_file, frac_samples = 0, averaged = F, use_fdr = T, specific_groups = NULL, defined_hierarchy = NULL, return_table = F) {

  ## build a list to test each hierarchical level
hierarchy <- list(phylum = c("kingdom","phylum"),
                  class = c("kingdom","phylum","class"),
                  order = c("kingdom","phylum","class","order"),
                  family = c("kingdom","phylum","class","order","family"),
                  genus = c("kingdom","phylum","class","order","family","genus"),
                  species = c("kingdom","phylum","class","order","family","genus","species"),
                  ASV = c("kingdom","phylum","class","order","family","genus","species","Feature.ID")
)

## build the value to cluster by

if(!is.null(defined_hierarchy)) {
  hierarchy = hierarchy[defined_hierarchy]
}

j=unlist(unname(hierarchy))
for(j in hierarchy) {
  
  taxonomy <- read.table(paste0('taxonomy-',taxonomy_extension,'/taxonomy_forR.csv'), fill = T, header = T, sep = ",") %>%
    unite(cluster,all_of(j),sep = "_", remove = F) %>% 
    dplyr::select(X.OTU.ID = Feature.ID, cluster)
  
  if(!is.null(specific_groups)) {
    taxonomy <- taxonomy %>% filter(cluster %in% specific_groups)
  }
  
  print(dim(taxonomy))
  ## read in OTU table, prep columns and row names
  otutable <- read.table(paste0('core-metrics-results-',file_path,"/rarefied_table.txt"), comment.char = "",skip = 1, fill = T, header = T, sep = "\t") %>%
    inner_join(taxonomy) %>%
    dplyr::select(-X.OTU.ID) %>%
    group_by(cluster) %>%
    summarize_all(sum)
  
  ## prep row names
  clusterIDs <- otutable$cluster
  
  ## prep column names
  sampleID <- otutable %>% dplyr::select(-cluster) %>% colnames()
  
  ## transpose table, including assigning sample names
  otutable2 <- otutable %>%
    dplyr::select(-cluster) %>%
    t() %>%
    data.frame() %>%
    mutate(X.SampleID = sampleID)
  
  ## assign column names  
  colnames(otutable2)[1:(dim(otutable2)[2]-1)] <- clusterIDs
  colnames(otutable2)[which(colnames(otutable2)=="")] = "No designation"
  
  ## merge with metadata
  cortable <- read.table(mapper_file, comment.char = "", header = T, fill = T, sep = "\t") %>% 
    inner_join(otutable2)
  
  if(dim(cortable)[1] == 0) {
    cortable <- read.table(mapper_file, comment.char = "", header = T, fill = T, sep = "\t") %>% 
      mutate(X.SampleID = paste0("X",X.SampleID)) %>%
      inner_join(otutable2) %>%
      mutate(X.SampleID = gsub("X", replacement = "", X.SampleID))
  }
  
  colnames(cortable)[which(colnames(cortable)=="")] = "No designation"
  
  
  ## average, if applicable
  
#  averaged = c("orchard")
  if (averaged != F) {
    int_table <- read.table(mapper_file, comment.char = "", header = T, fill = T, sep = "\t") %>% tibble() %>% tidyr::unite(col = "newcol", all_of(averaged)) %>% data.frame() %>%
      inner_join(otutable2)
    
    colnames(int_table)[which(colnames(int_table)=="")] = "No designation"
    
    int_table <- int_table %>%
      group_by(newcol) %>%
      summarise(across(colnames(otutable2), mean))
    
    cortable <- read.table(mapper_file, comment.char = "", header = T, fill = T, sep = "\t") %>% tidyr::unite(col = "newcol", all_of(averaged)) %>% group_by(newcol) %>% dplyr::sample_n(size = 1) %>% data.frame() %>% 
    inner_join(int_table %>% dplyr::select(-X.SampleID), by = "newcol")
    
  
  }
  
    ## a marker to determine where to stop variables to correlate and where taxonomic abundances begin  
  num_tests <- dim(cortable)[2] - dim(otutable2)[2] + 1
  
  sig <- data.frame(cluster = character(), pval = numeric(), rho = numeric(), sum = numeric(), fdr = numeric(), test = character(), taxonomy = character())
  
  #print(num_tests)
  #ggplot(cortable, aes(x = time, y = ` f__Lactobacillaceae`))+ geom_point() + geom_smooth(method = "lm") + theme_classic()
#cortable[,10]
  
#  print(str(cortable[,1:num_tests]))
  for (i in 1:num_tests) {
   # cat(paste0("testing ",i,"..."))
    if (class(cortable[,i])%in%(c("numeric", "integer"))) {
      #print(colnames(cortable[,i]))
      if(!is.na(sum(cortable[,i] >0, na.rm=T)) & sum(cortable[,i], na.rm=T)>0) {
        dir.create(paste0('core-metrics-results-',file_path,"/",colnames(cortable)[i]))
        out_df <- data.frame(cluster = character(), pval = numeric(), rho = numeric(), sum = numeric())
        for (k in (num_tests + 1) : dim(cortable)[2]) {
          if(sum(cortable[,k]>0)/length(cortable[,k]) > frac_samples) {
            new_df <- data.frame(cluster = as.character(colnames(cortable)[k]),
                                 pval = as.numeric(cor.test(cortable[,i],cortable[,k],method = "spearman", exact = F)$p.value),
                                 rho = as.numeric(cor.test(cortable[,i],cortable[,k],method = "spearman", exact = F)$estimate),
                                 sum = as.numeric(sum(cortable[,k]))
            )
            out_df <- rbind(out_df, new_df)
          }
        }
        out_df <- out_df %>%
          mutate(fdr = p.adjust(pval, method = "fdr")) %>%
          arrange(pval)
        write.csv(out_df, paste0('core-metrics-results-',file_path,"/",colnames(cortable)[i],"/",tail(j,1),".csv"))
        if(use_fdr == T) {
          sig <- rbind(sig, out_df %>%
                         filter(fdr < 0.05) %>%
                         mutate(test = as.character(colnames(cortable)[i]),
                                taxonomy = as.character(tail(j,1)))
          )
        } else {
          sig <- rbind(sig, out_df %>%
                         filter(pval < 0.05) %>%
                         mutate(test = as.character(colnames(cortable)[i]),
                                taxonomy = as.character(tail(j,1)))
          )
        }        
      }
    } else {
     # print(paste0("Can't test column",colnames(cortable[,i])))
    }
  }
  
  
  write.csv(sig, paste0('core-metrics-results-',file_path,"/all_sig_",tail(j,1),".csv"))
  if(return_table == T) {
    return(cortable)
  }
}
}


