## diets over time

mantel_distance_plot <- function(mantel_output, x_title = "microbiota distance (Bray-Curtis)", y_title = y_title, point_size = 1) {
  
  aab <- melt(as.matrix(mantel_output[[1]]))[melt(upper.tri(as.matrix(mantel_output[[1]])))$value,] %>% mutate(Var12 = paste0(Var1,Var2))
  aac <- melt(as.matrix(mantel_output[[2]]))[melt(upper.tri(as.matrix(mantel_output[[2]])))$value,] %>% mutate(Var1 = gsub("X","",Var1), Var2 = gsub("X","",Var2)) %>% mutate(Var12 = paste0(Var1,Var2))
  
  aad <- aab %>% inner_join(aac, by = "Var12")
  
  utmantel <- ggplot(aad, aes(x = value.x, y = value.y)) + 
    geom_jitter(size=point_size, width = max(aad$value.x)/75, alpha=0.5) + 
    stat_smooth(method = "lm", se=F) +
    theme_cowplot() + theme_pubclean() + labs(x = x_title, y = y_title )
  
  return(utmantel)
}


## fruits, starved, soils
check_grouping <- function (dataframe) {
  cat("Wolbachia groups: ")
  cat(names(table(dataframe$wolbachia)))
  print("")
  cat("Starved groups: ")
  cat(names(table(dataframe$starved)))
  print("")
  cat("Orchard subsites with more than one group:")
  cat(sum(table(dataframe$orchard_subsite)>1))
  #  print(table(dataframe$orchard_subsite)[table(dataframe$orchard_subsite)>1])
}

mantel_2groups <- function (file_path, file_path2, mapper_file, metric) {
  ## unweighted
  flies_unweighted <- read.table(paste('core-metrics-results-',file_path,'/',metric,'_distance_matrix/distance-matrix.tsv',sep=""), header=T, sep="\t") %>% mutate(X=as.character(X))
  flies_unifrac_table <- flies_unweighted %>% dplyr::select(X) %>% left_join(read.table(mapper_file,comment.char = "", header=T, fill=T, sep="\t"), by=c("X"="X.SampleID")) %>% droplevels()
  table(list(flies_unifrac_table$orchard_subsite)) # no duplicates
  colnames(flies_unweighted) <- c("X",flies_unifrac_table$orchard_subsite)
  flies_unweighted_dm <- as.dist(flies_unweighted[,2:dim(flies_unweighted)[2]])
  
  
  
  flies_unweighted2 <- read.table(paste('core-metrics-results-',file_path2,'/',metric,'_distance_matrix/distance-matrix.tsv',sep=""), header=T, sep="\t") %>% mutate(X=as.character(X))
  flies_unifrac_table2 <- flies_unweighted2 %>% dplyr::select(X) %>% left_join(read.table(mapper_file,comment.char = "", header=T, fill=T, sep="\t"), by=c("X"="X.SampleID")) %>% droplevels()
  table(list(flies_unifrac_table2$orchard_subsite)) # no duplicates
  colnames(flies_unweighted2) <- c("X",flies_unifrac_table2$orchard_subsite)
  flies_unweighted2 <- flies_unweighted2[,colnames(flies_unweighted)]
  flies_unweighted_dm2 <- as.dist(flies_unweighted2[,2:dim(flies_unweighted2)[2]])
  
  
  set.seed(43)
  mantel1 <- mantel(flies_unweighted_dm,flies_unweighted_dm2, method = "pearson", parallel = 7)
  
  cat(dim(flies_unweighted_dm))
  print(paste0("rho: ",round(mantel1$statistic,3),"; p: ",mantel1$signif))
  plot(flies_unweighted_dm,flies_unweighted_dm2)
  
}
mantel_2groups_all3 <- function (file_path, file_path2, mapper_file) {
  for(i in c("unweighted_unifrac","weighted_unifrac","bray_curtis")) {
    mantel_2groups(file_path,file_path2,mapper_file,i)
  }
}

## venn diagrams
get_core_features <- function(file_path = file_path, map_path = map_path, mapper_file = mapper_file, taxonomic_level = taxonomic_level, similarity_cutoff = 0.5, var1 = var1, var2 = var2, var3 = var3, var4 = var4, var5 = NULL, Group = Group, the_id = "X.SampleID") {
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
      dplyr::select(all_of(c(the_id,var1,var2,var3,var4, Group)))
  } else {
    map2 <- read.table(mapper_file,comment.char = "", header=T, fill=T,sep="\t") %>% 
      dplyr::select(all_of(c(the_id,var1,var2,var3,var4, var5,Group)))
  }
  #print(head(map2))
  
  phyl2 <- otu_table %>% 
    inner_join(otu2, by = "X.OTU.ID") %>%
    group_by(clustered_taxonomy) %>% 
    summarize_at(ref_names,sum, na.rm=T)
  rownames(phyl2) <- as.character(unlist(phyl2$clustered_taxonomy))
  
  phyl3 <- phyl2[,-1]
  
  hist(rowSums(phyl3>0)/dim(phyl3)[2], breaks=10)
  # max(rowSums(phyl3>0)/dim(phyl3)[2])
  # min(rowSums(phyl3>0)/dim(phyl3)[2])
  
  return(phyl2$clustered_taxonomy[rowSums(phyl3>0)/dim(phyl3)[2]>similarity_cutoff])
}
get_single_core_features <- function(file_path = file_path, map_path = map_path, mapper_file = mapper_file, taxonomic_level = taxonomic_level, var1 = var1, var2 = var2, var3 = var3, var4 = var4, var5 = NULL, Group = Group, the_id = "X.SampleID") {
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
  
  map2 <- read.table(mapper_file,comment.char = "", header=T, fill=T,sep="\t") %>% dplyr::select(any_of(c(the_id,var1, var2, var3, var4, var5, Group)))
  
  phyl2 <- otu_table %>% 
    inner_join(otu2, by = "X.OTU.ID") %>%
    group_by(clustered_taxonomy) %>% 
    summarize_at(ref_names,sum, na.rm=T)
  rownames(phyl2) <- as.character(unlist(phyl2$clustered_taxonomy))
  
  return(phyl2)
}
tabulate_individual_venns <- function(file_path = file_path, taxonomic_level = taxonomic_level, otu_table = otu_table, threshhold = 1) {
  
  starved_fruit_venn <- get_single_core_features(file_path = file_path, map_path = map_path, mapper_file = mapper_file, taxonomic_level = taxonomic_level,var1=var1, var2=var2, var3=var3,var4=var4, Group = Group)
  starved_fruit_venn$clustered_taxonomy = gsub("[", replacement = "", starved_fruit_venn$clustered_taxonomy, perl = T, fixed=T)
  starved_fruit_venn$clustered_taxonomy = gsub("]", replacement = "", starved_fruit_venn$clustered_taxonomy, perl = T, fixed=T)
  
  orchard_site_names <- otu_table %>% filter(X.SampleID %in% colnames(starved_fruit_venn)) %>% dplyr::pull(orchard_subsite) %>% table() %>% data.frame()
  
  group1name = strsplit(file_path, split = "_", fixed = T)[[1]][2]
  group2name = strsplit(file_path, split = "_", fixed = T)[[1]][3]
  
  out_df <- data.frame(orchard_subsite = character(), group1 = character(), group2 = character(), unique1 = numeric, shared = numeric(), unique2 = numeric)
  
  i=orchard_site_names$.[1]
  for (i in orchard_site_names$.) {
    sample_names <- otu_table %>% filter(orchard_subsite == i) %>% dplyr::pull(X.SampleID)
    filtered_names <- sample_names[sample_names %in% colnames(starved_fruit_venn)]
    two_way <- data.frame(starved_fruit_venn) %>% dplyr::select(clustered_taxonomy, all_of(filtered_names))
    groupinl1 = (two_way %>% filter(.[[2]] > threshhold) %>% dplyr::pull(clustered_taxonomy))
    groupinl2 = two_way %>% filter(.[[3]] > threshhold) %>% dplyr::pull(clustered_taxonomy)
    try(two_way_venns <- VennDiagram::get.venn.partitions(x = list(group1 = groupinl1, group2 = groupinl2)))
    if(class(try(two_way_venns <- VennDiagram::get.venn.partitions(x = list(group1 = groupinl1, group2 = groupinl2)),F))!="try-error") {
      two_way_venns <- VennDiagram::get.venn.partitions(x = list(group1 = groupinl1, group2 = groupinl2))
      out_df <- rbind(out_df, data.frame(orchard_subsite = as.character(i), 
                                         group1 = as.character(group1name), 
                                         group2 = as.character(group2name), 
                                         unique1 = as.numeric(two_way_venns %>% filter(group1 == T, group2 == F) %>% dplyr::pull(..count..)), 
                                         shared = as.numeric(two_way_venns %>% filter(group1 == T, group2 == T) %>% dplyr::pull(..count..)),
                                         unique2 = as.numeric(two_way_venns %>% filter(group1 == F, group2 == T) %>% dplyr::pull(..count..)))
      )
      
    }
  }
  return(out_df)
}


## neutral modelling
make_taxon_plot_matrix <- function(file_path, file_path2="") {
  otumat = read.table(paste0('core-metrics-results-',file_path,'/rarefied_table',file_path2,'.txt'), comment.char = "", header=T, fill=T, sep = "\t", skip = 1)
  rownames(otumat) <- otumat$X.OTU.ID
  otumat <- otumat %>% dplyr::select(-X.OTU.ID)
  otumat <- as.matrix(otumat)
  otumat
}


sncm.fit <- function(spp, pool=NULL, stats=TRUE, taxon=NULL){
  require(minpack.lm)
  require(Hmisc)
  require(stats4)
  
  options(warn=-1)
  
  #Calculate the number of individuals per community
  N <- mean(apply(spp, 1, sum))
  
  #Calculate the average relative abundance of each taxa across communities
  if(is.null(pool)){
    p.m <- apply(spp, 2, mean)
    p.m <- p.m[p.m != 0]
    p <- p.m/N
  } else {
    p.m <- apply(pool, 2, mean)
    p.m <- p.m[p.m != 0]
    p <- p.m/N
  }
  
  #Calculate the occurrence frequency of each taxa across communities
  spp.bi <- 1*(spp>0)
  freq <- apply(spp.bi, 2, mean)
  freq <- freq[freq != 0]
  
  #Combine
  C <- merge(p, freq, by=0)
  C <- C[order(C[,2]),]
  C <- as.data.frame(C)
  C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),] #Removes rows with any zero (absent in either source pool or local communities)
  p <- C.0[,2]
  freq <- C.0[,3]
  names(p) <- C.0[,1]
  names(freq) <- C.0[,1]
  
  #Calculate the limit of detection
  d = 1/N
  
  ##Fit model parameter m (or Nm) using Non-linear least squares (NLS)
  m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1-p), lower.tail=FALSE), start=list(m=0.1))
  m.ci <- confint(m.fit, 'm', level=0.95)
  
  ##Fit neutral model parameter m (or Nm) using Maximum likelihood estimation (MLE)
  sncm.LL <- function(m, sigma){
    R = freq - pbeta(d, N*m*p, N*m*(1-p), lower.tail=FALSE)
    R = dnorm(R, 0, sigma)
    -sum(log(R))
  }
  m.mle <- mle(sncm.LL, start=list(m=0.1, sigma=0.1), nobs=length(p))
  
  ##Calculate Akaike's Information Criterion (AIC)
  aic.fit <- AIC(m.mle, k=2)
  bic.fit <- BIC(m.mle)
  
  ##Calculate goodness-of-fit (R-squared and Root Mean Squared Error)
  freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1-p), lower.tail=FALSE)
  Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
  RMSE <- sqrt(sum((freq-freq.pred)^2)/(length(freq)-1))
  
  pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
  
  ##Calculate AIC for binomial model
  bino.LL <- function(mu, sigma){
    R = freq - pbinom(d, N, p, lower.tail=FALSE)
    R = dnorm(R, mu, sigma)
    -sum(log(R))
  }
  bino.mle <- mle(bino.LL, start=list(mu=0, sigma=0.1), nobs=length(p))
  
  aic.bino <- AIC(bino.mle, k=2)
  bic.bino <- BIC(bino.mle)
  
  ##Goodness of fit for binomial model
  bino.pred <- pbinom(d, N, p, lower.tail=FALSE)
  Rsqr.bino <- 1 - (sum((freq - bino.pred)^2))/(sum((freq - mean(freq))^2))
  RMSE.bino <- sqrt(sum((freq - bino.pred)^2)/(length(freq) - 1))
  
  bino.pred.ci <- binconf(bino.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
  
  ##Calculate AIC for Poisson model
  pois.LL <- function(mu, sigma){
    R = freq - ppois(d, N*p, lower.tail=FALSE)
    R = dnorm(R, mu, sigma)
    -sum(log(R))
  }
  pois.mle <- mle(pois.LL, start=list(mu=0, sigma=0.1), nobs=length(p))
  
  aic.pois <- AIC(pois.mle, k=2)
  bic.pois <- BIC(pois.mle)
  
  ##Goodness of fit for Poisson model
  pois.pred <- ppois(d, N*p, lower.tail=FALSE)
  Rsqr.pois <- 1 - (sum((freq - pois.pred)^2))/(sum((freq - mean(freq))^2))
  RMSE.pois <- sqrt(sum((freq - pois.pred)^2)/(length(freq) - 1))
  
  pois.pred.ci <- binconf(pois.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
  
  ##Results
  if(stats==TRUE){
    fitstats <- data.frame(m=numeric(), m.ci=numeric(), m.mle=numeric(), maxLL=numeric(), binoLL=numeric(), poisLL=numeric(), Rsqr=numeric(), Rsqr.bino=numeric(), Rsqr.pois=numeric(), RMSE=numeric(), RMSE.bino=numeric(), RMSE.pois=numeric(), AIC=numeric(), BIC=numeric(), AIC.bino=numeric(), BIC.bino=numeric(), AIC.pois=numeric(), BIC.pois=numeric(), N=numeric(), Samples=numeric(), Richness=numeric(), Detect=numeric())
    fitstats[1,] <- c(coef(m.fit), coef(m.fit)-m.ci[1], m.mle@coef['m'], m.mle@details$value, bino.mle@details$value, pois.mle@details$value, Rsqr, Rsqr.bino, Rsqr.pois, RMSE, RMSE.bino, RMSE.pois, aic.fit, bic.fit, aic.bino, bic.bino, aic.pois, bic.pois, N, nrow(spp), length(p), d)
    return(fitstats)
  } else {
    A <- cbind(p, freq, freq.pred, pred.ci[,2:3], bino.pred, bino.pred.ci[,2:3])
    A <- as.data.frame(A)
    colnames(A) <- c('p', 'freq', 'freq.pred', 'pred.lwr', 'pred.upr', 'bino.pred', 'bino.lwr', 'bino.upr')
    if(is.null(taxon)){
      B <- A[order(A[,1]),]
    } else {
      B <- merge(A, taxon, by=0, all=TRUE)
      row.names(B) <- B[,1]
      B <- B[,-1]
      B <- B[order(B[,1]),]
    }
    return(B)
  }
}


## from JoeExperiment.Rmd
mantel_distance_plot <- function(mantel_output, x_title = "microbiota distance (Bray-Curtis)", y_title = y_title, point_size = 1) {
  
  aab <- melt(as.matrix(mantel_output[[1]]))[melt(upper.tri(as.matrix(mantel_output[[1]])))$value,] %>% mutate(Var12 = paste0(Var1,Var2))
  aac <- melt(as.matrix(mantel_output[[2]]))[melt(upper.tri(as.matrix(mantel_output[[2]])))$value,] %>% mutate(Var1 = gsub("X","",Var1), Var2 = gsub("X","",Var2)) %>% mutate(Var12 = paste0(Var1,Var2))
  
  aad <- aab %>% inner_join(aac, by = "Var12")
  
  utmantel <- ggplot(aad, aes(x = value.x, y = value.y)) + 
    geom_jitter(size=point_size, width = max(aad$value.x)/75, alpha=0.5) + 
    stat_smooth(method = "lm", se=F) +
    theme_cowplot() + theme_pubclean() + labs(x = x_title, y = y_title )
  
  return(utmantel)
}

make_taxonomy_forR <- function(taxon_file_path) {
  taxa_table <- read.table(paste0(taxon_file_path,"/taxonomy.tsv"), sep = "\t", header = T) %>% 
    tidyr::separate(col = Taxon, into = c("kingdom","phylum","class","order","family","genus","species"), remove = T, sep = ";") %>%
    dplyr::select(-Confidence)
  write.csv(taxa_table, file = paste0(taxon_file_path,"/taxonomy_forR.csv"), quote = F, sep = ",", row.names = F)
}
