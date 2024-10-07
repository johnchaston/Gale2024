make_taxon_plot_matrix <- function(file_path) {
  otumat = read.table(paste0('core-metrics-results-',filepath,'/rarefied_table.txt'), comment.char = "", header=T, fill=T, sep = "\t", skip = 1)
  rownames(otumat) <- otumat$X.OTU.ID
  otumat <- otumat %>% dplyr::select(-X.OTU.ID)
  otumat <- as.matrix(otumat)
  otumat
}

make_taxonomy_forR <- function(taxon_file_path) {
  taxa_table <- read.table(paste0(taxon_file_path,"/taxonomy.tsv"), sep = "\t", header = T) %>% 
    tidyr::separate(col = Taxon, into = c("kingdom","phylum","class","order","family","genus","species"), remove = T, sep = ";") %>%
    dplyr::select(-Confidence)
  write.csv(taxa_table, file = paste0(taxon_file_path,"/taxonomy_forR.csv"), quote = F, sep = ",", row.names = F)
}

calc_mantel_microbiome_fromCFUs <- function(file_path = file_path, mapper_file = mapper_file, mapper_cols = mapper_cols) {
  
  ## make the micro
  dm <- read.table(paste(file_path,'/core-metrics-results-',file_path,'/bray_curtis_distance_matrix/distance-matrix.tsv',sep=""), header=T, sep="\t") %>% mutate(X=as.character(X))
  dmdist <- as.dist(dm[,2:dim(dm)[2]])
  
  ## make the metadata matrix
  metadata <- dm %>% 
    dplyr::select(X) %>% 
    left_join(read.table(mapper_file,comment.char = "", header=T, fill=T, sep="\t"), by=c("X"="X.SampleID")) %>% 
    droplevels() %>% 
    dplyr::select(X,any_of(mapper_cols))

  fut <- as.matrix(metadata %>% dplyr::select(-X))
  rownames(fut) <- metadata$X
  futmat <- dist(fut)
  
  cat("dimensions: ")
  print(dim(futmat))
  set.seed(43)
  mantel1 <- mantel(futmat,dmdist, method = "spearman", parallel = 7)
  print(paste0("rho: ",round(mantel1$statistic,3),"; p: ",mantel1$signif))
  
  plot(futmat,dmdist)
  
  ## unrarefied
  dm <- read.table(paste(file_path,'/bray_curtis_dm_',file_path,'/distance-matrix.tsv',sep=""), header=T, sep="\t") %>% mutate(X=as.character(X))
  dmdist <- as.dist(dm[,2:dim(dm)[2]])
  
  ## make the metadata matrix
  metadata <- dm %>% 
    dplyr::select(X) %>% 
    left_join(read.table(mapper_file,comment.char = "", header=T, fill=T, sep="\t"), by=c("X"="X.SampleID")) %>% 
    droplevels() %>% 
    dplyr::select(X, any_of(mapper_cols))
  
  fut <- as.matrix(metadata %>% dplyr::select(-X))
  rownames(fut) <- metadata$X
  futmat <- dist(fut)
  cat("dimensions: ")
  print(dim(futmat))
  
  set.seed(43)
  mantel1 <- mantel(futmat,dmdist, method = "spearman", parallel = 7)
  print(paste0("rho: ",round(mantel1$statistic,3),"; p: ",mantel1$signif))
  
  plot(futmat,dmdist)
  
}

make_taxon_plot_matrix <- function(file_path) {
  otumat = read.table(paste0('core-metrics-results-',file_path,'/rarefied_table.txt'), comment.char = "", header=T, fill=T, sep = "\t", skip = 1)
  rownames(otumat) <- otumat$X.OTU.ID
  otumat <- otumat %>% dplyr::select(-X.OTU.ID)
  otumat <- as.matrix(otumat)
  otumat
}

make_taxonomy_forR <- function(taxon_file_path) {
  taxa_table <- read.table(paste0(taxon_file_path,"/taxonomy.tsv"), sep = "\t", header = T) %>% 
    tidyr::separate(col = Taxon, into = c("kingdom","phylum","class","order","family","genus","species"), remove = T, sep = ";") %>%
    dplyr::select(-Confidence)
  write.csv(taxa_table, file = paste0(taxon_file_path,"/taxonomy_forR.csv"), quote = F, sep = ",", row.names = F)
}

#setwd('~/Dropbox/sequencing/2022-09-henry/')
#mapper_cols <- c("AvgRH","latitude","date","MaxTemp", "MaxcDNI")
#mapper_file = "mantel_henrymapper2_JMC.txt"
#file_path = "flies_mantel2v2"
#metric="bray_curtis"

calc_mantel_microbiome <- function(file_path = file_path, metric = "bray_curtis", mapper_file = mapper_file, mapper_cols = mapper_cols, filter = F) {
  
  ## make the micro
  dm <- read.table(paste('core-metrics-results-',file_path,'/',metric,'_distance_matrix/distance-matrix.tsv',sep=""), header=T, sep="\t") %>% mutate(X=as.character(X))
  dmdist <- as.dist(dm[,2:dim(dm)[2]])
  
  ## make the metadata matrix
  metadata <- dm %>% 
    dplyr::select(X) %>% 
    inner_join(read.table(mapper_file,comment.char = "", header=T, fill=T, sep="\t"), by=c("X"="X.SampleID")) %>% 
    droplevels() 
  metadata$Latitude
  fut <- as.matrix(metadata %>% dplyr::select(all_of(mapper_cols)))
  fut2 <- (fut)
  rownames(fut2) <- metadata$X
  futmat <- dist(fut2)
  
  if (filter == T) {
    dmdist <- as.dist(dm %>% filter(X %in% metadata$X) %>% dplyr::select(all_of(metadata$X)))
    cat("filtering on ...")
  }
  
  cat("dimensions: ")
  print(dim(futmat))
  print(dim(dmdist))
  
  sum(is.na(futmat))
  sum(is.na(dmdist))
  
  set.seed(43)
  mantel1 <- mantel(log10(futmat+1),log10(dmdist+1), method = "spearman", parallel = 7)
  print(paste0("rho: ",round(mantel1$statistic,3),"; p: ",mantel1$signif))
  # set.seed(43)
  # mantel2 <- mantel(dmdist,futmat, method = "spearman", parallel = 7)
  # print(paste0("rho: ",round(mantel2$statistic,3),"; p: ",mantel2$signif))
  
  #if(plot ==T) {
  plot(futmat,dmdist)
  #}
  
  return(list(futmat,dmdist))
  
}

##calc_mantel_microbiome <- function(file_path = file_path, metric = "bray_curtis", mapper_file = mapper_file, mapper_cols = mapper_cols) {
# 
#  ## make the micro
#  dm <- read.table(paste('core-metrics-results-',file_path,'/',metric,'_distance_matrix/distance-matrix.tsv',sep=""), header=T, sep="\t") %>% mutate(X=as.character(X))
#  dmdist <- as.dist(dm[,2:dim(dm)[2]])
#  
#  ## make the metadata matrix
#  metadata <- dm %>% 
#    dplyr::select(X) %>% 
#    left_join(read.table(mapper_file,comment.char = "", header=T, fill=T, sep="\t"), by=c("X"="X.SampleID")) %>% 
#    droplevels() 
#  fut <- as.matrix(metadata %>% dplyr::select(all_of(mapper_cols)))
#  fut2 <- (fut)
#  rownames(fut2) <- metadata$X
#  futmat <- dist(fut2)
#  cat("dimensions: ")
#  print(dim(futmat))
#  
#  set.seed(43)
#  mantel1 <- mantel(log10(futmat+1),log10(dmdist+1), method = "spearman", parallel = 7)
#  print(paste0("rho: ",round(mantel1$statistic,3),"; p: ",mantel1$signif))
#  # set.seed(43)
#  # mantel2 <- mantel(dmdist,futmat, method = "spearman", parallel = 7)
#  # print(paste0("rho: ",round(mantel2$statistic,3),"; p: ",mantel2$signif))
#  
#  #if(plot ==T) {
#  plot(futmat,dmdist)
#  #}
#  
#  return(list(futmat,dmdist))
#  
#}


# this was the original, i modified the version above to use only the specified mapper cols
# calc_mantel_microbiome <- function(file_path = file_path, metric = "bray_curtis", mapper_file = mapper_file, mapper_cols = mapper_cols) {
#   
#   ## make the micro
#   dm <- read.table(paste('core-metrics-results-',file_path,'/',metric,'_distance_matrix/distance-matrix.tsv',sep=""), header=T, sep="\t") %>% mutate(X=as.character(X))
#   dmdist <- as.dist(dm[,2:dim(dm)[2]])
#   
#   ## make the metadata matrix
#   metadata <- dm %>% 
#     dplyr::select(X) %>% 
#     left_join(read.table(mapper_file,comment.char = "", header=T, fill=T, sep="\t"), by=c("X"="X.SampleID")) %>% 
#     droplevels() 
#   
#   
#   metadata
#   fut <- as.matrix(metadata %>% dplyr::select(-X))
#   fut2 <- (fut)
#   rownames(fut2) <- metadata$X
#   futmat <- dist(fut2)
#   cat("dimensions: ")
#   print(dim(futmat))
#   
#   set.seed(43)
#   mantel1 <- mantel(futmat,dmdist, method = "pearson", parallel = 7)
#   print(paste0("rho: ",round(mantel1$statistic,3),"; p: ",mantel1$signif))
#   # set.seed(43)
#   # mantel2 <- mantel(dmdist,futmat, method = "spearman", parallel = 7)
#   # print(paste0("rho: ",round(mantel2$statistic,3),"; p: ",mantel2$signif))
#   
#   #if(plot ==T) {
#   plot(futmat,dmdist)
#   #}
#   
#   return(list(futmat,dmdist))
#   
# }


calc_mantel_microbiome_original <- function(file_path = file_path, metric = "bray_curtis", mapper_file = mapper_file, mapper_cols = mapper_cols) {
  
  ## make the micro
  dm <- read.table(paste('core-metrics-results-',file_path,'/',metric,'_distance_matrix/distance-matrix.tsv',sep=""), header=T, sep="\t") %>% mutate(X=as.character(X))
  dmdist <- as.dist(dm[,2:dim(dm)[2]])
  
  ## make the metadata matrix
  metadata <- dm %>% 
    dplyr::select(X) %>% 
    left_join(read.table(mapper_file,comment.char = "", header=T, fill=T, sep="\t"), by=c("X"="X.SampleID")) %>% 
    droplevels() %>% 
    dplyr::select(X,any_of(mapper_cols))
  
  fut <- as.matrix(metadata %>% dplyr::select(-X))
  rownames(fut) <- metadata$X
  futmat <- dist(fut)
  cat("dimensions: ")
  print(dim(futmat))
  
  set.seed(43)
  mantel1 <- mantel(futmat,dmdist, method = "spearman", parallel = 7)
  print(paste0("rho: ",round(mantel1$statistic,3),"; p: ",mantel1$signif))
  # set.seed(43)
  # mantel2 <- mantel(dmdist,futmat, method = "spearman", parallel = 7)
  # print(paste0("rho: ",round(mantel2$statistic,3),"; p: ",mantel2$signif))
  
  #if(plot ==T) {
    plot(futmat,dmdist)
  #}
  
  return(list(futmat,dmdist))
  
}
calc_mantel_microbiome_matrix <- function(file_path = file_path, metric = "bray_curtis", mapper_file = mapper_file, mapper_cols = mapper_cols) {
  
  ## make the micro
  dm <- read.table(paste('core-metrics-results-',file_path,'/',metric,'_distance_matrix/distance-matrix.tsv',sep=""), header=T, sep="\t") %>% mutate(X=as.character(X))
  dmdist <- as.dist(dm[,2:dim(dm)[2]])
  
  ## make the metadata matrix
  metadata <- dm %>% 
    dplyr::select(X) %>% 
    left_join(read.table(mapper_file,comment.char = "", header=T, fill=T, sep="\t"), by=c("X"="X.SampleID")) %>% 
    droplevels() %>% 
    dplyr::select(X,any_of(mapper_cols))
  
  fut <- as.matrix(metadata %>% dplyr::select(-X))
  rownames(fut) <- metadata$X
  futmat <- dist(fut)
  cat("dimensions: ")
  print(dim(futmat))
  
  set.seed(43)
  mantel1 <- mantel(futmat,dmdist, method = "spearman", parallel = 7)
  print(paste0("rho: ",round(mantel1$statistic,3),"; p: ",mantel1$signif))
  # set.seed(43)
  # mantel2 <- mantel(dmdist,futmat, method = "spearman", parallel = 7)
  # print(paste0("rho: ",round(mantel2$statistic,3),"; p: ",mantel2$signif))
  
  #if(plot ==T) {
  plot(futmat,dmdist)
  #}
  
  return(list(futmat,dmdist))
  
}
fst_file = "orch17.fst.csv"
  file_path = "mantel_2017"
  metric = "bray_curtis"
mantel_2017 <- function(file_path=file_path, metric=metric, fst_file = fst_file) {  
  ## make the micro
  dm <- read.table(paste('core-metrics-results-',file_path,'/',metric,'_distance_matrix/distance-matrix.tsv',sep=""), header=T, sep="\t") %>% mutate(X=as.character(X))
  dmdist <- as.dist(dm[,2:dim(dm)[2]])
  cat("dimensions: ")
  print(dim(dmdist))
  
  ## make the metadata matrix
  metadata <- dm %>% 
    dplyr::select(X) %>% 
    inner_join(read.table(fst_file,comment.char = "", header=T, fill=T, sep=","), by=c("X")) %>% 
    droplevels() 
  
  fut <- as.matrix(metadata %>% dplyr::select(-X))
  fut2 <- (fut)
  rownames(fut2) <- metadata$X
  futmat <- dist(fut2)
  cat("dimensions: ")
  print(dim(futmat))
  
  set.seed(43)
  mantel1 <- mantel(futmat,dmdist, method = "pearson", parallel = 7)
  print(paste0("rho: ",round(mantel1$statistic,3),"; p: ",mantel1$signif))
  
  plot(futmat,dmdist)
  
  return(list(futmat,dmdist))
  
}

mantel_2016 <- function(file_path = file_path, metric = metric, fst_file = fst_file) {
  dm <- read.table(paste('core-metrics-results-',file_path,'/',metric,'_distance_matrix/distance-matrix.tsv',sep=""), header=T, sep="\t") %>% mutate(X=as.character(X)) %>% inner_join(read.table("metadata_2016.txt", comment.char = "", header=T, fill=T, sep = "\t"), by=c("X"="X.SampleID"))
  colnames(dm)[2:(dim(dm)[2]-1)] <- dm$sampID
  dmdist <- as.dist(dm[,2:(dim(dm)[2]-1)])
  
  ## make the metadata matrix
  fst <- read.table(fst_file, sep = ",", header=T) 
  colnames(fst)[2:(dim(fst)[2])] <- fst$X
  fst2 <- fst %>% 
    dplyr::select((colnames(dm %>% dplyr::select(-sampID)))) %>% 
    filter(X %in% (colnames(dm %>% dplyr::select(-sampID)))) %>%
    arrange(match(X, (colnames(dm %>% dplyr::select(-sampID))))) 
  fst_dm <- fst2 %>% dplyr::select(-X) %>% as.matrix()
  fst_dm2 <- as.dist(fst_dm)
  
  cat("dimensions: ")
  print(dim(fst_dm2))
  print(dim(dmdist))
  sum(names(fst_dm2) == names(dmdist))
  
  set.seed(43)
  mantel1 <- mantel(fst_dm2,dmdist, method = "pearson", parallel = 7)
  print(paste0("rho: ",round(mantel1$statistic,3),"; p: ",mantel1$signif))
  plot(fst_dm2,dmdist)
  return(list(futmat,dmdist))
  
}



calc_nrel <- function(file_path, month_val, day_val, specify_location) {
  df1 <- read.table(file_path, header = T, skip = 2, fill = T, sep = ",") %>%
    filter(Month == month_val, Day == day_val) %>%
    dplyr::summarize(MaxTemp = max(Temperature), MinTemp = min(Temperature), AvgTemp = mean(Temperature), 
                     MaxcDHI = max(Clearsky.DHI), AvgcDHI = mean(Clearsky.DHI), 
                     MaxcDNI = max(Clearsky.DNI), AvgcDNI = mean(Clearsky.DNI), 
                     MaxcGHI = max(Clearsky.GHI), AvgcGHI = mean(Clearsky.GHI),  
                     MaxDP = max(Dew.Point), MinDP = min(Dew.Point), AvgDP = mean(Dew.Point), 
                     MaxDHI = max(DHI), AvgDHI = mean(DHI), 
                     MaxcDNI = max(DNI), AvgcDNI = mean(DNI), 
                     MaxFF = max(Fill.Flag), AvgFF = mean(Fill.Flag), 
                     MaxcGHI = max(GHI), AvgcGHI = mean(GHI),  
                     MaxRH = max(Relative.Humidity), MinRH = min(Relative.Humidity), AvgRH = mean(Relative.Humidity), 
                     MaxSZA = max(Solar.Zenith.Angle), MinSZA = min(Solar.Zenith.Angle), AvgSZA = mean(Solar.Zenith.Angle), 
                     MaxSA = max(Surface.Albedo), MinSA = min(Surface.Albedo), AvgSA = mean(Surface.Albedo), 
                     MaxPress = max(Pressure), MinPress = min(Pressure), AvgPres = mean(Pressure), 
                     MaxPW = max(Precipitable.Water), MinPW = min(Precipitable.Water), AvgPW = mean(Precipitable.Water), 
                     MaxWD = max(Wind.Direction), MinWD = min(Wind.Direction), AvgWD = mean(Wind.Direction), 
                     MaxWS = max(Wind.Speed), MinWS = min(Wind.Speed), AvgWS = mean(Wind.Speed), 
                     MaxUVO = max(Global.Horizontal.UV.Irradiance..280.400nm.), AvgUVO = mean(Global.Horizontal.UV.Irradiance..280.400nm.), 
                     MaxUVI = max(Global.Horizontal.UV.Irradiance..295.385nm.), AvgUVI = mean(Global.Horizontal.UV.Irradiance..295.385nm.)
    ) %>%
    mutate(location = specify_location) %>%
    relocate(location)
  return(df1)
}
calc_nrel_daylight <- function(file_path, month_val, day_val, specify_location,uncorrected_daylight,atmospheric_corrected_daylight) {
  df1 <- read.table(file_path, header = T, skip = 2, fill = T, sep = ",") %>%
    filter(Month == month_val, Day == day_val) %>%
    dplyr::summarize(MaxTemp = max(Temperature), MinTemp = min(Temperature), AvgTemp = mean(Temperature), 
                     MaxcDHI = max(Clearsky.DHI), AvgcDHI = mean(Clearsky.DHI), 
                     MaxcDNI = max(Clearsky.DNI), AvgcDNI = mean(Clearsky.DNI), 
                     MaxcGHI = max(Clearsky.GHI), AvgcGHI = mean(Clearsky.GHI),  
                     MaxDP = max(Dew.Point), MinDP = min(Dew.Point), AvgDP = mean(Dew.Point), 
                     MaxDHI = max(DHI), AvgDHI = mean(DHI), 
                     MaxcDNI = max(DNI), AvgcDNI = mean(DNI), 
                     MaxFF = max(Fill.Flag), AvgFF = mean(Fill.Flag), 
                     MaxcGHI = max(GHI), AvgcGHI = mean(GHI),  
                     MaxRH = max(Relative.Humidity), MinRH = min(Relative.Humidity), AvgRH = mean(Relative.Humidity), 
                     MaxSZA = max(Solar.Zenith.Angle), MinSZA = min(Solar.Zenith.Angle), AvgSZA = mean(Solar.Zenith.Angle), 
                     MaxSA = max(Surface.Albedo), MinSA = min(Surface.Albedo), AvgSA = mean(Surface.Albedo), 
                     MaxPress = max(Pressure), MinPress = min(Pressure), AvgPres = mean(Pressure), 
                     MaxPW = max(Precipitable.Water), MinPW = min(Precipitable.Water), AvgPW = mean(Precipitable.Water), 
                     MaxWD = max(Wind.Direction), MinWD = min(Wind.Direction), AvgWD = mean(Wind.Direction), 
                     MaxWS = max(Wind.Speed), MinWS = min(Wind.Speed), AvgWS = mean(Wind.Speed), 
                     MaxUVO = max(Global.Horizontal.UV.Irradiance..280.400nm.), AvgUVO = mean(Global.Horizontal.UV.Irradiance..280.400nm.), 
                     MaxUVI = max(Global.Horizontal.UV.Irradiance..295.385nm.), AvgUVI = mean(Global.Horizontal.UV.Irradiance..295.385nm.)
    ) %>%
    mutate(location = specify_location, uncorrected_daylight = uncorrected_daylight, atmospheric_corrected_daylight = atmospheric_corrected_daylight) %>%
    relocate(location)
  return(df1)
}

calc_nrel <- function(file_path, month_val, day_val, specify_location) {
  df1 <- read.table(file_path, header = T, skip = 2, fill = T, sep = ",") %>%
    filter(Month == month_val, Day == day_val) %>%
    dplyr::summarize(MaxTemp = max(Temperature), MinTemp = min(Temperature), AvgTemp = mean(Temperature), 
                     MaxcDHI = max(Clearsky.DHI), AvgcDHI = mean(Clearsky.DHI), 
                     MaxcDNI = max(Clearsky.DNI), AvgcDNI = mean(Clearsky.DNI), 
                     MaxcGHI = max(Clearsky.GHI), AvgcGHI = mean(Clearsky.GHI),  
                     MaxDP = max(Dew.Point), MinDP = min(Dew.Point), AvgDP = mean(Dew.Point), 
                     MaxDHI = max(DHI), AvgDHI = mean(DHI), 
                     MaxcDNI = max(DNI), AvgcDNI = mean(DNI), 
                     MaxFF = max(Fill.Flag), AvgFF = mean(Fill.Flag), 
                     MaxcGHI = max(GHI), AvgcGHI = mean(GHI),  
                     MaxRH = max(Relative.Humidity), MinRH = min(Relative.Humidity), AvgRH = mean(Relative.Humidity), 
                     MaxSZA = max(Solar.Zenith.Angle), MinSZA = min(Solar.Zenith.Angle), AvgSZA = mean(Solar.Zenith.Angle), 
                     MaxSA = max(Surface.Albedo), MinSA = min(Surface.Albedo), AvgSA = mean(Surface.Albedo), 
                     MaxPress = max(Pressure), MinPress = min(Pressure), AvgPres = mean(Pressure), 
                     MaxPW = max(Precipitable.Water), MinPW = min(Precipitable.Water), AvgPW = mean(Precipitable.Water), 
                     MaxWD = max(Wind.Direction), MinWD = min(Wind.Direction), AvgWD = mean(Wind.Direction), 
                     MaxWS = max(Wind.Speed), MinWS = min(Wind.Speed), AvgWS = mean(Wind.Speed), 
                     MaxUVO = max(Global.Horizontal.UV.Irradiance..280.400nm.), AvgUVO = mean(Global.Horizontal.UV.Irradiance..280.400nm.), 
                     MaxUVI = max(Global.Horizontal.UV.Irradiance..295.385nm.), AvgUVI = mean(Global.Horizontal.UV.Irradiance..295.385nm.)
    ) %>%
    mutate(location = specify_location) %>%
    relocate(location)
  return(df1)
}

calc_nrel_eu <- function(file_path, month_val, day_val, specify_location) {
  df1 <- read.table(file_path, header = T, skip = 2, fill = T, sep = ",") %>% filter(Month == month_val, Day == day_val) %>%
    dplyr::summarize(MaxTemp = max(Temperature), MinTemp = min(Temperature), AvgTemp = mean(Temperature), 
                      MaxcDHI = max(Clearsky.DHI), AvgcDHI = mean(Clearsky.DHI), 
                      MaxcDNI = max(Clearsky.DNI), AvgcDNI = mean(Clearsky.DNI), 
                      MaxcGHI = max(Clearsky.GHI), AvgcGHI = mean(Clearsky.GHI),  
                      MaxDP = max(Dew.Point), MinDP = min(Dew.Point), AvgDP = mean(Dew.Point), 
                      MaxDHI = max(DHI), AvgDHI = mean(DHI), 
                      MaxcDNI = max(DNI), AvgcDNI = mean(DNI), 
                      MaxFF = max(Fill.Flag), AvgFF = mean(Fill.Flag), 
                      MaxcGHI = max(GHI), AvgcGHI = mean(GHI),  
                      MaxRH = max(Relative.Humidity), MinRH = min(Relative.Humidity), AvgRH = mean(Relative.Humidity), 
                      MaxSZA = max(Solar.Zenith.Angle), MinSZA = min(Solar.Zenith.Angle), AvgSZA = mean(Solar.Zenith.Angle), 
                      MaxSA = max(Surface.Albedo), MinSA = min(Surface.Albedo), AvgSA = mean(Surface.Albedo), 
                      MaxPress = max(Pressure), MinPress = min(Pressure), AvgPres = mean(Pressure), 
                      MaxPW = max(Precipitable.Water), MinPW = min(Precipitable.Water), AvgPW = mean(Precipitable.Water), 
                      MaxWD = max(Wind.Direction), MinWD = min(Wind.Direction), AvgWD = mean(Wind.Direction), 
                      MaxWS = max(Wind.Speed), MinWS = min(Wind.Speed), AvgWS = mean(Wind.Speed)
#                      MaxUVO = max(Global.Horizontal.UV.Irradiance..280.400nm.), AvgUVO = mean(Global.Horizontal.UV.Irradiance..280.400nm.), 
#                      MaxUVI = max(Global.Horizontal.UV.Irradiance..295.385nm.), AvgUVI = mean(Global.Horizontal.UV.Irradiance..295.385nm.)
     ) %>%
     mutate(MaxUVO=0, AvgUVO=0, MaxUVI=0, AvgUVI=0, location = specify_location) %>%
     relocate(location)
  return(df1)
}

calc_nrel_daylight <- function(file_path, month_val, day_val, specify_location,uncorrected_daylight,atmospheric_corrected_daylight) {
  df1 <- read.table(file_path, header = T, skip = 2, fill = T, sep = ",") %>%
    filter(Month == month_val, Day == day_val) %>%
    dplyr::summarize(MaxTemp = max(Temperature), MinTemp = min(Temperature), AvgTemp = mean(Temperature), 
                     MaxcDHI = max(Clearsky.DHI), AvgcDHI = mean(Clearsky.DHI), 
                     MaxcDNI = max(Clearsky.DNI), AvgcDNI = mean(Clearsky.DNI), 
                     MaxcGHI = max(Clearsky.GHI), AvgcGHI = mean(Clearsky.GHI),  
                     MaxDP = max(Dew.Point), MinDP = min(Dew.Point), AvgDP = mean(Dew.Point), 
                     MaxDHI = max(DHI), AvgDHI = mean(DHI), 
                     MaxcDNI = max(DNI), AvgcDNI = mean(DNI), 
                     MaxFF = max(Fill.Flag), AvgFF = mean(Fill.Flag), 
                     MaxcGHI = max(GHI), AvgcGHI = mean(GHI),  
                     MaxRH = max(Relative.Humidity), MinRH = min(Relative.Humidity), AvgRH = mean(Relative.Humidity), 
                     MaxSZA = max(Solar.Zenith.Angle), MinSZA = min(Solar.Zenith.Angle), AvgSZA = mean(Solar.Zenith.Angle), 
                     MaxSA = max(Surface.Albedo), MinSA = min(Surface.Albedo), AvgSA = mean(Surface.Albedo), 
                     MaxPress = max(Pressure), MinPress = min(Pressure), AvgPres = mean(Pressure), 
                     MaxPW = max(Precipitable.Water), MinPW = min(Precipitable.Water), AvgPW = mean(Precipitable.Water), 
                     MaxWD = max(Wind.Direction), MinWD = min(Wind.Direction), AvgWD = mean(Wind.Direction), 
                     MaxWS = max(Wind.Speed), MinWS = min(Wind.Speed), AvgWS = mean(Wind.Speed), 
                     MaxUVO = max(Global.Horizontal.UV.Irradiance..280.400nm.), AvgUVO = mean(Global.Horizontal.UV.Irradiance..280.400nm.), 
                     MaxUVI = max(Global.Horizontal.UV.Irradiance..295.385nm.), AvgUVI = mean(Global.Horizontal.UV.Irradiance..295.385nm.)
    ) %>%
    mutate(location = specify_location, uncorrected_daylight = uncorrected_daylight, atmospheric_corrected_daylight = atmospheric_corrected_daylight) %>%
    relocate(location)
  return(df1)
}

mantel_EC21 <- function(mapper_name, bc_name, metadata_name) {
  test1 <- read.table(mapper_name, sep = "\t", comment.char = "", header = T) %>%
    mutate(Location = as.character(recode_factor(Geography, D = "Durham, NH", E = "Etna, ME", V = "Charlottesville, VA", MD = "Churchville, MD", PA = "Media, PA", L = "Middlefield, CT", H = "Harvard, MA", B = "Bowdoin, ME"))) %>% 
    mutate(latitude = as.character(recode_factor(Geography, D = "43.08745", E = "44.82082", V = "38.07474", MD = "39.539", PA = "39.90924", L = "41.49383", H = "42.50867", B = "44.0253"))) %>% 
    mutate(longitude = as.character(recode_factor(Geography, D = "-71.04666", E = "-69.11144", V = "-78.38155", MD = "-76.2399", PA = "-75.40529", L = "-72.71145", H = "-71.57071", B = "-69.94379"))) %>% 
    mutate(date = as.character(recode_factor(Geography, D = "44469.00", E = "44468.00", V = "44489.00", MD = "44490.00", PA = "44497.00", L = "44470.00", H = "44469.00", B = "44468.00"))) %>% 
    inner_join(read.table(paste0('~/Dropbox/east_coast_flies_manuscript/mantel/',metadata_name), sep = "\t", header=T, comment.char = ""), by = c("Location" = "location")) 
  
  write.table(test1, paste0("mantel_",mapper_name), row.names = F, quote = F, sep = "\t")
  
  mapper_file = paste0("mantel_",mapper_name)
  
  mapper_cols <- c("latitude","longitude","date",colnames(test1)[(dim(test1)[2]-40):dim(test1)[2]])
  print(mapper_cols)
  calc_mantel_microbiome_fromCFUs(bc_name,mapper_file,mapper_cols)
  # [1] "rho: -0.014; p: 0.676"
  # [1] "rho: 0.008; p: 0.321"
  
  mapper_cols <- c("latitude","longitude","date")
  print(mapper_cols)
  calc_mantel_microbiome_fromCFUs(bc_name,mapper_file,mapper_cols)
  # [1] "rho: 0.032; p: 0.057"
  # [1] "rho: -0.021; p: 0.949"
  
  mapper_cols <- c("latitude","longitude","date","MaxTemp","MinTemp","AvgTemp","MaxRH","MinRH","AvgRH","MaxUVO","AvgUVO","MaxUVI","AvgUVI")
  print(mapper_cols)
  calc_mantel_microbiome_fromCFUs(bc_name,mapper_file,mapper_cols)
  # [1] "rho: 0.052; p: 0.024"
  # [1] "rho: -0.034; p: 0.971"
  
  mapper_cols <- c("latitude")
  print(mapper_cols)
  calc_mantel_microbiome_fromCFUs(bc_name,mapper_file,mapper_cols)
  # [1] "rho: 0.052; p: 0.024"
  # [1] "rho: -0.034; p: 0.971"
}

mantel_EC21_M <- function(mapper_name, bc_name, metadata_name) {
  test1 <- read.table(mapper_name, sep = "\t", comment.char = "", header = T) %>%
    mutate(Location = as.character(recode_factor(Geography, D = "Durham, NH", E = "Etna, ME", V = "Charlottesville, VA", M = "Churchville, MD", PA = "Media, PA", L = "Middlefield, CT", H = "Harvard, MA", B = "Bowdoin, ME"))) %>% 
    mutate(latitude = as.character(recode_factor(Geography, D = "43.08745", E = "44.82082", V = "38.07474", M = "39.539", PA = "39.90924", L = "41.49383", H = "42.50867", B = "44.0253"))) %>% 
    mutate(longitude = as.character(recode_factor(Geography, D = "-71.04666", E = "-69.11144", V = "-78.38155", M = "-76.2399", PA = "-75.40529", L = "-72.71145", H = "-71.57071", B = "-69.94379"))) %>% 
    mutate(date = as.character(recode_factor(Geography, D = "44469.00", E = "44468.00", V = "44489.00", M = "44490.00", PA = "44497.00", L = "44470.00", H = "44469.00", B = "44468.00"))) %>% 
    inner_join(read.table(paste0('~/Dropbox/east_coast_flies_manuscript/mantel/',metadata_name), sep = "\t", header=T, comment.char = ""), by = c("Location" = "location")) 
  
  write.table(test1, paste0("mantel_",mapper_name), row.names = F, quote = F, sep = "\t")
  
  mapper_file = paste0("mantel_",mapper_name)
  
  mapper_cols <- c("latitude","longitude","date",colnames(test1)[(dim(test1)[2]-40):dim(test1)[2]])
  print(mapper_cols)
  calc_mantel_microbiome_fromCFUs(bc_name,mapper_file,mapper_cols)
  # [1] "rho: -0.014; p: 0.676"
  # [1] "rho: 0.008; p: 0.321"
  
  mapper_cols <- c("latitude","longitude","date")
  print(mapper_cols)
  calc_mantel_microbiome_fromCFUs(bc_name,mapper_file,mapper_cols)
  # [1] "rho: 0.032; p: 0.057"
  # [1] "rho: -0.021; p: 0.949"
  
  mapper_cols <- c("latitude","longitude","date","MaxTemp","MinTemp","AvgTemp","MaxRH","MinRH","AvgRH","MaxUVO","AvgUVO","MaxUVI","AvgUVI")
  print(mapper_cols)
  calc_mantel_microbiome_fromCFUs(bc_name,mapper_file,mapper_cols)
  # [1] "rho: 0.052; p: 0.024"
  # [1] "rho: -0.034; p: 0.971"
  
  mapper_cols <- c("latitude")
  print(mapper_cols)
  calc_mantel_microbiome_fromCFUs(bc_name,mapper_file,mapper_cols)
  # [1] "rho: 0.052; p: 0.024"
  # [1] "rho: -0.034; p: 0.971"
}

mantel_EC18 <- function(mapper_name, file_path, metadata_name, cols_to_test = NULL, metric_to_test = "bray_curtis") {
  test1 <- read.table(mapper_name, sep = "\t", comment.char = "", header = T) %>%
    inner_join(read.table(paste0(metadata_name), sep = "\t", header=T, comment.char = ""), by = c("location")) 
  
  write.table(test1, paste0("mantel_",mapper_name), row.names = F, quote = F, sep = "\t")
  
  mapper_file = paste0("mantel_",mapper_name)
  
  if (is.null(cols_to_test)) {
    mapper_cols <- c("latitude","longitude","date",colnames(test1)[(dim(test1)[2]-40):dim(test1)[2]])
    print(mapper_cols)
    for (i in c("bray_curtis","unweighted_unifrac","weighted_unifrac")) {
      cat(i)
      calc_mantel_microbiome(file_path,i,mapper_file,mapper_cols)
    }
    
    mapper_cols <- c("latitude","longitude","date")
    print(mapper_cols)
    for (i in c("bray_curtis","unweighted_unifrac","weighted_unifrac")) {
      cat(i)
      calc_mantel_microbiome(file_path,i,mapper_file,mapper_cols)
    }
    
    mapper_cols <- c("latitude","longitude","date","MaxTemp","MinTemp","AvgTemp","MaxRH","MinRH","AvgRH","MaxUVO","AvgUVO","MaxUVI","AvgUVI")
    print(mapper_cols)
    for (i in c("bray_curtis","unweighted_unifrac","weighted_unifrac")) {
      cat(i)
      calc_mantel_microbiome(file_path,i,mapper_file,mapper_cols)
    }
    
    mapper_cols <- c("latitude")
    print(mapper_cols)
    for (i in c("bray_curtis","unweighted_unifrac","weighted_unifrac")) {
      cat(i)
      calc_mantel_microbiome(file_path,i,mapper_file,mapper_cols)
    }
  } else {
    for (j in 1:length(cols_to_test)) {
      mapper_cols <- cols_to_test[[j]]
      print(mapper_cols)
      cat(metric_to_test)
      return(calc_mantel_microbiome(file_path,metric_to_test,mapper_file,mapper_cols))
    }
  }
}

mantel_UT20 <- function(mapper_name, file_path, metadata_name, cols_to_test = NULL, metric_to_test = "bray_curtis") {
  
  test1 <- read.table(mapper_name, sep = "\t", comment.char = "", header = T) %>%
    mutate(latitude = as.character(recode_factor(location, `Alpine, UT` = "40.4457589", `Santaquin, UT` = "39.966865", `Lindon1, UT` = "40.3327707", `Logan, UT` = "41.72573", `Ogden, UT` = "41.326252", `Lindon2, UT` = "40.3327707", `Lindon3, UT` = "40.3327707"))) %>% 
    mutate(longitude = as.character(recode_factor(location, `Alpine, UT` = "-111.7818936", `Santaquin, UT` = "-111.791562", `Lindon1, UT` = "-111.7219739", `Logan, UT` = "-111.809883", `Ogden, UT` = "-112.011375", `Lindon2, UT` = "-111.7219739", `Lindon3, UT` = "-111.7219739"))) %>%
    mutate(date = as.character(recode_factor(location, `Alpine, UT` = "44113", `Santaquin, UT` = "44114", `Lindon1, UT` = "44103", `Logan, UT` = "44112", `Ogden, UT` = "44113", `Lindon2, UT` = "44110", `Lindon3, UT` = "44114"))) %>%
    inner_join(read.table(paste0(metadata_name), sep = "\t", header=T, comment.char = ""), by = c("location")) 
  
  write.table(test1, paste0("mantel_",mapper_name), row.names = F, quote = F, sep = "\t")
  
  mapper_file = paste0("mantel_",mapper_name)
  
  if (is.null(cols_to_test)) {
    mapper_cols <- c("latitude","longitude","date",colnames(test1)[(dim(test1)[2]-40):dim(test1)[2]])
    print(mapper_cols)
    for (i in c("bray_curtis","unweighted_unifrac","weighted_unifrac")) {
      cat(i)
      calc_mantel_microbiome(file_path,i,mapper_file,mapper_cols)
    }
    
    mapper_cols <- c("latitude","longitude","date")
    print(mapper_cols)
    for (i in c("bray_curtis","unweighted_unifrac","weighted_unifrac")) {
      cat(i)
      calc_mantel_microbiome(file_path,i,mapper_file,mapper_cols)
    }
    
    mapper_cols <- c("latitude","longitude","date","MaxTemp","AvgRH","AvgUVO","AvgUVI")
    print(mapper_cols)
    for (i in c("bray_curtis","unweighted_unifrac","weighted_unifrac")) {
      cat(i)
      calc_mantel_microbiome(file_path,i,mapper_file,mapper_cols)
    }
    
    mapper_cols <- c("latitude")
    print(mapper_cols)
    for (i in c("bray_curtis","unweighted_unifrac","weighted_unifrac")) {
      cat(i)
      calc_mantel_microbiome(file_path,i,mapper_file,mapper_cols)
    }

  } else {
    
    for (j in 1:length(cols_to_test)) {
      mapper_cols <- cols_to_test[[j]]
      print(mapper_cols)
        cat(metric_to_test)
        return(calc_mantel_microbiome(file_path,metric_to_test,mapper_file,mapper_cols))
    }
  }
}

mantel_EC21seq <- function(mapper_name, file_path, metadata_name, cols_to_test = NULL, metric_to_test = "bray_curtis") {
  
  test1 <- read.table(mapper_name, sep = "\t", comment.char = "", header = T) %>%
    mutate(location = site) %>% 
    mutate(latitude = as.character(recode_factor(site, `Durham, NH` = "43.08745", `Etna, ME` = "44.82082", `Charlottesville, VA` = "38.07474", `Churchville, MD` = "39.539", `Media, PA` = "39.90924", `Middlefield, CT` = "41.49383", `Harvard, MA` = "42.50867", `Bowdoin, ME` = "44.0253"))) %>% 
    mutate(longitude = as.character(recode_factor(site, `Durham, NH` = "-71.04666", `Etna, ME` = "-69.11144", `Charlottesville, VA` = "-78.38155", `Churchville, MD` = "-76.2399", `Media, PA` = "-75.40529", `Middlefield, CT` = "-72.71145", `Harvard, MA` = "-71.57071", `Bowdoin, ME` = "-69.94379"))) %>% 
    mutate(longitude = as.character(recode_factor(site, `Durham, NH` = "44469.00", `Etna, ME` = "44468.00", `Charlottesville, VA` = "44489.00", `Churchville, MD` = "44490.00", `Media, PA` = "44497.00", `Middlefield, CT` = "44470.00", `Harvard, MA` = "44469.00", `Bowdoin, ME` = "44468.00"))) %>% 
    inner_join(read.table(paste0(metadata_name), sep = "\t", header=T, comment.char = ""), by = c("location")) 
  
  write.table(test1, paste0("mantel_",mapper_name), row.names = F, quote = F, sep = "\t")
  
  mapper_file = paste0("mantel_",mapper_name)
  if (is.null(cols_to_test)) {
    
  mapper_cols <- c("latitude","longitude","date",colnames(test1)[(dim(test1)[2]-40):dim(test1)[2]])
  print(mapper_cols)
  for (i in c("bray_curtis","unweighted_unifrac","weighted_unifrac")) {
    cat(i)
    calc_mantel_microbiome(file_path,i,mapper_file,mapper_cols)
  }
  
  mapper_cols <- c("latitude","longitude","date")
  print(mapper_cols)
  for (i in c("bray_curtis","unweighted_unifrac","weighted_unifrac")) {
    cat(i)
    calc_mantel_microbiome(file_path,i,mapper_file,mapper_cols)
  }
  
  mapper_cols <- c("latitude","longitude","date","MaxTemp","MinTemp","AvgTemp","MaxRH","MinRH","AvgRH","MaxUVO","AvgUVO","MaxUVI","AvgUVI")
  print(mapper_cols)
  for (i in c("bray_curtis","unweighted_unifrac","weighted_unifrac")) {
    cat(i)
    calc_mantel_microbiome(file_path,i,mapper_file,mapper_cols)
  }
  
  mapper_cols <- c("latitude")
  print(mapper_cols)
  for (i in c("bray_curtis","unweighted_unifrac","weighted_unifrac")) {
    cat(i)
    calc_mantel_microbiome(file_path,i,mapper_file,mapper_cols)
  }
    
  } else {
    for (j in 1:length(cols_to_test)) {
      mapper_cols <- cols_to_test[[j]]
      print(mapper_cols)
      cat(metric_to_test)
      return(calc_mantel_microbiome(file_path,metric_to_test,mapper_file,mapper_cols))
    }
  }
}

preliminary_2017_mantel <- function(file_path, meta2017=meta2017) {
  ## make the fst matrix
  fst <- read.table("orch17.fst_jmc.csv", sep = ",", header=T) %>%
    inner_join(meta2017, by = c("X" = "sampID")) %>%
    dplyr::select(-colnames(meta2017)[-c(2,7)]) %>%
    t() %>% 
    data.frame() 
  colnames(fst) <- fst[dim(fst)[1],]
  
  ## make the microbiota matrix
  dm <- read.table(paste('core-metrics-results-',file_path,'/',metric,'_distance_matrix/distance-matrix.tsv',sep=""), header=T, sep="\t") %>% 
    mutate(X=as.character(X)) %>% 
    inner_join(read.table("metadata_2017b.txt", comment.char = "", header=T, fill=T, sep = "\t"), by=c("X"="X.SampleID"))
  colnames(dm) = gsub(".","-",colnames(dm), fixed=T)
  
  names_in_both <- dm$X[dm$X %in% colnames(fst)]
  
  cat("Sum names not in the microbiota matrix: ")
  sum(!names_in_both %in% dm$X)
  cat("Sum names not in the fst matrix: ")
  sum(!names_in_both %in% colnames(fst))
  
  ## filter the matrices
  dm2 <- dm %>%
    filter(X %in% names_in_both) %>%
    arrange(match(X, names_in_both))%>%
    dplyr::select(names_in_both)
  dm3 <- as.dist(dm2)
  
  
  fst$sampID = rownames(fst)
  fst$sampID =gsub(".", "-",fst$sampID, fixed=T)
    
  fst2 <- fst %>%
    relocate(sampID) %>%
    inner_join(meta2017, by = c("sampID")) %>%
    dplyr::select(-colnames(meta2017)[-c(2,7)])
  
  rownames(fst2) <- fst2$X.SampleID   
  
  sum(table(fst2$X.SampleID)>1)
  
  fst3 <- fst2 %>%
    filter(X.SampleID %in% names_in_both) %>%
    arrange(match(X.SampleID, names_in_both)) %>%
    dplyr::select(names_in_both) 
  
  fst4 <- as.dist(fst3)
  
  cat("dimensions fst: ")
  print(dim(fst4))
  cat("dimenstions microbiota: ") 
  print(dim(dm3))
  cat("# of times names don't match for microbiota and fst")
  sum(names(dm3) == names(fst4))
  
  set.seed(43)
  print(mantel(fst4,dm3, method = "spearman", parallel = 7))
  plot(fst4,dm3)
}  

preliminary_2016_mantel <- function(file_path, metric, spring_or_fall) {
  dm <- read.table(paste('core-metrics-results-',file_path,'/',metric,'_distance_matrix/distance-matrix.tsv',sep=""), header=T, sep="\t") %>% mutate(X=as.character(X)) 
  colnames(dm)[2:(dim(dm)[2])] <- dm$X
  dmdist <- as.dist(dm[,2:(dim(dm)[2])])
  
  ## make the metadata matrix
  fst <- read.table(spring_or_fall, sep = ",", header=T) 
  colnames(fst)[2:(dim(fst)[2])] <- fst$X
  fst2 <- fst %>% 
    filter(X %in% (colnames(dm %>% dplyr::select(-X)))) %>%
    arrange(match(X, (colnames(dm %>% dplyr::select(-X))))) %>%
    dplyr::select((colnames(dm %>% dplyr::select(-X)))) 
  fst_dm <- fst2 %>% as.matrix()
  fst_dm2 <- as.dist(fst_dm)
  
  cat("dimensions: ")
  print(dim(fst_dm2))
  print(dim(dmdist))
  sum(names(fst_dm2) == names(dmdist))
  
  set.seed(43)
  mantel1 <- mantel(fst_dm2,dmdist, method = "spearman", parallel = 7)
  cat(metric)
  print(paste0("rho: ",round(mantel1$statistic,3),"; p: ",mantel1$signif))
  
  plot(fst_dm2,dmdist)
}


preliminary_2017_mantel_noreplication <- function(file_path, metric, spring_or_fall) {
  ## make the fst matrix
  fst <- read.table(spring_or_fall, sep = ",", header=T)
  fst_names <- fst$X
  
  dm <- read.table(paste('core-metrics-results-',file_path,'/',metric,'_distance_matrix/distance-matrix.tsv',sep=""), header=T, sep="\t") %>% mutate(X=as.character(X)) 
  dm_names <- dm$X
  
  names_in_both <- fst_names[fst_names %in% dm_names]
  names_to_drop <- dm_names[!dm_names %in% fst_names]
  
  if(length(names_to_drop)>0) {
    cat("go back and drop these names from your microbiota matrix: ")
    print(names_to_drop)
  } else {
    preliminary_2016_mantel(file_path,metric,spring_or_fall)
  }
}

# mantel_output <- henryorchard3
# point_size = 0.5
# y_title = "hi"
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

mantel_distance_plot_2017 <- function(mantel_output, x_title = "microbiota distance (Bray-Curtis)", y_title = y_title, point_size = 1) {
  
aab <- melt(as.matrix(mantel_output[[1]]))[melt(upper.tri(as.matrix(mantel_output[[1]])))$value,] %>% mutate(Var12 = paste0(Var1,Var2))
aac <- melt(as.matrix(mantel_output[[2]]))[melt(upper.tri(as.matrix(mantel_output[[2]])))$value,] %>% mutate(Var1 = gsub("X","",Var1), Var2 = gsub("X","",Var2),Var1 = gsub("\\.","-",Var1), Var2 = gsub("\\.","-",Var2)) %>% mutate(Var12 = paste0(Var1,Var2))

aad <- aab %>% inner_join(aac, by = "Var12")

utmantel <- ggplot(aad, aes(x = value.x, y = value.y)) + 
  geom_jitter(size=point_size, width = max(aad$value.x)/75, alpha=0.5) + 
  stat_smooth(method = "lm", se=F) +
  theme_cowplot() + theme_pubclean() + labs(x = x_title, y = y_title )

return(utmantel)
}



get_core_features <- function(file_path = file_path, map_path = map_path, taxonomic_level = taxonomic_level, similarity_cutoff = 0.5, var5 = NULL) {
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
  
  phyl3 <- phyl2[,-1]
  
  hist(rowSums(phyl3>0)/dim(phyl3)[2], breaks=10)
  # max(rowSums(phyl3>0)/dim(phyl3)[2])
  # min(rowSums(phyl3>0)/dim(phyl3)[2])
  
  return(phyl2$clustered_taxonomy[rowSums(phyl3>0)/dim(phyl3)[2]>similarity_cutoff])
}

get_single_core_features <- function(file_path = file_path, map_path = map_path, taxonomic_level = taxonomic_level, var5 = NULL) {
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
  
  return(phyl2)
}

tabulate_individual_venns <- function(file_path = file_path, taxonomic_level = taxonomic_level, otu_table = venns_otu_table, threshhold = 1) {
  starved_fruit_venn <- get_single_core_features(file_path = file_path, map_path = map_path, taxonomic_level = taxonomic_level)
  starved_fruit_venn$clustered_taxonomy = gsub("[", replacement = "", starved_fruit_venn$clustered_taxonomy, perl = T, fixed=T)
  starved_fruit_venn$clustered_taxonomy = gsub("]", replacement = "", starved_fruit_venn$clustered_taxonomy, perl = T, fixed=T)
  
  orchard_site_names <- otu_table %>% filter(X.SampleID %in% colnames(starved_fruit_venn)) %>% dplyr::pull(orchard_subsite) %>% table() %>% data.frame()
  
  group1name = strsplit(file_path, split = "_", fixed = T)[[1]][2]
  group2name = strsplit(file_path, split = "_", fixed = T)[[1]][3]
  
  out_df <- data.frame(orchard_subsite = character(), group1 = character(), group2 = character(), unique1 = numeric, shared = numeric(), unique2 = numeric)
  
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

tabulate_overall_venns <- function(file_path = file_path, taxonomic_level = taxonomic_level, otu_table = venns_otu_table, threshhold = 1) {
  starved_fruit_venn <- get_single_core_features(file_path = file_path, map_path = map_path, taxonomic_level = taxonomic_level)
  starved_fruit_venn$clustered_taxonomy = gsub("[", replacement = "", starved_fruit_venn$clustered_taxonomy, perl = T, fixed=T)
  starved_fruit_venn$clustered_taxonomy = gsub("]", replacement = "", starved_fruit_venn$clustered_taxonomy, perl = T, fixed=T)
  
  orchard_site_names <- otu_table %>% filter(X.SampleID %in% colnames(starved_fruit_venn)) %>% dplyr::pull(orchard_subsite) %>% table() %>% data.frame()
  
  group1name = strsplit(file_path, split = "_", fixed = T)[[1]][2]
  group2name = strsplit(file_path, split = "_", fixed = T)[[1]][3]
  
  out_df <- data.frame(orchard_subsite = character(), group1 = character(), group2 = character(), unique1 = numeric, shared = numeric(), unique2 = numeric)
  
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
