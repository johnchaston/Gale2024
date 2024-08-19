make_permanova_tables <- function(file_path = file_path, mapper_file = mapper_file, flies_unifrac_table = flies_unifrac_table, pform = pform, beta_div_tests = c("unweighted_unifrac","weighted_unifrac","bray_curtis"), seed_to_set = 42) {

write.table("Permanova tables", file = paste0('core-metrics-results-',file_path,"/",file_path,"_permanova_table.txt"), append = F, quote = F, sep = "\t", row.names = F, col.names = F)

i="unweighted_unifrac"
  for(i in beta_div_tests) {
  assign(x = paste0("flies_",i), value = read.table(paste('core-metrics-results-',file_path,'/',i,'_distance_matrix/distance-matrix.tsv',sep=""), header=T, sep="\t") %>% mutate(X=as.character(X)))
  assign(x = paste0("flies_",i,"_dm"), value = as.dist(get(paste0("flies_",i))[,2:dim(get(paste0("flies_",i)))[2]]))
  flies_unifrac_table <- get(paste0("flies_",i)) %>% dplyr::select(X) %>% left_join(read.table(mapper_file,comment.char = "", header=T, fill=T, sep="\t"), by=c("X"="X.SampleID")) %>% droplevels()
  set.seed(seed_to_set)
  assign(x = paste0("f_",i,"_permanova"), adonis2(as.formula(paste0("flies_",i,"_dm ",pform)), flies_unifrac_table, permutations=1000))
  assign(x = paste0("fig_",i), value = tableGrob(data.frame(Df = round(get(paste0("f_",i,"_permanova"))$Df, 2), SS = round(get(paste0("f_",i,"_permanova"))$SumOfSqs, 2), R2 = round(get(paste0("f_",i,"_permanova"))$R2, 2), Fval = round(get(paste0("f_",i,"_permanova"))$`F`, 2), p = round(get(paste0("f_",i,"_permanova"))$`Pr(>F)`, 2)),theme=ttheme_minimal()))
  plot(get(paste0("fig_",i)))
  write.table(i, file = paste0('core-metrics-results-',file_path,"/",file_path,"_permanova_table.txt"), append = T, quote = F, sep = "\t", row.names = F, col.names = F)
  write.table(round(get(paste0("f_",i,"_permanova")),2), file = paste0('core-metrics-results-',file_path,"/",file_path,"_permanova_table.txt"), append = T, quote = F, sep = "\t", row.names = T, col.names = T)
  
}
  
}