#make_permanova_tables <- function(file_path = file_path, mapper_file = mapper_file, flies_unifrac_table = flies_unifrac_table, pform = pform, beta_div_tests = c("unweighted_unifrac","weighted_unifrac","bray_curtis"), seed_to_set = 42, byval = "terms") {
#write.table("Permanova tables", file = paste0('core-metrics-results-',file_path,"/",file_path,"_permanova_table.txt"), append = F, quote = F, sep = "\t", row.names = F, col.names = F)
#i="unweighted_unifrac"
#  for(i in beta_div_tests) {
#  assign(x = paste0("flies_",i), value = read.table(paste('core-metrics-results-',file_path,'/',i,'_distance_matrix/distance-matrix.tsv',sep=""), header=T, sep="\t") %>% mutate(X=as.character(X)))
#  assign(x = paste0("flies_",i,"_dm"), value = as.dist(get(paste0("flies_",i))[,2:dim(get(paste0("flies_",i)))[2]]))
#  flies_unifrac_table <- get(paste0("flies_",i)) %>% dplyr::select(X) %>% left_join(read.table(mapper_file,comment.char = "", header=T, fill=T, sep="\t", quote = '"'), by=c("X"="X.SampleID")) %>% droplevels()
#  set.seed(seed_to_set)
#    if(byval == "margin") {
#  assign(x = paste0("f_",i,"_permanova"), adonis2(as.formula(paste0("flies_",i,"_dm ",pform)), flies_unifrac_table, permutations=1000, by = "margin"))
#      } else { 
#  assign(x = paste0("f_",i,"_permanova"), adonis2(as.formula(paste0("flies_",i,"_dm ",pform)), flies_unifrac_table, permutations=1000, by = "terms"))
#      }
#  assign(x = paste0("fig_",i), value = tableGrob(data.frame(Df = round(get(paste0("f_",i,"_permanova"))$Df, 2), SS = round(get(paste0("f_",i,"_permanova"))$SumOfSqs, 2), R2 = round(get(paste0("f_",i,"_permanova"))$R2, 2), Fval = round(get(paste0("f_",i,"_permanova"))$`F`, 2), p = round(get(paste0("f_",i,"_permanova"))$`Pr(>F)`, 2)),theme=ttheme_minimal()))
#  plot(get(paste0("fig_",i)))
#  write.table(i, file = paste0('core-metrics-results-',file_path,"/",file_path,"_permanova_table.txt"), append = T, quote = F, sep = "\t", row.names = F, col.names = F)
#  write.table(round(get(paste0("f_",i,"_permanova")),2), file = paste0('core-metrics-results-',file_path,"/",file_path,"_permanova_table.txt"), append = T, quote = F, sep = "\t", row.names = T, col.names = T)
#}
#}


make_permanova_tables <- function(file_path = file_path, mapper_file = mapper_file, flies_unifrac_table = flies_unifrac_table, pform = pform, beta_div_tests = c("unweighted_unifrac","weighted_unifrac","bray_curtis"), seed_to_set = 42, byval = "terms", stratavar=NULL, converttable =NULL) {

  write.table("Permanova tables", file = paste0('core-metrics-results-',file_path,"/",file_path,"_permanova_table.txt"), append = F, quote = F, sep = "\t", row.names = F, col.names = F)

   
  i="unweighted_unifrac"
  for(i in beta_div_tests) {
  assign(x = paste0("flies_",i), value = read.table(paste('core-metrics-results-',file_path,'/',i,'_distance_matrix/distance-matrix.tsv',sep=""), header=T, sep="\t") %>% mutate(X=as.character(X)))
  assign(x = paste0("flies_",i,"_dm"), value = as.dist(get(paste0("flies_",i))[,2:dim(get(paste0("flies_",i)))[2]]))
  flies_unifrac_table <- get(paste0("flies_",i)) %>% dplyr::select(X) %>% left_join(read.table(mapper_file,comment.char = "", header=T, fill=T, sep="\t", quote = '"'), by=c("X"="X.SampleID")) %>% droplevels() %>% mutate(time_fac = factor(time_num, ordered = T))
  if (byval == "margin" & is.null(stratavar)) {
    set.seed(seed_to_set)
    assign(x = paste0("f_",i,"_permanova"), adonis2(as.formula(paste0("flies_",i,"_dm ",pform)), flies_unifrac_table, permutations=1000, by = "margin"))
    print("margin")
  } else if (byval == "terms" & is.null(stratavar)) {
    set.seed(seed_to_set)
    assign(x = paste0("f_",i,"_permanova"), adonis2(as.formula(paste0("flies_",i,"_dm ",pform)), flies_unifrac_table, permutations=1000, by = "terms"))
    print("terms")
  } else if (byval == "margin" & !is.null(stratavar)) {
    set.seed(seed_to_set)
    assign(x = paste0("f_",i,"_permanova"), adonis2(as.formula(paste0("flies_",i,"_dm ",pform)), data = flies_unifrac_table, permutations=1000, by = "margin", strata = flies_unifrac_table[,paste0(stratavar)]))
    print("marginstrata")
  } else if (byval == "terms" & !is.null(stratavar)) {
    set.seed(seed_to_set)
    assign(x = paste0("f_",i,"_permanova"), adonis2(as.formula(paste0("flies_",i,"_dm ",pform)), data = flies_unifrac_table, permutations=1000, by = "terms", strata = flies_unifrac_table[,paste0(stratavar)]))
    print("termsstrata")
  }
  
  
  assign(x = paste0("fig_",i), value = tableGrob(data.frame(Df = round(get(paste0("f_",i,"_permanova"))$Df, 2), SS = round(get(paste0("f_",i,"_permanova"))$SumOfSqs, 2), R2 = round(get(paste0("f_",i,"_permanova"))$R2, 2), Fval = round(get(paste0("f_",i,"_permanova"))$`F`, 2), p = round(get(paste0("f_",i,"_permanova"))$`Pr(>F)`, 4)),theme=ttheme_minimal()))
  plot(get(paste0("fig_",i)))
  write.table(i, file = paste0('core-metrics-results-',file_path,"/",file_path,"_permanova_table.txt"), append = T, quote = F, sep = "\t", row.names = F, col.names = F)
  write.table(round(get(paste0("f_",i,"_permanova")),2), file = paste0('core-metrics-results-',file_path,"/",file_path,"_permanova_table.txt"), append = T, quote = F, sep = "\t", row.names = T, col.names = T)
  }
  if(!is.null(converttable)) {
    convert_permanova_to_word(file_path = file_path, table_row_names = converttable)
  }
}
  
make_permanova_tables2fac <- function(file_path = file_path, mapper_file = mapper_file, flies_unifrac_table = flies_unifrac_table, pform = pform, beta_div_tests = c("unweighted_unifrac","weighted_unifrac","bray_curtis"), seed_to_set = 42, byval = "terms", stratavar=NULL, converttable =NULL) {

  write.table("Permanova tables", file = paste0('core-metrics-results-',file_path,"/",file_path,"_permanova_table.txt"), append = F, quote = F, sep = "\t", row.names = F, col.names = F)

   
  i="unweighted_unifrac"
  for(i in beta_div_tests) {
  assign(x = paste0("flies_",i), value = read.table(paste('core-metrics-results-',file_path,'/',i,'_distance_matrix/distance-matrix.tsv',sep=""), header=T, sep="\t") %>% mutate(X=as.character(X)))
  assign(x = paste0("flies_",i,"_dm"), value = as.dist(get(paste0("flies_",i))[,2:dim(get(paste0("flies_",i)))[2]]))
  flies_unifrac_table <- get(paste0("flies_",i)) %>% dplyr::select(X) %>% left_join(read.table(mapper_file,comment.char = "", header=T, fill=T, sep="\t", quote = '"'), by=c("X"="X.SampleID")) %>% droplevels() %>% mutate(time_fac = factor(time_num, ordered = T))
  if (byval == "margin" & is.null(stratavar)) {
    set.seed(seed_to_set)
    assign(x = paste0("f_",i,"_permanova"), adonis2(as.formula(paste0("flies_",i,"_dm ",pform)), flies_unifrac_table, permutations=1000, by = "margin"))
    print("margin")
  } else if (byval == "terms" & is.null(stratavar)) {
    set.seed(seed_to_set)
    assign(x = paste0("f_",i,"_permanova"), adonis2(as.formula(paste0("flies_",i,"_dm ",pform)), flies_unifrac_table, permutations=1000, by = "terms"))
    print("terms")
  } else if (byval == "margin" & !is.null(stratavar)) {
    set.seed(seed_to_set)
    assign(x = paste0("f_",i,"_permanova"), adonis2(as.formula(paste0("flies_",i,"_dm ",pform)), data = flies_unifrac_table, permutations=1000, by = "margin", strata = flies_unifrac_table[,paste0(stratavar)]))
    print("marginstrata")
  } else if (byval == "terms" & !is.null(stratavar)) {
    set.seed(seed_to_set)
    assign(x = paste0("f_",i,"_permanova"), adonis2(as.formula(paste0("flies_",i,"_dm ",pform)), data = flies_unifrac_table, permutations=1000, by = "terms", strata = flies_unifrac_table[,paste0(stratavar)]))
    print("termsstrata")
  }
  
  
  assign(x = paste0("fig_",i), value = tableGrob(data.frame(Df = round(get(paste0("f_",i,"_permanova"))$Df, 2), SS = round(get(paste0("f_",i,"_permanova"))$SumOfSqs, 2), R2 = round(get(paste0("f_",i,"_permanova"))$R2, 2), Fval = round(get(paste0("f_",i,"_permanova"))$`F`, 2), p = round(get(paste0("f_",i,"_permanova"))$`Pr(>F)`, 4)),theme=ttheme_minimal()))
  plot(get(paste0("fig_",i)))
  write.table(i, file = paste0('core-metrics-results-',file_path,"/",file_path,"_permanova_table.txt"), append = T, quote = F, sep = "\t", row.names = F, col.names = F)
  write.table(round(get(paste0("f_",i,"_permanova")),2), file = paste0('core-metrics-results-',file_path,"/",file_path,"_permanova_table.txt"), append = T, quote = F, sep = "\t", row.names = T, col.names = T)
  }
  if(!is.null(converttable)) {
    convert_permanova_to_word(file_path = file_path, table_row_names = converttable)
  }
}


format_sig_digits <- function(x) {
  ## written by chat
  # Calculate how many significant digits the number naturally has
  # by removing the decimal point and leading/trailing zeros
  raw_digits <- nchar(gsub("(^0\\.|\\.0+$|\\.)", "", as.character(abs(x))))
  
  # Use 2 significant digits as a minimum, otherwise use the natural count
  digits_to_use <- pmax(2, raw_digits)
  
  # Apply formatting individually to each number
  sapply(seq_along(x), function(i) {
    signif(x[i], digits = digits_to_use[i])
  })
}

convert_permanova_to_word <- function(file_path = file_path, table_row_names = table_row_names) {


  
  perm_table <- read.table(paste0("core-metrics-results-",file_path,"/",file_path,"_permanova_table.txt"), fill = T, sep = "\t", skip = 2, header = F, row.names = NULL)
  
  num_rows = (which(perm_table$V1 == "weighted_unifrac")-1)
  
  unweighted_unifrac = perm_table[1:num_rows,3:6]
  weighted_unifrac = perm_table[(num_rows+2):(num_rows+1+num_rows),3:6]
  bray_curtis = perm_table[(num_rows*2+3):(num_rows*3+2),1:6]
  
  new_df <- cbind(bray_curtis, ` ` = "", weighted_unifrac, ` ` = "", unweighted_unifrac)
  new_df2 <- rbind(c("","Bray-Curtis","Bray-Curtis","Bray-Curtis","Bray-Curtis","Bray-Curtis","","Weighted Unifrac","Weighted Unifrac","Weighted Unifrac","Weighted Unifrac","","Unweighted Unifrac","Unweighted Unifrac","Unweighted Unifrac","Unweighted Unifrac"),c("","Df","SS","R2","F","p","","SS","R2","F","p","","SS","R2","F","p"),new_df)[-3,]
  colnames(new_df2) <- paste0("V",1:16)
  new_df2$V1 <- c("","",table_row_names,"Residual","Total")
  new_df2$V3[3:(num_rows+1)] <- as.character(sprintf("%.1f", as.numeric(new_df2$V3[3:(num_rows+1)])))
  new_df2$V4[3:(num_rows+1)] <- as.character(sprintf("%.2f", as.numeric(new_df2$V4[3:(num_rows+1)])))
  new_df2$V5[3:(num_rows+1)] <- as.character(sprintf("%.1f", as.numeric(new_df2$V5[3:(num_rows+1)])))
  new_df2$V8[3:(num_rows+1)] <- as.character(sprintf("%.1f", as.numeric(new_df2$V8[3:(num_rows+1)])))
  new_df2$V9[3:(num_rows+1)] <- as.character(sprintf("%.2f", as.numeric(new_df2$V9[3:(num_rows+1)])))
  new_df2$V10[3:(num_rows+1)] <- as.character(sprintf("%.1f", as.numeric(new_df2$V10[3:(num_rows+1)])))
  new_df2$V13[3:(num_rows+1)] <- as.character(sprintf("%.1f", as.numeric(new_df2$V13[3:(num_rows+1)])))
  new_df2$V14[3:(num_rows+1)] <- as.character(sprintf("%.2f", as.numeric(new_df2$V14[3:(num_rows+1)])))
  new_df2$V15[3:(num_rows+1)] <- as.character(sprintf("%.1f", as.numeric(new_df2$V15[3:(num_rows+1)])))
  new_df2$V5[(num_rows):(num_rows+1)] <- ""
  new_df2$V10[(num_rows):(num_rows+1)] <- ""
  new_df2$V15[(num_rows):(num_rows+1)] <- ""
  new_df2$V6[3:(num_rows-1)] <- ifelse(new_df2$V6[3:(num_rows-1)] == 0, 0.001,new_df2$V6[3:(num_rows-1)])
  new_df2$V11[3:(num_rows-1)] <- ifelse(new_df2$V11[3:(num_rows-1)] == 0, 0.001,new_df2$V11[3:(num_rows-1)])
  new_df2$V16[3:(num_rows-1)] <- ifelse(new_df2$V16[3:(num_rows-1)] == 0, 0.001,new_df2$V16[3:(num_rows-1)])
  col6replace <- which(as.character(new_df2$V6) == "0.001")
  col11replace <- which(as.character(new_df2$V11) == "0.001")
  col16replace <- which(as.character(new_df2$V16) == "0.001")
  
  format_sig_digits(as.numeric(new_df2$V6[3:(num_rows-1)]))
  
  ## working with flextable
  new_df3 <- flextable(new_df2)
  new_df3 <- compose(
    x = new_df3, 
    i = 2,
    j = c("V4","V9","V14"),
    value = as_paragraph("R", as_sup("2"))
  )
  
  new_df3 <- compose(x = new_df3, i = col6replace, j = "V6", value = as_paragraph("<10", as_sup("-3")))
  new_df3 <- compose(x = new_df3, i = col11replace, j = "V11", value = as_paragraph("<10", as_sup("-3")))
  new_df3 <- compose(x = new_df3, i = col16replace, j = "V16", value = as_paragraph("<10", as_sup("-3")))
  
  new_df3 <- merge_h(new_df3, i = 1)    # Merge identical vertical cells in Department
  
  new_df3 <- align(new_df3, align = "center", part = "all") # Center-align all text
  new_df3 <- align(new_df3, align = "right", j = "V1") # Vertically center merged text
  new_df3 <- autofit(new_df3)                    # Adjust column widths automatically
  
    new_df3 <- font(new_df3, fontname = "Cambria", part = "all")
  new_df3 <- fontsize(new_df3, size = 6, part = "all")
  
  doc <- read_docx()
  doc <- body_add_flextable(doc, value = new_df3)
  print(doc, target = paste0('core-metrics-results-',file_path,'/',file_path,'_permanova_table.docx'))
  print(new_df3)

}
