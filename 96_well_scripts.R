latitude_correlation_plot <- function(dataframe, x_axis_column, y_axis_column, log10_transformed = T, point_color = "black") {
  test1 <- cor.test(dataframe[,x_axis_column], dataframe[,y_axis_column], method = "spearman", exact = F)
  if(log10_transformed == T){
    return(ggplot(dataframe, aes(get(x_axis_column), log10(get(y_axis_column)+1))) + geom_jitter(width = 0.1, color = point_color) + geom_smooth(method = "lm") + theme_cowplot() + labs(x = "Latitude",y = bquote(Log[10]~"+1 CFU fly"^-1), title = bquote("S = "~.(test1$statistic)~"   "~rho~" = "~.(round(test1$estimate,2))~"   p = "~.(round(test1$p.value,4)))) + theme(plot.title = element_text(hjust = 1)))
  } else {
    return(ggplot(dataframe, aes(get(x_axis_column), get(y_axis_column))) + geom_jitter(width = 0.1, color = point_color) + geom_smooth(method = "lm") + theme_cowplot() + labs(x = "Latitude", y = bquote(Log[10]~"+1 CFU fly"^-1), title = bquote("S = "~.(test1$statistic)~"   "~rho~" = "~.(round(test1$estimate,2))~"   p = "~.(round(test1$p.value,4)))) + theme(plot.title = element_text(hjust = 1)))
  }
}


latitude_correlation_multiplot <- function(dataframe, x_axis_column, y_axis_column, log10_transformed = T, plot_colors = rep(c("black"),4), filename) {
  plot1 <- latitude_correlation_plot(dataframe, x_axis_column, y_axis_column[1], log10_transformed = log10_transformed, point_color = plot_colors[1])
  plot2 <- latitude_correlation_plot(dataframe, x_axis_column, y_axis_column[2], log10_transformed = log10_transformed, point_color = plot_colors[2])
  plot3 <- latitude_correlation_plot(dataframe, x_axis_column, y_axis_column[3], log10_transformed = log10_transformed, point_color = plot_colors[3])
  plot4 <- latitude_correlation_plot(dataframe, x_axis_column, y_axis_column[4], log10_transformed = F, point_color = plot_colors[4])
  jpeg(filename, h = 1000, w = 1000)
  grid.arrange(plot1, plot2, plot3, plot4, heights = c(2,2), widths = c(2,2))
  dev.off()
}


# make_chart_stats_96 <- function(bac, column, plotout, yaxis = NULL, factorspecs = NULL, aceto_bottom =T, label_adjust = 10000, bac) {
#   bac$testcol <- unlist(unname(bac[,column]))
#   bac <- bac %>%
#     mutate(aperc = cfuA/(cfuA+cfuL), lperc = cfuL/(cfuA+cfuL)) %>%
#     droplevels()
#   if(!is.null(factorspecs)) {
#     bac$testcol = factor(bac$testcol, levels = factorspecs)
#   } else {
#     bac$testcol = factor(bac$testcol)
#     print(levels(factor(bac$testcol)))
#   }
#   bac3 <- bac %>% 
#     mutate(total = cfuA + cfuL) %>%
#     group_by(repexp, exp,testcol) %>%
#     dplyr::summarize(cfuA = median(cfuA, na.rm=T), cfuL = median(cfuL, na.rm=T), total = median(total, na.rm=T)) %>%
#     ungroup() %>%
#     group_by(exp,testcol) %>%
#     dplyr::summarize(cfuA = mean(cfuA, na.rm=T), cfuL = mean(cfuL, na.rm=T), total = mean(total, na.rm=T)) %>%
#     ungroup() %>%
#     group_by(testcol) %>%
#     dplyr::summarize(cfuA = mean(cfuA, na.rm=T), cfuL = mean(cfuL, na.rm=T), total = mean(total, na.rm=T)) %>%
#     ungroup() %>%
#     mutate(total = cfuA + cfuL) %>%
#     reshape2::melt(measure.vars = c("cfuA","cfuL")) %>%
#     mutate(tempvar = paste0(testcol,"_",variable)) %>%
#     inner_join(bac %>% 
#                  mutate(total = cfuA + cfuL) %>%
#                  group_by(repexp, exp,testcol) %>%
#                  dplyr::summarize(cfuA = median(cfuA, na.rm=T), cfuL = median(cfuL, na.rm=T), total = median(total, na.rm=T)) %>%
#                  ungroup() %>%
#                  group_by(exp,testcol) %>%
#                  dplyr::summarize(cfuA = mean(cfuA, na.rm=T), cfuL = mean(cfuL, na.rm=T), total = mean(total, na.rm=T)) %>%
#                  ungroup() %>%
#                  group_by(testcol) %>%
#                  dplyr::summarize(semcfuA = sd(cfuA, na.rm=T)/sqrt(dplyr::n()), semcfuL = sd(cfuL, na.rm=T)/sqrt(dplyr::n()), semtotal = sd(total, na.rm=T)/sqrt(dplyr::n())) %>%
#                  ungroup() %>%
#                  reshape2::melt(measure.vars = c("semcfuA","semcfuL"))  %>%
#                  mutate(mergecol = gsub(pattern = "sem", replacement = "", x = variable)) %>%
#                  mutate(tempvar = paste0(testcol,"_",mergecol)), by = "tempvar") 
#   
#   bac3$addcol <- bac3$value.x
#   numvals <- length(table(list(bac$testcol %>% droplevels)))
#   
#   if (aceto_bottom == T) {
#     bac3$variable.x = factor(bac3$variable.x, levels = c("cfuL","cfuA"))
#     for (i in (numvals+1):(numvals*2)) {bac3$addcol[i] = bac3$total[i]}
#   } else {
#     bac3$variable.x = factor(bac3$variable.x, levels = c("cfuA","cfuL"))
#     for (i in 1:numvals) {bac3$addcol[i] = bac3$total[i]}
#   }
#   
#   if(kruskal.test(bac$cfuA ~ bac$testcol)$p.value < 0.05) {
#     if (length(table(bac$testcol))>2) {
#       efg <- dunn.test(bac$cfuA,bac$testcol, method = "bh", table = F, kw = T)
#       ghi <- cldList(comparison = efg$comparisons, p.value = efg$P.adjusted, threshold = 0.05)
#       ghi$Group2 <- factor(ghi$Group, levels = levels(bac$testcol))
#       cfuAclds <- as.character(ghi %>% arrange(Group2) %>% dplyr::select(Letter) %>% unlist() %>% unname())
#     } else {
#       cfuAclds <- c("a","b")
#     }
#   } else {
#     cfuAclds <- c(rep("a", numvals))
#   }
#   
#   if(kruskal.test(bac$cfuL ~ bac$testcol)$p.value < 0.05) {
#     if (length(table(bac$testcol))>2) {
#       efg <- dunn.test(bac$cfuL,bac$testcol, method = "bh", table = F, kw = F)
#       ghi <- cldList(comparison = efg$comparisons, p.value = efg$P.adjusted, threshold = 0.05)
#       ghi$Group2 <- factor(ghi$Group, levels = levels(bac$testcol))
#       cfuLclds <- as.character(ghi %>% arrange(Group2) %>% dplyr::select(Letter) %>% unlist() %>% unname())
#     } else {
#       cfuLclds <- c("a","b")
#     }
#   } else {
#     cfuLclds <- c(rep("a", numvals))
#   }
#   
#   cat("absolute abundance stats")
#   print(kruskal.test(bac$cfuA ~ bac$testcol))
#   print(kruskal.test(bac$cfuL ~ bac$testcol))
#   print(dim(bac))
#   
#   if(!is.null(yaxis)) {    
#     if(aceto_bottom != T) {
#       bac3plot <- ggplot(bac3, aes(x=testcol.x, y=value.x, fill = variable.x)) +
#         geom_bar(stat = "identity") + 
#         geom_errorbar(aes(ymin =addcol-value.y, ymax = addcol+value.y), width = 0.5) +
#         scale_fill_manual(values = c("red","blue"))+
#         theme(legend.position = "bottom") + 
#         ylim(yaxis)+
#         geom_text(x = c(1:numvals,1:numvals), aes(y = 0, label = c(rep("",numvals),cfuLclds)), vjust = 1)+
#         geom_text(x = c(1:numvals,1:numvals), aes(y = addcol+label_adjust+value.y, label = c(cfuAclds,rep("",numvals))))+
#         theme_cowplot()
#     } else {
#       bac3plot <- ggplot(bac3, aes(x=testcol.x, y=value.x, fill = variable.x)) +
#         geom_bar(stat = "identity") + 
#         geom_errorbar(aes(ymin =addcol-value.y, ymax = addcol+value.y), width = 0.5) +
#         scale_fill_manual(values = c("blue","red"))+
#         theme(legend.position = "bottom") +  
#         ylim(yaxis)+
#         geom_text(x = c(1:numvals,1:numvals), aes(y = 0, label = c(cfuAclds,rep("",numvals))), vjust = 1)+
#         geom_text(x = c(1:numvals,1:numvals), aes(y = addcol+label_adjust+value.y, label = c(rep("",numvals),cfuLclds)))+
#         theme_cowplot()
#     }
#   } else {
#     if(aceto_bottom != T) {
#       bac3plot <- ggplot(bac3, aes(x=testcol.x, y=value.x, fill = variable.x)) +
#         geom_bar(stat = "identity") + 
#         geom_errorbar(aes(ymin =addcol-value.y, ymax = addcol+value.y), width = 0.5) +
#         scale_fill_manual(values = c("red","blue"))+
#         theme(legend.position = "bottom") + 
#         geom_text(x = c(1:numvals,1:numvals), aes(y = 0, label = c(rep("",numvals),cfuLclds)), vjust = 1)+
#         geom_text(x = c(1:numvals,1:numvals), aes(y = addcol+label_adjust+value.y, label = c(cfuAclds,rep("",numvals))))+
#         theme_cowplot()
#     } else {
#       bac3plot <- ggplot(bac3, aes(x=testcol.x, y=value.x, fill = variable.x)) +
#         geom_bar(stat = "identity") + 
#         geom_errorbar(aes(ymin =addcol-value.y, ymax = addcol+value.y), width = 0.5) +
#         scale_fill_manual(values = c("blue","red"))+
#         theme(legend.position = "bottom") +  
#         geom_text(x = c(1:numvals,1:numvals), aes(y = 0, label = c(cfuAclds,rep("",numvals))), vjust = 1)+
#         geom_text(x = c(1:numvals,1:numvals), aes(y = addcol+label_adjust+value.y, label = c(rep("",numvals),cfuLclds)))+
#         theme_cowplot()
#     }
#   }
#   
#   # old
#   bac4 <- bac %>% 
#     mutate(total = cfuA + cfuL) %>%
#     group_by(repexp, exp,testcol) %>%
#     dplyr::summarize(cfuA = median(cfuA/total, na.rm=T), cfuL = median(cfuL/total, na.rm=T),  total = median(total/total, na.rm=T)) %>%
#     ungroup() %>%
#     group_by(exp,testcol) %>%
#     dplyr::summarize(cfuA = mean(cfuA, na.rm=T), cfuL = mean(cfuL, na.rm=T), total = mean(total, na.rm=T)) %>%
#     ungroup() %>%
#     group_by(testcol) %>%
#     dplyr::summarize(cfuA = mean(cfuA, na.rm=T), cfuL = mean(cfuL, na.rm=T), total = mean(total, na.rm=T)) %>%
#     ungroup() %>%
#     rowwise() %>%
#     mutate(totalcol = cfuA+cfuL) %>%
#     mutate(aperc = cfuA/totalcol, lperc = cfuL/totalcol) %>%
#     ungroup() %>%
#     reshape2::melt(measure.vars = c("aperc","lperc"))
#   
#   bac4$addcol <- bac4$value
#   if (aceto_bottom == T) {
#     bac4$variable = factor(bac4$variable, levels = c("lperc","aperc"))
#     for (i in (numvals+1):(numvals*2)) {bac4$addcol[i] = bac4$value[i-numvals]}
#   } else {
#     bac4$variable = factor(bac4$variable, levels = c("aperc","lperc"))
#     for (i in 1:numvals) {bac4$addcol[i] = bac4$value[i+numvals]}
#   }
#   
#   length(table(bac$testcol))
#   if(kruskal.test(bac$aperc ~ bac$testcol)$p.value < 0.05) {
#     if (length(table(bac$testcol))>2) {
#       efg <- dunn.test(bac$aperc,bac$testcol, method = "bh", table = F, kw = F)
#       ghi <- cldList(comparison = efg$comparisons, p.value = efg$P.adjusted, threshold = 0.05) 
#       ghi$Group2 <- factor(ghi$Group, levels = levels(bac$testcol))
#       apercclds <- as.character(ghi %>% arrange(Group2) %>% dplyr::select(Letter) %>% unlist() %>% unname())
#     } else {
#       apercclds <- c("a","b")
#     }
#   } else {
#     apercclds <- c(rep("a",numvals))
#   }
#   
#   if(kruskal.test(bac$lperc ~ bac$testcol)$p.value < 0.05) {
#     if (length(table(bac$testcol))>2) {
#       efg <- dunn.test(bac$lperc,bac$testcol, method = "bh", table = F, kw = F)
#       ghi <- cldList(comparison = efg$comparisons, p.value = efg$P.adjusted, threshold = 0.05) 
#       ghi$Group2 <- factor(ghi$Group, levels = levels(bac$testcol))
#       lpercclds <- as.character(ghi %>% arrange(Group2) %>% dplyr::select(Letter) %>% unlist() %>% unname())
#     } else {
#       lpercclds <- c("a","b")
#     }
#   } else {
#     lpercclds <-c(rep("a",numvals))
#   }
#   
#   cat("relative abundance stats")
#   print(kruskal.test(bac$aperc ~ bac$testcol))
#   print(kruskal.test(bac$lperc ~ bac$testcol))
#   
#   
#   if(aceto_bottom != T) {
#     
#     
#     bac4$variable <- factor(bac4$variable, levels = c("aperc", "lperc"))
#     bac4plot <- ggplot(bac4, aes(x=testcol, y=value, fill = variable)) +
#       geom_bar(stat = "identity", position = "stack") + theme(legend.position = "bottom") +
#       scale_fill_manual(values = c("red","blue"))+
#       geom_text(x = c(1:numvals,1:numvals), aes(y = 0, label = c(rep("",numvals),cfuLclds)), vjust = 1)+
#       geom_text(x = c(1:numvals,1:numvals), aes(y = addcol+0.05, label = c(cfuAclds,rep("",numvals))))+
#       theme_cowplot()
#     
#   } else {
#     bac4$variable <- factor(bac4$variable, levels = c("lperc", "aperc"))
#     
#     bac4plot <-1 #
#       ggplot(bac4, aes(x=testcol, y=value, fill = variable)) +
#       geom_bar(position = "dodge") + theme(legend.position = "bottom") +
#       scale_fill_manual(values = c("blue","red"))+
#       geom_text(x = c(1:numvals,1:numvals), aes(y = 0, label = c(cfuAclds,rep("",numvals))), vjust = 1)+
#       geom_text(x = c(1:numvals,1:numvals), aes(y = addcol+0.05, label = c(rep("",numvals),cfuLclds)))+
#       theme_cowplot()
#   }
#   
#   bac5plot <- make_chart_stats_96_log(bac = bac, column = column, yaxis = yaxis, factorspecs = factorspecs, aceto_bottom = aceto_bottom)
#   ifelse(plotout == 1, return(bac3plot),
#          ifelse(plotout == 2, return(bac5plot),
#                 ifelse(plotout == 3, return(bac4plot),
#                        ifelse(plotout == 4, return(grid.arrange(bac3plot, bac5plot, heights = c(1,1))),
#                               ifelse(plotout == 5, return(grid.arrange(bac3plot, bac4plot, heights = c(1,1))),
#                                      ifelse(plotout == 6, return(grid.arrange(bac4plot, bac5plot, heights = c(1,1))),
#                                             ifelse(plotout == 7, return(grid.arrange(bac3plot, bac5plot, bac4plot, heights = c(1,1,1))),
#                                                    ifelse(plotout == 8, return(list(bac, bac3, bac4))))))))))
# }

# make_chart_stats_96_log <- function(bac, column, yaxis = NULL, factorspecs = NULL, aceto_bottom =T, label_adjust = 0.2) {
#   #bac$testcol <- bac[,column]
#   bac$testcol <- unlist((bac[,column]))
#   bac <- bac %>%
#     mutate(aperc = cfuA/(cfuA+cfuL), lperc = cfuL/(cfuA+cfuL)) %>%
#     droplevels()
#   if(!is.null(factorspecs)) {
#     bac$testcol = factor(bac$testcol, levels = factorspecs)
#   } else {
#     bac$testcol = factor(bac$testcol)
#     print(levels(factor(bac$testcol)))
#   }
#   bac3 <- bac %>% 
#     group_by(repexp, exp,testcol) %>%
#     dplyr::summarize(cfuA = mean(log10(cfuA+1), na.rm=T), cfuL = mean(log10(cfuL+1), na.rm=T), total = mean(log10(total+1), na.rm=T)) %>%
#     ungroup() %>%
#     group_by(exp,testcol) %>%
#     dplyr::summarize(cfuA = mean(cfuA, na.rm=T), cfuL = mean(cfuL, na.rm=T), total = mean(total, na.rm=T)) %>%
#     ungroup() %>%
#     group_by(testcol) %>%
#     dplyr::summarize(cfuA = mean(cfuA, na.rm=T), cfuL = mean(cfuL, na.rm=T), total = mean(total, na.rm=T)) %>%
#     ungroup() %>%
#     reshape2::melt(measure.vars = c("cfuA","cfuL")) %>%
#     mutate(tempvar = paste0(testcol,"_",variable)) %>%
#     inner_join(bac %>% 
#                  mutate(total = cfuA + cfuL) %>%
#                  group_by(repexp, exp,testcol) %>%
#                  dplyr::summarize(cfuA = mean(log10(cfuA+1), na.rm=T), cfuL = mean(log10(cfuL+1), na.rm=T), total = mean(log10(total+1), na.rm=T)) %>%
#                  ungroup() %>%
#                  group_by(exp,testcol) %>%
#                  dplyr::summarize(cfuA = mean(cfuA, na.rm=T), cfuL = mean(cfuL, na.rm=T), total = mean(total, na.rm=T)) %>%
#                  ungroup() %>%
#                  group_by(testcol) %>%
#                  dplyr::summarize(semcfuA = sd(cfuA, na.rm=T)/sqrt(dplyr::n()), semcfuL = sd(cfuL, na.rm=T)/sqrt(dplyr::n()), semtotal = sd(total, na.rm=T)/sqrt(dplyr::n())) %>%
#                  ungroup() %>%
#                  reshape2::melt(measure.vars = c("semcfuA","semcfuL"))  %>%
#                  mutate(mergecol = gsub(pattern = "sem", replacement = "", x = variable)) %>%
#                  mutate(tempvar = paste0(testcol,"_",mergecol)), by = "tempvar") 
#   
#   bac3$addcol <- bac3$value.x
#   numvals <- length(table(list(bac$testcol %>% droplevels)))
#   
#   if (aceto_bottom == T) {
#     bac3$variable.x = factor(bac3$variable.x, levels = c("cfuL","cfuA"))
#     for (i in (numvals+1):(numvals*2)) {
#       bac3$addcol[i] = bac3$total[i] # sets the SEM center values
#       bac3$value.x[i] = bac3$total[i] - bac3$value.x[i-numvals] # adjusts the log-transformed values - not sure what to do about SEMs, though.
#     }
#   } else {
#     bac3$variable.x = factor(bac3$variable.x, levels = c("cfuA","cfuL"))
#     for (i in 1:numvals) {
#       bac3$addcol[i] = bac3$total[i]
#       bac3$value.x[i] = bac3$total[i]- bac3$value.x[i+numvals]# adjusts the log-transformed values - not sure what to do about SEMs, though.
#     }
#   }
#   
#   if(kruskal.test(bac$cfuA ~ bac$testcol)$p.value < 0.05) {
#     if (length(table(bac$testcol))>2) {
#       efg <- dunn.test(bac$cfuA,bac$testcol, method = "bh", table = F, kw = F)
#       ghi <- cldList(comparison = efg$comparisons, p.value = efg$P.adjusted, threshold = 0.05)
#       ghi$Group2 <- factor(ghi$Group, levels = levels(bac$testcol))
#       cfuAclds <- as.character(ghi %>% arrange(Group2) %>% dplyr::select(Letter) %>% unlist() %>% unname())
#     } else {
#       cfuAclds <- c("a","b")
#     }
#   } else {
#     cfuAclds <- c(rep("a", numvals))
#   }
#   
#   if(kruskal.test(bac$cfuL ~ bac$testcol)$p.value < 0.05) {
#     if (length(table(bac$testcol))>2) {
#       efg <- dunn.test(bac$cfuL,bac$testcol, method = "bh", table = F, kw = F)
#       ghi <- cldList(comparison = efg$comparisons, p.value = efg$P.adjusted, threshold = 0.05)
#       ghi$Group2 <- factor(ghi$Group, levels = levels(bac$testcol))
#       cfuLclds <- as.character(ghi %>% arrange(Group2) %>% dplyr::select(Letter) %>% unlist() %>% unname())
#     } else {
#       cfuLclds <- c("a","b")
#     }
#   } else {
#     cfuLclds <- c(rep("a", numvals))
#   }
#   
#   if(!is.null(yaxis)) {    
#     if(aceto_bottom != T) {
#       bac3plot <- ggplot(bac3, aes(x=testcol.x, y=value.x, fill = variable.x)) +
#         geom_bar(stat = "identity") + 
#         geom_errorbar(aes(ymin =addcol-value.y, ymax = addcol+value.y), width = 0.5) +
#         scale_fill_manual(values = c("red","blue"))+
#         theme(legend.position = "bottom") + 
#         ylim(yaxis)+
#         geom_text(x = c(1:numvals,1:numvals), aes(y = addcol-label_adjust-value.y, label = c(rep("",numvals),cfuLclds)))+
#         geom_text(x = c(1:numvals,1:numvals), aes(y = addcol+label_adjust+value.y, label = c(cfuAclds,rep("",numvals))))+
#         theme_cowplot()
#     } else {
#       bac3plot <- ggplot(bac3, aes(x=testcol.x, y=value.x, fill = variable.x)) +
#         geom_bar(stat = "identity") + 
#         geom_errorbar(aes(ymin =addcol-value.y, ymax = addcol+value.y), width = 0.5) +
#         scale_fill_manual(values = c("blue","red"))+
#         theme(legend.position = "bottom") +  
#         ylim(yaxis)+
#         geom_text(x = c(1:numvals,1:numvals), aes(y = addcol-label_adjust-value.y, label = c(cfuAclds,rep("",numvals))))+
#         geom_text(x = c(1:numvals,1:numvals), aes(y = addcol+label_adjust+value.y, label = c(rep("",numvals),cfuLclds)))+
#         theme_cowplot()
#     }
#   } else {
#     if(aceto_bottom != T) {
#       bac3plot <- ggplot(bac3, aes(x=testcol.x, y=value.x, fill = variable.x)) +
#         geom_bar(stat = "identity") + 
#         geom_errorbar(aes(ymin =addcol-value.y, ymax = addcol+value.y), width = 0.5) +
#         scale_fill_manual(values = c("red","blue"))+
#         theme(legend.position = "bottom") + 
#         geom_text(x = c(1:numvals,1:numvals), aes(y = addcol-label_adjust-value.y, label = c(rep("",numvals),cfuLclds)))+
#         geom_text(x = c(1:numvals,1:numvals), aes(y = addcol+label_adjust+value.y, label = c(cfuAclds,rep("",numvals))))+
#         theme_cowplot()
#     } else {
#       bac3plot <- ggplot(bac3, aes(x=testcol.x, y=value.x, fill = variable.x)) +
#         geom_bar(stat = "identity") + 
#         geom_errorbar(aes(ymin =addcol-value.y, ymax = addcol+value.y), width = 0.5) +
#         scale_fill_manual(values = c("blue","red"))+
#         theme(legend.position = "bottom") +  
#         geom_text(x = c(1:numvals,1:numvals), aes(y = addcol-label_adjust-value.y, label = c(cfuAclds,rep("",numvals))))+
#         geom_text(x = c(1:numvals,1:numvals), aes(y = addcol+label_adjust+value.y, label = c(rep("",numvals),cfuLclds)))+
#         theme_cowplot()
#     }
#   }
# }

# bac = unstarved_plates_converted %>% filter(Sex == "M") %>% filter(total < 1000000& total > 2000& !is.na(cfuA) & !is.na(cfuL)& Geography != "MD"& Geography != "PA") %>% droplevels()
# colum = "GT"
# plotout = 3
# yaxis = NULL
# factorspecs = NULL
# aceto_bottom =F 
# label_adjust = 10000

make_chart_stats_96_perclog <- function(bac, column, plotout, yaxis = NULL, factorspecs = NULL, aceto_bottom =T, label_adjust = 10000) {
  
  bac$testcol <- unlist((bac[,column]))
  bac <- bac %>%
    mutate(aperc = cfuA/(cfuA+cfuL), lperc = cfuL/(cfuA+cfuL)) %>%
    droplevels()
  if(!is.null(factorspecs)) {
    bac$testcol = factor(bac$testcol, levels = factorspecs)
  } else {
    bac$testcol = factor(bac$testcol)
    print(levels(factor(bac$testcol)))
  }
  # old
  
  numvals <- length(table(list(bac$testcol %>% droplevels)))
  
  bac4 <- bac %>% 
    mutate(total = cfuA + cfuL) %>%
    group_by(repexp, exp,testcol) %>%
    dplyr::summarize(cfuA = median(cfuA/total, na.rm=T), cfuL = median(cfuL/total, na.rm=T), total = median(total/total, na.rm=T)) %>%
    ungroup() %>%
    group_by(exp,testcol) %>%
    dplyr::summarize(cfuA = mean(cfuA, na.rm=T), cfuL = mean(cfuL, na.rm=T), total = mean(total, na.rm=T)) %>%
    ungroup() %>%
    group_by(testcol) %>%
    dplyr::summarize(aperc = mean(cfuA, na.rm=T), lperc = mean(cfuL, na.rm=T), total = mean(total, na.rm=T)) %>%
    ungroup() %>%
    rowwise() %>%
  #  mutate(totalcol = cfuA+cfuL) %>%
  #  mutate(aperc = cfuA/totalcol, lperc = cfuL/totalcol) %>%
    ungroup() %>%
    reshape2::melt(measure.vars = c("aperc","lperc"))
  
  bac4$addcol <- bac4$value
  if (aceto_bottom == T) {
    bac4$variable = factor(bac4$variable, levels = c("lperc","aperc"))
    for (i in (numvals+1):(numvals*2)) {bac4$addcol[i] = bac4$value[i-numvals]}
  } else {
    bac4$variable = factor(bac4$variable, levels = c("aperc","lperc"))
    for (i in 1:numvals) {bac4$addcol[i] = bac4$value[i+numvals]}
  }
  
  length(table(bac$testcol))
  if(kruskal.test(bac$aperc ~ bac$testcol)$p.value < 0.05) {
    if (length(table(bac$testcol))>2) {
      efg <- dunn.test(bac$aperc,bac$testcol, method = "bh", table = F, kw = F)
      ghi <- cldList(comparison = efg$comparisons, p.value = efg$P.adjusted, threshold = 0.05) 
      ghi$Group2 <- factor(ghi$Group, levels = levels(bac$testcol))
      apercclds <- as.character(ghi %>% arrange(Group2) %>% dplyr::select(Letter) %>% unlist() %>% unname())
    } else {
      apercclds <- c("a","b")
    }
  } else {
    apercclds <- c(rep("a",numvals))
  }
  
  if(kruskal.test(bac$lperc ~ bac$testcol)$p.value < 0.05) {
    if (length(table(bac$testcol))>2) {
      efg <- dunn.test(bac$lperc,bac$testcol, method = "bh", table = F, kw = F)
      ghi <- cldList(comparison = efg$comparisons, p.value = efg$P.adjusted, threshold = 0.05) 
      ghi$Group2 <- factor(ghi$Group, levels = levels(bac$testcol))
      lpercclds <- as.character(ghi %>% arrange(Group2) %>% dplyr::select(Letter) %>% unlist() %>% unname())
    } else {
      lpercclds <- c("a","b")
    }
  } else {
    lpercclds <-c(rep("a",numvals))
  }  
  if(aceto_bottom != T) {
    
    log10(bac4$value+1)
    
    bac4plot <- ggplot(bac4, aes(x=testcol, y=value, fill = variable)) +
      geom_bar(stat = "identity", position = "stack") + theme(legend.position = "bottom") +
      scale_fill_manual(values = c("red","blue"))+
      geom_text(x = c(1:numvals,1:numvals), aes(y = 0, label = c(rep("",numvals),cfuLclds)), vjust = 1)+
      geom_text(x = c(1:numvals,1:numvals), aes(y = addcol+0.05, label = c(cfuAclds,rep("",numvals))))+
      theme_cowplot() + 
    scale_y_continuous(trans = "log10") + coord_cartesian(ylim = c(0.001,1))
      coord_trans(y="log2")
  bac4plot  
  } else {
    bac4plot <- ggplot(bac4, aes(x=testcol, y=value, fill = variable)) +
      geom_bar(stat = "identity", position = "stack") + theme(legend.position = "bottom") +
      scale_fill_manual(values = c("blue","red"))+
      geom_text(x = c(1:numvals,1:numvals), aes(y = 0, label = c(cfuAclds,rep("",numvals))), vjust = 1)+
      geom_text(x = c(1:numvals,1:numvals), aes(y = addcol+0.05, label = c(rep("",numvals),cfuLclds)))+
      theme_cowplot()
  }
  
  
  return(bac4plot)
  
}

# cor_table <- cor_df
# vec <- NULL
# cor_column <- c("cfuA","latitude")
# logplot = F


# cortable <- readin7 %>%  filter(Treatment == "A" & Sex == "F") %>%  mutate(total = cfuA + cfuL)
# vec = vec
# cor_column =  c("perc","latitude")
# logplot = F


cortests_96well <- function(cortable, vec, cor_column, logplot = F) {
  readinA <- cortable %>% data.frame()
  
  if(!is.null(vec)){
    for(i in 1:length(vec)) {
      vec_subset <- c(vec[i:length(vec)], cor_column[2])
      readinA <- readinA %>%
        group_by(across(all_of(vec_subset))) %>%
        dplyr::summarize(cfuA = median(cfuA, na.rm=T), cfuL = median(cfuL, na.rm=T), total = median(total, na.rm=T), perc = median(perc, na.rm=T)) %>%
        #      dplyr::summarize(cfuA = mean(cfuA, na.rm=T), cfuL = mean(cfuL, na.rm=T), total = mean(total, na.rm=T)) %>%
        ungroup() %>%
        data.frame()
    }
  }
  
  print(cor.test(readinA[,paste0(cor_column[1])], readinA[,paste0(cor_column[2])], method = "spearman", exact=F))
  
  if(logplot == T) {
    ggplot(readinA, aes(x = get(cor_column[2]), y = log10(get(cor_column[1])+1))) + 
      geom_jitter(width = 0.15, height = 0.05) + 
      stat_smooth(method = "lm")
  } else {
    ggplot(readinA, aes(x = get(cor_column[2]), y = get(cor_column[1]))) + 
      geom_jitter(width = 0.15, height = 0.05) + 
      stat_smooth(method = "lm")
  }
}

# file_path = "maggie"
# mapper_file = "cfu_mapper.tsv"
# formula = "~ Treatment * Geography * starved * Sex + exp/repexp"
# permutations = 50
# bc_analysis <- function(file_path, mapper_file, formula, permutations = 999) {
#   
#   ## weighted
#   flies_unweighted <- read.table(paste0(file_path,'/bray_curtis_dm_',file_path,'/distance-matrix.tsv'), header=T, sep="\t", fill = T) %>% mutate(X=as.character(X))
#   flies_unweighted_dm <- as.dist(flies_unweighted[,2:dim(flies_unweighted)[2]])
#   flies_unifrac_table <- flies_unweighted %>% 
#     dplyr::select(X) %>% 
#     left_join(read.table(mapper_file,comment.char = "", header=T, fill=T, sep="\t"), by=c("X"="X.SampleID")) %>%
#     droplevels()
#   
#   set.seed(42)
#   weighted_permanova <- adonis2(as.formula(paste0("flies_unweighted_dm ",formula)), flies_unifrac_table, permutations=permutations)
#   
#   weighted_permanova$Df
#   w_df <- data.frame(variables = rownames(weighted_permanova), Df = round(weighted_permanova$Df, 2), SS = round(weighted_permanova$SumOfSqs, 2), R2 = round(weighted_permanova$R2, 2), Fval = round(weighted_permanova$`F`, 2), p = round(weighted_permanova$`Pr(>F)`, 2))
#   wtable <- tableGrob(data.frame(variables = rownames(weighted_permanova), Df = round(weighted_permanova$Df, 2), SS = round(weighted_permanova$SumOfSqs, 2), R2 = round(weighted_permanova$R2, 2), Fval = round(weighted_permanova$`F`, 2), p = round(weighted_permanova$`Pr(>F)`, 2)),theme=ttheme_minimal())
#   
#   ## unweighted
#   flies_unweighted <- read.table(paste0(file_path,'/core-metrics-results-',file_path,'/bray_curtis_distance_matrix/distance-matrix.tsv'), header=T, sep="\t") %>% mutate(X=as.character(X))
#   flies_unweighted_dm <- as.dist(flies_unweighted[,2:dim(flies_unweighted)[2]])
#   flies_unifrac_table <- flies_unweighted %>% 
#     dplyr::select(X) %>% 
#     left_join(read.table(mapper_file,comment.char = "", header=T, fill=T, sep="\t"), by=c("X"="X.SampleID")) %>%
#     droplevels()
#   
#   set.seed(42)
#   unweighted_permanova <- adonis2(as.formula(paste0("flies_unweighted_dm ",formula)), flies_unifrac_table, permutations=permutations)
#   uw_df <- data.frame(variables = rownames(unweighted_permanova), Df = round(unweighted_permanova$Df, 2), SS = round(unweighted_permanova$SumOfSqs, 2), R2 = round(unweighted_permanova$R2, 2), Fval = round(unweighted_permanova$`F`, 2), p = round(unweighted_permanova$`Pr(>F)`, 2))
#   uwtable <- tableGrob(data.frame(variables = rownames(unweighted_permanova), Df = round(unweighted_permanova$Df, 2), SS = round(unweighted_permanova$SumOfSqs, 2), R2 = round(unweighted_permanova$R2, 2), Fval = round(unweighted_permanova$`F`, 2), p = round(unweighted_permanova$`Pr(>F)`, 2)),theme=ttheme_minimal())
#   
#   # print(unweighted_permanova)
#   #  print(weighted_permanova)
#   
#   return(list(grid.arrange(wtable , uwtable), w_df, uw_df))
#   
# }


# bac = all_counts %>% filter(total < 1000000& total > 2000& !is.na(cfuA) & !is.na(cfuL), starved == "N") %>% droplevels()
# column = "TS"
# plotout = 3
# plot_aceto = F
#column = "GT"
make_chart_stats_96_jitter <- function(bac, column, plotout, yaxis = NULL, factorspecs = NULL, plot_aceto =T, label_adjust = 10000, jitter_width = 0.15, jitter_height = 0, meancol = "black", meanwidth = .25, baradjust = 1.02, plot_color = "black", aceto_bottom = F) {
  bac$testcol <- unlist(unname(bac[,column]))
  bac <- bac %>%
    mutate(aperc = cfuA/(cfuA+cfuL), lperc = cfuL/(cfuA+cfuL)) %>%
    droplevels()
  if(!is.null(factorspecs)) {
    bac$testcol = factor(bac$testcol, levels = factorspecs)
  } else {
    bac$testcol = factor(bac$testcol)
    print(levels(factor(bac$testcol)))
  }
  bac3 <- bac %>% 
    mutate(total = cfuA + cfuL) %>%
    group_by(repexp, exp,testcol) %>%
    dplyr::summarize(cfuA = median(cfuA, na.rm=T), cfuL = median(cfuL, na.rm=T), total = median(total, na.rm=T)) %>%
    ungroup() %>%
    group_by(exp,testcol) %>%
    dplyr::summarize(cfuA = mean(cfuA, na.rm=T), cfuL = mean(cfuL, na.rm=T), total = mean(total, na.rm=T)) %>%
    ungroup() %>%
    group_by(testcol) %>%
    dplyr::summarize(cfuA = mean(cfuA, na.rm=T), cfuL = mean(cfuL, na.rm=T), total = mean(total, na.rm=T)) %>%
    ungroup() %>%
    mutate(total = cfuA + cfuL) %>%
    reshape2::melt(measure.vars = c("cfuA","cfuL")) %>%
    mutate(tempvar = paste0(testcol,"_",variable)) %>%
    inner_join(bac %>% 
                 mutate(total = cfuA + cfuL) %>%
                 group_by(repexp, exp,testcol) %>%
                 dplyr::summarize(cfuA = median(cfuA, na.rm=T), cfuL = median(cfuL, na.rm=T), total = median(total, na.rm=T)) %>%
                 ungroup() %>%
                 group_by(exp,testcol) %>%
                 dplyr::summarize(cfuA = mean(cfuA, na.rm=T), cfuL = mean(cfuL, na.rm=T), total = mean(total, na.rm=T)) %>%
                 ungroup() %>%
                 group_by(testcol) %>%
                 dplyr::summarize(semcfuA = sd(cfuA, na.rm=T)/sqrt(dplyr::n()), semcfuL = sd(cfuL, na.rm=T)/sqrt(dplyr::n()), semtotal = sd(total, na.rm=T)/sqrt(dplyr::n())) %>%
                 ungroup() %>%
                 reshape2::melt(measure.vars = c("semcfuA","semcfuL"))  %>%
                 mutate(mergecol = gsub(pattern = "sem", replacement = "", x = variable)) %>%
                 mutate(tempvar = paste0(testcol,"_",mergecol)), by = "tempvar") 
  
  bac3$addcol <- bac3$value.x
  numvals <- length(table(list(bac$testcol %>% droplevels)))
  
  if (aceto_bottom == T) {
    bac3$variable.x = factor(bac3$variable.x, levels = c("cfuL","cfuA"))
    for (i in (numvals+1):(numvals*2)) {bac3$addcol[i] = bac3$total[i]}
  } else {
    bac3$variable.x = factor(bac3$variable.x, levels = c("cfuA","cfuL"))
    for (i in 1:numvals) {bac3$addcol[i] = bac3$total[i]}
  }
  
  if(kruskal.test(bac$cfuA ~ bac$testcol)$p.value < 0.05) {
    if (length(table(bac$testcol))>2) {
      efg <- dunn.test(bac$cfuA,bac$testcol, method = "bh", table = F, kw = T)
      ghi <- cldList(comparison = efg$comparisons, p.value = efg$P.adjusted, threshold = 0.05)
      ghi$Group2 <- factor(ghi$Group, levels = levels(bac$testcol))
      cfuAclds <- as.character(ghi %>% arrange(Group2) %>% dplyr::select(Letter) %>% unlist() %>% unname())
    } else {
      cfuAclds <- c("a","b")
    }
  } else {
    cfuAclds <- c(rep("a", numvals))
  }
  
  if(kruskal.test(bac$cfuL ~ bac$testcol)$p.value < 0.05) {
    if (length(table(bac$testcol))>2) {
      efg <- dunn.test(bac$cfuL,bac$testcol, method = "bh", table = F, kw = F)
      ghi <- cldList(comparison = efg$comparisons, p.value = efg$P.adjusted, threshold = 0.05)
      ghi$Group2 <- factor(ghi$Group, levels = levels(bac$testcol))
      cfuLclds <- as.character(ghi %>% arrange(Group2) %>% dplyr::select(Letter) %>% unlist() %>% unname())
    } else {
      cfuLclds <- c("a","b")
    }
  } else {
    cfuLclds <- c(rep("a", numvals))
  }
  
  cat("absolute abundance stats")
  print(kruskal.test(bac$cfuA ~ bac$testcol))
  print(kruskal.test(bac$cfuL ~ bac$testcol))
  print(dim(bac))
  
  if(!is.null(yaxis)) {    
    if(plot_aceto == T) {
        bac3plot <- ggplot(bac, aes(x=testcol, y=log10(cfuA+1))) +
          geom_jitter(width = jitter_width, height = jitter_height, col = plot_color) +
          geom_crossbar(data=bac3 %>% filter(variable.x == "cfuA"), aes(x = testcol.x, y = log10(value.x+1), ymin = log10(value.x+1), ymax = log10(value.x+1)), size=.5, col=meancol, width = meanwidth) + 
          geom_text(data = bac3 %>% filter(variable.x == "cfuA"), aes(x = testcol.x, y = max(log10(bac$cfuA +1))*baradjust), label = c(cfuAclds), vjust = 1, inherit.aes = T)+
          theme_cowplot()+ 
          ylim(yaxis)
    } else {
        bac3plot <- ggplot(bac, aes(x=testcol, y=log10(cfuL+1))) +
          geom_jitter(width = jitter_width, height = jitter_height, col = plot_color) +
          geom_crossbar(data=bac3 %>% filter(variable.x == "cfuL"), aes(x = testcol.x, y = log10(value.x+1), ymin = log10(value.x+1), ymax = log10(value.x+1)), size=.5, col=meancol, width = meanwidth) + 
          geom_text(data = bac3 %>% filter(variable.x == "cfuL"), aes(x = testcol.x, y = max(log10(bac$cfuL +1))*baradjust), label = c(cfuLclds), vjust = 1, inherit.aes = T)+
          theme_cowplot()+ 
          ylim(yaxis)
    }
  } else {
      if(plot_aceto == T) {
        bac3plot <-ggplot(bac, aes(x=testcol, y=log10(cfuA+1))) +
            geom_jitter(width = jitter_width, height = jitter_height, col = plot_color) +
            geom_crossbar(data=bac3 %>% filter(variable.x == "cfuA"), aes(x = testcol.x, y = log10(value.x+1), ymin = log10(value.x+1), ymax = log10(value.x+1)), size=.5, col=meancol, width = meanwidth) + 
            geom_text(data = bac3 %>% filter(variable.x == "cfuA"), aes(x = testcol.x, y = max(log10(bac$cfuA +1))*baradjust), label = c(cfuAclds), vjust = 1, inherit.aes = T)+
            theme_cowplot()
      } else {
          bac3plot <- ggplot(bac, aes(x=testcol, y=log10(cfuL+1))) +
            geom_jitter(width = jitter_width, height = jitter_height, col = plot_color) +
            geom_crossbar(data=bac3 %>% filter(variable.x == "cfuL"), aes(x = testcol.x, y = log10(value.x+1), ymin = log10(value.x+1), ymax = log10(value.x+1)), size=.5, col=meancol, width = meanwidth) + 
            geom_text(data = bac3 %>% filter(variable.x == "cfuL"), aes(x = testcol.x, y = max(log10(bac$cfuL +1))*baradjust), label = c(cfuLclds), vjust = 1, inherit.aes = T)+
            theme_cowplot()
      }
    }
  
  
  if(!is.null(yaxis)) {    
    if(plot_aceto == T) {
        bac5plot <- ggplot(bac, aes(x=testcol, y=cfuA)) +
          geom_jitter(width = jitter_width, height=jitter_height, col = plot_color) +
          geom_crossbar(data=bac3 %>% filter(variable.x == "cfuA"), aes(x = testcol.x, y = value.x, ymin = value.x, ymax = value.x), size=.5, col=meancol, width = meanwidth) + 
          geom_text(data = bac3 %>% filter(variable.x == "cfuA"), aes(x = testcol.x, y = max(bac$cfuA)*baradjust), label = c(cfuAclds), vjust = 1, inherit.aes = T)+
          theme_cowplot()+ 
          ylim(yaxis)
    } else {
      bac5plot <- ggplot(bac, aes(x=testcol, y=cfuL)) +
          geom_jitter(width = jitter_width, height=jitter_height, col = plot_color) +
          geom_crossbar(data=bac3 %>% filter(variable.x == "cfuL"), aes(x = testcol.x, y = (value.x), ymin = (value.x), ymax = (value.x)), size=.5, col=meancol, width = meanwidth) + 
        geom_text(data = bac3 %>% filter(variable.x == "cfuL"), aes(x = testcol.x, y = max(bac$cfuL)*baradjust), label = c(cfuLclds), vjust = 1, inherit.aes = T)+
        theme_cowplot()+ 
          ylim(yaxis)
    }
  } else {
    if(plot_aceto == T) {
      bac5plot <- ggplot(bac, aes(x=testcol, y=cfuA)) +
          geom_jitter(width = jitter_width, height=jitter_height, col = plot_color) +
          geom_crossbar(data=bac3 %>% filter(variable.x == "cfuA"), aes(x = testcol.x, y = (value.x), ymin = (value.x), ymax = (value.x)), size=.5, col=meancol, width = meanwidth) + 
        geom_text(data = bac3 %>% filter(variable.x == "cfuA"), aes(x = testcol.x, y = max(bac$cfuA)*baradjust), label = c(cfuAclds), vjust = 1, inherit.aes = T)+
        theme_cowplot()
    } else {
      bac5plot <- ggplot(bac, aes(x=testcol, y=cfuL)) +
          geom_jitter(width = jitter_width, height=jitter_height, col = plot_color) +
          geom_crossbar(data=bac3 %>% filter(variable.x == "cfuL"), aes(x = testcol.x, y = (value.x), ymin = (value.x), ymax = (value.x)), size=.5, col=meancol, width = meanwidth) + 
        geom_text(data = bac3 %>% filter(variable.x == "cfuL"), aes(x = testcol.x, y = max(bac$cfuL)*baradjust), label = c(cfuLclds), vjust = 1, inherit.aes = T)+
        theme_cowplot()
    }
  }
  
  # old
  bac4 <- bac %>% 
    mutate(total = cfuA + cfuL) %>%
    group_by(repexp, exp,testcol) %>%
    dplyr::summarize(cfuA = median(cfuA/total, na.rm=T), cfuL = median(cfuL/total, na.rm=T),  total = median(total/total, na.rm=T)) %>%
    ungroup() %>%
    group_by(exp,testcol) %>%
    dplyr::summarize(cfuA = mean(cfuA, na.rm=T), cfuL = mean(cfuL, na.rm=T), total = mean(total, na.rm=T)) %>%
    ungroup() %>%
    group_by(testcol) %>%
    dplyr::summarize(cfuA = mean(cfuA, na.rm=T), cfuL = mean(cfuL, na.rm=T), total = mean(total, na.rm=T)) %>%
    ungroup() %>%
    rowwise() %>%
    mutate(totalcol = cfuA+cfuL) %>%
    mutate(aperc = cfuA/totalcol, lperc = cfuL/totalcol) %>%
    ungroup() %>%
    reshape2::melt(measure.vars = c("aperc","lperc"))
  
  bac4A <- bac %>% 
    mutate(total = cfuA + cfuL) %>%
    mutate(percA = cfuA/total, percL = cfuL/total)
  
  bac4$addcol <- bac4$value
  if (aceto_bottom == T) {
    bac4$variable = factor(bac4$variable, levels = c("lperc","aperc"))
    for (i in (numvals+1):(numvals*2)) {bac4$addcol[i] = bac4$value[i-numvals]}
  } else {
    bac4$variable = factor(bac4$variable, levels = c("aperc","lperc"))
    for (i in 1:numvals) {bac4$addcol[i] = bac4$value[i+numvals]}
  }
  
  length(table(bac$testcol))
  if(kruskal.test(bac$aperc ~ bac$testcol)$p.value < 0.05) {
    if (length(table(bac$testcol))>2) {
      efg <- dunn.test(bac$aperc,bac$testcol, method = "bh", table = F, kw = F)
      ghi <- cldList(comparison = efg$comparisons, p.value = efg$P.adjusted, threshold = 0.05) 
      ghi$Group2 <- factor(ghi$Group, levels = levels(bac$testcol))
      apercclds <- as.character(ghi %>% arrange(Group2) %>% dplyr::select(Letter) %>% unlist() %>% unname())
      print(kruskal.test(bac$aperc ~ bac$testcol))
    } else {
      apercclds <- c("a","b")
    }
  } else {
    apercclds <- c(rep("a",numvals))
  }
  
  if(kruskal.test(bac$lperc ~ bac$testcol)$p.value < 0.05) {
    if (length(table(bac$testcol))>2) {
      efg <- dunn.test(bac$lperc,bac$testcol, method = "bh", table = F, kw = F)
      ghi <- cldList(comparison = efg$comparisons, p.value = efg$P.adjusted, threshold = 0.05) 
      ghi$Group2 <- factor(ghi$Group, levels = levels(bac$testcol))
      lpercclds <- as.character(ghi %>% arrange(Group2) %>% dplyr::select(Letter) %>% unlist() %>% unname())
      print(kruskal.test(bac$lperc ~ bac$testcol))
    } else {
      lpercclds <- c("a","b")
    }
  } else {
    lpercclds <-c(rep("a",numvals))
  }

  cat("relative abundance stats")
  print(kruskal.test(bac$aperc ~ bac$testcol))
  print(kruskal.test(bac$lperc ~ bac$testcol))
  
    if(aceto_bottom == T) {
      # if(logplot == T) {
      #   bac4plot <- ggplot(bac4A, aes(x=testcol, y=log10(cfuA+1))) +
      #     geom_jitter(width = jitter_width, height = jitter_height, col = plot_color) +
      #     geom_crossbar(data=bac4 %>% filter(variable == "aperc"), aes(x = testcol, y = log10(value+1), ymin = log10(value+1), ymax = log10(value+1)), size=.5, col="black", width = .25) + 
      #     theme_cowplot()
      # } else {
        bac4plot <- ggplot(bac4A, aes(x=testcol, y=aperc)) +
          geom_jitter(width = jitter_width, height=jitter_height, col = plot_color) +
          geom_crossbar(data=bac4 %>% filter(variable == "aperc"), aes(x = testcol, y = value, ymin = value, ymax = value), size=.5, col=meancol, width = meanwidth) + 
          geom_text(data = bac4 %>% filter(variable == "aperc"), aes(x = testcol, y = max(bac$aperc)*baradjust), label = c(apercclds), vjust = 1, inherit.aes = T)+
          theme_cowplot()
    #  }
    } else {
      # if(logplot == T) {
      #   bac4plot <- ggplot(bac4A, aes(x=testcol, y=log10(cfuL+1))) +
      #     geom_jitter(width = jitter_width, height = jitter_height, col = plot_color) +
      #     geom_crossbar(data=bac4 %>% filter(variable == "lperc"), aes(x = testcol, y = log10(value+1), ymin = log10(value+1), ymax = log10(value+1)), size=.5, col="black", width = .25) + 
      #     theme_cowplot()
      # } else {
        bac4plot <- ggplot(bac4A, aes(x=testcol, y=lperc)) +
          geom_jitter(width = jitter_width, height=jitter_height, col = plot_color) +
          geom_crossbar(data=bac4 %>% filter(variable == "lperc"), aes(x = testcol, y = value, ymin = value, ymax = value), size=.5, col=meancol, width = meanwidth) + 
          geom_text(data = bac4 %>% filter(variable == "lperc"), aes(x = testcol, y = max(bac$lperc)*baradjust), label = c(lpercclds), vjust = 1, inherit.aes = T)+
          theme_cowplot()
    #  }
    }

  ifelse(plotout == 1, return(bac3plot),
         ifelse(plotout == 2, return(bac5plot),
                ifelse(plotout == 3, return(bac4plot),
                       ifelse(plotout == 4, return(grid.arrange(bac3plot, bac5plot, heights = c(1,1))),
                              ifelse(plotout == 5, return(grid.arrange(bac3plot, bac4plot, heights = c(1,1))),
                                     ifelse(plotout == 6, return(grid.arrange(bac4plot, bac5plot, heights = c(1,1))),
                                            ifelse(plotout == 7, return(grid.arrange(bac3plot, bac5plot, bac4plot, heights = c(1,1,1))))))))))
}


# file_path = "figure2g"
# mapper_file = "cfu_mapper_gh.tsv"
# formula = "~ Geography * Sex"
# 
# permutations = 1000
# PA Control M 1 NA 2286 4 Flies 1 NA_1_ 
# PA Control M 1 NA 2286 4 Flies 1 NA_1_   
# dim(flies_prerarefied_dm)

# Geography * Sex * starved + exp/PLATE/code2
# 
# file_path = "figureoptimize1data"
# mapper_file = "cfu_mapper_optimize1.tsv"
# formula = "~ Treatment"
# permutations = 1000

bc_analysis <- function(file_path, mapper_file, formula, permutations = 999) {
  
  try(rm(flies_prerarefied),T)
  try(rm(flies_rarefied),T)
  
  
  ## prerarefied
  flies_prerarefied <- read.table(paste0(file_path,'/bray_curtis_dm_',file_path,'/distance-matrix.tsv'), header=T, sep="\t") %>% mutate(X=as.character(X))
  flies_prerarefied_dm <- as.dist(flies_prerarefied[,2:dim(flies_prerarefied)[2]])
  flies_unifrac_table <- flies_prerarefied %>% 
    dplyr::select(X) %>% 
    left_join(read.table(mapper_file,comment.char = "", header=T, fill=T, sep="\t"), by=c("X"="X.SampleID")) %>%
    droplevels()
  
  set.seed(42)
  prerarefied_permanova <- vegan::adonis2(as.formula(paste0("flies_prerarefied_dm ",formula)), flies_unifrac_table, permutations=permutations)
  
#  prerarefied_permanova <- vegan::adonis2("flies_prerarefied_dm ~ Geography * Sex", flies_unifrac_table, permutations=permutations)
  
  # sum(table(flies_unifrac_table$Geography))
  # sum(table(flies_unifrac_table$Sex))
  # dim(flies_prerarefied_dm)
  # dim(flies_unifrac_table)
  prerarefied_permanova$Df
  w_df <- data.frame(variables = rownames(prerarefied_permanova), Df = round(prerarefied_permanova$Df, 2), SS = round(prerarefied_permanova$SumOfSqs, 2), R2 = round(prerarefied_permanova$R2, 2), Fval = round(prerarefied_permanova$`F`, 2), p = round(prerarefied_permanova$`Pr(>F)`, 2))
  wtable <- tableGrob(data.frame(variables = rownames(prerarefied_permanova), Df = round(prerarefied_permanova$Df, 2), SS = round(prerarefied_permanova$SumOfSqs, 2), R2 = round(prerarefied_permanova$R2, 2), Fval = round(prerarefied_permanova$`F`, 2), p = round(prerarefied_permanova$`Pr(>F)`, 2)),theme=ttheme_minimal())
  
  ## rarefied
  flies_rarefied <- read.table(paste0(file_path,'/core-metrics-results-',file_path,'/bray_curtis_distance_matrix/distance-matrix.tsv'), header=T, sep="\t") %>% mutate(X=as.character(X))
  flies_rarefied_dm <- as.dist(flies_rarefied[,2:dim(flies_rarefied)[2]])
  flies_unifrac_table <- flies_rarefied %>% 
    dplyr::select(X) %>% 
    left_join(read.table(mapper_file,comment.char = "", header=T, fill=T, sep="\t"), by=c("X"="X.SampleID")) %>%
    droplevels()
  
  set.seed(42)
  rarefied_permanova <- adonis2(as.formula(paste0("flies_rarefied_dm ",formula)), flies_unifrac_table, permutations=permutations)
  uw_df <- data.frame(variables = rownames(rarefied_permanova), Df = round(rarefied_permanova$Df, 2), SS = round(rarefied_permanova$SumOfSqs, 2), R2 = round(rarefied_permanova$R2, 2), Fval = round(rarefied_permanova$`F`, 2), p = round(rarefied_permanova$`Pr(>F)`, 2))
  uwtable <- tableGrob(data.frame(variables = rownames(rarefied_permanova), Df = round(rarefied_permanova$Df, 2), SS = round(rarefied_permanova$SumOfSqs, 2), R2 = round(rarefied_permanova$R2, 2), Fval = round(rarefied_permanova$`F`, 2), p = round(rarefied_permanova$`Pr(>F)`, 2)),theme=ttheme_minimal())
  
  print(prerarefied_permanova)
  print(rarefied_permanova)
  
  #    return(list(grid.arrange(wtable , uwtable), w_df, uw_df))
  
}


## oroginal
# make_chart_stats_96 <- function(bac, column, plotout, yaxis = NULL, factorspecs = NULL, aceto_bottom =T, label_adjust = 10000, logplot = logplot) {
#   #bac$testcol <- bac[,column]
#   bac$testcol <- unlist((bac[,column]))
#   bac <- bac %>%
#     mutate(aperc = cfuA/(cfuA+cfuL), lperc = cfuL/(cfuA+cfuL)) %>%
#     droplevels()
#   if(!is.null(factorspecs)) {
#     bac$testcol = factor(bac$testcol, levels = factorspecs)
#   } else {
#     bac$testcol = factor(bac$testcol)
#     print(levels(factor(bac$testcol)))
#   }
#   if(logplot == F) {
#     bac3 <- bac %>% 
#       mutate(total = cfuA + cfuL) %>%
#       group_by(repexp, exp,testcol) %>%
#       dplyr::summarize(cfuA = median(cfuA, na.rm=T), cfuL = median(cfuL, na.rm=T), total = median(total, na.rm=T)) %>%
#       ungroup() %>%
#       group_by(exp,testcol) %>%
#       dplyr::summarize(cfuA = mean(cfuA, na.rm=T), cfuL = mean(cfuL, na.rm=T), total = mean(total, na.rm=T)) %>%
#       ungroup() %>%
#       group_by(testcol) %>%
#       dplyr::summarize(cfuA = mean(cfuA, na.rm=T), cfuL = mean(cfuL, na.rm=T), total = mean(total, na.rm=T)) %>%
#       ungroup() %>%
#       reshape2::melt(measure.vars = c("cfuA","cfuL")) %>%
#       mutate(tempvar = paste0(testcol,"_",variable)) %>%
#       inner_join(bac %>% 
#                    mutate(total = cfuA + cfuL) %>%
#                    group_by(repexp, exp,testcol) %>%
#                    dplyr::summarize(cfuA = median(cfuA, na.rm=T), cfuL = median(cfuL, na.rm=T), total = median(total, na.rm=T)) %>%
#                    ungroup() %>%
#                    group_by(exp,testcol) %>%
#                    dplyr::summarize(cfuA = mean(cfuA, na.rm=T), cfuL = mean(cfuL, na.rm=T), total = mean(total, na.rm=T)) %>%
#                    ungroup() %>%
#                    group_by(testcol) %>%
#                    dplyr::summarize(semcfuA = sd(cfuA, na.rm=T)/sqrt(dplyr::n()), semcfuL = sd(cfuL, na.rm=T)/sqrt(dplyr::n()), semtotal = sd(total, na.rm=T)/sqrt(dplyr::n())) %>%
#                    ungroup() %>%
#                    reshape2::melt(measure.vars = c("semcfuA","semcfuL"))  %>%
#                    mutate(mergecol = gsub(pattern = "sem", replacement = "", x = variable)) %>%
#                    mutate(tempvar = paste0(testcol,"_",mergecol)), by = "tempvar") 
#   } else {
#     bac3 <- bac %>% 
#       mutate(total = cfuA + cfuL) %>%
#       group_by(repexp, exp,testcol) %>%
#       dplyr::summarize(cfuA = mean(log10(cfuA+1), na.rm=T), cfuL = mean(log10(cfuL+1), na.rm=T), total = mean(log10(total+1), na.rm=T)) %>%
#       ungroup() %>%
#       group_by(exp,testcol) %>%
#       dplyr::summarize(cfuA = mean(cfuA, na.rm=T), cfuL = mean(cfuL, na.rm=T), total = mean(total, na.rm=T)) %>%
#       ungroup() %>%
#       group_by(testcol) %>%
#       dplyr::summarize(cfuA = mean(cfuA, na.rm=T), cfuL = mean(cfuL, na.rm=T), total = mean(total, na.rm=T)) %>%
#       ungroup() %>%
#       reshape2::melt(measure.vars = c("cfuA","cfuL")) %>%
#       mutate(tempvar = paste0(testcol,"_",variable)) %>%
#       inner_join(bac %>% 
#                    mutate(total = cfuA + cfuL) %>%
#                    group_by(repexp, exp,testcol) %>%
#                    dplyr::summarize(cfuA = mean(log10(cfuA+1), na.rm=T), cfuL = mean(log10(cfuL+1), na.rm=T), total = mean(log10(total+1), na.rm=T)) %>%
#                    ungroup() %>%
#                    group_by(exp,testcol) %>%
#                    dplyr::summarize(cfuA = mean(cfuA, na.rm=T), cfuL = mean(cfuL, na.rm=T), total = mean(total, na.rm=T)) %>%
#                    ungroup() %>%
#                    group_by(testcol) %>%
#                    dplyr::summarize(semcfuA = sd(cfuA, na.rm=T)/sqrt(dplyr::n()), semcfuL = sd(cfuL, na.rm=T)/sqrt(dplyr::n()), semtotal = sd(total, na.rm=T)/sqrt(dplyr::n())) %>%
#                    ungroup() %>%
#                    reshape2::melt(measure.vars = c("semcfuA","semcfuL"))  %>%
#                    mutate(mergecol = gsub(pattern = "sem", replacement = "", x = variable)) %>%
#                    mutate(tempvar = paste0(testcol,"_",mergecol)), by = "tempvar") 
#   }
#   
#   bac3$addcol <- bac3$value.x
#   numvals <- length(table(list(bac$testcol %>% droplevels)))
#   
#   if (aceto_bottom == T) {
#     bac3$variable.x = factor(bac3$variable.x, levels = c("cfuL","cfuA"))
#     for (i in (numvals+1):(numvals*2)) {bac3$addcol[i] = bac3$total[i]}
#   } else {
#     bac3$variable.x = factor(bac3$variable.x, levels = c("cfuA","cfuL"))
#     for (i in 1:numvals) {bac3$addcol[i] = bac3$total[i]}
#   }
#   
#   
#   if(kruskal.test(bac$cfuA ~ bac$testcol)$p.value < 0.05) {
#     efg <- dunn.test(bac$cfuA,bac$testcol, method = "bh", table = F, kw = F)
#     ghi <- cldList(comparison = efg$comparisons, p.value = efg$P.adjusted, threshold = 0.05)
#     ghi$Group2 <- factor(ghi$Group, levels = levels(bac$testcol))
#     cfuAclds <- as.character(ghi %>% arrange(Group2) %>% dplyr::select(Letter) %>% unlist() %>% unname())
#   } else {
#     cfuAclds <- c(rep("a", numvals))
#   }
#   
#   if(kruskal.test(bac$cfuL ~ bac$testcol)$p.value < 0.05) {
#     efg <- dunn.test(bac$cfuL,bac$testcol, method = "bh", table = F, kw = F)
#     ghi <- cldList(comparison = efg$comparisons, p.value = efg$P.adjusted, threshold = 0.05)
#     ghi$Group2 <- factor(ghi$Group, levels = levels(bac$testcol))
#     cfuLclds <- as.character(ghi %>% arrange(Group2) %>% dplyr::select(Letter) %>% unlist() %>% unname())
#   } else {
#     cfuLclds <- c(rep("a", numvals))
#   }
#   
#   if(!is.null(yaxis)) {    
#     if(aceto_bottom != T) {
#       bac3plot <- ggplot(bac3, aes(x=testcol.x, y=value.x, fill = variable.x)) +
#         geom_bar(stat = "identity") + 
#         geom_errorbar(aes(ymin =addcol-value.y, ymax = addcol+value.y), width = 0.5) +
#         scale_fill_manual(values = c("red","blue"))+
#         theme(legend.position = "bottom") + 
#         ylim(yaxis)+
#         geom_text(x = c(1:numvals,1:numvals), aes(y = addcol-label_adjust-value.y, label = c(rep("",numvals),cfuAclds)))+
#         geom_text(x = c(1:numvals,1:numvals), aes(y = addcol+label_adjust+value.y, label = c(cfuLclds,rep("",numvals))))+
#         theme_cowplot()
#     } else {
#       bac3plot <- ggplot(bac3, aes(x=testcol.x, y=value.x, fill = variable.x)) +
#         geom_bar(stat = "identity") + 
#         geom_errorbar(aes(ymin =addcol-value.y, ymax = addcol+value.y), width = 0.5) +
#         scale_fill_manual(values = c("blue","red"))+
#         theme(legend.position = "bottom") +  
#         ylim(yaxis)+
#         geom_text(x = c(1:numvals,1:numvals), aes(y = addcol-label_adjust-value.y, label = c(cfuAclds,rep("",numvals))))+
#         geom_text(x = c(1:numvals,1:numvals), aes(y = addcol+label_adjust+value.y, label = c(rep("",numvals),cfuLclds)))+
#         theme_cowplot()
#     }
#   } else {
#     if(aceto_bottom != T) {
#       bac3plot <- ggplot(bac3, aes(x=testcol.x, y=value.x, fill = variable.x)) +
#         geom_bar(stat = "identity") + 
#         geom_errorbar(aes(ymin =addcol-value.y, ymax = addcol+value.y), width = 0.5) +
#         scale_fill_manual(values = c("red","blue"))+
#         theme(legend.position = "bottom") + 
#         geom_text(x = c(1:numvals,1:numvals), aes(y = addcol-label_adjust-value.y, label = c(rep("",numvals),cfuAclds)))+
#         geom_text(x = c(1:numvals,1:numvals), aes(y = addcol+label_adjust+value.y, label = c(cfuLclds,rep("",numvals))))+
#         theme_cowplot()
#     } else {
#       bac3plot <- ggplot(bac3, aes(x=testcol.x, y=value.x, fill = variable.x)) +
#         geom_bar(stat = "identity") + 
#         geom_errorbar(aes(ymin =addcol-value.y, ymax = addcol+value.y), width = 0.5) +
#         scale_fill_manual(values = c("blue","red"))+
#         theme(legend.position = "bottom") +  
#         geom_text(x = c(1:numvals,1:numvals), aes(y = addcol-label_adjust-value.y, label = c(cfuAclds,rep("",numvals))))+
#         geom_text(x = c(1:numvals,1:numvals), aes(y = addcol+label_adjust+value.y, label = c(rep("",numvals),cfuLclds)))+
#         theme_cowplot()
#     }
#   }
#   
#   # old
#   bac4 <- bac %>%
#     mutate(total = cfuA + cfuL) %>%
#     group_by(repexp, exp,testcol) %>%
#     dplyr::summarize(cfuA = median(cfuA, na.rm=T), cfuL = median(cfuL, na.rm=T),  total = median(total, na.rm=T)) %>%
#     ungroup() %>%
#     group_by(exp,testcol) %>%
#     dplyr::summarize(cfuA = mean(cfuA, na.rm=T), cfuL = mean(cfuL, na.rm=T), total = mean(total, na.rm=T)) %>%
#     ungroup() %>%
#     group_by(testcol) %>%
#     dplyr::summarize(cfuA = mean(cfuA, na.rm=T), cfuL = mean(cfuL, na.rm=T), total = mean(total, na.rm=T)) %>%
#     ungroup() %>%
#     rowwise() %>%
#     mutate(totalcol = cfuA+cfuL) %>%
#     mutate(aperc = cfuA/totalcol, lperc = cfuL/totalcol) %>%
#     ungroup() %>%
#     reshape2::melt(measure.vars = c("aperc","lperc"))
#   
#   bac4 <- bac %>% 
#     mutate(total = cfuA + cfuL) %>%
#     group_by(repexp, exp,testcol) %>%
#     dplyr::summarize(cfuA = median(cfuA/total, na.rm=T), cfuL = median(cfuL/total, na.rm=T),  total = median(total/total, na.rm=T)) %>%
#     ungroup() %>%
#     group_by(exp,testcol) %>%
#     dplyr::summarize(cfuA = mean(cfuA, na.rm=T), cfuL = mean(cfuL, na.rm=T), total = mean(total, na.rm=T)) %>%
#     ungroup() %>%
#     group_by(testcol) %>%
#     dplyr::summarize(cfuA = mean(cfuA, na.rm=T), cfuL = mean(cfuL, na.rm=T), total = mean(total, na.rm=T)) %>%
#     ungroup() %>%
#     rowwise() %>%
#     mutate(totalcol = cfuA+cfuL) %>%
#     mutate(aperc = cfuA/totalcol, lperc = cfuL/totalcol) %>%
#     ungroup() %>%
#     reshape2::melt(measure.vars = c("aperc","lperc"))
#   
#   bac4$addcol <- bac4$value
#   if (aceto_bottom == T) {
#     bac4$variable = factor(bac4$variable, levels = c("lperc","aperc"))
#     for (i in (numvals+1):(numvals*2)) {bac4$addcol[i] = bac4$value[i-numvals]}
#   } else {
#     bac4$variable = factor(bac4$variable, levels = c("aperc","lperc"))
#     for (i in 1:numvals) {bac4$addcol[i] = bac4$value[i+numvals]}
#   }
#   
#   length(table(bac$testcol))
#   if(kruskal.test(bac$aperc ~ bac$testcol)$p.value < 0.05) {
#     efg <- dunn.test(bac$aperc,bac$testcol, method = "bh", table = F, kw = F)
#     ghi <- cldList(comparison = efg$comparisons, p.value = efg$P.adjusted, threshold = 0.05) 
#     ghi$Group2 <- factor(ghi$Group, levels = levels(bac$testcol))
#     apercclds <- as.character(ghi %>% arrange(Group2) %>% dplyr::select(Letter) %>% unlist() %>% unname())
#   } else {
#     apercclds <- c(rep("a",numvals))
#   }
#   
#   if(kruskal.test(bac$lperc ~ bac$testcol)$p.value < 0.05) {
#     efg <- dunn.test(bac$lperc,bac$testcol, method = "bh", table = F, kw = F)
#     ghi <- cldList(comparison = efg$comparisons, p.value = efg$P.adjusted, threshold = 0.05) 
#     ghi$Group2 <- factor(ghi$Group, levels = levels(bac$testcol))
#     lpercclds <- as.character(ghi %>% arrange(Group2) %>% dplyr::select(Letter) %>% unlist() %>% unname())
#   } else {
#     lpercclds <-c(rep("a",numvals))
#   }
#   
#   if(aceto_bottom != T) {
#     
#     bac4plot <- ggplot(bac4, aes(x=testcol, y=value, fill = variable)) +
#       geom_bar(stat = "identity", position = "stack") + theme(legend.position = "bottom") +
#       scale_fill_manual(values = c("red","blue"))+
#       geom_text(x = c(1:numvals,1:numvals), aes(y = addcol-0.05, label = c(rep("",numvals),cfuAclds)))+
#       geom_text(x = c(1:numvals,1:numvals), aes(y = addcol+0.05, label = c(cfuLclds,rep("",numvals))))+
#       theme_cowplot()
#     
#   } else {
#     bac4plot <- ggplot(bac4, aes(x=testcol, y=value, fill = variable)) +
#       geom_bar(stat = "identity", position = "stack") + theme(legend.position = "bottom") +
#       scale_fill_manual(values = c("blue","red"))+
#       geom_text(x = c(1:numvals,1:numvals), aes(y = addcol-0.05, label = c(cfuAclds,rep("",numvals))))+
#       geom_text(x = c(1:numvals,1:numvals), aes(y = addcol+0.05, label = c(rep("",numvals),cfuLclds)))+
#       theme_cowplot()
#   }
#   
#   ifelse(plotout == 1, return(bac3plot), ifelse(plotout == 2, return(bac4plot), ifelse(plotout==3, return(grid.arrange(bac3plot, bac4plot, heights = c(1,1))))))
# }


## revised 02/28/23 - blocks printing the stats becuase the stats get labelled wrong.
make_chart_stats_96_bigmean <- function(bac, column, plotout, yaxis = NULL, factorspecs = NULL, aceto_bottom =T, label_adjust = 10000, logplot = logplot) {
  #bac$testcol <- bac[,column]
  bac$testcol <- unlist((bac[,column]))
  bac <- bac %>%
    mutate(aperc = cfuA/(cfuA+cfuL), lperc = cfuL/(cfuA+cfuL)) %>%
    droplevels()
  if(!is.null(factorspecs)) {
    bac$testcol = factor(bac$testcol, levels = factorspecs)
  } else {
    bac$testcol = factor(bac$testcol)
    print(levels(factor(bac$testcol)))
  }
  if(logplot == F) {
    bac3 <- bac %>% 
      mutate(total = cfuA + cfuL) %>%
      group_by(testcol) %>%
      dplyr::summarize(cfuA = mean(cfuA, na.rm=T), cfuL = mean(cfuL, na.rm=T), total = mean(total, na.rm=T)) %>%
      ungroup() %>%
      reshape2::melt(measure.vars = c("cfuA","cfuL")) %>%
      mutate(tempvar = paste0(testcol,"_",variable)) %>%
      inner_join(bac %>% 
                   mutate(total = cfuA + cfuL) %>%
                   group_by(testcol) %>%
                   dplyr::summarize(semcfuA = sd(cfuA, na.rm=T)/sqrt(dplyr::n()), semcfuL = sd(cfuL, na.rm=T)/sqrt(dplyr::n()), semtotal = sd(total, na.rm=T)/sqrt(dplyr::n())) %>%
                   ungroup() %>%
                   reshape2::melt(measure.vars = c("semcfuA","semcfuL"))  %>%
                   mutate(mergecol = gsub(pattern = "sem", replacement = "", x = variable)) %>%
                   mutate(tempvar = paste0(testcol,"_",mergecol)), by = "tempvar") 
  } else {
    bac3 <- bac %>% 
      mutate(total = cfuA + cfuL) %>%
      group_by(repexp, exp,testcol) %>%
      dplyr::summarize(cfuA = mean(log10(cfuA+1), na.rm=T), cfuL = mean(log10(cfuL+1), na.rm=T), total = mean(log10(total+1), na.rm=T)) %>%
      ungroup() %>%
      reshape2::melt(measure.vars = c("cfuA","cfuL")) %>%
      mutate(tempvar = paste0(testcol,"_",variable)) %>%
      inner_join(bac %>% 
                   mutate(total = cfuA + cfuL) %>%
                   group_by(testcol) %>%
                   dplyr::summarize(semcfuA = sd(log10(cfuA+1), na.rm=T)/sqrt(dplyr::n()), semcfuL = sd(log10(cfuL+1), na.rm=T)/sqrt(dplyr::n()), semtotal = sd(log10(total+1), na.rm=T)/sqrt(dplyr::n())) %>%
                   ungroup() %>%
                   reshape2::melt(measure.vars = c("semcfuA","semcfuL"))  %>%
                   mutate(mergecol = gsub(pattern = "sem", replacement = "", x = variable)) %>%
                   mutate(tempvar = paste0(testcol,"_",mergecol)), by = "tempvar") 
  }
  
  print(paste0("N:",dim(bac)[1]))
  bac3$addcol <- bac3$value.x
  numvals <- length(table(list(bac$testcol)))
  
  if (aceto_bottom == T) {
    bac3$variable.x = factor(bac3$variable.x, levels = c("cfuL","cfuA"))
    for (i in (numvals+1):(numvals*2)) {bac3$addcol[i] = bac3$total[i]}
  } else {
    bac3$variable.x = factor(bac3$variable.x, levels = c("cfuA","cfuL"))
    for (i in 1:numvals) {bac3$addcol[i] = bac3$total[i]}
  }
  
  cat("Stats on AAB abundances")
  print(kruskal.test(bac$cfuA ~ bac$testcol))
  if(kruskal.test(bac$cfuA ~ bac$testcol)$p.value < 0.05) {
    efg <- dunn.test(bac$cfuA,bac$testcol, method = "bh", table = F, kw = F)
    ghi <- cldList(comparison = efg$comparisons, p.value = efg$P.adjusted, threshold = 0.05)
    ghi$Group2 <- factor(ghi$Group, levels = levels(bac$testcol))
    cfuAclds <- as.character(ghi %>% arrange(Group2) %>% dplyr::select(Letter) %>% unlist() %>% unname())
    print("")
    cat("AAB Plot_values:")
    cat(ghi$Group)
    cat(" ")
    cat(cfuAclds)
  } else {
    cfuAclds <- c(rep("a", numvals))
    print("")
    cat("AAB Plot_values:")
    cat(cfuAclds)
  }
  
  cat("Stats on LAB abundances")
  print(kruskal.test(bac$cfuL ~ bac$testcol))
  
  if(kruskal.test(bac$cfuL ~ bac$testcol)$p.value < 0.05) {
    efg <- dunn.test(bac$cfuL,bac$testcol, method = "bh", table = F, kw = F)
    ghi <- cldList(comparison = efg$comparisons, p.value = efg$P.adjusted, threshold = 0.05)
    ghi$Group2 <- factor(ghi$Group, levels = levels(bac$testcol))
    cfuLclds <- as.character(ghi %>% arrange(Group2) %>% dplyr::select(Letter) %>% unlist() %>% unname())
    print("")
    cat("LAB Plot_values:")
    cat(ghi$Group)
    cat(" ")
    cat(cfuLclds)
  } else {
    cfuLclds <- c(rep("a", numvals))
    print("")
    cat("LAB Plot_values:")
    cat(cfuLclds)
  }
  
  if(!is.null(yaxis)) {    
    if(aceto_bottom != T) {
      bac3plot <- ggplot(bac3, aes(x=testcol.x, y=value.x, fill = variable.x)) +
        geom_bar(stat = "identity") + 
        geom_errorbar(aes(ymin =addcol-value.y, ymax = addcol+value.y), width = 0.5) +
        scale_fill_manual(values = c("red","blue"))+
        theme(legend.position = "bottom") + 
        ylim(yaxis)+
        #     geom_text(x = c(1:numvals,1:numvals), aes(y = addcol-label_adjust-value.y, label = c(rep("",numvals),cfuAclds)))+
        #    geom_text(x = c(1:numvals,1:numvals), aes(y = addcol+label_adjust+value.y, label = c(cfuLclds,rep("",numvals))))+
        theme_cowplot()
    } else {
      bac3plot <- ggplot(bac3, aes(x=testcol.x, y=value.x, fill = variable.x)) +
        geom_bar(stat = "identity") + 
        geom_errorbar(aes(ymin =addcol-value.y, ymax = addcol+value.y), width = 0.5) +
        scale_fill_manual(values = c("blue","red"))+
        theme(legend.position = "bottom") +  
        ylim(yaxis)+
        #   geom_text(x = c(1:numvals,1:numvals), aes(y = addcol-label_adjust-value.y, label = c(cfuAclds,rep("",numvals))))+
        #  geom_text(x = c(1:numvals,1:numvals), aes(y = addcol+label_adjust+value.y, label = c(rep("",numvals),cfuLclds)))+
        theme_cowplot()
    }
  } else {
    if(aceto_bottom != T) {
      bac3plot <- ggplot(bac3, aes(x=testcol.x, y=value.x, fill = variable.x)) +
        geom_bar(stat = "identity") + 
        geom_errorbar(aes(ymin =addcol-value.y, ymax = addcol+value.y), width = 0.5) +
        scale_fill_manual(values = c("red","blue"))+
        theme(legend.position = "bottom") + 
        #     geom_text(x = c(1:numvals,1:numvals), aes(y = addcol-label_adjust-value.y, label = c(rep("",numvals),cfuAclds)))+
        #    geom_text(x = c(1:numvals,1:numvals), aes(y = addcol+label_adjust+value.y, label = c(cfuLclds,rep("",numvals))))+
        theme_cowplot()
    } else {
      bac3plot <- ggplot(bac3, aes(x=testcol.x, y=value.x, fill = variable.x)) +
        geom_bar(stat = "identity") + 
        geom_errorbar(aes(ymin =addcol-value.y, ymax = addcol+value.y), width = 0.5) +
        scale_fill_manual(values = c("blue","red"))+
        theme(legend.position = "bottom") +  
        #   geom_text(x = c(1:numvals,1:numvals), aes(y = addcol-label_adjust-value.y, label = c(cfuAclds,rep("",numvals))))+
        #  geom_text(x = c(1:numvals,1:numvals), aes(y = addcol+label_adjust+value.y, label = c(rep("",numvals),cfuLclds)))+
        theme_cowplot()
    }
  }
  
  # old
  bac4 <- bac %>% 
    mutate(total = cfuA + cfuL) %>%
    group_by(testcol) %>%
    dplyr::summarize(cfuA = mean(cfuA, na.rm=T), cfuL = mean(cfuL, na.rm=T), total = mean(total, na.rm=T)) %>%
    ungroup() %>%
    rowwise() %>%
    mutate(totalcol = cfuA+cfuL) %>%
    mutate(aperc = cfuA/totalcol, lperc = cfuL/totalcol) %>%
    ungroup() %>%
    reshape2::melt(measure.vars = c("aperc","lperc"))
  
  bac4$addcol <- bac4$value
  if (aceto_bottom == T) {
    bac4$variable = factor(bac4$variable, levels = c("lperc","aperc"))
    for (i in (numvals+1):(numvals*2)) {bac4$addcol[i] = bac4$value[i-numvals]}
  } else {
    bac4$variable = factor(bac4$variable, levels = c("aperc","lperc"))
    for (i in 1:numvals) {bac4$addcol[i] = bac4$value[i+numvals]}
  }
  
  
  cat("Stats on abundances")
  print(kruskal.test(bac$aperc ~ bac$testcol))
  if(kruskal.test(bac$aperc ~ bac$testcol)$p.value < 0.05) {
    efg <- dunn.test(bac$aperc,bac$testcol, method = "bh", table = F, kw = F)
    ghi <- cldList(comparison = efg$comparisons, p.value = efg$P.adjusted, threshold = 0.05) 
    ghi$Group2 <- factor(ghi$Group, levels = levels(bac$testcol))
    apercclds <- as.character(ghi %>% arrange(Group2) %>% dplyr::select(Letter) %>% unlist() %>% unname())
    print("")
    cat("AAB Plot_values:")
    cat(ghi$Group)
    cat(" ")
    cat(apercclds)
  } else {
    apercclds <- c(rep("a",numvals))
    print("")
    cat("AAB Plot_values:")
    cat(apercclds)
  }
  
  if(kruskal.test(bac$lperc ~ bac$testcol)$p.value < 0.05) {
    efg <- dunn.test(bac$lperc,bac$testcol, method = "bh", table = F, kw = F)
    ghi <- cldList(comparison = efg$comparisons, p.value = efg$P.adjusted, threshold = 0.05) 
    ghi$Group2 <- factor(ghi$Group, levels = levels(bac$testcol))
    lpercclds <- as.character(ghi %>% arrange(Group2) %>% dplyr::select(Letter) %>% unlist() %>% unname())
    print("")
    cat("LAB Plot_values: ")
    cat(ghi$Group)
    cat(" ")
    cat(lpercclds)
  } else {
    lpercclds <-c(rep("a",numvals))
    print("")
    cat("LAB Plot_values:")
    cat(lpercclds)
  }
  
  if(aceto_bottom != T) {
    
    bac4plot <- ggplot(bac4, aes(x=testcol, y=value, fill = variable)) +
      geom_bar(stat = "identity", position = "stack") + theme(legend.position = "bottom") +
      scale_fill_manual(values = c("red","blue"))+
      #      geom_text(x = c(1:numvals,1:numvals), aes(y = addcol-0.05, label = c(rep("",numvals),cfuAclds)))+
      #    geom_text(x = c(1:numvals,1:numvals), aes(y = addcol+0.05, label = c(cfuLclds,rep("",numvals))))+
      theme_cowplot()
    
  } else {
    bac4plot <- ggplot(bac4, aes(x=testcol, y=value, fill = variable)) +
      geom_bar(stat = "identity", position = "stack") + theme(legend.position = "bottom") +
      scale_fill_manual(values = c("blue","red"))+
      # geom_text(x = c(1:numvals,1:numvals), aes(y = addcol-0.05, label = c(cfuAclds,rep("",numvals))))+
      #  geom_text(x = c(1:numvals,1:numvals), aes(y = addcol+0.05, label = c(rep("",numvals),cfuLclds)))+
      theme_cowplot()
  }
  
  ifelse(plotout == 1, return(bac3plot), ifelse(plotout == 2, return(bac4plot), ifelse(plotout==3, return(grid.arrange(bac3plot, bac4plot, heights = c(1,1))))))
}


make_chart_stats_96 <- function(bac, column, plotout, yaxis = NULL, factorspecs = NULL, aceto_bottom =T, label_adjust = 10000, logplot = logplot) {
  #bac$testcol <- bac[,column]
  bac$testcol <- unlist((bac[,column]))
  bac <- bac %>%
    mutate(aperc = cfuA/(cfuA+cfuL), lperc = cfuL/(cfuA+cfuL)) %>%
    droplevels()
  if(!is.null(factorspecs)) {
    bac$testcol = factor(bac$testcol, levels = factorspecs)
  } else {
    bac$testcol = factor(bac$testcol)
    print(levels(factor(bac$testcol)))
  }
  if(logplot == F) {
    bac3 <- bac %>% 
      mutate(total = cfuA + cfuL) %>%
      group_by(repexp, exp,testcol) %>%
      dplyr::summarize(cfuA = median(cfuA, na.rm=T), cfuL = median(cfuL, na.rm=T), total = median(total, na.rm=T)) %>%
      ungroup() %>%
      group_by(exp,testcol) %>%
      dplyr::summarize(cfuA = mean(cfuA, na.rm=T), cfuL = mean(cfuL, na.rm=T), total = mean(total, na.rm=T)) %>%
      ungroup() %>%
      group_by(testcol) %>%
      dplyr::summarize(cfuA = mean(cfuA, na.rm=T), cfuL = mean(cfuL, na.rm=T), total = mean(total, na.rm=T)) %>%
      ungroup() %>%
      reshape2::melt(measure.vars = c("cfuA","cfuL")) %>%
      mutate(tempvar = paste0(testcol,"_",variable)) %>%
      inner_join(bac %>% 
                   mutate(total = cfuA + cfuL) %>%
                   group_by(repexp, exp,testcol) %>%
                   dplyr::summarize(cfuA = median(cfuA, na.rm=T), cfuL = median(cfuL, na.rm=T), total = median(total, na.rm=T)) %>%
                   ungroup() %>%
                   group_by(exp,testcol) %>%
                   dplyr::summarize(cfuA = mean(cfuA, na.rm=T), cfuL = mean(cfuL, na.rm=T), total = mean(total, na.rm=T)) %>%
                   ungroup() %>%
                   group_by(testcol) %>%
                   dplyr::summarize(semcfuA = sd(cfuA, na.rm=T)/sqrt(dplyr::n()), semcfuL = sd(cfuL, na.rm=T)/sqrt(dplyr::n()), semtotal = sd(total, na.rm=T)/sqrt(dplyr::n())) %>%
                   ungroup() %>%
                   reshape2::melt(measure.vars = c("semcfuA","semcfuL"))  %>%
                   mutate(mergecol = gsub(pattern = "sem", replacement = "", x = variable)) %>%
                   mutate(tempvar = paste0(testcol,"_",mergecol)), by = "tempvar") 
  } else {
    bac3 <- bac %>% 
      mutate(total = cfuA + cfuL) %>%
      group_by(repexp, exp,testcol) %>%
      dplyr::summarize(cfuA = mean(log10(cfuA+1), na.rm=T), cfuL = mean(log10(cfuL+1), na.rm=T), total = mean(log10(total+1), na.rm=T)) %>%
      ungroup() %>%
      group_by(exp,testcol) %>%
      dplyr::summarize(cfuA = mean(cfuA, na.rm=T), cfuL = mean(cfuL, na.rm=T), total = mean(total, na.rm=T)) %>%
      ungroup() %>%
      group_by(testcol) %>%
      dplyr::summarize(cfuA = mean(cfuA, na.rm=T), cfuL = mean(cfuL, na.rm=T), total = mean(total, na.rm=T)) %>%
      ungroup() %>%
      reshape2::melt(measure.vars = c("cfuA","cfuL")) %>%
      mutate(tempvar = paste0(testcol,"_",variable)) %>%
      inner_join(bac %>% 
                   mutate(total = cfuA + cfuL) %>%
                   group_by(repexp, exp,testcol) %>%
                   dplyr::summarize(cfuA = mean(log10(cfuA+1), na.rm=T), cfuL = mean(log10(cfuL+1), na.rm=T), total = mean(log10(total+1), na.rm=T)) %>%
                   ungroup() %>%
                   group_by(exp,testcol) %>%
                   dplyr::summarize(cfuA = mean(cfuA, na.rm=T), cfuL = mean(cfuL, na.rm=T), total = mean(total, na.rm=T)) %>%
                   ungroup() %>%
                   group_by(testcol) %>%
                   dplyr::summarize(semcfuA = sd(cfuA, na.rm=T)/sqrt(dplyr::n()), semcfuL = sd(cfuL, na.rm=T)/sqrt(dplyr::n()), semtotal = sd(total, na.rm=T)/sqrt(dplyr::n())) %>%
                   ungroup() %>%
                   reshape2::melt(measure.vars = c("semcfuA","semcfuL"))  %>%
                   mutate(mergecol = gsub(pattern = "sem", replacement = "", x = variable)) %>%
                   mutate(tempvar = paste0(testcol,"_",mergecol)), by = "tempvar") 
  }
  
  
  print(paste0("N:",dim(bac)[1]))
  bac3$addcol <- bac3$value.x
  numvals <- length(table(list(bac$testcol)))
  
  if (aceto_bottom == T) {
    bac3$variable.x = factor(bac3$variable.x, levels = c("cfuL","cfuA"))
    for (i in (numvals+1):(numvals*2)) {bac3$addcol[i] = bac3$total[i]}
  } else {
    bac3$variable.x = factor(bac3$variable.x, levels = c("cfuA","cfuL"))
    for (i in 1:numvals) {bac3$addcol[i] = bac3$total[i]}
  }
  
  cat("Stats on AAB abundances")
  print(kruskal.test(bac$cfuA ~ bac$testcol))
  if(kruskal.test(bac$cfuA ~ bac$testcol)$p.value < 0.05) {
    efg <- dunn.test(bac$cfuA,bac$testcol, method = "bh", table = F, kw = F)
    ghi <- cldList(comparison = efg$comparisons, p.value = efg$P.adjusted, threshold = 0.05)
    ghi$Group2 <- factor(ghi$Group, levels = levels(bac$testcol))
    cfuAclds <- as.character(ghi %>% arrange(Group2) %>% dplyr::select(Letter) %>% unlist() %>% unname())
    print("")
    cat("AAB Plot_values:")
    cat(ghi$Group)
    cat(" ")
    cat(cfuAclds)
  } else {
    cfuAclds <- c(rep("a", numvals))
    print("")
    cat("AAB Plot_values:")
    cat(cfuAclds)
  }
  
  cat("Stats on LAB abundances")
  print(kruskal.test(bac$cfuL ~ bac$testcol))
  
  if(kruskal.test(bac$cfuL ~ bac$testcol)$p.value < 0.05) {
    efg <- dunn.test(bac$cfuL,bac$testcol, method = "bh", table = F, kw = F)
    ghi <- cldList(comparison = efg$comparisons, p.value = efg$P.adjusted, threshold = 0.05)
    ghi$Group2 <- factor(ghi$Group, levels = levels(bac$testcol))
    cfuLclds <- as.character(ghi %>% arrange(Group2) %>% dplyr::select(Letter) %>% unlist() %>% unname())
    print("")
    cat("LAB Plot_values:")
    cat(ghi$Group)
    cat(" ")
    cat(cfuLclds)
  } else {
    cfuLclds <- c(rep("a", numvals))
    print("")
    cat("LAB Plot_values:")
    cat(cfuLclds)
  }
  
  if(!is.null(yaxis)) {    
    if(aceto_bottom != T) {
      bac3plot <- ggplot(bac3, aes(x=testcol.x, y=value.x, fill = variable.x)) +
        geom_bar(stat = "identity") + 
        geom_errorbar(aes(ymin =addcol-value.y, ymax = addcol+value.y), width = 0.5) +
        scale_fill_manual(values = c("red","blue"))+
        theme(legend.position = "bottom") + 
        ylim(yaxis)+
   #     geom_text(x = c(1:numvals,1:numvals), aes(y = addcol-label_adjust-value.y, label = c(rep("",numvals),cfuAclds)))+
    #    geom_text(x = c(1:numvals,1:numvals), aes(y = addcol+label_adjust+value.y, label = c(cfuLclds,rep("",numvals))))+
        theme_cowplot()
    } else {
      bac3plot <- ggplot(bac3, aes(x=testcol.x, y=value.x, fill = variable.x)) +
        geom_bar(stat = "identity") + 
        geom_errorbar(aes(ymin =addcol-value.y, ymax = addcol+value.y), width = 0.5) +
        scale_fill_manual(values = c("blue","red"))+
        theme(legend.position = "bottom") +  
        ylim(yaxis)+
     #   geom_text(x = c(1:numvals,1:numvals), aes(y = addcol-label_adjust-value.y, label = c(cfuAclds,rep("",numvals))))+
      #  geom_text(x = c(1:numvals,1:numvals), aes(y = addcol+label_adjust+value.y, label = c(rep("",numvals),cfuLclds)))+
        theme_cowplot()
    }
  } else {
    if(aceto_bottom != T) {
      bac3plot <- ggplot(bac3, aes(x=testcol.x, y=value.x, fill = variable.x)) +
        geom_bar(stat = "identity") + 
        geom_errorbar(aes(ymin =addcol-value.y, ymax = addcol+value.y), width = 0.5) +
        scale_fill_manual(values = c("red","blue"))+
        theme(legend.position = "bottom") + 
   #     geom_text(x = c(1:numvals,1:numvals), aes(y = addcol-label_adjust-value.y, label = c(rep("",numvals),cfuAclds)))+
    #    geom_text(x = c(1:numvals,1:numvals), aes(y = addcol+label_adjust+value.y, label = c(cfuLclds,rep("",numvals))))+
        theme_cowplot()
    } else {
      bac3plot <- ggplot(bac3, aes(x=testcol.x, y=value.x, fill = variable.x)) +
        geom_bar(stat = "identity") + 
        geom_errorbar(aes(ymin =addcol-value.y, ymax = addcol+value.y), width = 0.5) +
        scale_fill_manual(values = c("blue","red"))+
        theme(legend.position = "bottom") +  
     #   geom_text(x = c(1:numvals,1:numvals), aes(y = addcol-label_adjust-value.y, label = c(cfuAclds,rep("",numvals))))+
      #  geom_text(x = c(1:numvals,1:numvals), aes(y = addcol+label_adjust+value.y, label = c(rep("",numvals),cfuLclds)))+
        theme_cowplot()
    }
  }
  
  # old
  bac4 <- bac %>%
    mutate(total = cfuA + cfuL) %>%
    group_by(repexp, exp,testcol) %>%
    dplyr::summarize(cfuA = median(cfuA, na.rm=T), cfuL = median(cfuL, na.rm=T),  total = median(total, na.rm=T)) %>%
    ungroup() %>%
    group_by(exp,testcol) %>%
    dplyr::summarize(cfuA = mean(cfuA, na.rm=T), cfuL = mean(cfuL, na.rm=T), total = mean(total, na.rm=T)) %>%
    ungroup() %>%
    group_by(testcol) %>%
    dplyr::summarize(cfuA = mean(cfuA, na.rm=T), cfuL = mean(cfuL, na.rm=T), total = mean(total, na.rm=T)) %>%
    ungroup() %>%
    rowwise() %>%
    mutate(totalcol = cfuA+cfuL) %>%
    mutate(aperc = cfuA/totalcol, lperc = cfuL/totalcol) %>%
    ungroup() %>%
    reshape2::melt(measure.vars = c("aperc","lperc"))
  
  bac4 <- bac %>% 
    mutate(total = cfuA + cfuL) %>%
    group_by(repexp, exp,testcol) %>%
    dplyr::summarize(cfuA = median(cfuA/total, na.rm=T), cfuL = median(cfuL/total, na.rm=T),  total = median(total/total, na.rm=T)) %>%
    ungroup() %>%
    group_by(exp,testcol) %>%
    dplyr::summarize(cfuA = mean(cfuA, na.rm=T), cfuL = mean(cfuL, na.rm=T), total = mean(total, na.rm=T)) %>%
    ungroup() %>%
    group_by(testcol) %>%
    dplyr::summarize(cfuA = mean(cfuA, na.rm=T), cfuL = mean(cfuL, na.rm=T), total = mean(total, na.rm=T)) %>%
    ungroup() %>%
    rowwise() %>%
    mutate(totalcol = cfuA+cfuL) %>%
    mutate(aperc = cfuA/totalcol, lperc = cfuL/totalcol) %>%
    ungroup() %>%
    reshape2::melt(measure.vars = c("aperc","lperc"))
  
  bac4$addcol <- bac4$value
  if (aceto_bottom == T) {
    bac4$variable = factor(bac4$variable, levels = c("lperc","aperc"))
    for (i in (numvals+1):(numvals*2)) {bac4$addcol[i] = bac4$value[i-numvals]}
  } else {
    bac4$variable = factor(bac4$variable, levels = c("aperc","lperc"))
    for (i in 1:numvals) {bac4$addcol[i] = bac4$value[i+numvals]}
  }
  
  
  cat("Stats on abundances")
  print(kruskal.test(bac$aperc ~ bac$testcol))
  if(kruskal.test(bac$aperc ~ bac$testcol)$p.value < 0.05) {
    efg <- dunn.test(bac$aperc,bac$testcol, method = "bh", table = F, kw = F)
    ghi <- cldList(comparison = efg$comparisons, p.value = efg$P.adjusted, threshold = 0.05) 
    ghi$Group2 <- factor(ghi$Group, levels = levels(bac$testcol))
    apercclds <- as.character(ghi %>% arrange(Group2) %>% dplyr::select(Letter) %>% unlist() %>% unname())
    print("")
    cat("AAB Plot_values:")
    cat(ghi$Group)
    cat(" ")
    cat(apercclds)
  } else {
    apercclds <- c(rep("a",numvals))
    print("")
    cat("AAB Plot_values:")
    cat(apercclds)
  }
  
  if(kruskal.test(bac$lperc ~ bac$testcol)$p.value < 0.05) {
    efg <- dunn.test(bac$lperc,bac$testcol, method = "bh", table = F, kw = F)
    ghi <- cldList(comparison = efg$comparisons, p.value = efg$P.adjusted, threshold = 0.05) 
    ghi$Group2 <- factor(ghi$Group, levels = levels(bac$testcol))
    lpercclds <- as.character(ghi %>% arrange(Group2) %>% dplyr::select(Letter) %>% unlist() %>% unname())
    print("")
    cat("LAB Plot_values: ")
    cat(ghi$Group)
    cat(" ")
    cat(lpercclds)
  } else {
    lpercclds <-c(rep("a",numvals))
    print("")
    cat("LAB Plot_values:")
    cat(lpercclds)
  }
  
  if(aceto_bottom != T) {
    
    bac4plot <- ggplot(bac4, aes(x=testcol, y=value, fill = variable)) +
      geom_bar(stat = "identity", position = "stack") + theme(legend.position = "bottom") +
      scale_fill_manual(values = c("red","blue"))+
      #      geom_text(x = c(1:numvals,1:numvals), aes(y = addcol-0.05, label = c(rep("",numvals),cfuAclds)))+
      #    geom_text(x = c(1:numvals,1:numvals), aes(y = addcol+0.05, label = c(cfuLclds,rep("",numvals))))+
      theme_cowplot()
    
  } else {
    bac4plot <- ggplot(bac4, aes(x=testcol, y=value, fill = variable)) +
      geom_bar(stat = "identity", position = "stack") + theme(legend.position = "bottom") +
      scale_fill_manual(values = c("blue","red"))+
      # geom_text(x = c(1:numvals,1:numvals), aes(y = addcol-0.05, label = c(cfuAclds,rep("",numvals))))+
      #  geom_text(x = c(1:numvals,1:numvals), aes(y = addcol+0.05, label = c(rep("",numvals),cfuLclds)))+
      theme_cowplot()
  }
  
  ifelse(plotout == 1, return(bac3plot), ifelse(plotout == 2, return(bac4plot), ifelse(plotout==3, return(grid.arrange(bac3plot, bac4plot, heights = c(1,1))))))
}





make_chart_stats_96_1mean <- function(bac, column, plotout, yaxis = NULL, factorspecs = NULL, aceto_bottom =T, label_adjust = 10000, logplot = logplot) {
  bac$testcol <- unlist((bac[,column]))
  bac <- bac %>%
    mutate(aperc = cfuA/(cfuA+cfuL), lperc = cfuL/(cfuA+cfuL)) %>%
    droplevels()
  if(!is.null(factorspecs)) {
    bac$testcol = factor(bac$testcol, levels = factorspecs)
  } else {
    bac$testcol = factor(bac$testcol)
    print(levels(factor(bac$testcol)))
  }
  if(logplot == F) {
    bac3 <- bac %>% 
      mutate(total = cfuA + cfuL) %>%
      group_by(testcol) %>%
      dplyr::summarize(cfuA = mean(cfuA, na.rm=T), cfuL = mean(cfuL, na.rm=T), total = mean(total, na.rm=T)) %>%
      ungroup() %>%
      reshape2::melt(measure.vars = c("cfuA","cfuL")) %>%
      mutate(tempvar = paste0(testcol,"_",variable)) %>%
      inner_join(bac %>% 
                   mutate(total = cfuA + cfuL) %>%
                   group_by(testcol) %>%
                   dplyr::summarize(semcfuA = sd(cfuA, na.rm=T)/sqrt(dplyr::n()), semcfuL = sd(cfuL, na.rm=T)/sqrt(dplyr::n()), semtotal = sd(total, na.rm=T)/sqrt(dplyr::n())) %>%
                   ungroup() %>%
                   reshape2::melt(measure.vars = c("semcfuA","semcfuL"))  %>%
                   mutate(mergecol = gsub(pattern = "sem", replacement = "", x = variable)) %>%
                   mutate(tempvar = paste0(testcol,"_",mergecol)), by = "tempvar") 
  } else {
    bac3 <- bac %>% 
      mutate(total = cfuA + cfuL) %>%
      group_by(testcol) %>%
      dplyr::summarize(cfuA = mean(cfuA, na.rm=T), cfuL = mean(cfuL, na.rm=T), total = mean(total, na.rm=T)) %>%
      ungroup() %>%
      reshape2::melt(measure.vars = c("cfuA","cfuL")) %>%
      mutate(tempvar = paste0(testcol,"_",variable)) %>%
      inner_join(bac %>% 
                   mutate(total = cfuA + cfuL) %>%
                   group_by(testcol) %>%
                   dplyr::summarize(semcfuA = sd(cfuA, na.rm=T)/sqrt(dplyr::n()), semcfuL = sd(cfuL, na.rm=T)/sqrt(dplyr::n()), semtotal = sd(total, na.rm=T)/sqrt(dplyr::n())) %>%
                   ungroup() %>%
                   reshape2::melt(measure.vars = c("semcfuA","semcfuL"))  %>%
                   mutate(mergecol = gsub(pattern = "sem", replacement = "", x = variable)) %>%
                   mutate(tempvar = paste0(testcol,"_",mergecol)), by = "tempvar") 
  }
  
  
  print(paste0("N:",dim(bac)[1]))
  bac3$addcol <- bac3$value.x
  numvals <- length(table(list(bac$testcol)))
  
  if (aceto_bottom == T) {
    bac3$variable.x = factor(bac3$variable.x, levels = c("cfuL","cfuA"))
    for (i in (numvals+1):(numvals*2)) {bac3$addcol[i] = bac3$total[i]}
  } else {
    bac3$variable.x = factor(bac3$variable.x, levels = c("cfuA","cfuL"))
    for (i in 1:numvals) {bac3$addcol[i] = bac3$total[i]}
  }
  
  cat("Stats on AAB abundances")
  if(length(table(bac$testcol))<3) {
    print(wilcox.test(bac$cfuA ~ bac$testcol))
  } else {
    print(kruskal.test(bac$cfuA ~ bac$testcol))
    if(kruskal.test(bac$cfuA ~ bac$testcol)$p.value < 0.05) {
      efg <- dunn.test(bac$cfuA,bac$testcol, method = "bh", table = F, kw = F)
      ghi <- cldList(comparison = efg$comparisons, p.value = efg$P.adjusted, threshold = 0.05)
      ghi$Group2 <- factor(ghi$Group, levels = levels(bac$testcol))
      cfuAclds <- as.character(ghi %>% arrange(Group2) %>% dplyr::select(Letter) %>% unlist() %>% unname())
      print("")
      cat("AAB Plot_values:")
      cat(ghi$Group)
      cat(" ")
      cat(cfuAclds)
    } else {
      cfuAclds <- c(rep("a", numvals))
      print("")
      cat("AAB Plot_values:")
      cat(cfuAclds)
    }
  }
  
  cat("Stats on LAB abundances")
  
  if(length(table(bac$testcol))<3) {
    print(wilcox.test(bac$cfuL ~ bac$testcol))
  } else {
    print(kruskal.test(bac$cfuL ~ bac$testcol))
    if(kruskal.test(bac$cfuL ~ bac$testcol)$p.value < 0.05) {
    efg <- dunn.test(bac$cfuL,bac$testcol, method = "bh", table = F, kw = F)
    ghi <- cldList(comparison = efg$comparisons, p.value = efg$P.adjusted, threshold = 0.05)
    ghi$Group2 <- factor(ghi$Group, levels = levels(bac$testcol))
    cfuLclds <- as.character(ghi %>% arrange(Group2) %>% dplyr::select(Letter) %>% unlist() %>% unname())
    print("")
    cat("LAB Plot_values:")
    cat(ghi$Group)
    cat(" ")
    cat(cfuLclds)
  } else {
    cfuLclds <- c(rep("a", numvals))
    print("")
    cat("LAB Plot_values:")
    cat(cfuLclds)
  }
  }
  
  if(!is.null(yaxis)) {    
    if(aceto_bottom != T) {
      bac3plot <- ggplot(bac3, aes(x=testcol.x, y=value.x, fill = variable.x)) +
        geom_bar(stat = "identity") + 
        geom_errorbar(aes(ymin =addcol-value.y, ymax = addcol+value.y), width = 0.5) +
        scale_fill_manual(values = c("red","blue"))+
        theme(legend.position = "bottom") + 
        ylim(yaxis)+
        #     geom_text(x = c(1:numvals,1:numvals), aes(y = addcol-label_adjust-value.y, label = c(rep("",numvals),cfuAclds)))+
        #    geom_text(x = c(1:numvals,1:numvals), aes(y = addcol+label_adjust+value.y, label = c(cfuLclds,rep("",numvals))))+
        theme_cowplot()
    } else {
      bac3plot <- ggplot(bac3, aes(x=testcol.x, y=value.x, fill = variable.x)) +
        geom_bar(stat = "identity") + 
        geom_errorbar(aes(ymin =addcol-value.y, ymax = addcol+value.y), width = 0.5) +
        scale_fill_manual(values = c("blue","red"))+
        theme(legend.position = "bottom") +  
        ylim(yaxis)+
        #   geom_text(x = c(1:numvals,1:numvals), aes(y = addcol-label_adjust-value.y, label = c(cfuAclds,rep("",numvals))))+
        #  geom_text(x = c(1:numvals,1:numvals), aes(y = addcol+label_adjust+value.y, label = c(rep("",numvals),cfuLclds)))+
        theme_cowplot()
    }
  } else {
    if(aceto_bottom != T) {
      bac3plot <- ggplot(bac3, aes(x=testcol.x, y=value.x, fill = variable.x)) +
        geom_bar(stat = "identity") + 
        geom_errorbar(aes(ymin =addcol-value.y, ymax = addcol+value.y), width = 0.5) +
        scale_fill_manual(values = c("red","blue"))+
        theme(legend.position = "bottom") + 
        #     geom_text(x = c(1:numvals,1:numvals), aes(y = addcol-label_adjust-value.y, label = c(rep("",numvals),cfuAclds)))+
        #    geom_text(x = c(1:numvals,1:numvals), aes(y = addcol+label_adjust+value.y, label = c(cfuLclds,rep("",numvals))))+
        theme_cowplot()
    } else {
      bac3plot <- ggplot(bac3, aes(x=testcol.x, y=value.x, fill = variable.x)) +
        geom_bar(stat = "identity") + 
        geom_errorbar(aes(ymin =addcol-value.y, ymax = addcol+value.y), width = 0.5) +
        scale_fill_manual(values = c("blue","red"))+
        theme(legend.position = "bottom") +  
        #   geom_text(x = c(1:numvals,1:numvals), aes(y = addcol-label_adjust-value.y, label = c(cfuAclds,rep("",numvals))))+
        #  geom_text(x = c(1:numvals,1:numvals), aes(y = addcol+label_adjust+value.y, label = c(rep("",numvals),cfuLclds)))+
        theme_cowplot()
    }
  }
  
  # old
  bac4 <- bac %>%
    mutate(total = cfuA + cfuL) %>%
    group_by(testcol) %>%
    dplyr::summarize(cfuA = mean(cfuA, na.rm=T), cfuL = mean(cfuL, na.rm=T), total = mean(total, na.rm=T)) %>%
    ungroup() %>%
    rowwise() %>%
    mutate(totalcol = cfuA+cfuL) %>%
    mutate(aperc = cfuA/totalcol, lperc = cfuL/totalcol) %>%
    ungroup() %>%
    reshape2::melt(measure.vars = c("aperc","lperc"))
  
  bac4 <- bac %>% 
    mutate(total = cfuA + cfuL) %>%
    group_by(testcol) %>%
    dplyr::summarize(cfuA = mean(cfuA, na.rm=T), cfuL = mean(cfuL, na.rm=T), total = mean(total, na.rm=T)) %>%
    ungroup() %>%
    rowwise() %>%
    mutate(totalcol = cfuA+cfuL) %>%
    mutate(aperc = cfuA/totalcol, lperc = cfuL/totalcol) %>%
    ungroup() %>%
    reshape2::melt(measure.vars = c("aperc","lperc"))
  
  bac4$addcol <- bac4$value
  if (aceto_bottom == T) {
    bac4$variable = factor(bac4$variable, levels = c("lperc","aperc"))
    for (i in (numvals+1):(numvals*2)) {bac4$addcol[i] = bac4$value[i-numvals]}
  } else {
    bac4$variable = factor(bac4$variable, levels = c("aperc","lperc"))
    for (i in 1:numvals) {bac4$addcol[i] = bac4$value[i+numvals]}
  }
  
  
  cat("Stats on abundances")
  if(length(table(bac$testcol))<3) {
    print(wilcox.test(bac$aperc ~ bac$testcol))
  } else {
    print(kruskal.test(bac$aperc ~ bac$testcol))
    if(kruskal.test(bac$aperc ~ bac$testcol)$p.value < 0.05) {
    efg <- dunn.test(bac$aperc,bac$testcol, method = "bh", table = F, kw = F)
    ghi <- cldList(comparison = efg$comparisons, p.value = efg$P.adjusted, threshold = 0.05) 
    ghi$Group2 <- factor(ghi$Group, levels = levels(bac$testcol))
    apercclds <- as.character(ghi %>% arrange(Group2) %>% dplyr::select(Letter) %>% unlist() %>% unname())
    print("")
    cat("AAB Plot_values:")
    cat(ghi$Group)
    cat(" ")
    cat(apercclds)
  } else {
    apercclds <- c(rep("a",numvals))
    print("")
    cat("AAB Plot_values:")
    cat(apercclds)
  }
  }
  
  if(length(table(bac$testcol))<3) {
  } else {
    if(kruskal.test(bac$lperc ~ bac$testcol)$p.value < 0.05) {
    efg <- dunn.test(bac$lperc,bac$testcol, method = "bh", table = F, kw = F)
    ghi <- cldList(comparison = efg$comparisons, p.value = efg$P.adjusted, threshold = 0.05) 
    ghi$Group2 <- factor(ghi$Group, levels = levels(bac$testcol))
    lpercclds <- as.character(ghi %>% arrange(Group2) %>% dplyr::select(Letter) %>% unlist() %>% unname())
    print("")
    cat("LAB Plot_values: ")
    cat(ghi$Group)
    cat(" ")
    cat(lpercclds)
  } else {
    lpercclds <-c(rep("a",numvals))
    print("")
    cat("LAB Plot_values:")
    cat(lpercclds)
  }
  }
  
  if(aceto_bottom != T) {
    
    bac4plot <- ggplot(bac4, aes(x=testcol, y=value, fill = variable)) +
      geom_bar(stat = "identity", position = "stack") + theme(legend.position = "bottom") +
      scale_fill_manual(values = c("red","blue"))+
      #      geom_text(x = c(1:numvals,1:numvals), aes(y = addcol-0.05, label = c(rep("",numvals),cfuAclds)))+
      #    geom_text(x = c(1:numvals,1:numvals), aes(y = addcol+0.05, label = c(cfuLclds,rep("",numvals))))+
      theme_cowplot()
    
  } else {
    bac4plot <- ggplot(bac4, aes(x=testcol, y=value, fill = variable)) +
      geom_bar(stat = "identity", position = "stack") + theme(legend.position = "bottom") +
      scale_fill_manual(values = c("blue","red"))+
      # geom_text(x = c(1:numvals,1:numvals), aes(y = addcol-0.05, label = c(cfuAclds,rep("",numvals))))+
      #  geom_text(x = c(1:numvals,1:numvals), aes(y = addcol+0.05, label = c(rep("",numvals),cfuLclds)))+
      theme_cowplot()
  }
  
  ifelse(plotout == 1, return(bac3plot), ifelse(plotout == 2, return(bac4plot), ifelse(plotout==3, return(grid.arrange(bac3plot, bac4plot, heights = c(1,1))))))
}



convert_96_to_df <- function(plate1, plate_name) {
  ## from maggie_june2022_meta.Rmd
  new_df <- data.frame(plate=character(), well_position = character(), aab = character(), lb = character(), lp = character())
  for(i in 1:dim(plate1)[1]) {
    for (j in 1:dim(plate1)[2]) {
      vec <- as.character(plate1[i,j])
      vec2 <- strsplit(x = vec, split = ",")[[1]]
      if(length(vec2) == 1) {
        new_df_row <- data.frame(plate=as.character(plate_name), well_position = as.character(paste0(j, LETTERS[i])), aab = as.character(vec2), lb = as.character(0), lp = as.character(0))
      }
      if(length(vec2) == 2) {
        new_df_row <- data.frame(plate=as.character(plate_name), well_position = as.character(paste0(j, LETTERS[i])), aab = as.character(vec2[1]), lb = as.character(vec2[2]), lp = as.character(0))
      }
      if(length(vec2) == 3) {
        new_df_row <- data.frame(plate=as.character(plate_name), well_position = as.character(paste0(j, LETTERS[i])), aab = as.character(vec2[1]), lb = as.character(vec2[3]), lp = as.character(vec2[2]))
      }
      new_df <- rbind(new_df, new_df_row)
    }
  }
  return(new_df)
}
prep_plates <- function(plate1_vals, plate1_names, name_of_plate) {
  ## from maggie_june2022_meta.Rmd
  return(convert_96_to_df(plate1_vals, name_of_plate) %>% 
           inner_join(convert_96_to_df(plate1_names, name_of_plate) %>% dplyr::select(well_position, name = aab), by = "well_position") %>%
           filter(aab != "X") %>%
           tidyr::separate(aab, c("aab_dil", "aab_count"), ":", convert = T) %>%
           tidyr::separate(lb, c("lb_dil", "lb_count"), ":", convert = T) %>%
           tidyr::separate(lp, c("lp_dil", "lp_count"), ":", convert = T) %>%
           tidyr::separate(name, c("Genotype", "Treatment", "Sex","Replicate"), " ", convert = T, remove = F) %>% dplyr::select(-Replicate) %>% dplyr::rename(Replicate = name) %>%
           rowwise() %>%
           mutate(aab_dil = ifelse(aab_count == 0,ifelse(is.na(aab_dil),0, aab_dil), aab_dil)) %>%
           mutate(lb_dil = ifelse(lb_count == 0,ifelse(is.na(lb_dil),0, lb_dil), lb_dil)) %>%
           mutate(lp_dil = ifelse(lp_count == 0,ifelse(is.na(lp_dil),0, lp_dil), lp_dil)) %>%
           ungroup() %>%
           mutate(aab_dil = as.numeric(as.character(aab_dil)), lb_dil = as.numeric(as.character(lb_dil)), lp_dil = as.numeric(as.character(lp_dil)), aab_count = as.numeric(as.character(aab_count)), lb_count = as.numeric(as.character(lb_count)), lp_count = as.numeric(as.character(lp_count))) %>%
           mutate(aab = aab_count * (8^(aab_dil-1)) * 50) %>%
           mutate(lb = lb_count * (8^(lb_dil-1)) * 50) %>%
           mutate(lp = lp_count * (8^(lp_dil-1)) * 50)
  )
}

# datain = figure4adata_forplot
# humidityvar = c("C")
# sexvar = c("F","M")
# # datain = f4ce
# # humidityvar = "B"
# # sexvar = "F"


calculate_figure4means_pearson <- function(datain, humidityvar, sexvar) {
  
  ##
  # perc of means
  ##
  
  figure4asam_means = datain %>% 
    filter(Humidity%in%humidityvar, Sex %in% sexvar) %>%
    group_by(Geography) %>% 
    dplyr::summarize(meanA = mean(cfuA), 
                     semA = sd(cfuA)/sqrt(dplyr::n()), 
                     meanL = mean(cfuL), 
                     semL = sd(cfuL)/sqrt(dplyr::n()), 
                     count = dplyr::n()
    ) %>% 
    ungroup() %>% 
    mutate(percA = meanA/(meanA+meanL), percL = meanL/(meanA+meanL), total = meanA+meanL) %>% 
    droplevels()
  
  figure4asam_means <- figure4asam_means %>%
    mutate(latitude = c(1:length(table(figure4asam_means$Geography)))) 
  
  cat("\nperc of means\n")
  abc <- cor.test(figure4asam_means$latitude, figure4asam_means$meanA, method = "pearson")
  print(paste0("AAB: S = ",round(abc$statistic,2),", rho = ", round(abc$estimate,2),", p = ", round(abc$p.value,2)))
  abc <- cor.test(figure4asam_means$latitude, figure4asam_means$total, method = "pearson")
  print(paste0("Total: S = ",round(abc$statistic,2),", rho = ", round(abc$estimate,2),", p = ", round(abc$p.value,2)))
  abc <- cor.test(figure4asam_means$latitude, figure4asam_means$meanL, method = "pearson")
  print(paste0("LAB = ",round(abc$statistic,2),", rho = ", round(abc$estimate,2),", p = ", round(abc$p.value,2)))
  abc <- cor.test(figure4asam_means$latitude, figure4asam_means$percL, method = "pearson")
  print(paste0("percL = ",round(abc$statistic,2),", rho = ", round(abc$estimate,2),", p = ", round(abc$p.value,2)))
  
  plot(figure4asam_means$latitude, figure4asam_means$percL)  
  
  ##
  # mean of percs
  ## 
  
  figure4asam_means = datain %>% 
    filter(Humidity%in%humidityvar, Sex %in% sexvar) %>%
    mutate(percA = cfuA/(cfuL+cfuA), percL = cfuL/(cfuA+cfuL), total = 1) %>% 
    group_by(Geography) %>% 
    dplyr::summarize(meanA = mean(percA), meanL = mean(percL), total=1) %>% 
    ungroup() 
  
  print(table(figure4asam_means$Geography))
  figure4asam_means <- figure4asam_means %>%
    mutate(latitude = c(1:sum(table(figure4asam_means$Geography)!=0))) 
  print(figure4asam_means)
  
  cat("\nmean of percs\n")
  abc <- cor.test(figure4asam_means$latitude, figure4asam_means$meanA, method = "pearson")
  print(paste0("AAB: S = ",round(abc$statistic,2),", rho = ", round(abc$estimate,2),", p = ", round(abc$p.value,2)))
  abc <- cor.test(figure4asam_means$latitude, figure4asam_means$meanL, method = "pearson")
  print(paste0("LAB = ",round(abc$statistic,2),", rho = ", round(abc$estimate,2),", p = ", round(abc$p.value,2)))
  plot(figure4asam_means$latitude, figure4asam_means$meanL)  
  
  ##
  # cor of percs
  ##
  
  figure4asam_means = datain %>%
    filter(Humidity%in%humidityvar, Sex %in% sexvar) %>%
    mutate(meanA = cfuA, meanL = cfuL, total = cfuA+cfuL, percL = cfuL/total)
  print(figure4asam_means)
  cat("\nperc of raw data\n")
  abc <- cor.test(figure4asam_means$latitude, figure4asam_means$meanA, method = "pearson", exact=F)
  #  plot(figure4asam_means$latitude, log10(figure4asam_means$meanL+1))
  print(paste0("AAB: S = ",round(abc$statistic,2),", rho = ", round(abc$estimate,2),", p = ", round(abc$p.value,2)))
  abc <- cor.test(figure4asam_means$latitude, figure4asam_means$total, method = "pearson", exact=F)
  print(paste0("Total: S = ",round(abc$statistic,2),", rho = ", round(abc$estimate,2),", p = ", round(abc$p.value,2)))
  abc <- cor.test(figure4asam_means$latitude, figure4asam_means$meanL, method = "pearson", exact=F)
  print(paste0("LAB = ",round(abc$statistic,2),", rho = ", round(abc$estimate,2),", p = ", round(abc$p.value,2)))
  abc <- cor.test(figure4asam_means$latitude, figure4asam_means$percL, method = "pearson", exact=F)
  print(paste0("percL = ",round(abc$statistic,2),", rho = ", round(abc$estimate,2),", p = ", round(abc$p.value,2)))
  
  plot(figure4asam_means$latitude, figure4asam_means$percL)  
  
  ##
  # cor of perc means
  ## 
}

calculate_figure4means <- function(datain, humidityvar, sexvar) {
  
  ##
  # perc of means
  ##
  
  figure4asam_means = datain %>% 
    filter(Humidity%in%humidityvar, Sex %in% sexvar) %>%
    group_by(Geography) %>% 
    dplyr::summarize(meanA = mean(cfuA), 
                     semA = sd(cfuA)/sqrt(dplyr::n()), 
                     meanL = mean(cfuL), 
                     semL = sd(cfuL)/sqrt(dplyr::n()), 
                     count = dplyr::n()
    ) %>% 
    ungroup() %>% 
    mutate(percA = meanA/(meanA+meanL), percL = meanL/(meanA+meanL), total = meanA+meanL) %>% 
    droplevels()
  
  figure4asam_means <- figure4asam_means %>%
    mutate(latitude = c(1:length(table(figure4asam_means$Geography)))) 

  cat("\nperc of means\n")
  abc <- cor.test(figure4asam_means$latitude, figure4asam_means$meanA, method = "spearman")
  print(paste0("AAB: S = ",round(abc$statistic,2),", rho = ", round(abc$estimate,2),", p = ", round(abc$p.value,2)))
  abc <- cor.test(figure4asam_means$latitude, figure4asam_means$total, method = "spearman")
  print(paste0("Total: S = ",round(abc$statistic,2),", rho = ", round(abc$estimate,2),", p = ", round(abc$p.value,2)))
  abc <- cor.test(figure4asam_means$latitude, figure4asam_means$meanL, method = "spearman")
  print(paste0("LAB = ",round(abc$statistic,2),", rho = ", round(abc$estimate,2),", p = ", round(abc$p.value,2)))
  abc <- cor.test(figure4asam_means$latitude, figure4asam_means$percL, method = "spearman")
  print(paste0("percL = ",round(abc$statistic,2),", rho = ", round(abc$estimate,2),", p = ", round(abc$p.value,2)))
 
plot(figure4asam_means$latitude, figure4asam_means$percL)  
  
  ##
  # mean of percs
  ## 
  
  figure4asam_means = datain %>% 
    filter(Humidity%in%humidityvar, Sex %in% sexvar) %>%
    mutate(percA = cfuA/(cfuL+cfuA), percL = cfuL/(cfuA+cfuL), total = 1) %>% 
    group_by(Geography) %>% 
    dplyr::summarize(meanA = mean(percA), meanL = mean(percL), total=1) %>% 
    ungroup() 
  
  print(table(figure4asam_means$Geography))
  figure4asam_means <- figure4asam_means %>%
    mutate(latitude = c(1:sum(table(figure4asam_means$Geography)!=0))) 
  
  cat("\nmean of percs\n")
  abc <- cor.test(figure4asam_means$latitude, figure4asam_means$meanA, method = "spearman")
  print(paste0("AAB: S = ",round(abc$statistic,2),", rho = ", round(abc$estimate,2),", p = ", round(abc$p.value,2)))
  abc <- cor.test(figure4asam_means$latitude, figure4asam_means$meanL, method = "spearman")
  print(paste0("LAB = ",round(abc$statistic,2),", rho = ", round(abc$estimate,2),", p = ", round(abc$p.value,2)))

  ##
  # cor of percs
  ##
  
  figure4asam_means = datain %>%
    filter(Humidity%in%humidityvar, Sex %in% sexvar) %>%
    mutate(meanA = cfuA, meanL = cfuL, total = cfuA+cfuL, percL = cfuL/total)

  cat("\nperc of raw data\n")
  abc <- cor.test(figure4asam_means$latitude, figure4asam_means$meanA, method = "spearman", exact=F)
#  plot(figure4asam_means$latitude, log10(figure4asam_means$meanL+1))
  print(paste0("AAB: S = ",round(abc$statistic,2),", rho = ", round(abc$estimate,2),", p = ", round(abc$p.value,2)))
  abc <- cor.test(figure4asam_means$latitude, figure4asam_means$total, method = "spearman", exact=F)
  print(paste0("Total: S = ",round(abc$statistic,2),", rho = ", round(abc$estimate,2),", p = ", round(abc$p.value,2)))
  abc <- cor.test(figure4asam_means$latitude, figure4asam_means$meanL, method = "spearman", exact=F)
  print(paste0("LAB = ",round(abc$statistic,2),", rho = ", round(abc$estimate,2),", p = ", round(abc$p.value,2)))
  abc <- cor.test(figure4asam_means$latitude, figure4asam_means$percL, method = "spearman", exact=F)
  print(paste0("percL = ",round(abc$statistic,2),", rho = ", round(abc$estimate,2),", p = ", round(abc$p.value,2)))
  
   
}

calculate_figure4medians <- function(datain, humidityvar, sexvar) {
  
  ##
  # perc of medians
  ##
  
  figure4asam_medians = datain %>% 
    filter(Humidity%in%humidityvar, Sex %in% sexvar) %>%
    group_by(Geography) %>% 
    dplyr::summarize(medianA = median(cfuA), 
                     semA = sd(cfuA)/sqrt(dplyr::n()), 
                     medianL = median(cfuL), 
                     semL = sd(cfuL)/sqrt(dplyr::n()), 
                     count = dplyr::n()
    ) %>% 
    ungroup() %>% 
    mutate(percA = medianA/(medianA+medianL), percL = medianL/(medianA+medianL), total = medianA+medianL) %>% 
    droplevels()
  
  figure4asam_medians <- figure4asam_medians %>%
    mutate(latitude = c(1:length(table(figure4asam_medians$Geography)))) 
  
  cat("\nperc of medians\n")
  abc <- cor.test(figure4asam_medians$latitude, figure4asam_medians$medianA, method = "spearman")
  print(paste0("AAB: S = ",round(abc$statistic,2),", rho = ", round(abc$estimate,2),", p = ", round(abc$p.value,2)))
  abc <- cor.test(figure4asam_medians$latitude, figure4asam_medians$total, method = "spearman")
  print(paste0("Total: S = ",round(abc$statistic,2),", rho = ", round(abc$estimate,2),", p = ", round(abc$p.value,2)))
  abc <- cor.test(figure4asam_medians$latitude, figure4asam_medians$medianL, method = "spearman")
  print(paste0("LAB = ",round(abc$statistic,2),", rho = ", round(abc$estimate,2),", p = ", round(abc$p.value,2)))
  abc <- cor.test(figure4asam_medians$latitude, figure4asam_medians$percL, method = "spearman")
  print(paste0("percL = ",round(abc$statistic,2),", rho = ", round(abc$estimate,2),", p = ", round(abc$p.value,2)))
  
  ##
  # median of percs
  ## 
  
  figure4asam_medians = datain %>% 
    filter(Humidity%in%humidityvar, Sex %in% sexvar) %>%
    mutate(percA = cfuA/(cfuL+cfuA), percL = cfuL/(cfuA+cfuL), total = 1) %>% 
    group_by(Geography) %>% 
    dplyr::summarize(medianA = median(percA), medianL = median(percL), total=1) %>% 
    ungroup() 
  
  figure4asam_medians <- figure4asam_medians %>%
    mutate(latitude = c(1:length(table(figure4asam_medians$Geography)))) 
  
  cat("\nmedian of percs\n")
  print(figure4asam_medians)
  abc <- cor.test(figure4asam_medians$latitude, figure4asam_medians$medianA, method = "spearman")
  print(paste0("AAB: S = ",round(abc$statistic,2),", rho = ", round(abc$estimate,2),", p = ", round(abc$p.value,2)))
  abc <- cor.test(figure4asam_medians$latitude, figure4asam_medians$medianL, method = "spearman")
  print(paste0("LAB = ",round(abc$statistic,2),", rho = ", round(abc$estimate,2),", p = ", round(abc$p.value,2)))
  
  ##
  # cor of percs
  ##
  
  figure4asam_medians = datain %>%
    filter(Humidity%in%humidityvar, Sex %in% sexvar) %>%
    mutate(medianA = cfuA, medianL = cfuL, total = cfuA+cfuL)
  
  cat("\nperc of raw data\n")
  abc <- cor.test(figure4asam_medians$latitude, figure4asam_medians$medianA, method = "spearman", exact=F)
  #  plot(figure4asam_medians$latitude, log10(figure4asam_medians$medianL+1))
  print(paste0("AAB: S = ",round(abc$statistic,2),", rho = ", round(abc$estimate,2),", p = ", round(abc$p.value,2)))
  abc <- cor.test(figure4asam_medians$latitude, figure4asam_medians$total, method = "spearman", exact=F)
  print(paste0("Total: S = ",round(abc$statistic,2),", rho = ", round(abc$estimate,2),", p = ", round(abc$p.value,2)))
  abc <- cor.test(figure4asam_medians$latitude, figure4asam_medians$medianL, method = "spearman", exact=F)
  print(paste0("LAB = ",round(abc$statistic,2),", rho = ", round(abc$estimate,2),", p = ", round(abc$p.value,2)))
  abc <- cor.test(figure4asam_medians$latitude, figure4asam_medians$percL, method = "spearman", exact=F)
  print(paste0("percL = ",round(abc$statistic,2),", rho = ", round(abc$estimate,2),", p = ", round(abc$p.value,2)))
  
  
}


calculate_figure4meansv2 <- function(datain, humidityvar, sexvar) {

  figure4asam_means = datain %>% 
    filter(Humidity%in%humidityvar, Sex %in% sexvar) %>%
    mutate(percA = cfuA/(cfuL+cfuA), percL = cfuL/(cfuA+cfuL), total = 1) %>% 
    group_by(Geography) %>% 
    dplyr::summarize(meanA = mean(percA), meanL = mean(percL), total=1) %>% 
    ungroup() 

figure4asam_means <- figure4asam_means %>%
    mutate(latitude = c(1:length(table(figure4asam_means$Geography)))) 
  
  print(figure4asam_means)
  abc <- cor.test(1:max(figure4asam_means$latitude), figure4asam_means$meanA, method = "spearman")
  print(paste0("AAB: S = ",round(abc$statistic,2),", rho = ", round(abc$estimate,2),", p = ", round(abc$p.value,2)))
  abc <- cor.test(1:max(figure4asam_means$latitude), figure4asam_means$total, method = "spearman")
  print(paste0("Total: S = ",round(abc$statistic,2),", rho = ", round(abc$estimate,2),", p = ", round(abc$p.value,2)))
  abc <- cor.test(1:max(figure4asam_means$latitude), figure4asam_means$meanL, method = "spearman")
  print(paste0("LAB = ",round(abc$statistic,2),", rho = ", round(abc$estimate,2),", p = ", round(abc$p.value,2)))
  abc <- cor.test(1:max(figure4asam_means$latitude), figure4asam_means$percL, method = "spearman")
  print(paste0("percL = ",round(abc$statistic,2),", rho = ", round(abc$estimate,2),", p = ", round(abc$p.value,2)))
  
}


#datain = figure2bdata
#bacteriavar = "CV"
#sexvar = c("F","M")
calculate_figure2means <- function(datain, bacteriavar, sexvar) {
  
  figure4asam_means = datain %>% 
    filter(Bacteria%in%bacteriavar, Sex %in% sexvar) %>%
    group_by(Geography) %>% 
    dplyr::summarize(meanA = mean(cfuA), 
                     semA = sd(cfuA)/sqrt(dplyr::n()), 
                     meanL = mean(cfuL), 
                     semL = sd(cfuL)/sqrt(dplyr::n()), 
                     count = dplyr::n()
    ) %>% 
    ungroup() %>% 
    mutate(percA = meanA/(meanA+meanL), percL = meanL/(meanA+meanL), total = meanA+meanL) %>% 
    droplevels()
  
  
  
  figure4asam_means <- figure4asam_means %>%
    mutate(latitude = c(1:length(table(figure4asam_means$Geography)))) 
  
  abc <- cor.test(1:max(figure4asam_means$latitude), figure4asam_means$meanA, method = "spearman")
  print(paste0("AAB: S = ",round(abc$statistic,2),", rho = ", round(abc$estimate,2),", p = ", round(abc$p.value,2)))
  abc <- cor.test(1:max(figure4asam_means$latitude), figure4asam_means$total, method = "spearman")
  print(paste0("Total: S = ",round(abc$statistic,2),", rho = ", round(abc$estimate,2),", p = ", round(abc$p.value,2)))
  abc <- cor.test(1:max(figure4asam_means$latitude), figure4asam_means$meanL, method = "spearman")
  print(paste0("LAB = ",round(abc$statistic,2),", rho = ", round(abc$estimate,2),", p = ", round(abc$p.value,2)))
  abc <- cor.test(1:max(figure4asam_means$latitude), figure4asam_means$percL, method = "spearman")
  print(paste0("percL = ",round(abc$statistic,2),", rho = ", round(abc$estimate,2),", p = ", round(abc$p.value,2)))
  
}

make_chart_stats_96_jitter_1mean_orig <- function(bac, column, plotout, yaxis = NULL, factorspecs = NULL, logplot = F, label_adjust = 10000, baradjust = 0, baradjust_bottom = 0, loglabeladjust = 0, aceto_bottom = F, yaxis_perc = NULL) {
  bac$testcol <- unlist(unname(bac[,column]))
  bac <- bac %>%
    mutate(aperc = cfuA/(cfuA+cfuL), lperc = cfuL/(cfuA+cfuL)) %>%
    droplevels()
  if(!is.null(factorspecs)) {
    bac$testcol = factor(bac$testcol, levels = factorspecs)
  } else {
    bac$testcol = factor(bac$testcol)
    print(levels(factor(bac$testcol)))
  }
  
  if(logplot == F) {
    bac3 <- bac %>% 
      mutate(total = cfuA + cfuL) %>%
      group_by(testcol) %>%
      dplyr::summarize(cfuA = mean(cfuA, na.rm=T), cfuL = mean(cfuL, na.rm=T), total = mean(total, na.rm=T)) %>%
      ungroup() %>%
      reshape2::melt(measure.vars = c("cfuA","cfuL")) %>%
      mutate(tempvar = paste0(testcol,"_",variable)) %>%
      inner_join(bac %>% 
                   mutate(total = cfuA + cfuL) %>%
                   group_by(testcol) %>%
                   dplyr::summarize(semcfuA = sd(cfuA, na.rm=T)/sqrt(dplyr::n()), semcfuL = sd(cfuL, na.rm=T)/sqrt(dplyr::n())) %>%
                   ungroup() %>%
                   reshape2::melt(measure.vars = c("semcfuA","semcfuL"))  %>%
                   mutate(mergecol = gsub(pattern = "sem", replacement = "", x = variable)) %>%
                   mutate(tempvar = paste0(testcol,"_",mergecol)), by = "tempvar") 
  } else {
    if(aceto_bottom==F) {
      bac3 <- 1
      bac3 <- bac %>% 
        group_by(testcol) %>%
        dplyr::summarize(cfuLmean = mean(log10(cfuL+1), na.rm=T),  totalmean = mean(log10(cfuA+cfuL+1), na.rm=T)) %>%
        ungroup() %>%
        mutate(cfuAmean = totalmean-cfuLmean) %>%
        dplyr::select(testcol, cfuL = cfuLmean, total = totalmean, cfuA = cfuAmean) %>%
       reshape2::melt(measure.vars = c("cfuA","cfuL")) %>%
        mutate(tempvar = paste0(testcol,"_",variable)) %>%
        inner_join(bac %>% 
                     group_by(testcol) %>%
                     dplyr::summarize(semcfuA = sd(log10(cfuA+1), na.rm=T)/sqrt(dplyr::n()), semcfuL = sd(log10(cfuL+1), na.rm=T)/sqrt(dplyr::n())) %>%
                     ungroup() %>%
                     reshape2::melt(measure.vars = c("semcfuA","semcfuL"))  %>%
                     mutate(mergecol = gsub(pattern = "sem", replacement = "", x = variable)) %>%
                     mutate(tempvar = paste0(testcol,"_",mergecol)), by = "tempvar") 
    } else {
      print("hi")
      bac3 <- bac %>% 
        group_by(testcol) %>%
        dplyr::summarize(cfuAmean = mean(log10(cfuA+1), na.rm=T),  totalmean = mean(log10(cfuA+cfuL+1), na.rm=T)) %>%
        ungroup() %>%
        mutate(cfuLmean = totalmean-cfuAmean) %>%
        dplyr::select(testcol, cfuL = cfuLmean, total = totalmean, cfuA = cfuAmean) %>%
        reshape2::melt(measure.vars = c("cfuA","cfuL")) %>%
        mutate(tempvar = paste0(testcol,"_",variable)) %>%
        inner_join(bac %>% 
                     group_by(testcol) %>%
                     dplyr::summarize(semcfuA = sd(log10(cfuA+1), na.rm=T)/sqrt(dplyr::n()), semcfuL = sd(log10(cfuL+1), na.rm=T)/sqrt(dplyr::n())) %>%
                     ungroup() %>%
                     reshape2::melt(measure.vars = c("semcfuA","semcfuL"))  %>%
                     mutate(mergecol = gsub(pattern = "sem", replacement = "", x = variable)) %>%
                     mutate(tempvar = paste0(testcol,"_",mergecol)), by = "tempvar") 
    }
  }
  
bac3$value.y = ifelse(is.na(bac3$value.y),0,bac3$value.y) ## added 2023-08-12 - maybe remove if things get broken
  
  print(paste0("N:",dim(bac)[1]))
  bac3$addcol <- bac3$value.x
  numvals <- length(table(list(bac$testcol %>% droplevels)))
  if (aceto_bottom == T) {
    bac3$variable.x = factor(bac3$variable.x, levels = c("cfuL","cfuA"))
    for (i in (numvals+1):(numvals*2)) {bac3$addcol[i] = bac3$total[i]}
  } else {
    bac3$variable.x = factor(bac3$variable.x, levels = c("cfuA","cfuL"))
    for (i in 1:numvals) {bac3$addcol[i] = bac3$total[i]}
  }
  

  try(rm(ghi, mid_test),T)
  cat("Stats on AAB abundances")
  print(kruskal.test(bac$cfuA ~ bac$testcol))
  mid_test <- kruskal.test(bac$cfuA ~ bac$testcol)
  print(paste0("chi2 ",mid_test$parameter,", ",dim(bac)[1]," = ",round(mid_test$statistic,2),", p = ",round(mid_test$p.value,4)))
  if(kruskal.test(bac$cfuA ~ bac$testcol)$p.value < 0.05) {
    if (length(table(bac$testcol))>2) {
      efg <- dunn.test(bac$cfuA,bac$testcol, method = "bh", table = F, kw = T)
      ghi <- cldList(comparison = efg$comparisons, p.value = efg$P.adjusted, threshold = 0.05)
      ghi$Group2 <- factor(ghi$Group, levels = levels(bac$testcol))
      cfuAclds <- as.character(ghi %>% arrange(Group2) %>% dplyr::select(Letter) %>% unlist() %>% unname())
      print(ghi)
    } else {
      cfuAclds <- c("a","b")
    }
    cat("AAB Plot_values:")
  } else {
    cfuAclds <- c(rep("a", numvals))
    # cat("AAB Plot_values:")
    # cat(" ")
    # cat(cfuAclds)
  }
  
   try(rm(ghi, mid_test),T)
   cat("\n\nStats on LAB abundances")
  print(kruskal.test(bac$cfuL ~ bac$testcol))
  mid_test <- kruskal.test(bac$cfuL ~ bac$testcol)
  print(paste0("chi2 ",mid_test$parameter,", ",dim(bac)[1]," = ",round(mid_test$statistic,2),", p = ",round(mid_test$p.value,4)))
  if(kruskal.test(bac$cfuL ~ bac$testcol)$p.value < 0.05) {
    if (length(table(bac$testcol))>2) {
      efg <- dunn.test(bac$cfuL,bac$testcol, method = "bh", table = F, kw = F)
      ghi <- cldList(comparison = efg$comparisons, p.value = efg$P.adjusted, threshold = 0.05)
      ghi$Group2 <- factor(ghi$Group, levels = levels(bac$testcol))
      cfuLclds <- as.character(ghi %>% arrange(Group2) %>% dplyr::select(Letter) %>% unlist() %>% unname())
      print(ghi)
    } else {
      cfuLclds <- c("a","b")
    }
    cat("LAB Plot_values:")
  } else {
    cfuLclds <- c(rep("a", numvals))
    # cat("LAB Plot_values:")
    # cat(" ")
    # cat(cfuLclds)
  }
  
  # cat("absolute abundance stats")
  # print(kruskal.test(bac$cfuA ~ bac$testcol))
  # print(kruskal.test(bac$cfuL ~ bac$testcol))
  # print(dim(bac))
  
    if(!is.null(yaxis)) {    
      if(aceto_bottom != T) {
        bac3plot <-1 
        bac3plot <- ggplot(bac3, aes(x=testcol.x, y=value.x, fill = variable.x)) + geom_bar(stat = "identity") + 
          geom_errorbar(aes(ymin =addcol-value.y, ymax = addcol+value.y), width = 0.5) +
          scale_fill_manual(values = c("red","blue"))+
          theme(legend.position = "bottom") + 
          ylim(yaxis)+
          geom_text(data = bac3 %>% filter(variable.x == "cfuA"), aes(x = testcol.x, y = ((total+value.y)+baradjust), label = c(cfuAclds), vjust = 1)) + #, inherit.aes = T)) +
          geom_text(data = bac3 %>% filter(variable.x == "cfuL"), aes(x = testcol.x, y = 0+baradjust_bottom, label = c(cfuLclds), vjust = 1)) + #, inherit.aes = T)+
          theme_cowplot()
      } else {
        bac3plot <- ggplot(bac3, aes(x=testcol.x, y=value.x, fill = variable.x)) +
          geom_bar(stat = "identity") + 
          geom_errorbar(aes(ymin =addcol-value.y, ymax = addcol+value.y), width = 0.5) +
          scale_fill_manual(values = c("blue","red"))+
          theme(legend.position = "bottom") +  
          ylim(yaxis)+
          geom_text(data = bac3 %>% filter(variable.x == "cfuL"), aes(x = testcol.x, y = ((total+value.y)+baradjust), label = c(cfuLclds), vjust = 1)) +
          geom_text(data = bac3 %>% filter(variable.x == "cfuA"), aes(x = testcol.x, y = 0+baradjust_bottom, label = c(cfuAclds), vjust = 1)) + 
          theme_cowplot()
      }
    } else {
      if(aceto_bottom != T) {
        bac3plot <- ggplot(bac3, aes(x=testcol.x, y=value.x, fill = variable.x)) + geom_bar(stat = "identity") + 
          geom_errorbar(aes(ymin =addcol-value.y, ymax = addcol+value.y), width = 0.5) +
          scale_fill_manual(values = c("red","blue"))+
          theme(legend.position = "bottom") + 
          geom_text(data = bac3 %>% filter(variable.x == "cfuA"), aes(x = testcol.x, y = ((total+value.y)+baradjust), label = c(cfuAclds), vjust = 1)) +
          geom_text(data = bac3 %>% filter(variable.x == "cfuL"), aes(x = testcol.x, y = 0+baradjust_bottom, label = c(cfuLclds), vjust = 1)) +
          theme_cowplot()
      } else {
        bac3plot <- ggplot(bac3, aes(x=testcol.x, y=value.x, fill = variable.x)) +
          geom_bar(stat = "identity") + 
          geom_errorbar(aes(ymin =addcol-value.y, ymax = addcol+value.y), width = 0.5) +
          scale_fill_manual(values = c("blue","red"))+
          theme(legend.position = "bottom") +  
          geom_text(data = bac3 %>% filter(variable.x == "cfuL"), aes(x = testcol.x, y = ((total+value.y)+baradjust), label = c(cfuLclds), vjust = 1)) +
          geom_text(data = bac3 %>% filter(variable.x == "cfuA"), aes(x = testcol.x, y = 0+baradjust_bottom, label = c(cfuAclds), vjust = 1)) + 
          theme_cowplot()
      }
    }

  
#  if(logplot==F) { 
#   if(!is.null(yaxis)) {    
#     if(aceto_bottom != T) {
#       bac3plot <-1 
#       bac3plot <- ggplot(bac3, aes(x=testcol.x, y=value.x, fill = variable.x)) + geom_bar(stat = "identity") + 
#         geom_errorbar(aes(ymin =addcol-value.y, ymax = addcol+value.y), width = 0.5) +
#         scale_fill_manual(values = c("red","blue"))+
#         theme(legend.position = "bottom") + 
#         ylim(yaxis)+
#         geom_text(data = bac3 %>% filter(variable.x == "cfuA"), aes(x = testcol.x, y = 0, label = c(cfuAclds), vjust = 1)) + #, inherit.aes = T)) +
#         geom_text(data = bac3 %>% filter(variable.x == "cfuL"), aes(x = testcol.x, y = (bac$cfuL + bac$cfuA+1)*baradjust), label = c(cfuLclds), vjust = 1) + #, inherit.aes = T)+
#         theme_cowplot()
#     } else {
#       bac3plot <- ggplot(bac3, aes(x=testcol.x, y=value.x, fill = variable.x)) +
#         geom_bar(stat = "identity") + 
#         geom_errorbar(aes(ymin =addcol-value.y, ymax = addcol+value.y), width = 0.5) +
#         scale_fill_manual(values = c("blue","red"))+
#         theme(legend.position = "bottom") +  
#         ylim(yaxis)+
#         geom_text(data = bac3 %>% filter(variable.x == "cfuL"), aes(x = testcol.x, y = 0, label = c(cfuLclds), vjust = 1)) +
#         geom_text(data = bac3 %>% filter(variable.x == "cfuA"), aes(x = testcol.x, y = (bac$cfuL + bac$cfuA+1)*baradjust), label = c(cfuAclds), vjust = 1) + 
#         theme_cowplot()
#     }
#   } else {
#     if(aceto_bottom != T) {
#       bac3plot <- ggplot(bac3, aes(x=testcol.x, y=value.x, fill = variable.x)) + geom_bar(stat = "identity") + 
#         geom_errorbar(aes(ymin =addcol-value.y, ymax = addcol+value.y), width = 0.5) +
#         scale_fill_manual(values = c("red","blue"))+
#         theme(legend.position = "bottom") + 
#         geom_text(data = bac3 %>% filter(variable.x == "cfuA"), aes(x = length(table(bac3$testcol.x)), y = 0, label = c(cfuAclds), vjust = 1)) +
#         geom_text(data = bac3 %>% filter(variable.x == "cfuL"), aes(x = testcol.x, y = (bac$cfuL + bac$cfuA+1)*baradjust), label = c(cfuLclds), vjust = 1) +
#         theme_cowplot()
#     } else {
#       bac3plot <- ggplot(bac3, aes(x=testcol.x, y=value.x, fill = variable.x)) +
#         geom_bar(stat = "identity") + 
#         geom_errorbar(aes(ymin =addcol-value.y, ymax = addcol+value.y), width = 0.5) +
#         scale_fill_manual(values = c("blue","red"))+
#         theme(legend.position = "bottom") +  
#         geom_text(data = bac3 %>% filter(variable.x == "cfuL"), aes(x = testcol.x, y = 0, label = c(cfuLclds), vjust = 1)) +
#         geom_text(data = bac3 %>% filter(variable.x == "cfuA"), aes(x = testcol.x, y = (bac$cfuL + bac$cfuA+1)*baradjust), label = c(cfuAclds), vjust = 1) + 
#         theme_cowplot()
#     }
#   }
#  } else {
#    if(!is.null(yaxis)) {    
#      if(aceto_bottom != T) {
#        bac3plot <-1 
#        bac3plot <- ggplot(bac3, aes(x=testcol.x, y=value.x, fill = variable.x)) + geom_bar(stat = "identity") + 
#          geom_errorbar(aes(ymin =addcol-value.y, ymax = addcol+value.y), width = 0.5) +
#          scale_fill_manual(values = c("red","blue"))+
#          theme(legend.position = "bottom") + 
#          ylim(yaxis)+
#          geom_text(data = bac3 %>% filter(variable.x == "cfuA"), aes(x = testcol.x, y = 0, label = c(cfuAclds), vjust = 1)) + #, inherit.aes = T)) +
#          geom_text(data = bac3 %>% filter(variable.x == "cfuL"), aes(x = testcol.x, y = (log10(bac$cfuL+1) + log10(bac$cfuA+1)+1)*baradjust), label = c(cfuLclds), vjust = 1) + #, inherit.aes = T)+
#          theme_cowplot()
#      } else {
#        bac3plot <- ggplot(bac3, aes(x=testcol.x, y=value.x, fill = variable.x)) +
#          geom_bar(stat = "identity") + 
#          geom_errorbar(aes(ymin =addcol-value.y, ymax = addcol+value.y), width = 0.5) +
#          scale_fill_manual(values = c("blue","red"))+
#          theme(legend.position = "bottom") +  
#          ylim(yaxis)+
#          geom_text(data = bac3 %>% filter(variable.x == "cfuL"), aes(x = testcol.x, y = 0, label = c(cfuLclds), vjust = 1)) +
#          geom_text(data = bac3 %>% filter(variable.x == "cfuA"), aes(x = testcol.x, y = (log10(bac$cfuL+1) + log10(bac$cfuA+1)+1)*baradjust), label = c(cfuAclds), vjust = 1) + 
#          theme_cowplot()
#      }
#    } else {
#      if(aceto_bottom != T) {
#        bac3plot <- ggplot(bac3, aes(x=testcol.x, y=value.x, fill = variable.x)) + geom_bar(stat = "identity") + 
#          geom_errorbar(aes(ymin =addcol-value.y, ymax = addcol+value.y), width = 0.5) +
#          scale_fill_manual(values = c("red","blue"))+
#          theme(legend.position = "bottom") + 
#          geom_text(data = bac3 %>% filter(variable.x == "cfuA"), aes(x = testcol.x, y = 0, label = c(cfuAclds), vjust = 1)) +
#          geom_text(data = bac3 %>% filter(variable.x == "cfuL"), aes(x = testcol.x, y = (log10(bac$cfuL+1) + log10(bac$cfuA+1)+1)*baradjust), label = c(cfuLclds), vjust = 1) +
#          theme_cowplot()
#      } else {
#        ggplot(bac3, aes(x=testcol.x, y=value.x, fill = variable.x)) +
#          geom_bar(stat = "identity") + 
#          geom_errorbar(aes(ymin =addcol-value.y, ymax = addcol+value.y), width = 0.5) +
#          scale_fill_manual(values = c("blue","red"))+
#          theme(legend.position = "bottom") +  
#          geom_text(data = bac3 %>% filter(variable.x == "cfuL"), aes(x = testcol.x, y = 0, label = c(cfuLclds), vjust = 1)) +
#          geom_text(data = bac3 %>% filter(variable.x == "cfuA"), aes(x = testcol.x, y = (*baradjust), label = c(cfuAclds), vjust = 1) +
#          theme_cowplot()
#      }
#    }
# }
  
  
  # old
  bac4 <- bac %>% 
    mutate(total = cfuA + cfuL) %>%
    group_by(testcol) %>%
    dplyr::summarize(cfuA = mean(cfuA, na.rm=T), cfuL = mean(cfuL, na.rm=T), total = mean(total, na.rm=T)) %>%
    ungroup() %>%
    rowwise() %>%
    mutate(totalcol = cfuA+cfuL) %>%
    mutate(aperc = cfuA/totalcol, lperc = cfuL/totalcol) %>%
    ungroup() %>%
    reshape2::melt(measure.vars = c("aperc","lperc"))
  
  bac4A <- bac %>% 
    mutate(total = cfuA + cfuL) %>%
    mutate(percA = cfuA/total, percL = cfuL/total)
  
  bac4$addcol <- bac4$value
  if (aceto_bottom == T) {
    bac4$variable = factor(bac4$variable, levels = c("lperc","aperc"))
    for (i in (numvals+1):(numvals*2)) {bac4$addcol[i] = bac4$value[i-numvals]}
  } else {
    bac4$variable = factor(bac4$variable, levels = c("aperc","lperc"))
    for (i in 1:numvals) {bac4$addcol[i] = bac4$value[i+numvals]}
  }
  
  try(rm(ghi, mid_test),T)
  cat("\n\nStats on relative abundances")
  print(kruskal.test(bac$aperc ~ bac$testcol))
  mid_test <- kruskal.test(bac$aperc ~ bac$testcol)
  print(paste0("chi2 ",mid_test$parameter,", ",dim(bac)[1]," = ",round(mid_test$statistic,2),", p = ",round(mid_test$p.value,4)))
  if(kruskal.test(bac$aperc ~ bac$testcol)$p.value < 0.05) {
    if (length(table(bac$testcol))>2) {
      efg <- dunn.test(bac$aperc,bac$testcol, method = "bh", table = F, kw = F)
      ghi <- cldList(comparison = efg$comparisons, p.value = efg$P.adjusted, threshold = 0.05) 
      ghi$Group2 <- factor(ghi$Group, levels = levels(bac$testcol))
      apercclds <- as.character(ghi %>% arrange(Group2) %>% dplyr::select(Letter) %>% unlist() %>% unname())
      print(kruskal.test(bac$aperc ~ bac$testcol))
      cat("AAB Plot_values:")
      print(ghi)
    } else {
      apercclds <- c("a","b")
    }
  } else {
    apercclds <- c(rep("a",numvals))
    # print("")
    # cat("AAB Plot_values:")
    # cat(" ")
    # cat(apercclds)
  }
  
  # cat("relative abundance stats")
  # print(kruskal.test(bac$aperc ~ bac$testcol))
  # print(kruskal.test(bac$lperc ~ bac$testcol))
  # 
  
  if(!is.null(yaxis_perc)) {    
    if(aceto_bottom == F) {
      bac4g <- bac4 %>% filter(variable == "lperc") %>% rbind(bac4 %>% filter(variable == "aperc") %>% mutate(value = value - (1-yaxis_perc[2])))
      bac4plot <- ggplot(bac4g, aes(x=testcol, y=value, fill = variable)) + geom_bar(stat = "identity") +
        scale_fill_manual(values = c("red","blue"))+
        theme(legend.position = "bottom") + 
        geom_text(data = bac4 %>% filter(variable == "aperc"), aes(x = testcol, y = yaxis_perc[2]+loglabeladjust, label = c(apercclds), vjust = 0))+
        ylim(yaxis_perc*1.001)+
        theme_cowplot()
    } else {
      bac4g <- bac4 %>% filter(variable == "aperc") %>% rbind(bac4 %>% filter(variable == "lperc") %>% mutate(value = value - (1-yaxis_perc[2])))
      bac4g$variablerev = factor(bac4g$variable, levels = rev(levels(bac4g$variable)))
      bac4plot <- ggplot(bac4g, aes(x=testcol, y=value, fill = variablerev)) + geom_bar(stat = "identity") +
        scale_fill_manual(values = c("blue","red"))+
        theme(legend.position = "bottom") + 
        geom_text(data = bac4g %>% filter(variable == "aperc"), aes(x = testcol, y = yaxis_perc[2]+loglabeladjust, label = c(apercclds), vjust = 0))+
        ylim(yaxis_perc*1.001) +
        theme_cowplot()
    }
  } else {
    if(aceto_bottom == F) {
      bac4plot <- ggplot(bac4, aes(x=testcol, y=value, fill = variable)) + geom_bar(stat = "identity") +
        scale_fill_manual(values = c("red","blue"))+
        theme(legend.position = "bottom") + 
        geom_text(data = bac4 %>% filter(variable == "aperc"), aes(x = testcol, y = 1+loglabeladjust, label = c(apercclds), vjust = 0))+
        theme_cowplot()
    } else {
      bac4c <- bac4 %>% filter(variable == "lperc") %>% rbind(bac4 %>% filter(variable == "aperc"))
      bac4c$variablerev = factor(bac4c$variable, levels = rev(levels(bac4c$variable)))
      bac4plot <- ggplot(bac4c, aes(x=testcol, y=value, fill = variablerev)) + geom_bar(stat = "identity") +
        scale_fill_manual(values = c("blue","red"))+
        theme(legend.position = "bottom") + 
        geom_text(data = bac4c %>% filter(variable == "aperc"), aes(x = testcol, y = 1+loglabeladjust, label = c(apercclds), vjust = 0))+
        theme_cowplot()
    }
  }
  
  ifelse(plotout == 1, 
         return(bac3plot),
         ifelse(plotout == 2, 
                return(bac4plot),
                ifelse(plotout == 3, 
                       return(grid.arrange(bac3plot, bac4plot, heights = c(1,1))))))
}



# bac = optimize1 %>% filter(!is.na(cfuA) & !is.na(cfuL) & starved == "N")
# column = "Sex"
# plotout = 2
# aceto_bottom = F
# logplot = T
# yaxis_perc = c(0,0.005)
# loglabeladjust = 10.05
#   

make_chart_stats_96_jitter_1median <- function(bac, column, plotout, yaxis = NULL, factorspecs = NULL, logplot = F, label_adjust = 10000, baradjust = 0, baradjust_bottom = 0, loglabeladjust = 0, aceto_bottom = F, yaxis_perc = NULL, perc_errorbars = NULL) {
  bac$testcol <- unlist(unname(bac[,column]))
  bac <- bac %>%
    mutate(aperc = cfuA/(cfuA+cfuL), lperc = cfuL/(cfuA+cfuL)) %>%
    droplevels()
  if(!is.null(factorspecs)) {
    bac$testcol = factor(bac$testcol, levels = factorspecs)
  } else {
    bac$testcol = factor(bac$testcol)
    print(levels(factor(bac$testcol)))
  }
  
  if(logplot == F) {
    bac3 <- bac %>% 
      mutate(total = cfuA + cfuL) %>%
      group_by(testcol) %>%
      dplyr::summarize(cfuA = median(cfuA, na.rm=T), cfuL = median(cfuL, na.rm=T), total = median(total, na.rm=T)) %>%
      ungroup() %>%
      reshape2::melt(measure.vars = c("cfuA","cfuL")) %>%
      mutate(tempvar = paste0(testcol,"_",variable)) %>%
      inner_join(bac %>% 
                   mutate(total = cfuA + cfuL) %>%
                   group_by(testcol) %>%
                   dplyr::summarize(semcfuA = sd(cfuA, na.rm=T)/sqrt(dplyr::n()), semcfuL = sd(cfuL, na.rm=T)/sqrt(dplyr::n())) %>%
                   ungroup() %>%
                   reshape2::melt(measure.vars = c("semcfuA","semcfuL"))  %>%
                   mutate(mergecol = gsub(pattern = "sem", replacement = "", x = variable)) %>%
                   mutate(tempvar = paste0(testcol,"_",mergecol)), by = "tempvar") 
  } else {
    if(aceto_bottom==F) {
      bac3 <- 1
      bac3 <- bac %>% 
        group_by(testcol) %>%
        dplyr::summarize(cfuLmedian = median(log10(cfuL+1), na.rm=T),  totalmedian = median(log10(cfuA+cfuL+1), na.rm=T)) %>%
        ungroup() %>%
        mutate(cfuAmedian = totalmedian-cfuLmedian) %>%
        dplyr::select(testcol, cfuL = cfuLmedian, total = totalmedian, cfuA = cfuAmedian) %>%
        reshape2::melt(measure.vars = c("cfuA","cfuL")) %>%
        mutate(tempvar = paste0(testcol,"_",variable)) %>%
        inner_join(bac %>% 
                     group_by(testcol) %>%
                     dplyr::summarize(semcfuA = sd(log10(cfuA+1), na.rm=T)/sqrt(dplyr::n()), semcfuL = sd(log10(cfuL+1), na.rm=T)/sqrt(dplyr::n())) %>%
                     ungroup() %>%
                     reshape2::melt(measure.vars = c("semcfuA","semcfuL"))  %>%
                     mutate(mergecol = gsub(pattern = "sem", replacement = "", x = variable)) %>%
                     mutate(tempvar = paste0(testcol,"_",mergecol)), by = "tempvar") 
    } else {
      print("hi")
      bac3 <- bac %>% 
        group_by(testcol) %>%
        dplyr::summarize(cfuAmedian = median(log10(cfuA+1), na.rm=T),  totalmedian = median(log10(cfuA+cfuL+1), na.rm=T)) %>%
        ungroup() %>%
        mutate(cfuLmedian = totalmedian-cfuAmedian) %>%
        dplyr::select(testcol, cfuL = cfuLmedian, total = totalmedian, cfuA = cfuAmedian) %>%
        reshape2::melt(measure.vars = c("cfuA","cfuL")) %>%
        mutate(tempvar = paste0(testcol,"_",variable)) %>%
        inner_join(bac %>% 
                     group_by(testcol) %>%
                     dplyr::summarize(semcfuA = sd(log10(cfuA+1), na.rm=T)/sqrt(dplyr::n()), semcfuL = sd(log10(cfuL+1), na.rm=T)/sqrt(dplyr::n())) %>%
                     ungroup() %>%
                     reshape2::melt(measure.vars = c("semcfuA","semcfuL"))  %>%
                     mutate(mergecol = gsub(pattern = "sem", replacement = "", x = variable)) %>%
                     mutate(tempvar = paste0(testcol,"_",mergecol)), by = "tempvar") 
    }
  }
  
  bac3$value.y = ifelse(is.na(bac3$value.y),0,bac3$value.y) ## added 2023-08-12 - maybe remove if things get broken
  
  print(paste0("N:",dim(bac)[1]))
  bac3$addcol <- bac3$value.x
  numvals <- length(table(list(bac$testcol %>% droplevels)))
  if (aceto_bottom == T) {
    bac3$variable.x = factor(bac3$variable.x, levels = c("cfuL","cfuA"))
    for (i in (numvals+1):(numvals*2)) {bac3$addcol[i] = bac3$total[i]}
  } else {
    bac3$variable.x = factor(bac3$variable.x, levels = c("cfuA","cfuL"))
    for (i in 1:numvals) {bac3$addcol[i] = bac3$total[i]}
  }
  
  
  try(rm(ghi, mid_test),T)
  cat("Stats on AAB abundances")
  print(kruskal.test(bac$cfuA ~ bac$testcol))
  mid_test <- kruskal.test(bac$cfuA ~ bac$testcol)
  print(paste0("chi2 ",mid_test$parameter,", ",dim(bac)[1]," = ",round(mid_test$statistic,2),", p = ",round(mid_test$p.value,4)))
  if(kruskal.test(bac$cfuA ~ bac$testcol)$p.value < 0.05) {
    if (length(table(bac$testcol))>2) {
      efg <- dunn.test(bac$cfuA,bac$testcol, method = "bh", table = F, kw = T)
      ghi <- cldList(comparison = efg$comparisons, p.value = efg$P.adjusted, threshold = 0.05)
      ghi$Group2 <- factor(ghi$Group, levels = levels(bac$testcol))
      cfuAclds <- as.character(ghi %>% arrange(Group2) %>% dplyr::select(Letter) %>% unlist() %>% unname())
      print(ghi)
    } else {
      cfuAclds <- c("a","b")
    }
    cat("AAB Plot_values:")
  } else {
    cfuAclds <- c(rep("a", numvals))
  }
  
  try(rm(ghi, mid_test),T)
  cat("\n\nStats on LAB abundances")
  print(kruskal.test(bac$cfuL ~ bac$testcol))
  mid_test <- kruskal.test(bac$cfuL ~ bac$testcol)
  print(paste0("chi2 ",mid_test$parameter,", ",dim(bac)[1]," = ",round(mid_test$statistic,2),", p = ",round(mid_test$p.value,4)))
  if(kruskal.test(bac$cfuL ~ bac$testcol)$p.value < 0.05) {
    if (length(table(bac$testcol))>2) {
      efg <- dunn.test(bac$cfuL,bac$testcol, method = "bh", table = F, kw = F)
      ghi <- cldList(comparison = efg$comparisons, p.value = efg$P.adjusted, threshold = 0.05)
      ghi$Group2 <- factor(ghi$Group, levels = levels(bac$testcol))
      cfuLclds <- as.character(ghi %>% arrange(Group2) %>% dplyr::select(Letter) %>% unlist() %>% unname())
      print(ghi)
    } else {
      cfuLclds <- c("a","b")
    }
    cat("LAB Plot_values:")
  } else {
    cfuLclds <- c(rep("a", numvals))
  }
  
  print(bac3)
  if(aceto_bottom != T) {
    bac3plot <-1 
    bac3plot <- ggplot(bac3, aes(x=testcol.x, y=value.x, fill = variable.x)) + geom_bar(stat = "identity") + 
      geom_errorbar(aes(ymin =addcol-value.y, ymax = addcol+value.y), width = 0.5) +
      scale_fill_manual(values = c("red","blue"))+
      theme(legend.position = "bottom") + 
      geom_text(data = bac3 %>% filter(variable.x == "cfuA"), aes(x = testcol.x, y = ((total+value.y)+baradjust), label = c(cfuAclds), vjust = 1)) + #, inherit.aes = T)) +
      geom_text(data = bac3 %>% filter(variable.x == "cfuL"), aes(x = testcol.x, y = 0+baradjust_bottom, label = c(cfuLclds), vjust = 1)) + #, inherit.aes = T)+
      theme_cowplot()
  } else {
    bac3plot <- ggplot(bac3, aes(x=testcol.x, y=value.x, fill = variable.x)) +
      geom_bar(stat = "identity") + 
      geom_errorbar(aes(ymin =addcol-value.y, ymax = addcol+value.y), width = 0.5) +
      scale_fill_manual(values = c("blue","red"))+
      theme(legend.position = "bottom") +  
      geom_text(data = bac3 %>% filter(variable.x == "cfuL"), aes(x = testcol.x, y = ((total+value.y)+baradjust), label = c(cfuLclds), vjust = 1)) +
      geom_text(data = bac3 %>% filter(variable.x == "cfuA"), aes(x = testcol.x, y = 0+baradjust_bottom, label = c(cfuAclds), vjust = 1)) + 
      theme_cowplot()
  }
  
  if(!is.null(yaxis)) {    
    bac3plot <- bac3plot + ylim(yaxis)
  }  
  
  # bac4 <- bac %>% 
  #   mutate(total = cfuA + cfuL) %>%
  #   group_by(testcol) %>%
  #   dplyr::summarize(cfuA = median(cfuA, na.rm=T), cfuL = median(cfuL, na.rm=T), total = median(total, na.rm=T)) %>%
  #   ungroup() %>%
  #   rowwise() %>%
  #   mutate(totalcol = cfuA+cfuL) %>%
  #   mutate(aperc = cfuA/totalcol, lperc = cfuL/totalcol) %>%
  #   ungroup() %>%
  #   reshape2::melt(measure.vars = c("aperc","lperc"))
  # 
  bac4 <- bac %>% 
    mutate(total = cfuA + cfuL) %>%
    mutate(percA = cfuA/total, percL = cfuL/total) %>%
    group_by(testcol) %>%
    dplyr::summarize(medianA = median(percA, na.rm=T), medianL = median(percL, na.rm=T), semA = sd(percA)/sqrt(dplyr::n()), semL = sd(percL)/sqrt(dplyr::n())) %>%
    ungroup() %>%
    mutate(aperc = medianA, lperc = medianL) %>%
    reshape2::melt(measure.vars = c("aperc","lperc"))
  
  if (aceto_bottom == T) {
    bac4 <- bac4 %>%
      rowwise() %>%
      mutate(addcol = ifelse(variable == "aperc", medianA, 1)) %>%
      ungroup() 
  } else {
    bac4 <- bac4 %>%
      rowwise() %>%
      mutate(addcol = ifelse(variable == "aperc",1,medianL)) %>%
      ungroup() 
  }
  
  try(rm(ghi, mid_test),T)
  cat("\n\nStats on relative abundances")
  print(kruskal.test(bac$aperc ~ bac$testcol))
  mid_test <- kruskal.test(bac$aperc ~ bac$testcol)
  print(paste0("chi2 ",mid_test$parameter,", ",dim(bac)[1]," = ",round(mid_test$statistic,2),", p = ",round(mid_test$p.value,4)))
  if(kruskal.test(bac$aperc ~ bac$testcol)$p.value < 0.05) {
    if (length(table(bac$testcol))>2) {
      efg <- dunn.test(bac$aperc,bac$testcol, method = "bh", table = F, kw = F)
      ghi <- cldList(comparison = efg$comparisons, p.value = efg$P.adjusted, threshold = 0.05) 
      ghi$Group2 <- factor(ghi$Group, levels = levels(bac$testcol))
      apercclds <- as.character(ghi %>% arrange(Group2) %>% dplyr::select(Letter) %>% unlist() %>% unname())
      print(kruskal.test(bac$aperc ~ bac$testcol))
      cat("AAB Plot_values:")
      print(ghi)
    } else {
      apercclds <- c("a","b")
    }
  } else {
    apercclds <- c(rep("a",numvals))
  }
  
  print(bac4)
  if(!is.null(yaxis_perc)) {    
    if(aceto_bottom == F) {
      bac4g <- bac4 %>% filter(variable == "lperc") %>% rbind(bac4 %>% filter(variable == "aperc") %>% mutate(value = value - (1-yaxis_perc[2])))
      bac4plot <- ggplot(bac4g, aes(x=testcol, y=value, fill = variable)) + geom_bar(stat = "identity") +
        scale_fill_manual(values = c("red","blue"))+
        theme(legend.position = "bottom") + 
        geom_text(data = bac4 %>% filter(variable == "aperc"), aes(x = testcol, y = yaxis_perc[2]+loglabeladjust, label = c(apercclds), vjust = 0))+
        ylim(yaxis_perc*1.001)+
        theme_cowplot()
    } else {
      bac4g <- bac4 %>% filter(variable == "aperc") %>% rbind(bac4 %>% filter(variable == "lperc") %>% mutate(value = value - (1-yaxis_perc[2])))
      bac4g$variablerev = factor(bac4g$variable, levels = rev(levels(bac4g$variable)))
      bac4plot <- ggplot(bac4g, aes(x=testcol, y=value, fill = variablerev)) + geom_bar(stat = "identity") +
        scale_fill_manual(values = c("blue","red"))+
        theme(legend.position = "bottom") + 
        geom_text(data = bac4g %>% filter(variable == "aperc"), aes(x = testcol, y = yaxis_perc[2]+loglabeladjust, label = c(apercclds), vjust = 0))+
        ylim(yaxis_perc*1.001) +
        theme_cowplot()
    }
  } else {
    if(aceto_bottom == F) {
      bac4plot <- ggplot(bac4, aes(x=testcol, y=value, fill = variable)) + geom_bar(stat = "identity") +
        scale_fill_manual(values = c("red","blue"))+
        theme(legend.position = "bottom") + 
        geom_text(data = bac4 %>% filter(variable == "aperc"), aes(x = testcol, y = 1+loglabeladjust, label = c(apercclds), vjust = 0))+
        #       geom_errorbar(aes(ymin = addcol - semA , ymax = addcol+semA), width = 0.5) +
        theme_cowplot()
    } else {
      bac4c <- bac4 %>% filter(variable == "lperc") %>% rbind(bac4 %>% filter(variable == "aperc"))
      bac4c$variablerev = factor(bac4c$variable, levels = rev(levels(bac4c$variable)))
      bac4plot <- ggplot(bac4c, aes(x=testcol, y=value, fill = variablerev)) + geom_bar(stat = "identity") +
        scale_fill_manual(values = c("blue","red"))+
        theme(legend.position = "bottom") + 
        geom_text(data = bac4c %>% filter(variable == "aperc"), aes(x = testcol, y = 1+loglabeladjust, label = c(apercclds), vjust = 0))+
        #        geom_errorbar(aes(ymin = addcol - semA , ymax = addcol+semA), width = 0.5) +
        theme_cowplot()
    }
  }
  
  if(!is.null(perc_errorbars)) {    
    bac4plot <- bac4plot + geom_errorbar(aes(ymin = addcol - semA , ymax = addcol+semA), width = 0.5)
  }
  
  ifelse(plotout == 1, 
         return(bac3plot),
         ifelse(plotout == 2, 
                return(bac4plot),
                ifelse(plotout == 3, 
                       return(grid.arrange(bac3plot, bac4plot, heights = c(1,1))))))
}

## change name to the version 7 name to use on version 6
make_chart_stats_96_jitter_1mean <- function(bac, column, plotout, yaxis = NULL, factorspecs = NULL, logplot = F, label_adjust = 10000, baradjust = 0, baradjust_bottom = 0, loglabeladjust = 0, aceto_bottom = F, yaxis_perc = NULL, perc_errorbars = NULL) {
  bac$testcol <- unlist(unname(bac[,column]))
  bac <- bac %>%
    mutate(aperc = cfuA/(cfuA+cfuL), lperc = cfuL/(cfuA+cfuL)) %>%
    droplevels()
  if(!is.null(factorspecs)) {
    bac$testcol = factor(bac$testcol, levels = factorspecs)
  } else {
    bac$testcol = factor(bac$testcol)
    print(levels(factor(bac$testcol)))
  }
  
  if(logplot == F) {
    bac3 <- bac %>% 
      mutate(total = cfuA + cfuL) %>%
      group_by(testcol) %>%
      dplyr::summarize(cfuA = mean(cfuA, na.rm=T), cfuL = mean(cfuL, na.rm=T), total = mean(total, na.rm=T)) %>%
      ungroup() %>%
      reshape2::melt(measure.vars = c("cfuA","cfuL")) %>%
      mutate(tempvar = paste0(testcol,"_",variable)) %>%
      inner_join(bac %>% 
                   mutate(total = cfuA + cfuL) %>%
                   group_by(testcol) %>%
                   dplyr::summarize(semcfuA = sd(cfuA, na.rm=T)/sqrt(dplyr::n()), semcfuL = sd(cfuL, na.rm=T)/sqrt(dplyr::n())) %>%
                   ungroup() %>%
                   reshape2::melt(measure.vars = c("semcfuA","semcfuL"))  %>%
                   mutate(mergecol = gsub(pattern = "sem", replacement = "", x = variable)) %>%
                   mutate(tempvar = paste0(testcol,"_",mergecol)), by = "tempvar") 
  } else {
    if(aceto_bottom==F) {
      bac3 <- 1
      bac3 <- bac %>% 
        group_by(testcol) %>%
        dplyr::summarize(cfuLmean = mean(log10(cfuL+1), na.rm=T),  totalmean = mean(log10(cfuA+cfuL+1), na.rm=T)) %>%
        ungroup() %>%
        mutate(cfuAmean = totalmean-cfuLmean) %>%
        dplyr::select(testcol, cfuL = cfuLmean, total = totalmean, cfuA = cfuAmean) %>%
        reshape2::melt(measure.vars = c("cfuA","cfuL")) %>%
        mutate(tempvar = paste0(testcol,"_",variable)) %>%
        inner_join(bac %>% 
                     group_by(testcol) %>%
                     dplyr::summarize(semcfuA = sd(log10(cfuA+1), na.rm=T)/sqrt(dplyr::n()), semcfuL = sd(log10(cfuL+1), na.rm=T)/sqrt(dplyr::n())) %>%
                     ungroup() %>%
                     reshape2::melt(measure.vars = c("semcfuA","semcfuL"))  %>%
                     mutate(mergecol = gsub(pattern = "sem", replacement = "", x = variable)) %>%
                     mutate(tempvar = paste0(testcol,"_",mergecol)), by = "tempvar") 
    } else {
      print("hi")
      bac3 <- bac %>% 
        group_by(testcol) %>%
        dplyr::summarize(cfuAmean = mean(log10(cfuA+1), na.rm=T),  totalmean = mean(log10(cfuA+cfuL+1), na.rm=T)) %>%
        ungroup() %>%
        mutate(cfuLmean = totalmean-cfuAmean) %>%
        dplyr::select(testcol, cfuL = cfuLmean, total = totalmean, cfuA = cfuAmean) %>%
        reshape2::melt(measure.vars = c("cfuA","cfuL")) %>%
        mutate(tempvar = paste0(testcol,"_",variable)) %>%
        inner_join(bac %>% 
                     group_by(testcol) %>%
                     dplyr::summarize(semcfuA = sd(log10(cfuA+1), na.rm=T)/sqrt(dplyr::n()), semcfuL = sd(log10(cfuL+1), na.rm=T)/sqrt(dplyr::n())) %>%
                     ungroup() %>%
                     reshape2::melt(measure.vars = c("semcfuA","semcfuL"))  %>%
                     mutate(mergecol = gsub(pattern = "sem", replacement = "", x = variable)) %>%
                     mutate(tempvar = paste0(testcol,"_",mergecol)), by = "tempvar") 
    }
  }
  
  bac3$value.y = ifelse(is.na(bac3$value.y),0,bac3$value.y) ## added 2023-08-12 - maybe remove if things get broken
  
  print(paste0("N:",dim(bac)[1]))
  bac3$addcol <- bac3$value.x
  numvals <- length(table(list(bac$testcol %>% droplevels)))
  if (aceto_bottom == T) {
    bac3$variable.x = factor(bac3$variable.x, levels = c("cfuL","cfuA"))
    for (i in (numvals+1):(numvals*2)) {bac3$addcol[i] = bac3$total[i]}
  } else {
    bac3$variable.x = factor(bac3$variable.x, levels = c("cfuA","cfuL"))
    for (i in 1:numvals) {bac3$addcol[i] = bac3$total[i]}
  }
  
  
  try(rm(ghi, mid_test),T)
  cat("Stats on AAB abundances")
  print(kruskal.test(bac$cfuA ~ bac$testcol))
  mid_test <- kruskal.test(bac$cfuA ~ bac$testcol)
  print(paste0("chi2 ",mid_test$parameter,", ",dim(bac)[1]," = ",round(mid_test$statistic,2),", p = ",round(mid_test$p.value,4)))
  if(kruskal.test(bac$cfuA ~ bac$testcol)$p.value < 0.05) {
    if (length(table(bac$testcol))>2) {
      efg <- dunn.test(bac$cfuA,bac$testcol, method = "bh", table = F, kw = T)
      ghi <- cldList(comparison = efg$comparisons, p.value = efg$P.adjusted, threshold = 0.05)
      ghi$Group2 <- factor(ghi$Group, levels = levels(bac$testcol))
      cfuAclds <- as.character(ghi %>% arrange(Group2) %>% dplyr::select(Letter) %>% unlist() %>% unname())
      print(ghi)
    } else {
      cfuAclds <- c("a","b")
    }
    cat("AAB Plot_values:")
  } else {
    cfuAclds <- c(rep("a", numvals))
  }
  
  try(rm(ghi, mid_test),T)
  cat("\n\nStats on LAB abundances")
  print(kruskal.test(bac$cfuL ~ bac$testcol))
  mid_test <- kruskal.test(bac$cfuL ~ bac$testcol)
  print(paste0("chi2 ",mid_test$parameter,", ",dim(bac)[1]," = ",round(mid_test$statistic,2),", p = ",round(mid_test$p.value,4)))
  if(kruskal.test(bac$cfuL ~ bac$testcol)$p.value < 0.05) {
    if (length(table(bac$testcol))>2) {
      efg <- dunn.test(bac$cfuL,bac$testcol, method = "bh", table = F, kw = F)
      ghi <- cldList(comparison = efg$comparisons, p.value = efg$P.adjusted, threshold = 0.05)
      ghi$Group2 <- factor(ghi$Group, levels = levels(bac$testcol))
      cfuLclds <- as.character(ghi %>% arrange(Group2) %>% dplyr::select(Letter) %>% unlist() %>% unname())
      print(ghi)
    } else {
      cfuLclds <- c("a","b")
    }
    cat("LAB Plot_values:")
  } else {
    cfuLclds <- c(rep("a", numvals))
  }
  
  print(bac3)
  if(aceto_bottom != T) {
    bac3plot <-1 
    bac3plot <- ggplot(bac3, aes(x=testcol.x, y=value.x, fill = variable.x)) + geom_bar(stat = "identity") + 
      geom_errorbar(aes(ymin =addcol-value.y, ymax = addcol+value.y), width = 0.5) +
      scale_fill_manual(values = c("red","blue"))+
      theme(legend.position = "bottom") + 
      geom_text(data = bac3 %>% filter(variable.x == "cfuA"), aes(x = testcol.x, y = ((total+value.y)+baradjust), label = c(cfuAclds), vjust = 1)) + #, inherit.aes = T)) +
      geom_text(data = bac3 %>% filter(variable.x == "cfuL"), aes(x = testcol.x, y = 0+baradjust_bottom, label = c(cfuLclds), vjust = 1)) + #, inherit.aes = T)+
      theme_cowplot()
  } else {
    bac3plot <- ggplot(bac3, aes(x=testcol.x, y=value.x, fill = variable.x)) +
      geom_bar(stat = "identity") + 
      geom_errorbar(aes(ymin =addcol-value.y, ymax = addcol+value.y), width = 0.5) +
      scale_fill_manual(values = c("blue","red"))+
      theme(legend.position = "bottom") +  
      geom_text(data = bac3 %>% filter(variable.x == "cfuL"), aes(x = testcol.x, y = ((total+value.y)+baradjust), label = c(cfuLclds), vjust = 1)) +
      geom_text(data = bac3 %>% filter(variable.x == "cfuA"), aes(x = testcol.x, y = 0+baradjust_bottom, label = c(cfuAclds), vjust = 1)) + 
      theme_cowplot()
  }
  
  if(!is.null(yaxis)) {    
    bac3plot <- bac3plot + ylim(yaxis)
  }  

   bac4 <- bac %>% 
     mutate(total = cfuA + cfuL) %>%
     mutate(percA = cfuA/total, percL = cfuL/total) %>%
     group_by(testcol) %>%
     dplyr::summarize(meanA = mean(percA, na.rm=T), meanL = mean(percL, na.rm=T), semA = sd(percA)/sqrt(dplyr::n()), semL = sd(percL)/sqrt(dplyr::n())) %>%
     ungroup() %>%
     mutate(aperc = meanA, lperc = meanL) %>%
     reshape2::melt(measure.vars = c("aperc","lperc"))
   
  if (aceto_bottom == T) {
    bac4 <- bac4 %>%
      rowwise() %>%
      mutate(addcol = ifelse(variable == "aperc", meanA, 1)) %>%
      ungroup() 
  } else {
    bac4 <- bac4 %>%
      rowwise() %>%
      mutate(addcol = ifelse(variable == "aperc",1,meanL)) %>%
      ungroup() 
  }
 
try(rm(ghi, mid_test),T)
  cat("\n\nStats on relative abundances")
  print(kruskal.test(bac$aperc ~ bac$testcol))
  mid_test <- kruskal.test(bac$aperc ~ bac$testcol)
  print(paste0("chi2 ",mid_test$parameter,", ",dim(bac)[1]," = ",round(mid_test$statistic,2),", p = ",round(mid_test$p.value,4)))
  if(kruskal.test(bac$aperc ~ bac$testcol)$p.value < 0.05) {
    if (length(table(bac$testcol))>2) {
      efg <- dunn.test(bac$aperc,bac$testcol, method = "bh", table = F, kw = F)
      ghi <- cldList(comparison = efg$comparisons, p.value = efg$P.adjusted, threshold = 0.05) 
      ghi$Group2 <- factor(ghi$Group, levels = levels(bac$testcol))
      apercclds <- as.character(ghi %>% arrange(Group2) %>% dplyr::select(Letter) %>% unlist() %>% unname())
      print(kruskal.test(bac$aperc ~ bac$testcol))
      cat("AAB Plot_values:")
      print(ghi)
    } else {
      apercclds <- c("a","b")
    }
  } else {
    apercclds <- c(rep("a",numvals))
  }
  
  print(bac4)
  if(!is.null(yaxis_perc)) {    
    if(aceto_bottom == F) {
      bac4g <- bac4 %>% filter(variable == "lperc") %>% rbind(bac4 %>% filter(variable == "aperc") %>% mutate(value = value - (1-yaxis_perc[2])))
      bac4plot <- ggplot(bac4g, aes(x=testcol, y=value, fill = variable)) + geom_bar(stat = "identity") +
        scale_fill_manual(values = c("red","blue"))+
        theme(legend.position = "bottom") + 
        geom_text(data = bac4 %>% filter(variable == "aperc"), aes(x = testcol, y = yaxis_perc[2]+loglabeladjust, label = c(apercclds), vjust = 0))+
        ylim(yaxis_perc*1.001)+
        theme_cowplot()
      print(1)
    } else {
      bac4g <- bac4 %>% filter(variable == "aperc") %>% rbind(bac4 %>% filter(variable == "lperc") %>% mutate(value = value - (1-yaxis_perc[2])))
      bac4g$variablerev = factor(bac4g$variable, levels = rev(levels(bac4g$variable)))
      print(2)
      bac4plot <- ggplot(bac4g, aes(x=testcol, y=value, fill = variablerev)) + geom_bar(stat = "identity") +
        scale_fill_manual(values = c("blue","red"))+
        theme(legend.position = "bottom") + 
        geom_text(data = bac4g %>% filter(variable == "aperc"), aes(x = testcol, y = yaxis_perc[2]+loglabeladjust, label = c(apercclds), vjust = 0))+
        ylim(yaxis_perc*1.001) +
        theme_cowplot()
    }
  } else {
    if(aceto_bottom == F) {
      print(3)
      bac4plot <- ggplot(bac4, aes(x=testcol, y=value, fill = variable)) + geom_bar(stat = "identity") +
        scale_fill_manual(values = c("red","blue"))+
        theme(legend.position = "bottom") + 
        geom_text(data = bac4 %>% filter(variable == "aperc"), aes(x = testcol, y = 1+loglabeladjust, label = c(apercclds), vjust = 0))+
 #       geom_errorbar(aes(ymin = addcol - semA , ymax = addcol+semA), width = 0.5) +
        theme_cowplot()
    } else {
      print(4)
      bac4c <- bac4 %>% filter(variable == "lperc") %>% rbind(bac4 %>% filter(variable == "aperc"))
      bac4c$variablerev = factor(bac4c$variable, levels = rev(levels(bac4c$variable)))
      bac4plot <- ggplot(bac4c, aes(x=testcol, y=value, fill = variablerev)) + geom_bar(stat = "identity") +
        scale_fill_manual(values = c("blue","red"))+
        theme(legend.position = "bottom") + 
        geom_text(data = bac4c %>% filter(variable == "aperc"), aes(x = testcol, y = 1+loglabeladjust, label = c(apercclds), vjust = 0))+
#        geom_errorbar(aes(ymin = addcol - semA , ymax = addcol+semA), width = 0.5) +
        theme_cowplot()
    }
  }
  
  if(!is.null(perc_errorbars)) {    
    bac4plot <- bac4plot + geom_errorbar(aes(ymin = addcol - semA , ymax = addcol+semA), width = 0.5)
  }
  
  ifelse(plotout == 1, 
         return(bac3plot),
         ifelse(plotout == 2, 
                return(bac4plot),
                ifelse(plotout == 3, 
                       return(grid.arrange(bac3plot, bac4plot, heights = c(1,1))))))
}

  




make_chart_stats_96_bar_jitter_violin <- function(bac, column, plotout, yaxis = NULL, factorspecs = NULL, logplot = F, label_adjust = 10000, baradjust = 0, baradjust_bottom = 0, loglabeladjust = 0, aceto_bottom = F, yaxis_perc = NULL, perc_errorbars = NULL, add_jitter = NULL, add_violin = NULL) {

bac$testcol <- unlist(unname(bac[,column]))
bac <- bac %>%
  mutate(aperc = cfuA/(cfuA+cfuL), lperc = cfuL/(cfuA+cfuL)) %>%
  droplevels()
if(!is.null(factorspecs)) {
  bac$testcol = factor(bac$testcol, levels = factorspecs)
} else {
  bac$testcol = factor(bac$testcol)
  print(levels(factor(bac$testcol)))
}

if(logplot == F) {
  bac3 <- bac %>% 
    mutate(total = cfuA + cfuL) %>%
    group_by(testcol) %>%
    dplyr::summarize(cfuA = mean(cfuA, na.rm=T), cfuL = mean(cfuL, na.rm=T), total = mean(total, na.rm=T)) %>%
    ungroup() %>%
    reshape2::melt(measure.vars = c("cfuA","cfuL")) %>%
    mutate(tempvar = paste0(testcol,"_",variable)) %>%
    inner_join(bac %>% 
                 mutate(total = cfuA + cfuL) %>%
                 group_by(testcol) %>%
                 dplyr::summarize(semcfuA = sd(cfuA, na.rm=T)/sqrt(dplyr::n()), semcfuL = sd(cfuL, na.rm=T)/sqrt(dplyr::n())) %>%
                 ungroup() %>%
                 reshape2::melt(measure.vars = c("semcfuA","semcfuL"))  %>%
                 mutate(mergecol = gsub(pattern = "sem", replacement = "", x = variable)) %>%
                 mutate(tempvar = paste0(testcol,"_",mergecol)), by = "tempvar") 
} else {
  if(aceto_bottom==F) {
    bac3 <- 1
    bac3 <- bac %>% 
      group_by(testcol) %>%
      dplyr::summarize(cfuLmean = mean(log10(cfuL+1), na.rm=T),  totalmean = mean(log10(cfuA+cfuL+1), na.rm=T)) %>%
      ungroup() %>%
      mutate(cfuAmean = totalmean-cfuLmean) %>%
      dplyr::select(testcol, cfuL = cfuLmean, total = totalmean, cfuA = cfuAmean) %>%
      reshape2::melt(measure.vars = c("cfuA","cfuL")) %>%
      mutate(tempvar = paste0(testcol,"_",variable)) %>%
      inner_join(bac %>% 
                   group_by(testcol) %>%
                   dplyr::summarize(semcfuA = sd(log10(cfuA+1), na.rm=T)/sqrt(dplyr::n()), semcfuL = sd(log10(cfuL+1), na.rm=T)/sqrt(dplyr::n())) %>%
                   ungroup() %>%
                   reshape2::melt(measure.vars = c("semcfuA","semcfuL"))  %>%
                   mutate(mergecol = gsub(pattern = "sem", replacement = "", x = variable)) %>%
                   mutate(tempvar = paste0(testcol,"_",mergecol)), by = "tempvar") 
  } else {
    print("hi")
    bac3 <- bac %>% 
      group_by(testcol) %>%
      dplyr::summarize(cfuAmean = mean(log10(cfuA+1), na.rm=T),  totalmean = mean(log10(cfuA+cfuL+1), na.rm=T)) %>%
      ungroup() %>%
      mutate(cfuLmean = totalmean-cfuAmean) %>%
      dplyr::select(testcol, cfuL = cfuLmean, total = totalmean, cfuA = cfuAmean) %>%
      reshape2::melt(measure.vars = c("cfuA","cfuL")) %>%
      mutate(tempvar = paste0(testcol,"_",variable)) %>%
      inner_join(bac %>% 
                   group_by(testcol) %>%
                   dplyr::summarize(semcfuA = sd(log10(cfuA+1), na.rm=T)/sqrt(dplyr::n()), semcfuL = sd(log10(cfuL+1), na.rm=T)/sqrt(dplyr::n())) %>%
                   ungroup() %>%
                   reshape2::melt(measure.vars = c("semcfuA","semcfuL"))  %>%
                   mutate(mergecol = gsub(pattern = "sem", replacement = "", x = variable)) %>%
                   mutate(tempvar = paste0(testcol,"_",mergecol)), by = "tempvar") 
  }
}

bac3$value.y = ifelse(is.na(bac3$value.y),0,bac3$value.y) ## added 2023-08-12 - maybe remove if things get broken

print(paste0("N:",dim(bac)[1]))
bac3$addcol <- bac3$value.x
numvals <- length(table(list(bac$testcol %>% droplevels)))
if (aceto_bottom == T) {
  bac3$variable.x = factor(bac3$variable.x, levels = c("cfuL","cfuA"))
  for (i in (numvals+1):(numvals*2)) {bac3$addcol[i] = bac3$total[i]}
} else {
  bac3$variable.x = factor(bac3$variable.x, levels = c("cfuA","cfuL"))
  for (i in 1:numvals) {bac3$addcol[i] = bac3$total[i]}
}


try(rm(ghi, mid_test),T)
cat("Stats on AAB abundances")
#print(kruskal.test(bac$cfuA ~ bac$testcol))
mid_test <- kruskal.test(bac$cfuA ~ bac$testcol)
print(paste0("chi2 ",mid_test$parameter,", ",dim(bac)[1]," = ",round(mid_test$statistic,2),", p = ",round(mid_test$p.value,4)))
if(kruskal.test(bac$cfuA ~ bac$testcol)$p.value < 0.05) {
  if (length(table(bac$testcol))>2) {
    efg <- dunn.test(bac$cfuA,bac$testcol, method = "bh", table = F, kw = T)
    ghi <- cldList(comparison = efg$comparisons, p.value = efg$P.adjusted, threshold = 0.05)
    ghi$Group2 <- factor(ghi$Group, levels = levels(bac$testcol))
    cfuAclds <- as.character(ghi %>% arrange(Group2) %>% dplyr::select(Letter) %>% unlist() %>% unname())
  #  print(ghi)
  } else {
    cfuAclds <- c("a","b")
  }
 # cat("AAB Plot_values:")
} else {
  cfuAclds <- c(rep("a", numvals))
}

try(rm(ghi, mid_test),T)
cat("\n\nStats on LAB abundances")
#print(kruskal.test(bac$cfuL ~ bac$testcol))
mid_test <- kruskal.test(bac$cfuL ~ bac$testcol)
print(paste0("chi2 ",mid_test$parameter,", ",dim(bac)[1]," = ",round(mid_test$statistic,2),", p = ",round(mid_test$p.value,4)))
if(kruskal.test(bac$cfuL ~ bac$testcol)$p.value < 0.05) {
  if (length(table(bac$testcol))>2) {
    efg <- dunn.test(bac$cfuL,bac$testcol, method = "bh", table = F, kw = F)
    ghi <- cldList(comparison = efg$comparisons, p.value = efg$P.adjusted, threshold = 0.05)
    ghi$Group2 <- factor(ghi$Group, levels = levels(bac$testcol))
    cfuLclds <- as.character(ghi %>% arrange(Group2) %>% dplyr::select(Letter) %>% unlist() %>% unname())
#    print(ghi)
  } else {
    cfuLclds <- c("a","b")
  }
#  cat("LAB Plot_values:")
} else {
  cfuLclds <- c(rep("a", numvals))
}

#print(bac3)
if(aceto_bottom != T) {
  bac3plot <-1 
  bac3plot <- ggplot(bac3, aes(x=testcol.x, y=value.x, fill = variable.x)) + geom_bar(stat = "identity") + 
    geom_errorbar(aes(ymin =addcol-value.y, ymax = addcol+value.y), width = 0.5) +
    scale_fill_manual(values = c("red","blue"))+
    theme(legend.position = "bottom") + 
    geom_text(data = bac3 %>% filter(variable.x == "cfuA"), aes(x = testcol.x, y = ((total+value.y)+baradjust), label = c(cfuAclds), vjust = 1)) + #, inherit.aes = T)) +
    geom_text(data = bac3 %>% filter(variable.x == "cfuL"), aes(x = testcol.x, y = 0+baradjust_bottom, label = c(cfuLclds), vjust = 1)) + #, inherit.aes = T)+
    theme_cowplot()
} else {
  bac3plot <- ggplot(bac3, aes(x=testcol.x, y=value.x, fill = variable.x)) +
    geom_bar(stat = "identity") + 
    geom_errorbar(aes(ymin =addcol-value.y, ymax = addcol+value.y), width = 0.5) +
    scale_fill_manual(values = c("blue","red"))+
    theme(legend.position = "bottom") +  
    geom_text(data = bac3 %>% filter(variable.x == "cfuL"), aes(x = testcol.x, y = ((total+value.y)+baradjust), label = c(cfuLclds), vjust = 1)) +
    geom_text(data = bac3 %>% filter(variable.x == "cfuA"), aes(x = testcol.x, y = 0+baradjust_bottom, label = c(cfuAclds), vjust = 1)) + 
    theme_cowplot()
}

if(!is.null(yaxis)) {    
  bac3plot <- bac3plot + ylim(yaxis)
}  


  ## make the mean and raw data data frames
  bac4 <- bac %>% 
    mutate(total = cfuA + cfuL) %>%
    #mutate(percA = cfuA/total, percL = cfuL/total) %>%
    group_by(testcol) %>%
    dplyr::summarize(meanA = mean(cfuA, na.rm=T), meanL = mean(cfuL, na.rm=T), total = mean(total, na.rm=T)) %>%
    #    dplyr::summarize(meanA = mean(cfuA, na.rm=T), meanL = mean(percL, na.rm=T), semA = sd(percA)/sqrt(dplyr::n()), semL = sd(percL)/sqrt(dplyr::n())) %>%
    ungroup() %>%
    rowwise() %>%
    mutate(totalcol = meanA + meanL) %>%
    mutate(aperc = meanA/totalcol, lperc = meanL/totalcol) %>%
    ungroup() %>%
    reshape2::melt(measure.vars = c("aperc","lperc")) 
  bac4$addcol <- bac4$value
  if (aceto_bottom == T) {
    bac4$variable = factor(bac4$variable, levels = c("lperc","aperc"))
    for (i in (numvals+1):(numvals*2)) {bac4$addcol[i] = bac4$value[i-numvals]}
  } else {
    bac4$variable = factor(bac4$variable, levels = c("aperc","lperc"))
    for (i in 1:numvals) {bac4$addcol[i] = bac4$value[i+numvals]}
  }
  
  bac4A <- bac %>% 
    mutate(total = cfuA + cfuL) %>%
    mutate(percA = cfuA/total, percL = cfuL/total) %>%
    reshape2::melt(measure.vars = c("aperc","lperc"))
  bac4A$addcol <- bac4A$value
  if (aceto_bottom == T) {
    bac4A$variable = factor(bac4A$variable, levels = c("lperc","aperc"))
    for (i in (numvals+1):(numvals*2)) {bac4A$addcol[i] = bac4A$value[i-numvals]}
  } else {
    bac4A$variable = factor(bac4A$variable, levels = c("aperc","lperc"))
    for (i in 1:numvals) {bac4A$addcol[i] = bac4A$value[i+numvals]}
  }

  ## run the stats
  try(rm(ghi, mid_test),T)
  cat("\n\nStats on relative abundances")
 # print(kruskal.test(bac$aperc ~ bac$testcol))
  mid_test <- kruskal.test(bac$aperc ~ bac$testcol)
  print(paste0("chi2 ",mid_test$parameter,", ",dim(bac)[1]," = ",round(mid_test$statistic,2),", p = ",round(mid_test$p.value,4)))
  if(kruskal.test(bac$aperc ~ bac$testcol)$p.value < 0.05) {
    if (length(table(bac$testcol))>2) {
      efg <- dunn.test(bac$aperc,bac$testcol, method = "bh", table = F, kw = F)
      ghi <- cldList(comparison = efg$comparisons, p.value = efg$P.adjusted, threshold = 0.05) 
      ghi$Group2 <- factor(ghi$Group, levels = levels(bac$testcol))
      apercclds <- as.character(ghi %>% arrange(Group2) %>% dplyr::select(Letter) %>% unlist() %>% unname())
 #     print(kruskal.test(bac$aperc ~ bac$testcol))
      cat("AAB Plot_values:")
#      print(ghi)
    } else {
      apercclds <- c("a","b")
    }
  } else {
    apercclds <- c(rep("a",numvals))
  }



 #   print(bac4)
    if(!is.null(yaxis_perc)) {    
      if(aceto_bottom == F) {
        bac4g <- bac4 %>% filter(variable == "lperc") %>% rbind(bac4 %>% filter(variable == "aperc") %>% mutate(value = value - (1-yaxis_perc[2])))
 #       print(bac4g)
 #       print(1)
        bac4plot <- ggplot(data = NULL) + 
          geom_bar(data = bac4g, aes(x = testcol, y = value, fill = variable), stat="identity") + 
          scale_fill_manual(values = c("red","blue"))+
          theme(legend.position = "bottom") + 
          geom_text(data = bac4 %>% filter(variable == "aperc"), aes(x = testcol, y = yaxis_perc[2]+loglabeladjust, label = c(apercclds), vjust = 0))+
          ylim(yaxis_perc*1.001) +
          theme_cowplot() 
        } else {
   #       print(2)
          bac4g <- bac4 %>% filter(variable == "aperc") %>% rbind(bac4 %>% filter(variable == "lperc") %>% mutate(value = value - (1-yaxis_perc[2])))
        bac4g$variablerev = factor(bac4g$variable, levels = rev(levels(bac4g$variable)))
        bac4plot <- ggplot(data = NULL) + 
          geom_bar(data = bac4g, aes(x = testcol, y = value, fill = variablerev), stat="identity") + 
          scale_fill_manual(values = c("blue","red"))+
          theme(legend.position = "bottom") + 
          geom_text(data = bac4g %>% filter(variable == "aperc"), aes(x = testcol, y = yaxis_perc[2]+loglabeladjust, label = c(apercclds), vjust = 0))+
          ylim(yaxis_perc*1.001) +
          theme_cowplot() 
      }
    } else {
      bac4g <- bac4
      if(aceto_bottom == F) {
    #    print(3)
        bac4plot <- ggplot(data = NULL) + 
          geom_bar(data = bac4g, aes(x = testcol, y = value, fill = variable), stat="identity") + 
          scale_fill_manual(values = c("red","blue"))+
          theme(legend.position = "bottom") + 
          geom_text(data = bac4g %>% filter(variable == "aperc"), aes(x = testcol, y = yaxis_perc[2]+loglabeladjust, label = c(apercclds), vjust = 0))+
          theme_cowplot() 
      } else {
   #     print(4)
        bac4c <- bac4 %>% filter(variable == "lperc") %>% rbind(bac4 %>% filter(variable == "aperc"))
        bac4c$variablerev = factor(bac4c$variable, levels = rev(levels(bac4c$variable)))
        bac4plot <- ggplot(data = NULL) + 
          geom_bar(data = bac4g, aes(x = testcol, y = value, fill = variablerev), stat="identity") + 
          scale_fill_manual(values = c("red","blue"))+
          theme(legend.position = "bottom") + 
          geom_text(data = bac4g %>% filter(variable == "aperc"), aes(x = testcol, y = yaxis_perc[2]+loglabeladjust, label = c(apercclds), vjust = 0))+
          theme_cowplot() 
      }
    }
#print("gothere5")

if(!is.null(add_jitter)) {    
  set.seed(43)
  bac4plot <- bac4plot + geom_jitter(data = bac4A %>% filter(variable == "lperc"), aes(x = testcol, y=value, group = testcol), width = 0.1, alpha = .4, size = .85, color = "white")
}

if(!is.null(add_violin)) {    
  set.seed(43)
  bac4plot <- bac4plot + geom_violin(data = bac4A, aes(x = testcol, y=value), alpha = 0, width = 0.6)
}

#print(bac4plot)
ifelse(plotout == 1, 
       return(bac3plot),
       ifelse(plotout == 2, 
              return(bac4plot),
              ifelse(plotout == 3, 
                     return(grid.arrange(bac3plot, bac4plot, heights = c(1,1))))))



}



calculate_figuremeans <- function(datain, humidityvar, sexvar) {
  
  ##
  # perc of means
  ##
  
  figure4asam_means = datain %>% 
    filter(Humidity%in%humidityvar, Sex %in% sexvar) %>%
    group_by(Geography) %>% 
    dplyr::summarize(meanA = mean(cfuA), 
                     semA = sd(cfuA)/sqrt(dplyr::n()), 
                     meanL = mean(cfuL), 
                     semL = sd(cfuL)/sqrt(dplyr::n()), 
                     count = dplyr::n()
    ) %>% 
    ungroup() %>% 
    mutate(percA = meanA/(meanA+meanL), percL = meanL/(meanA+meanL), total = meanA+meanL) %>% 
    droplevels()
  
  figure4asam_means <- figure4asam_means %>%
    mutate(latitude = c(1:length(table(figure4asam_means$Geography)))) 
  
  cat("\nperc of means\n")
  abc <- cor.test(figure4asam_means$latitude, figure4asam_means$meanA, method = "spearman")
  print(paste0("AAB: S = ",round(abc$statistic,2),", rho = ", round(abc$estimate,2),", p = ", round(abc$p.value,2)))
  abc <- cor.test(figure4asam_means$latitude, figure4asam_means$total, method = "spearman")
  print(paste0("Total: S = ",round(abc$statistic,2),", rho = ", round(abc$estimate,2),", p = ", round(abc$p.value,2)))
  abc <- cor.test(figure4asam_means$latitude, figure4asam_means$meanL, method = "spearman")
  print(paste0("LAB = ",round(abc$statistic,2),", rho = ", round(abc$estimate,2),", p = ", round(abc$p.value,2)))
  abc <- cor.test(figure4asam_means$latitude, figure4asam_means$percL, method = "spearman")
  print(paste0("percL = ",round(abc$statistic,2),", rho = ", round(abc$estimate,2),", p = ", round(abc$p.value,2)))
  
  plot(figure4asam_means$latitude, figure4asam_means$percL)  
  
}
