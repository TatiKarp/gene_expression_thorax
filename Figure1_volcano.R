## This script will create a volcano plot of the differentially expressed genes between asthma and healthy in ATLANTIS 
library(ggplot2)
library(ggrepel)

setwd("/Users/tatiana/Work/RP2/ATLANTIS")
de.results <- read.csv("./Umi_dedup/Dif_expr/DE.genes.ATLANTIS.csv")

# define the y axis coordinate (-log10 P-value) for the FDR = 0.05 
fdr_yintersept <- -log10((de.results$PValue[sum(de.results$FDR < 0.05)] + de.results$PValue[sum(de.results$FDR < 0.05)+1])/2 )

p <- ggplot(data=de.results, aes(x=logFC, y=-log10(PValue))) +
  geom_point(aes(color= ifelse((FDR < 0.05)&(logFC>0),'firebrick', ifelse((FDR < 0.05)&(logFC<0),'blue','gray')))) + 
  geom_hline(yintercept = fdr_yintersept, col="black", linewidth = 0.3, linetype = 2) + #look at table
  geom_text_repel(aes(label=ifelse(((FDR < 0.05 )& (logFC>1))|((FDR < 0.05 )& (logFC<(-1))), (ifelse((!is.na(hgnc_symbol)), hgnc_symbol, '')), ''),
                      lineheight=0.5, hjust= 0.5, vjust= 0.4),
                  max.overlaps =100, 
                  direction = 'both',
                  alpha = 0.6, 
                  box.padding=0.3,
                  point.padding=0.5,
                  size = 5)+# segment.color = 'grey50'))
  xlab (expression (log[2]~fold~change))+
  ylab (expression(-log[10]~Pvalue))+
  scale_color_identity(name = '', breaks= c('firebrick', 'blue','gray') ,
                       labels = c('higher expressed','lower expressed', 'not significant'), guide = "legend")+
  theme(text = element_text(family = "Arial"), 
        axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="bold"),
        legend.position = "bottom",
        panel.grid.major =  element_line(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", linewidth =0.6),
        legend.text=element_text(size=15))+
  xlim(-4, 4)

png("./Umi_dedup/Dif_expr/Volcano_ERJ.png",
    width=2100, height=2100, res=300)
print(p)
dev.off()