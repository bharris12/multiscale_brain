library(scrattch)
library(tidyverse)
library(RColorBrewer)

metacell<-read.csv('~/biccn_paper/data/metacell_cluster_assignments.csv')
colnames(metacell)<-c("sample_id","cluster_id")
metacell$cluster_label<-paste0(metacell$cluster_id,'_metacell')
color_df<-as.data.frame(paste0(c(brewer.pal(9, "Set1"),brewer.pal(9, "Set3")),"60"))

metacell$cluster_color<-color_df[metacell$cluster_id]
metacell$cluster_color<-color_df[metacell$cluster_id,1]

bulk<-read.csv('~/biccn_paper/data/bulk_cluster_assignments.csv')
colnames(bulk)<-c("sample_id","cluster_id")
bulk$cluster_label<-paste0(bulk$cluster_id,'_bulk')

color_df2<-as.data.frame(rev(paste0(c(brewer.pal(9, "Set1"),brewer.pal(9, "Set3")),"60")))
bulk$cluster_color<-color_df2[bulk$cluster_id,1]

common.cells<-metacell$sample_id

df1 = metacell %>% filter(sample_id %in% common.cells) %>% select(cols)
df2 = bulk %>% filter(sample_id %in% common.cells) %>% select(cols)
colnames(df2)[-1] = paste0("bulk_",     colnames(df2)[-1])
df<-left_join(df1, df2)
scrattch.vis::build_river_plot(df,c('cluster','bulk_cluster'),pad=0, show_labels = F) + ggsave('~/biccn_paper/figures/river_plot_metacells_bulk.pdf',device = 'pdf',dpi='retina')
