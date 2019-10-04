library(tidyverse)
library(writexl)
library(readxl)

setwd("~/Bioinformatics/Jumbophages/")

# Read in all the data for the systems...
  typeIII<-read_xlsx("output/hits_typeIII.xlsx")
  typeIE<-read_xlsx("output/hits_typeIE.xlsx")
  typeIF<-read_xlsx("output/hits_typeIF.xlsx")
 
all.hits<-typeIII%>%bind_rows(typeIE,typeIF)  
  
# Remove duplicate spacer hits.
  all.hits<-all.hits%>%group_by(system,shuffled,spacer.seq)%>%mutate(spacer.hits=n())%>%
    filter(score==max(score))%>%filter(target.length==max(target.length))%>%
    arrange(target.database)%>%filter(row_number()==1)%>%
    filter(target.length==max(target.length))%>%
    group_by(system,shuffled,spacer.seq)%>%mutate(spacer.hits=n())   

# Now everything is unique...
  length.cutoff<-200000 # >200 kb are jumbophages

scoring.table.200kb<-all.hits%>%filter(target.length>=length.cutoff)%>%mutate(score=ifelse(score<=18,18,score))%>%group_by(system,shuffled,score)%>%
  summarise(score.freq=n())%>%unite("heading",system,shuffled)%>%spread(heading,score.freq)

scoring.table.sub200kb<-all.hits%>%filter(target.length<length.cutoff)%>%mutate(score=ifelse(score<=18,18,score))%>%group_by(system,shuffled,score)%>%
  summarise(score.freq=n())%>%unite("heading",system,shuffled)%>%spread(heading,score.freq)

write_xlsx(scoring.table.200kb,"figures/Scoring_table_all_grouped_200kb.xlsx")
write_xlsx(scoring.table.sub200kb,"figures/Scoring_table_all_grouped_sub200kb.xlsx")


scoring.table.all<-all.hits%>%mutate(score=ifelse(score<=18,18,score))%>%group_by(system,shuffled,score)%>%
  summarise(score.freq=n())%>%unite("heading",system,shuffled)%>%spread(heading,score.freq)

write_xlsx(scoring.table.all,"figures/Scoring_table_all_grouped.xlsx")

scoring.table.NF<-all.hits%>%filter(nucleoid.forming!="nonNF")%>%mutate(score=ifelse(score<=18,18,score))%>%group_by(system,shuffled,score)%>%
  summarise(score.freq=n())%>%unite("heading",system,shuffled)%>%spread(heading,score.freq)

write_xlsx(scoring.table.NF,"figures/Scoring_table_all_grouped_NF.xlsx")

#

tax.table<-read_delim("taxonomy/bac120_taxonomy_r89.tsv",delim="\t",col_names=c("GCF_tax","taxonomy"))
tax.table<-tax.table%>%mutate(GCF_tax=gsub("RS_","",GCF_tax))%>%mutate(GCF_tax=gsub("\\..*","",GCF_tax))

NF.hits<-all.hits%>%filter(nucleoid.forming!="nonNF")%>%filter(score>=25)%>%mutate(GCF_tax=gsub("\\..*","",GCF))%>%left_join(tax.table,by="GCF_tax")
write_xlsx(NF.hits,"figures/All_nucleoid_hits.xlsx")

# Processing top hits

hits.top<-read_xlsx("output/results_server/Final_analysis_processed_Genbank_10-09-19_IMGVR_150kb.xlsx")

hits.top<-hits.top%>%filter(score>=26 & nucleoid.forming!="nonNF")

aa<-hits.top%>%group_by(GCF,target.GCA)%>%mutate(host.cluster.hits=n(),targets.in.relationship=length(unique(target.ID)))

# Keep one representative example for each host-target relationship (i.e. if there are more than 1 spacers targeting the phage, just count once)
  unique.host.target.hits<-aa%>%group_by(GCF,target.GCA,target.ID)%>%arrange(desc(score))%>%mutate(scores.all=paste0(score,collapse=";"))%>%filter(score==max(score))%>%filter(row_number()==1)

# Summary
total.host.target.relationships<-nrow(unique.host.target.hits)
unique.hosts<-length(unique(unique.host.target.hits$GCF))  
unique.target.contigs<-length(unique(unique.host.target.hits$target.ID))  
unique.target.clusters<-length(unique(unique.host.target.hits$target.GCA))  

#####

