library(tidyverse)
library(writexl)
library(readxl)

setwd("~/Bioinformatics/Jumbophages/")

# Read in the outputs for > 200 kb refseq

Process_Matches<-function(filename,metaspacers,IMG.VR){

matches.in<-read_delim(paste("output/results_server/",filename,".txt",sep=""),delim="\t",col_names=c("spacer.ID","target.ID","strand","position","detials","mismatches","gaps","alignment.seq","alignment.length"))
  matches.in$alignment.length<-as.numeric(matches.in$alignment.length)
  matches.in$mismatches<-as.numeric(matches.in$mismatches)
  matches.in$gaps<-as.numeric(matches.in$gaps)

# Remake the ID data
  matches<-matches.in%>%separate(spacer.ID,c("GCF","CRISPR"),sep="-Cr")%>%separate(CRISPR,c("CRISPR","spacer"),sep="-Sp")%>%
    separate(spacer,c("spacer","shared.spacer"),sep="-dup-")%>%
    separate(target.ID,c("target.ID","shuffled"),sep="-")%>%mutate(shuffled=ifelse(is.na(shuffled)==F,shuffled,"non-shuffled"))

if (IMG.VR==F){  
  # Read in phage names
    phage.names<-read_delim("output/Phage_name_lookup.txt",delim="\t",col_names=T)  
  # Add target genome length data
    genbank.phage.lengths<-read_delim("output/genbank_phage_genome_lengths_lookup.txt",delim="\t",skip=1,col_names=c("target.ID","target.GCA","target.length"))%>%
      select(target.GCA,target.length)
  phage.names<-phage.names%>%left_join(genbank.phage.lengths,by="target.GCA")%>%rename(target.ID=target.id,phage=species)
}  

if (IMG.VR==T){  
  # Read in phage names
  phage.names<-read_delim("output/IMG_VR_metadata_all.txt",delim="\t",col_names=c("target.ID","phage","target.GCA","target.length","target.completeness"))
}  
  
if (metaspacers==F){  
    # Join the spacer/repeat assignment information
      spacer.details<-read_delim("output/typeIII_spacers.txt",col_names=c("GCF","CRISPR","spacer.seq","spacer","spacer.ID","CRISPR.score","CRISPR.type","repeat.seq","repeat.length","ID"),delim="\t")
  # Join
    matches.details<-matches%>%left_join(spacer.details,by=c("GCF","CRISPR","spacer"))%>%mutate(spacer.length=nchar(spacer.seq),notes="")%>%
      mutate(score=spacer.length-2*mismatches-gaps*3)
  # Read in strain names
    strain.names<-read_delim("output/GCF_strain_names_typeIII.txt",delim="\t",col_names=c("assembly","strain"))%>%
      mutate(assembly=substr(assembly,1,15))#%>%group_by(assembly)%>%mutate(count=n())
  # Join - need to reformat GCF part first    
      matches.details<-matches.details%>%mutate(assembly=substr(GCF,1,15))%>%left_join(strain.names,by="assembly")%>%left_join(phage.names,by="target.ID")
}  

# Format column order    
  matches.details<-matches.details%>%select(GCF,strain,target.ID,target.GCA,phage,target.length,shuffled,notes,spacer.length,alignment.length,mismatches,gaps,score,CRISPR.type,CRISPR.score,repeat.seq,repeat.length,spacer.seq,everything())%>%
    arrange(desc(score))

  filename<-gsub(".txt","",filename)
  filename<-gsub(".fasta","",filename)

# Some of the IMGVR alignments have NNNNNN - remove these...
  matches.details<-matches.details%>%filter(grepl("NN",alignment.seq)==F)  
  
# Output the data
  write_xlsx(matches.details,paste0("output/",filename,"_target-matches.xlsx"))  

if (IMG.VR==F){    
  # Read in the nucleoid-forming classification data...
    all.nucleoid<-read_xlsx("output/identified_nucleoid_phages_Genbank.xlsx")
  # Bind to the spacer hit result data...
    classified.results<-matches.details%>%left_join(all.nucleoid,by="target.GCA")%>%mutate(target.completeness=100)
}
  
if (IMG.VR==T){    
    # Read in the nucleoid-forming classification data...
      all.nucleoid<-read_xlsx("output/identified_nucleoid_phages_IMGVR_metadata.xlsx")
    # Remove duplicate columns
      all.nucleoid<-all.nucleoid%>%select(-target.length,-phage,-IMGVR.cluster,-cluster.count,-target.completeness)
    # Bind to the spacer hit result data...
      classified.results<-matches.details%>%left_join(all.nucleoid,by="target.ID")
}
  # Tidy up the nucleoid-forming column for NA entries # shouldn't need this anymore...
    classified.results$nucleoid.forming[is.na(classified.results$nucleoid.forming)]<-"nonNF"
  # Tidy column order etc..
    classified.results<-classified.results%>%select(GCF,target.GCA,strain,phage,nucleoid.forming,score,notes,everything())%>%filter(gaps==0)
  # Reduce redundancy - same host spacer targeting multiple phages...
    # Start by taking only the top socing hits
      classified.results<-classified.results%>%group_by(spacer.seq)%>%mutate(spacer.hits.total=n())%>%filter(score==max(score))%>%mutate(spacer.hits.topscoring=n())
    # Keep targets with the highest estimated completeness
      classified.results<-classified.results%>%group_by(spacer.seq)%>%filter(target.completeness==max(target.completeness))
    # If there a multiple spacers matching a contig, then keep the contigs with the most hits
      classified.results<-classified.results%>%group_by(GCF,target.ID)%>%mutate(host.target.hits=n())
      classified.results<-classified.results%>%group_by(spacer.seq)%>%filter(host.target.hits==max(host.target.hits))
      classified.results<-classified.results%>%mutate(spacer.hits.final=n())
  # Collapse target details and keep the first entry (lowest GCA number?)  
    classified.results<-classified.results%>%arrange(desc(target.length),target.ID,target.GCA)%>%mutate(target.GCA=paste0(unique(target.GCA),collapse=";"))%>%
      mutate(phage=paste0(unique(phage), collapse=";"))%>%
      mutate(nucleoid.forming=paste0(unique(nucleoid.forming), collapse=";"))%>%
      mutate(target.ID=paste0(unique(target.ID), collapse=";"))%>%
      mutate(strand=paste0(strand, collapse=";"))%>%
      mutate(position=paste0(position, collapse=";"))%>%
      mutate(alignment.seq=paste0(unique(alignment.seq), collapse=";"))%>%
      mutate(shell.ID=paste0(unique(shell.ID), collapse=";"))%>%
      mutate(tubulin.ID=paste0(unique(tubulin.ID), collapse=";"))%>%
      mutate(shell.evalue=paste0(unique(shell.evalue), collapse=";"))%>%
      mutate(tubulin.evalue=paste0(unique(tubulin.evalue), collapse=";"))%>%
      mutate(shell.description=paste0(unique(shell.description), collapse=";"))%>%
      mutate(tubulin.description=paste0(unique(tubulin.description), collapse=";"))%>%
      mutate(target.length=max(target.length))%>%
      unique()%>%group_by(spacer.seq)%>%mutate(spacer.hits.final=n())%>%arrange(desc(score)) 

  # Output
    write_xlsx(classified.results,paste0("output/Processed_III_",filename,".xlsx"))
    
return (classified.results)
}
  
# Data were first processed as <> 150 kb, but these are later merged and filtered into <> 200 kb.
hits_Genbank_gb150kb<-Process_Matches("search_results_typeIII_spacers_unique_Refseq_genbank_caudovirales_150kb.fasta",F,F)%>%
  mutate(system="TypeIII",target.database="gb>150kb")
hits_Genbank_gb150kb_shuffled<-Process_Matches("search_results_typeIII_spacers_unique_Refseq_caudovirales_150kb.fasta_shuffled",F,F)%>%
  mutate(system="TypeIII",target.database="gb>150kb")

hits_Genbank_gb_sub150kb<-Process_Matches("search_results_typeIII_spacers_unique_Refseq_caudovirales_sub150kb.fasta",F,F)%>%
  mutate(system="TypeIII",target.database="gb<150kb")
hits_Genbank_gb_sub150kb_shuffled<-Process_Matches("search_results_typeIII_spacers_unique_Refseq_caudovirales_sub150kb.fasta_shuffled",F,F)%>%
  mutate(system="TypeIII",target.database="gb<150kb")

hits_Genbank_IMGVR150kb<-Process_Matches("search_results_typeIII_spacers_unique_Refseq_IMGVR_150kb.fasta",F,T)%>%
  mutate(system="TypeIII",target.database="IMGVR>150kb")
hits_Genbank_IMGVR150kb_shuffled<-Process_Matches("search_results_typeIII_spacers_unique_Refseq_IMGVR_150kb.fasta_shuffled",F,T)%>%
  mutate(system="TypeIII",target.database="IMGVR>150kb")

hits_Genbank_IMGVRsub150kb<-Process_Matches("search_results_typeIII_spacers_unique_Refseq_IMGVR_sub150kb.fasta",F,T)%>%
  mutate(system="TypeIII",target.database="IMGVR<150kb")
hits_Genbank_IMGVRsub150kb_shuffled<-Process_Matches("search_results_typeIII_spacers_unique_Refseq_IMGVR_sub150kb.fasta_shuffled",F,T)%>%
  mutate(system="TypeIII",target.database="IMGVR<150kb")

# Bind all the data into one big table...

all_hits<-hits_Genbank_gb150kb%>%bind_rows(hits_Genbank_gb150kb_shuffled)%>%
  bind_rows(hits_Genbank_IMGVR150kb)%>%bind_rows(hits_Genbank_IMGVR150kb_shuffled)%>%
  bind_rows(hits_Genbank_gb_sub150kb)%>%bind_rows(hits_Genbank_gb_sub150kb_shuffled)%>%
  bind_rows(hits_Genbank_IMGVRsub150kb)%>%bind_rows(hits_Genbank_IMGVRsub150kb_shuffled)

# Output the data - use another script to analyse all types together...
  write_xlsx(all_hits,"output/hits_typeIII.xlsx")
