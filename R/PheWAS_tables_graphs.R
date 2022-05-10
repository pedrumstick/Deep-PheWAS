#!/usr/bin/env Rscript
'Creates result tables and graphs for the association analyses from either the GRS or plink methods. The two usages below represent inputting data from plink_results and from R_association_results.
By default with mandatory inputs the fubnction will save a single table in the folder /path_to_results_file/group, to produce graphs input the per_group_name_graph, per_snp_graph or R_association_graph depending on result source and transformation of the data required. There are many options to select related to filtering for MAC in plink results and general appearence of the graphs. By using the provided filter inputs it is possible to edit any indicidual graphs using the arary of options.

Usage: PheWAS_tables_graphs.R --results_file=<FILE> --analysis_name=<name> --plink_results --SNP_list=<FILE> [--save_folder=<FOLDER> --group_filter=<text> --PheWAS_ID_filter=<FILE> --PheWAS_manifest_overide=<FILE> --max_pheno=<number> --sig_FDR=<number> --no_save_table_all_results --no_graph_all --no_graph_sig --max_FDR_graph=<number> --SNP_filter=<FILE> --group_name_filter=<FILE> --save_raw_plink --MAC=<number> --MAC_case=<number> --MAC_control=<number> --per_group_name_graph --per_snp_graph --save_table_per_group_name --save_table_per_snp --sex_split]

       PheWAS_tables_graphs.R --results_file=<FILE> --analysis_name=<name> --R_association_results [--save_folder=<FOLDER> --group_filter=<text> --PheWAS_ID_filter=<FILE> --PheWAS_manifest_overide=<FILE> --max_pheno=<number> --sig_FDR=<number> --no_save_table_all_results --no_graph_all --no_graph_sig --max_FDR_graph=<number> --R_association_graph --sex_split]

Options:
    -h --help  Show this screen.
    -v --version  Show version.
    
    Mandatory with any input
    --results_file=<FILE>             Full file path of the results file RDS R list object.
    --analysis_name=<name>            Name for the analysis.
    
    Plink result input
    --plink_results                   Select if results are from PheWAS_association_PLINK.R
    --SNP_list=<FILE>                 Full file path to the SNP_list file, used for extracting variants with extracting_snps.R
    
    R association testing input
    --R_association_results           Select if results are from R_association_testing.R
    
    Options for both inputs
    --save_folder=<FOLDER>            Full file path to save folder location, if not used will default to folder that results_file originates from and
                                      files will be saved in the respective groups files. Example if results_file was found in
                                      /genetic_test/SNP_results/results_file and within that file were results from two groups, group_A and group_B,
                                      then tables and graphs would be produced per-group and stored in /genetic_test/SNP_results/group_A and
                                      /genetic_test/SNP_results/group_B. The analysis will be default have already created these folder structures.
    --group_filter=<text>             Comma separated text input, use if wanting to filter the group to which the table and grph functions are applied. 
                                      Inputted groups are the ones that are retained.
    --PheWAS_ID_filter=<FILE>         Full file path to file containing list of PheWAS_IDs, these will be retained and must include any sex split PheWAS_IDs 
                                      in full. Use if wanting to apply the table and/or graphing functions a subset of phenotypes.
    --PheWAS_manifest_overide=<FILE>  Full file path of the alternative PheWAS_manifest file.
    --max_pheno=<number>              Manual override for inputting maximum phenotypes analysed. Used for calculating FDR. The default used the largest number                                        of associations in the results_file.
    --sig_FDR=<number>                Value of FDR for which associations are reported as significant. Will alter the line of significance in both graph types
                                      and the reported phenotypes in the sognifgant_pheno graphs. [default: 0.01]
    --no_save_table_all_results       Select if not wanting to save a table of all the results.
    --no_graph_all                    Select if not wanting to produce a graph of all the results.
    --no_graph_sig                    Select if not wanting to produce a graph of the significant results by FDR.
    --max_FDR_graph=<number>          Value used in making graphs. Used when associations are 0, replaces this value with 10x-max_FDR_graph. This allows 
                                      plotting of the most significant results. [default: 300]
    --sex_split                       Select if wanting to split the graphs by sex. Used only with sex-specific phenotypes. Will create three versions of each                                        selected output, female, male and combined. Does not create split tables.
    
    Options for plink results
    --SNP_filter=<FILE>               Full file path to file containing list of rsids for SNPs. Use if wanting to apply the table and/or graphing functions
                                      a subset of results. Only works for plink_results.
    --group_name_filter=<FILE>        Full file path to file containing list of group_names that correspond with the SNP_list. Use if wanting to apply the 
                                      table and/or graphing functions a subset of results. Only works for plink_results.
    --save_raw_plink                  Select if wanting to save the raw Plink results unfiltered or altered. Is an option during PheWAS_association_PLINK.R,
                                      repeated here.
    --MAC=<number>                    MAC (minor allele count) filter applied to all associations, only applicable in results from PheWAS_association_PLINK.R.
                                      [default: 20]
    --MAC_case=<number>               MAC (minor allele count) filter applied only to cases, only applicable in results from PheWAS_association_PLINK.R
                                      and only in binary phenotypes. [default: 5]
    --MAC_control=<number>            MAC (minor allele count) filter applied only to controls, only applicable in results from PheWAS_association_PLINK.R
                                      and only in binary phenotypes. [default: 10]    
    --per_group_name_graph            Select if wanting to produce graphs per-group_name. This is used when looking to report the most significant finding 
                                      across several SNPs for a single construct, potential and gene or a sentinal SNP with a credible set. It is an 
                                      column in the SNP_list file.
    --per_snp_graph                   Select if wanting to produce graphs per SNP/variant.
    --save_table_per_group_name       Select if wanting to save a table of results per_group_name. Will extract the lowest FDR value for a  
                                      PheWAS_ID-group_name combination and report that. Only works for plink_results.
    --save_table_per_snp              Select if wanting to save a table of results per-SNP/genetic variant. Will be saved in a created folder named
                                      /analysis_name_group_per_SNP_tables. Example if the group was groupA and analysis_name top_SNPs the folder would be 
                                      /top_SNPs_groupA_per_SNP_tables. 

    Options R association results
    --R_association_graph             Select if wanting to produce graphs for R_association_results.
' -> doc

suppressMessages(library(docopt))

arguments <- docopt(doc, version = 'v0.1 PheWAS_tables_graphs.R')

suppressMessages(library(data.table))
suppressMessages(library(tidyverse))
suppressMessages(library(viridis))
suppressMessages(library(ggrepel))
suppressMessages(library(scales))
library(here)

# Functions ---------------------------------------------------------------
fdr_calc <- function(x) {
  
  fdr <- x %>% 
    distinct(PheWAS_ID, .keep_all = T) %>% 
    mutate(FDR=p.adjust(p = P, method = "fdr", n = max_pheno_tests)) %>% 
    arrange(FDR)
  
  return(fdr)
}
making_graphs <- function(a,b,c,d) {
  
  data <- a
  data_test <- data %>% 
    filter(FDR<=FDR_figure)
  
  if (nrow(data_test)<1) {
    return()
  } else {
    data <- data %>% 
      mutate(FDR=ifelse(FDR==0,max_FDR_graph,FDR))
    
    order <- data  %>% 
      group_by(phenotype_group) %>% 
      summarise(min_FDR=min(FDR)) %>% 
      arrange(min_FDR) %>% 
      mutate(group_number =seq(1:(nrow(.)))) %>% 
      select(phenotype_group,group_number)
    
    all_pheno_data <- data %>% 
      left_join(order,by="phenotype_group") %>% 
      arrange(FDR) %>% 
      arrange(group_number) %>% 
      mutate(annotate=ifelse(FDR<=FDR_figure,T,F),
             short_desc_T=paste0(short_desc," ",group),
             short_desc=factor(short_desc, levels = unique(.$short_desc)),
             seq=seq(nrow(.)),
             y_axis_info=-log10(FDR)) %>% 
      mutate(short_desc_T=factor(short_desc_T, levels = unique(.$short_desc_T)),
             seq=seq+150*(group_number-1)) %>% 
      mutate(size=case_when(is.na(OR) ~ sqrt(Beta^2),
                            is.na(Beta) ~ sqrt((log(OR))^2)))
    
    labels= summarize(group_by(all_pheno_data, group_number), tick=mean(unique(seq)),label=as.character(phenotype_group[1]))
    labels=labels[order(labels$tick),]
    
    round_any = function(x, accuracy, f=ceiling){f(x/ accuracy) * accuracy}
    
    shapes = c("positive" = 24, "negative" = 25)
    max.x<-max(all_pheno_data$seq)
    max.y <- round_any(max(all_pheno_data$y_axis_info),5)
    
    if(b=="per_group_name") {
      name_file <- paste0(unique(all_pheno_data$collective_name),d)
    } else if(b=="per_snp") {
      name_file <- paste0(unique(all_pheno_data$graph_save_name),"_",unique(all_pheno_data$group),d)
    } else if(b=="R") {
      name_file <- paste0(unique(all_pheno_data$name_group),d)
    }
    
    
    # making the graphs
    all_pheno_graph <- ggplot(all_pheno_data, aes(x = seq, y = -log10(FDR), color = phenotype_group, fill = phenotype_group, )) +
      geom_point(aes(shape=effect_direction)) +
      scale_shape_manual(values = shapes) +
      scale_colour_viridis(option = "H",discrete=T) +
      scale_fill_viridis(option= "H",discrete=T) +
      scale_y_continuous(limits=c(0,max.y)) +
      scale_x_continuous(name="Phenotype groups", limits=c(1,max.x), breaks=labels$tick, labels=labels$label, expand=c(.01,0)) +
      geom_hline(yintercept=-log10(FDR_figure),colour="red", alpha=I(1/3),size=1) +
      geom_text_repel(aes(label=short_desc),colour="black",data=all_pheno_data[all_pheno_data$annotate,],
                      size=2,angle=0,max.overlaps = 15) +
      guides(color=F,fill=F,shape=guide_legend("Direction of Effect")) +
      scale_size_area(name="Effect size",breaks=breaks_extended(3),max_size = 4) +
      theme(panel.background=element_blank(), 
            panel.grid.minor=element_blank(),
            axis.text.y=element_text(size=9, colour="black", hjust=1, vjust=.5), 
            axis.text.x=element_text(size=7, colour="black", angle=-45, hjust=0, vjust=0), 
            axis.line =element_line(colour="black"),
            axis.ticks=element_line(colour="black"),
            legend.key=element_blank()) 
    
    sig_pheno_data <- all_pheno_data %>% 
      filter(FDR<0.01)
    
    sig_results_all <- ggplot(sig_pheno_data, aes(x = short_desc, y = -log10(FDR), color = phenotype_group, fill = phenotype_group)) +
      geom_point(aes(shape=effect_direction)) +
      scale_shape_manual(values = shapes) +
      scale_x_discrete(name="Phenotypes") +
      scale_colour_viridis(option = "H",discrete=T) +
      scale_fill_viridis(option= "H",discrete=T) +
      geom_hline(aes(yintercept=-log10(FDR_figure)),colour="red", alpha=I(1/3),size=1) +
      scale_y_continuous(limits=c(0, max(-log10(sig_pheno_data$FDR)))) +
      guides(fill=guide_legend(title="Phenotypic Category", ncol=1, override.aes=list(shape=24,size=1)), 
             col=guide_legend(title="Phenotypic Category", ncol=1), 
             shape=guide_legend(title="Direction of Effect", ncol=1,override.aes=list(size=1))) +
      scale_size_area(name="Effect size",breaks=breaks_extended(3),max_size = 4) +
      theme(panel.background=element_blank(), 
            panel.grid.minor=element_blank(),
            axis.text.y=element_text(size=9, colour="black", hjust=1, vjust=.5), 
            axis.text.x=element_text(size=6, colour="black", angle=90, hjust=1, vjust=0), 
            axis.line =element_line(colour="black"),
            axis.ticks=element_line(colour="black"),
            legend.key=element_blank(),
            legend.key.height=unit(0.5,"line"))
    if(c=="both") {
      ggsave(filename = paste0(graph_save_location,"/",name_file,"_sig_pheno.png"),
             plot = sig_results_all, device = "png", dpi = 300, width = 9,height = 6,units = "in")
      ggsave(filename = paste0(graph_save_location,"/",name_file,"_all_pheno.png"),
             plot = all_pheno_graph, device = "png", dpi = 300, width = 9,height = 6,units = "in")
      
    } else if (c=="sig_only") {
      ggsave(filename = paste0(graph_save_location,"/",name_file,"_sig_pheno.png"),
             plot = sig_results_all, device = "png", dpi = 300, width = 9,height = 6,units = "in")
      
    } else if(c=="all_pheno") {
      ggsave(filename = paste0(graph_save_location,"/",name_file,"_all_pheno.png"),
             plot = all_pheno_graph, device = "png", dpi = 300, width = 9,height = 6,units = "in")
    }
  }
}
Deep_PheWAS_graphs_tables <- function(x,y){
  # define save folder
  if(is.null(arguments$save_folder)) {
    save_root <- paste0(dirname(arguments$results_file),"/",x)
  } else {
    save_root <- arguments$save_folder
  }
  dir.create(save_root,recursive = T)
  
  results <- y 
  
  if(is.null(arguments$PheWAS_ID_filter)){
    PheWAS_ID_filter <- unique(results$PheWAS_ID)
  } else {
    PheWAS_IDs <- fread(arguments$PheWAS_ID_filter) %>% 
      pull(1)
    PheWAS_ID_filter <- unique(PheWAS_IDs)
  }
  
  results_PheWAS_ID_filter <- results %>% 
    filter(PheWAS_ID %in% PheWAS_ID_filter)
  
  if(arguments$plink_results){
    if(arguments$save_raw_plink){
      save_name_plink <- paste0(arguments$analysis_name,"_",x,"_plink_results_raw.csv")
      fwrite(results,paste0(save_root,"/",save_name_plink))
    }
    
    main_table <- results_PheWAS_ID_filter %>% 
      left_join(SNP_list) %>% 
      filter(ID %in% SNP_list$ID,
             rsid %in% SNP_filter,
             group_name %in% group_name_filter) %>% 
      mutate(across(c(POS,A1_CT,ALLELE_CT,A1_CASE_CT,A1_CTRL_CT,
                      A1_FREQ,A1_CASE_FREQ,A1_CTRL_FREQ,MACH_R2,OBS_CT,
                      OR,`LOG(OR)_SE`,L95,U95,Z_STAT,P,BETA,SE,T_STAT), as.numeric),
             Beta=BETA,
             group=paste0(x),
             minor_allele=ifelse(A1_FREQ<0.5,A1,REF),
             MAF=ifelse(A1_FREQ<0.5,A1_FREQ,1-A1_FREQ),
             SE=ifelse(is.na(SE),`LOG(OR)_SE`,SE),
             alter_direction=case_when(coded_allele==A1 ~ 0,
                                       coded_allele!=A1 ~ 1),
             Beta=ifelse(is.na(Beta),NA,ifelse(alter_direction==1,-1*Beta,Beta)),
             OR=ifelse(is.na(OR),NA,ifelse(alter_direction==1,1/OR,OR)),
             N_L95=ifelse(is.na(OR),ifelse(alter_direction==1,-1*U95,L95),ifelse(alter_direction==1,1/U95,L95)),
             N_U95=ifelse(is.na(OR),ifelse(alter_direction==1,-1*L95,U95),ifelse(alter_direction==1,1/L95,U95)),
             effect_direction=ifelse(is.na(OR),ifelse(is.na(Beta),NA,ifelse(Beta>0,"positive","negative")),ifelse(OR>1,"positive","negative")),
             Z_STAT=ifelse(alter_direction==1,-1*Z_STAT,Z_STAT),
             T_STAT=ifelse(alter_direction==1,-1*T_STAT,T_STAT),
             Z_T_STAT=ifelse(is.na(Z_STAT),T_STAT,Z_STAT),
             Error_flag=ifelse(ERRCODE!=".",ERRCODE,NA),
             ID=as.factor(ID),
             MAC=ifelse(A1_FREQ<0.5,A1_CT,ALLELE_CT-A1_CT),
             MAC_cases=ifelse(A1_FREQ<0.5,A1_CASE_CT,(A1_CASE_CT/A1_CASE_FREQ)*(1-A1_CASE_FREQ)),
             MAC_controls=ifelse(A1_FREQ<0.5,A1_CTRL_CT,(A1_CTRL_CT/A1_CTRL_FREQ)*(1-A1_CTRL_FREQ)),
             expected_MAC_cases=ifelse(A1_FREQ<0.5,(A1_CASE_CT/A1_CASE_FREQ)*A1_CTRL_FREQ,(A1_CASE_CT/A1_CASE_FREQ)*(1-A1_CTRL_FREQ)),
             MAC_diff=MAC_cases-expected_MAC_cases,
             keep=ifelse(!is.na(Beta),2,ifelse(expected_MAC_cases<=2&MAC_cases>=3,1,ifelse(expected_MAC_cases>=7&(MAC_cases<=5&MAC_cases>=1),1,0))),
             ratio=1/((A1_CASE_CT/A1_CASE_FREQ)/(A1_CTRL_CT/A1_CTRL_FREQ)),
             sex_pheno=sub(".*_", "", PheWAS_ID),
             join_name = str_remove(PheWAS_ID,"_male|_female")) %>% 
      filter(MAC >=MAC_figure,
             case_when(keep==0 ~ MAC_cases >= MAC_cases_N & MAC_controls >= MAC_control_N,
                       keep==1 ~ MAC_cases>=1 & MAC_controls >= MAC_control_N,
                       keep==2 ~ MAC >= MAC_figure)) %>% 
      left_join(PheWAS_manifest,by=c("join_name"="PheWAS_ID")) %>% 
      mutate(short_desc=ifelse(sex_pheno!=PheWAS_ID,paste0(short_desc," (",sex_pheno,")"), short_desc)) %>% 
      select(group,collective_name=group_name,PheWAS_ID,category,description=phenotype,N_ID=OBS_CT,rsid,P,OR,Beta,L95=N_L95,U95=N_U95,coded_allele,non_coded_allele,minor_allele,MAF,MAC,MAC_cases,MAC_controls,chromosome=`#CHROM`,position=POS,Z_T_STAT,SE,effect_direction,category,phenotype_group=pheno_group,phenoytpe_group_narrow=group_narrow,short_desc,Info_score=MACH_R2,firth=`FIRTH?`,TEST,Error_flag,ID,graph_save_name)
    
    main_table_fdr_split <- main_table %>% 
      group_split(ID)
  }
  if(arguments$R_association_results){
    main_table <- results_PheWAS_ID_filter %>%  
      mutate(across(c(P,OR,Beta,SE,L95,U95), as.numeric)) %>% 
      select(-phenotype,-phenotype_group,-group_narrow,-short_desc) %>% 
      left_join(updated_manifest,by="PheWAS_ID")%>% 
      filter(PheWAS_ID %in% updated_manifest$PheWAS_ID) %>% 
      select(name,PheWAS_ID,category,phenotype,P,OR,Beta,L95,U95,SE,phenotype_group,group_narrow,short_desc,effect_direction,group,name_group)
    
    main_table_fdr_split <- as_tibble(main_table) %>% 
      group_split(name)
  }
  if(length(main_table_fdr_split)==0){
    return()
  }
  main_table_fdr <- lapply(main_table_fdr_split,fdr_calc) %>% 
    rbindlist(.) %>% 
    relocate(FDR,.before = P) %>% 
    arrange(FDR) %>% 
    mutate(sex_pheno_identifier="")
  
  if(arguments$sex_split){
    main_table_fdr <- main_table_fdr %>% 
      mutate(male=ifelse(str_detect(PheWAS_ID,"_male$"),1,0),
             female=ifelse(str_detect(PheWAS_ID,"_female$"),1,0),
             combined=ifelse(male==1 | female ==1,0,1),
             sex_pheno_identifier=case_when(male==1 ~ "male",
                                            female==1 ~ "female",
                                            combined==1 ~ "")) %>%
      select(-female,-male,-combined)
  }
  
  main_table_fdr_split <- main_table_fdr %>% 
    group_split(sex_pheno_identifier)
  
  if(isFALSE(arguments$no_save_table_all_results)) {
    fwrite(main_table_fdr,paste0(save_root,"/",arguments$analysis_name,"_",x,"_filtered_results.csv"))
  }
  
  if(arguments$plink_results){
    if(arguments$save_table_per_group_name) {
      per_group_name <- main_table_fdr %>% 
        group_by(collective_name,PheWAS_ID) %>% 
        slice_min(FDR) 
      per_group_name_split <- per_group_name %>% 
        ungroup() %>% 
        group_split(collective_name)
      
      grouping <- x
      per_group_name_folder <- paste0(save_root,"/",arguments$analysis_name,"_",x,"_results_per_collective_name/")
      dir.create(per_group_name_folder)
      fwrite(per_group_name,paste0(save_root,x,"_results_per_collective_name_all_",arguments$analysis_name,".csv"))
      lapply(per_group_name_split,function(x) fwrite(x,paste0(per_group_name_folder,unique(x$collective_name),"_",grouping,"_results.csv")))
    } 
    if(arguments$save_table_per_snp){
      per_snp_tables <- main_table_fdr %>% 
        group_split(ID)
      grouping <- x
      per_SNP_folder <- paste0(save_root,"/",arguments$analysis_name,"_",x,"_results_per_SNP/")
      dir.create(per_SNP_folder)
      lapply(per_snp_tables,function(x) fwrite(x,paste0(per_SNP_folder,unique(x$graph_save_name),"_",grouping,"_results.csv")))
    }
    if(arguments$per_group_name_graph) {
      per_group_function <- function(x,y) {
        per_group_name <- x %>% 
          group_by(collective_name,PheWAS_ID) %>% 
          slice_min(FDR) %>% 
          ungroup() %>% 
          group_split(collective_name)
        
        per_group_name_folder <- paste0(save_root,"/",arguments$analysis_name,"_",y,"_results_per_collective_name/")
        dir.create(per_group_name_folder)
        
        graph_save_location <<- per_group_name_folder
        
        sex_label <- unique(x$sex_pheno_identifier)
        
        mapply(making_graphs,per_group_name,MoreArgs = list(b="per_group_name",c=graph_choice,d=sex_label))
      }
      mapply(per_group_function,main_table_fdr_split,MoreArgs = list(y=x))
    }
    if(arguments$per_snp_graph) {
      
      per_snp_function <- function(x,y){
        per_SNP <- x %>% 
          group_split(ID)
        sex_label <- unique(x$sex_pheno_identifier)
        per_SNP_folder <- paste0(save_root,"/",arguments$analysis_name,"_",y,"_results_per_SNP/")
        dir.create(per_SNP_folder)
        graph_save_location <<- per_SNP_folder
        
        mapply(making_graphs,per_SNP,MoreArgs = list(b="per_snp",c=graph_choice,d=sex_label))
      }
      mapply(per_snp_function,main_table_fdr_split,MoreArgs = list(y=x))
    }
  }
  
  if(arguments$R_association_graph) {
    R_association_fuction <- function(x){
      R_association_graph <- x %>% 
        group_split(name)
      sex_label <- unique(x$sex_pheno_identifier)
      graph_save_location <<- save_root
      
      mapply(making_graphs,R_association_graph,MoreArgs = list(b="R",c=graph_choice,d=sex_label))
    }
    mapply(R_association_fuction,main_table_fdr_split)
  }
  
}
# Load in defining variables ----------------------------------------------
if(is.null(arguments$PheWAS_manifest_overide)){
  PheWAS_manifest <- fread(here("data","PheWAS_manifest.csv"))
} else {
  PheWAS_manifest <- fread(arguments$PheWAS_manifest_overide)
}
# create an updated manifests that includes sex-specific phenotypes, equates ot maximum possible phenotypes
PheWAS_manifest_edit <- PheWAS_manifest %>% 
  select(PheWAS_ID,sex,phenotype,category,pheno_group,group_narrow,short_desc,included_in_analysis) %>% 
  filter(included_in_analysis==1)
sex_names <- c("Male","Female")
PheWAS_manifest_male <- PheWAS_manifest_edit %>% 
  filter(!sex %in% sex_names) %>% 
  mutate(PheWAS_ID=paste0(PheWAS_ID,"_male"),
         short_desc=paste(short_desc,"_male"))
PheWAS_manifest_female <- PheWAS_manifest_edit %>% 
  filter(!sex %in% sex_names) %>% 
  mutate(PheWAS_ID=paste0(PheWAS_ID,"_female"),
         short_desc=paste(short_desc,"_female"))
updated_manifest <- PheWAS_manifest_edit %>% 
  bind_rows(PheWAS_manifest_male,PheWAS_manifest_female) %>% 
  select(PheWAS_ID,phenotype,category,phenotype_group=pheno_group,group_narrow,short_desc)
#read in all results
all_results <- readRDS(arguments$results_file)
# group_filter
if(!is.null(arguments$group_filter)) {
  group_filter <- unlist(strsplit(arguments$group_filter,","))
} else {
  group_filter <- names(all_results)
}
all_results <- all_results[group_filter]
group_names <- names(all_results)
# edit to plink results
if(arguments$plink_results){
  # SNP_list edit only for plink list
  SNP_list <- fread(arguments$SNP_list) %>% 
    mutate(group=group_name,
           coded_number=match(str_sub(coded_allele,1,1),LETTERS),
           non_coded_number=match(str_sub(non_coded_allele,1,1),LETTERS),
           multi=case_when(nchar(coded_allele) > 1 & nchar(non_coded_allele)== 1 ~ 1,
                           nchar(non_coded_allele) > 1 & nchar(coded_allele)== 1 ~ 2,
                           nchar(coded_allele) > 1 & nchar(non_coded_allele) <1 ~ 3),
           allele_order=ifelse(non_coded_number>coded_number,paste0(coded_allele,"_",non_coded_allele),
                               ifelse(coded_number>non_coded_number,paste0(non_coded_allele,"_",coded_allele),
                                      ifelse(coded_number==non_coded_number,
                                             ifelse(multi==1,paste0(non_coded_allele,"_",coded_allele),
                                                    ifelse(multi==2,paste0(coded_allele,"_",non_coded_allele),
                                                           ifelse(nchar(coded_allele)>nchar(non_coded_allele),paste0(non_coded_allele,"_",coded_allele),paste0(coded_allele,"_",non_coded_allele)))),NA))),
           ID=paste0(rsid,"_",allele_order)) %>% 
    select(chromosome,rsid,group_name,coded_allele,non_coded_allele,graph_save_name,ID)
  # create SNP_filter
  if(!is.null(arguments$SNP_filter)) {
    SNP_filter <- fread(arguments$SNP_filter) %>% 
      pull(rsid)
  } else {
    SNP_filter <- SNP_list %>% 
      pull(rsid)
  }
  # group_name filter
  if(!is.null(arguments$group_name_filter)) {
    group_name_filter <- fread(arguments$group_name_filter) %>% 
      pull(group_name)
  } else {
    group_name_filter <- SNP_list %>% 
      pull(group_name)
  }
}
# inputs into graphs and tables limits to only those that should have been analysed, but does not limit to those selected for graphs. Or is manually selected
if(is.null(arguments$max_pheno)) {
  max_pheno_tests <- max(unlist(lapply(all_results, function(x) nrow(data.frame(ids=unique(x$PheWAS_ID)) %>% filter(ids %in% updated_manifest$PheWAS_ID)))))
} else {
  max_pheno_tests <- as.numeric(arguments$max_pheno)
}
FDR_figure <- as.numeric(arguments$sig_FDR)
MAC_figure <- as.numeric(arguments$MAC)
MAC_cases_N <- as.numeric(arguments$MAC_case)
MAC_control_N <- as.numeric(arguments$MAC_control)
max_FDR_graph_parse <- as.numeric(arguments$max_FDR_graph)
max_FDR_graph <- 10^(-max_FDR_graph_parse)
# graph choices
if(arguments$no_graph_all & isFALSE(arguments$no_graph_sig)) {
  graph_choice <- "sig_only"
} else if(isFALSE(arguments$no_graph_all) & arguments$no_graph_sig) {
  graph_choice <- "all_pheno"
} else {
  graph_choice <- "both"
}
# Run functions -----------------------------------------------------------
mapply(Deep_PheWAS_graphs_tables,group_names,all_results)
