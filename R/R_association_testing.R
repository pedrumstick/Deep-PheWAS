#!/usr/bin/env Rscript
'
Runs regression models using R. Is slower than the plink method but is more flexible designed to be used when the input is not a standard SNP model or one that can easily be inputted into standard formats for genetic data such bed, bgen or VCF. Use cases for this method include genomic risk scores(GRS)/polygenic risk score (PRS) and multiallelic copy number variants (CNV). Whilst this could be used for SNP analysis it would be significantly slower than other provided methods and has no built in filter for MAC/MAF or correction for case-control imbalance.

User must select whether the analysis is to be done on a GRS/PRS or on non-GRS/PRS  data such as an alternative genetic measurement such as CNVs. 

GRS/PRS is inputted via a csv file with trait,group,phenotype_data,genetic_data,variable this acts as a map for the functions . See user guide for full description of this input but in brief, it allows the user to input any number of traits that may be analysed in any of the availble groupings (such as ancestry) and point to the location of the appropriate genetic and phenotype files. GRS/PRS data is analysed per-trait-group combination, results are saved as a R list object per trait that contain all trait-group combinations, this allows multiple traits to be analysed across multiple groups in a single input.

Non-GRS/PRS data is assumed to be an alternative genetic measurement with a column of IDs followed by 1-n columns of genetic data, example CNVs. Results are per column of the input file (minus ID column) with a single R object saved containing a list of results per-group such that each non-ID column will be tested for association with all available phenotypes in each of the inputted groups. Unlike the GRS/PRS input the full location of all the phenotype files is required as input.

If a covariate file is included then it must contain a column age, when selecting covariates note currently only participants with full covariate data will be analysed. Other options perform filtering on the phenotypes that are analysed, by default all phenotypes that are held in the respective phenotype tables created by pre_association_preparation.R are analysed. Minimum case numbers are options to allow the input of genetic variables that do not fully overlap with the samples used for making the phenotypes.


Usage:
    R_association_testing.R --analysis_folder=<folder> --GRS_input=<FILE> [--covariates=<FILE> --N_cores=<number> --PheWAS_manifest_overide=<FILE>] [--exclude_sex_phenotypes | --sex_include_file=<FILE> | --phenotype_overide=<FILE>]
    
    R_association_testing.R --analysis_folder=<folder> --non_GRS_data=<FILE> --phenotype_files=<FILE> --analysis_name=<name> [--covariates=<FILE> --group_name_overide=<text> --N_cores=<number> --PheWAS_manifest_overide=<FILE>] [--exclude_sex_phenotypes | --sex_include_file=<FILE> | --phenotype_overide=<FILE>]
    
    
Options:
    -h --help  Show this screen.
    -v --version  Show version.
    
    Mandatory Inputs
    --analysis_folder=<FOLDER>            Full file path of the root folder where the analysis results will be saved. If GRS data a trait folder
                                          containing a folder for each group will be created within the analysis folder.
    One of
    --GRS_input=<FILE>                    Full file path of the GRS_input csv file. See user guide for more information on format.             
    --non_GRS_data=<FILE>                 Full file path to the non_GRS genetic data. See user guide for more information on format.
    
    If non_GRS_data must include
    --phenotype_files=<FILE>              Comma separated full file paths of the phenotype data. For example 
                                          /home/phenotypes/EUR_pheno, /home/phenotypes/AFR_pheno 
                                          Only required for non_GRS_data.
    --analysis_name=<name>                Name for the analysis, is used later in saving tables, so should distignuish between other analyses. For
                                          GRS analysis this name is always the trait being analysed.
    Options
    --covariates=<FILE>                   Full file path for the covariates file see user guide for more information.
    --group_name_overide=<text>           Comma separated text input for group names. By default, group names are created by extracted all text 
                                          before the first underscore (_) from the phenotype file names (note not the full path but the files names) of 
                                          the inputs from phenotype_files. 
                                          Example /home/phenotypes/EUR_pheno
                                          Would extract EUR as the group name. The order for names should match the order of the input for phenotype_files.
    
    --N_cores=<number>                    Number of cores requested if wanting to use parallel computing.
    
    --exclude_sex_phenotypes              Select if wanting to exclude the sex specific phenotypes (if any) from the analysis.
    --sex_include_file=<FILE>             Full file path to list of desired included sex specific phenotypes, defaults to select all.
    --phenotype_overide=<FILE>            Full file path to a text file containing a single column of the desired PheWAS_IDs for inclusion in the analysis.
    
    --binary_Case_N=<number>              Number that represents the minimum number of cases for binary phenotype inclusion. 
                                          [default: 50]
    --quantitative_Case_N=<number>        Number that represents the minimum number of cases for quantitative phenotype inclusion. 
                                          [default: 100]
    --PheWAS_manifest_overide=<FILE>      Full file path of the alternative PheWAS_manifest file.
    
' -> doc

suppressMessages(library(docopt))

arguments <- docopt(doc, version = 'v0.1 R_association_testing.R')

suppressMessages(library(data.table))
suppressMessages(library(tidyverse))
suppressMessages(library(parallel))
library(here)

# Functions ---------------------------------------------------------------
GRS_association <- function(e,a,b,c,d) {
  
  results_name <- a  
  
  if(!is.null(arguments$GRS_input)){
    table_save_name <- results_name
  } else if (!is.null(arguments$non_GRS_data)){
    table_save_name <- analysis_save_name
  }
  
  variable_name <- d
  
  ## Load the GRS
  GRS <- fread(b) %>% 
    rename(eid=1)
  
  # combine GRS with covariates
  GRS_cov <- GRS %>% 
    left_join(covariates) %>% 
    drop_na()
  
  # link phenotypes with GRS
  phenotypes <- fread(c) %>% 
    select(1,any_of(all_phenos$PheWAS_ID))
  
  available_phenotypes <- colnames(phenotypes)
  
  if(!is.null(arguments$covariates)) {
    association_guide <- all_phenos %>% 
      filter(include==1) %>% 
      filter(PheWAS_ID %in% available_phenotypes) %>% 
      mutate(temp_name=str_remove(PheWAS_ID,"_male|_female"),
             age_col_complete = ifelse(age_col=="named",paste0(temp_name,"_age"),"age"),
             analysis_option = ifelse(analysis=="quant","gaussian","binomial")) %>% 
      select(PheWAS_ID,age_col_complete,analysis_option)
  } else {
    association_guide <- all_phenos %>% 
      filter(include==1) %>% 
      filter(PheWAS_ID %in% available_phenotypes) %>% 
      mutate(temp_name=str_remove(PheWAS_ID,"_male|_female"),
             age_col_complete = "",
             analysis_option = ifelse(analysis=="quant","gaussian","binomial")) %>% 
      select(PheWAS_ID,age_col_complete,analysis_option)
    
  }
  
  GRS_association_testing <- function (a,b,c,d) {
    extracted_cols <- c(a,b)
    
    temp_name <- str_remove(a,"_male|_female")
    
    phenotype_extracted <- phenotypes %>% 
      select(eid,any_of(extracted_cols)) %>% 
      left_join(GRS_cov) %>% 
      drop_na()
    
    phenotype_selected <- phenotype_extracted %>% 
      select({{a}}) %>% 
      rename(selected_phenotype=1)
    
    type <- PheWAS_manifest %>% 
      filter(PheWAS_ID=={{temp_name}}) %>% 
      select(analysis) %>% 
      pull()
    
    if (type=="quant"){
      if(nrow(phenotype_selected)<quant_min_cases){
        return()
      }
      
    } else if(type=="binary") {
      if(sum(phenotype_selected$selected_phenotype)<binary_min_cases)
        return()
    }
    
    if(b=="") {
      formula_covariates <- paste(colnames(covariate_col_names), collapse = " + ")
    } else {
      formula_covariates <- paste0(" + ", paste(c(b,colnames(covariate_col_names)), collapse = " + "))
    }
    
    if (c=="binomial") {
      my_form <- paste0("as.factor(",a,")~ ",d,formula_covariates)
      my_formula <- as.formula(my_form)
      
      my_glm <- glm(my_formula, family = c, data = phenotype_extracted, na.action = na.omit)
      
      p_value <- (coef(summary(my_glm))[,4])[d]
      SE <- unname((summary(my_glm)$coefficients[, 2])[d])
      P <- unname(p_value)
      if (is.na(P)) {
        P <- NA
        odds_df <- data.frame(list(OR=NA,L95=NA,U95=NA,SE=NA)) %>% 
          rename(OR=1,L95=2,U95=3,SE=4)
        
      } else if (P > 0.05) {
        OR_raw <- exp(coef(my_glm))[d]
        OR <- unname(OR_raw)
        if (is.na(OR)) {
          odds_df <- data.frame(list(OR=NA,L95=NA,U95=NA,SE=NA)) %>% 
            rename(OR=1,L95=2,U95=3,SE=NA)
        } else {
          odds <- exp(coef(my_glm)[d])
          odds_df <- data.frame(list(OR=unname(odds),L95=NA,U95=NA,SE=SE)) %>% 
            rename(OR=1,L95=2,U95=3,SE=4)
        }
      } else {
        OR_raw <- exp(coef(my_glm))[d]
        OR <- unname(OR_raw)
        if (is.na(OR)) {
          odds_df <- data.frame(list(OR=NA,L95=NA,U95=NA,SE=NA)) %>% 
            rename(OR=1,L95=2,U95=3,SE=4)
        } else {
          odds_CI <- exp(odds_CI <- c(Beta = coef(my_glm)[d], confint(my_glm,d),SE=SE))             
          odds_df <- data.frame(as.list(odds_CI)) %>% 
            rename(OR=1,L95=2,U95=3,SE=4)
        }
      }
      
      
    } else if (c=="gaussian") {
      my_form <- paste0(a," ~ ",d,formula_covariates)
      my_formula <- as.formula(my_form)
      
      my_glm <- glm(my_formula, family = c, data = phenotype_extracted, na.action = na.omit)
      
      p_value <- (coef(summary(my_glm))[,4])[d]
      SE <- unname((summary(my_glm)$coefficients[, 2])[d])
      P <- unname(p_value)
      if (is.na(P)) {
        P <- NA
        odds_df <- data.frame(list(Beta=NA,L95=NA,U95=NA,SE=NA)) %>% 
          rename(Beta=1,L95=2,U95=3,SE=4)
        
      } else if (P > 0.05) {
        OR_raw <- coef(my_glm)[d]
        OR <- unname(OR_raw)
        if (is.na(OR)) {
          odds_df <- data.frame(list(Beta=NA,L95=NA,U95=NA,SE=NA)) %>% 
            rename(Beta=1,L95=2,U95=3,SE=4)
        } else {
          odds <- coef(my_glm)[d]
          odds_df <- data.frame(list(Beta=unname(odds),L95=NA,U95=NA,SE=SE)) %>% 
            rename(Beta=1,L95=2,U95=3,SE=4)
        }
      } else {
        OR_raw <- coef(my_glm)[d]
        OR <- unname(OR_raw)
        if (is.na(OR)) {
          odds_df <- data.frame(list(Beta=NA,L95=NA,U95=NA,SE=NA)) %>% 
            rename(Beta=1,L95=2,U95=3,SE=4)
        } else {
          odds_CI <- c(Beta = coef(my_glm)[d], confint(my_glm,d),SE=SE)             
          odds_df <- data.frame(as.list(odds_CI)) %>% 
            rename(Beta=1,L95=2,U95=3,SE=4)
        }
      }
    }
    
    P_df <- data.frame(P)
    
    name_df <- data.frame(a) %>% 
      rename(PheWAS_ID=1)
    
    df_results <- name_df %>% 
      bind_cols(P_df,odds_df)
    
    return(df_results)
    
  }
  ## all association results have tables with Phewas_ID,P,OR,Beta
  association_results <- mapply(GRS_association_testing,
                                association_guide$PheWAS_ID,
                                association_guide$age_col_complete,
                                association_guide$analysis_option,
                                MoreArgs = list(d=variable_name),
                                SIMPLIFY = F) %>% 
    reduce(bind_rows)
  
  if ("Beta" %in% colnames(association_results) & "OR" %in% colnames(association_results)) {
    association_results_edit <- association_results %>% 
      rowwise() %>% 
      mutate(P=ifelse(P>0.05,signif(P, 2),signif(P, 3)),
             L95=ifelse(is.na(Beta),signif(L95, 5),signif(L95, 2)),
             U95=ifelse(is.na(Beta),signif(U95, 5),signif(U95, 2)),
             Beta=signif(Beta, 2),
             name=results_name,
             SE=signif(SE, 2),
             effect_direction=ifelse(is.na(OR),ifelse(is.na(Beta),NA,ifelse(Beta>0,"positive","negative")),ifelse(OR>1,"positive","negative")),
             OR=signif(OR, 5),
             group=e,
             name_group=paste0(name,"_",group),
             table_save=table_save_name,
             sex_pheno=sub(".*_", "", PheWAS_ID),
             join_name = str_remove(PheWAS_ID,"_male|_female")) %>%
      left_join(PheWAS_manifest,by=c("join_name"="PheWAS_ID")) %>% 
      mutate(short_desc=ifelse(sex_pheno!=PheWAS_ID,paste0(short_desc," (",sex_pheno,")"), short_desc)) %>% 
      select(name,PheWAS_ID,phenotype,P,OR,Beta,L95,U95,SE,phenotype_group=pheno_group,group_narrow,short_desc,effect_direction,group,name_group,table_save)
    
  } else if ("Beta" %in% colnames(association_results) & !"OR" %in% colnames(association_results)) {
    association_results_edit <- association_results %>% 
      rowwise() %>% 
      mutate(P=ifelse(P>0.05,signif(P, 2),signif(P, 3)),
             L95=signif(L95, 2),
             U95=signif(U95, 2),
             Beta=signif(Beta, 2),
             name=results_name,
             SE=signif(SE, 2),
             effect_direction=ifelse(is.na(Beta),NA,ifelse(Beta>0,"positive","negative")),
             group=e,
             name_group=paste0(name,"_",group),
             table_save=table_save_name,
             sex_pheno=sub(".*_", "", PheWAS_ID),
             join_name = str_remove(PheWAS_ID,"_male|_female")) %>%
      left_join(PheWAS_manifest,by=c("join_name"="PheWAS_ID")) %>% 
      mutate(short_desc=ifelse(sex_pheno!=PheWAS_ID,paste0(short_desc," (",sex_pheno,")"), short_desc)) %>%
      select(name,PheWAS_ID,phenotype,P,Beta,L95,U95,SE,phenotype_group=pheno_group,group_narrow,short_desc,effect_direction,group,name_group,table_save)
    
  } else if (!"Beta" %in% colnames(association_results) & "OR" %in% colnames(association_results)) {
    association_results_edit <- association_results %>% 
      rowwise() %>% 
      mutate(P=ifelse(P>0.05,signif(P, 2),signif(P, 3)),
             L95=signif(L95, 5),
             U95=signif(U95, 5),
             name=results_name,
             SE=signif(SE, 2),
             effect_direction=ifelse(is.na(OR),NA,ifelse(OR>1,"positive","negative")),
             OR=signif(OR, 5),
             group=e,
             name_group=paste0(name,"_",group),
             table_save=table_save_name,
             sex_pheno=sub(".*_", "", PheWAS_ID),
             join_name = str_remove(PheWAS_ID,"_male|_female")) %>%
      left_join(PheWAS_manifest,by=c("join_name"="PheWAS_ID")) %>% 
      mutate(short_desc=ifelse(sex_pheno!=PheWAS_ID,paste0(short_desc," (",sex_pheno,")"), short_desc)) %>% 
      select(name,PheWAS_ID,phenotype,P,OR,L95,U95,SE,phenotype_group=pheno_group,group_narrow,short_desc,effect_direction,group,name_group,table_save)
    
  }
  
  return(association_results_edit)
  
}

# Load in defining variables ----------------------------------------------
if(is.null(arguments$PheWAS_manifest_overide)){
  PheWAS_manifest <- fread(here("data","PheWAS_manifest.csv"))
} else {
  PheWAS_manifest <- fread(arguments$PheWAS_manifest_overide)
}

if(!is.null(arguments$N_cores)) {
  N_cores <- as.numeric(arguments$N_cores)
} else {
  N_cores <- NA
}

quant_min_cases <- as.numeric(arguments$quantitative_Case_N)
binary_min_cases <-as.numeric(arguments$binary_Case_N)

if(!is.null(arguments$covariates)) {
  covariates <- fread(arguments$covariates)
  covariate_col_names <- covariates %>% 
    select(-eid,-age)
} else {
  covariates <- data.frame(eid=NA)
  covariate_col_names <- covariates %>% 
    select(-eid)
}


if(!is.null(arguments$phenotype_overide)) {
  
  phenotype_overide <- fread(arguments$phenotype_overide) %>% 
    pull(1)
  
  age_phenos <- PheWAS_manifest %>% 
    filter(age_col!="age") %>% 
    filter(PheWAS_ID %in% phenotype_overide) %>% 
    mutate(age_col_complete = ifelse(age_col=="named",paste0(PheWAS_ID,"_age"),"age")) %>% 
    select(PheWAS_ID=age_col_complete)
  
  sex_phenos <- PheWAS_manifest %>% 
    filter(included_in_analysis==1) %>%
    mutate(male=paste0(PheWAS_ID,"_male"),
           female=paste0(PheWAS_ID,"_female"),
           male_analysis=analysis,female_analysis=analysis,
           male_age_col=age_col,female_age_col=age_col) %>% 
    select(PheWAS_ID,male,female,analysis,male_analysis,female_analysis,age_col,male_age_col,female_age_col)
  
  all_phenos_complete <- data.frame(PheWAS_ID=c(sex_phenos$PheWAS_ID,sex_phenos$male,sex_phenos$female), analysis=c(sex_phenos$analysis,sex_phenos$male_analysis,sex_phenos$female_analysis), age_col=c(sex_phenos$age_col,sex_phenos$male_age_col,sex_phenos$female_age_col))
  
  all_phenos <- all_phenos_complete %>%
    filter(PheWAS_ID %in% phenotype_overide) %>% 
    mutate(include=1) %>% 
    bind_rows(age_phenos)
  
} else if (arguments$exclude_sex_phenotypes) {
  
  age_phenos <- PheWAS_manifest %>% 
    filter(age_col!="age") %>% 
    filter(included_in_analysis==1) %>% 
    mutate(age_col_complete = ifelse(age_col=="named",paste0(PheWAS_ID,"_age"),"age")) %>% 
    select(PheWAS_ID=age_col_complete)
  
  all_phenos <- PheWAS_manifest %>% 
    filter(included_in_analysis==1) %>% 
    mutate(include=1) %>% 
    select(PheWAS_ID,analysis,age_col,include) %>% 
    bind_rows(age_phenos)
  
} else if(is.null(arguments$sex_include_file)) {
  
  age_phenos <- PheWAS_manifest %>% 
    filter(age_col!="age") %>% 
    filter(included_in_analysis==1) %>% 
    mutate(age_col_complete = ifelse(age_col=="named",paste0(PheWAS_ID,"_age"),"age")) %>% 
    select(PheWAS_ID=age_col_complete)
  
  sex_phenos <- PheWAS_manifest %>% 
    filter(included_in_analysis==1) %>%
    mutate(male=paste0(PheWAS_ID,"_male"),
           female=paste0(PheWAS_ID,"_female"),
           male_analysis=analysis,female_analysis=analysis,
           male_age_col=age_col,female_age_col=age_col) %>% 
    select(PheWAS_ID,male,female,analysis,male_analysis,female_analysis,age_col,male_age_col,female_age_col)
  
  all_phenos_complete <- data.frame(PheWAS_ID=c(sex_phenos$PheWAS_ID,sex_phenos$male,sex_phenos$female), analysis=c(sex_phenos$analysis,sex_phenos$male_analysis,sex_phenos$female_analysis), age_col=c(sex_phenos$age_col,sex_phenos$male_age_col,sex_phenos$female_age_col))
  
  all_phenos <- all_phenos_complete %>% 
    mutate(include=1) %>% 
    bind_rows(age_phenos)
  
} else {
  
  age_phenos <- PheWAS_manifest %>% 
    filter(age_col!="age") %>% 
    filter(included_in_analysis==1) %>% 
    mutate(age_col_complete = ifelse(age_col=="named",paste0(PheWAS_ID,"_age"),"age")) %>% 
    select(PheWAS_ID=age_col_complete)
  
  sex_phenos_list <- fread(arguments$sex_include_file) %>% 
    left_join(PheWAS_manifest) %>% 
    mutate(male=paste0(PheWAS_ID,"_male"),
           female=paste0(PheWAS_ID,"_female"),
           male_analysis=analysis,female_analysis=analysis,
           male_age_col=age_col,female_age_col=age_col) %>% 
    select(PheWAS_ID,male,female,analysis,male_analysis,female_analysis,age_col,male_age_col,female_age_col)
  
  sex_phenos <- data.frame(PheWAS_ID=c(sex_phenos_list$male,sex_phenos_list$female), analysis=c(sex_phenos_list$male_analysis,sex_phenos_list$female_analysis), age_col=c(sex_phenos_list$male_age_col,sex_phenos_list$female_age_col))
  
  all_phenos <- PheWAS_manifest %>% 
    filter(included_in_analysis==1) %>% 
    bind_rows(sex_phenos) %>% 
    mutate(include=1) %>% 
    select(PheWAS_ID,analysis,age_col,include) %>% 
    bind_rows(age_phenos)
  
}

# Run functions -----------------------------------------------------------
if(!is.null(arguments$GRS_input)) {
  
  GRS_map_edit <- fread(arguments$GRS_input) %>% 
    mutate(save_location=paste0(arguments$analysis_folder,"/",trait,"/",trait,"_","results_list.RDS")) %>% 
    group_split(trait)
  
  # create folders to host results
  making_GRS_folders <- function(x) {
    
    GRS_per_trait <- x
    trait_name <- unique(x$trait)
    
    lapply(GRS_per_trait$group, function(x) dir.create(paste0(arguments$analysis_folder,"/",trait_name,"/",x),recursive = T))
    
  }
  
  lapply(GRS_map_edit,making_GRS_folders)

  per_group_per_trait_GRS <- function(x) {
    
    selected_trait <- x
    
    save_name <- unique(selected_trait$save_location)
    
    GRS_results <- mapply(GRS_association,
                            selected_trait$group,
                            selected_trait$trait,
                            selected_trait$genetic_data,
                            selected_trait$phenotype_data,
                            selected_trait$column_name,
                            SIMPLIFY = F,
                            USE.NAMES = T)
    saveRDS(GRS_results,save_name)
  }
  
  if(!is.null(arguments$N_cores)) {
    mclapply(GRS_map_edit,per_group_per_trait_GRS,
             mc.cores = N_cores)
  } else {
    lapply(GRS_map_edit,per_group_per_trait_GRS)
  }
  
} else if(!is.null(arguments$non_GRS_data)){
  
  genetic_data <- fread(arguments$non_GRS_data) %>% 
    rename(eid=1)
  
  tested_variables <- genetic_data %>% 
    select(-eid) %>% 
    colnames(.)
  
  analysis_save_name <- arguments$analysis_name
  
  # groups
  if(is.null(arguments$group_name_overide)){
  phenotype_group_name <-gsub("_.*", "", basename(unlist(strsplit(arguments$phenotype_files,","))))
  } else {
    phenotype_group_name <- unlist(strsplit(arguments$group_name_overide,","))
  }
  phenotype_files <- unlist(strsplit(arguments$phenotype_files,","))
  
  # making folders
  lapply(phenotype_group_name, function(x) dir.create(paste0(arguments$analysis_folder,"/",x),recursive = T))
  
  per_group_none_GRS  <- function(x,y){
    GRS_results <- mapply(GRS_association,
                          a=tested_variables,
                          d=tested_variables,
                          MoreArgs = list(b=arguments$non_GRS_data,c=y,e=x),
                          SIMPLIFY = F,
                          USE.NAMES = T) %>% 
      reduce(bind_rows)
    return(GRS_results)
  }
  
  if(!is.null(arguments$N_cores)) {
    all_results <- mcmapply(per_group_none_GRS,
                            phenotype_group_name,
                            phenotype_files,
                            mc.cores = N_cores,
                            SIMPLIFY = F,
                            USE.NAMES = T)
    
  } else {
    all_results <-  mapply(per_group_none_GRS,
                           phenotype_group_name,
                           phenotype_files,
                           SIMPLIFY = F,
                           USE.NAMES = T)
  }
  
  saveRDS(all_results,paste0(arguments$analysis_folder,"/",arguments$analysis_name,"_all_group_results_list.RDS"))
  
}

