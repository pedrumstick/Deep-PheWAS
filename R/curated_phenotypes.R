#!/usr/bin/env Rscript
'Creates curated phenotypes and curated concepts using a central map file, also creates population_control ID lists that are used in creating the curated phenotypes. The phenotype creating script iterates over itself (default 5 times) as some curated phenotypes use other curated phenotypes in their definition.  

By default the script looks for files in the root Deep_PheWAS directory for curated_phenotype_map.csv (included as a starting file) and control_populations (created in the phenotype_data.R script see user guide for file format), the user can specify a folder for stored phenotypes if using the default setting on previous scripts these files will be concepts.RDS, phecodes.RDS, quantitative.RDS,range_ID.RDS, all are required.

Usage:
    curated_phenotypes.R [--curated_phenotype_map=<FILE> --phenotype_folder=<FOLDER> --control_populations=<FILE> --N_iterations=<number> --phenotype_save_file=<FILE> --update_list]
    
Options:
    -h --help                                 Show this screen.
    -v --version                              Show version.
    
    --curated_phenotype_map=<FILE>            Full file path of the curated_phenotyope_map file. [default: data/curated_phenotype_map.csv]
    --phenotype_folder=<FOLDER>               Full file path of the folder containing the phenotype data (if using default names will contain,
                                              concepts.RDS,phecodes.RDS,quantitative.RDS,range_ID.RDS). Function will search for files with .RDS and combine
                                              safest method is to store phenotype files is a seperate directory. [default: data/phenotypes/]
    
    --N_iterations=<number>                   Number of iterations curated phenotype script is to run through. [default: 5]
                                            
    --phenotype_save_file=<FILE>              Full file path for the save file for the generated concepts RDS used for phenotype creation. 
                                              [default: data/phenotypes/curated_phenotypes.RDS]
                                              
    --control_populations=<FILE>              Full file path for the control_populations file. Contains columns of IDs, columns names are used to create
                                              lists of IDs. Saved as a list called control_populations.RDS in the folder inputted in the
                                              phenotype_folder flag. The lists are used to define some of the curated phenotype control populations 
                                              as directed by the curated_phenotype_map file. The un-edited version of the curated_phenotype_map file 
                                              uses two populations, all_pop and primary_care_pop, which represent all IDs in the sample and all IDs with
                                              availble primary care data.
                                              To create curated phenotypes the names of these list must match the curated_phenotype_map file. 
                                              The phenotype_date.R script creates a control_populations file that is used by default. 
                                              [default: data/control_populations]
    --update_list                             Option to run code as an update to existing curated_phenotype list object will use phenotype_save_file to load
                                              existing list and then save over that file upon completion.
' -> doc

suppressMessages(library(docopt))

arguments <- docopt(doc, version = 'v0.1 curated_phenotypes.R')

suppressMessages(library(data.table))
suppressMessages(library(tidyverse))
suppressMessages(library(parallel))
suppressMessages(library(lubridate))
library(here)

# Functions ---------------------------------------------------------------
curated_phenotype_creation <- function(a,x) {
  
  pheno_name <- a
  message(pheno_name)
  all_concepts <- x %>% 
    select(PheWAS_ID,starts_with("C_")) %>% 
    pivot_longer(starts_with("C_"), names_to = "Concepts", values_drop_na = T) %>% 
    pull(value)
  
  all_concepts_edit <- unlist(lapply(seq(1:length(all_concepts)),function(x) unlist(strsplit(all_concepts[[x]],","))))
  
  (!all_concepts_edit %in% c(names(phenotypes),names(curated_pheno_list)))
  
  if(!all(all_concepts_edit %in% c(names(phenotypes),names(curated_pheno_list)))) {
    
    return()
    
  } else {
    
    case_input <- x %>% 
      filter(`Case/control`=="case") %>% 
      group_split(Case_N)
    
    control_input <- x %>% 
      filter(`Case/control`=="control") %>% 
      group_split(Control_N)
    
    per_line_edit <- function(x) {
      #removing the fixed cols and then removing cols with NA from the columns requiring conversion
      fixed_details <- x %>% 
        select(PheWAS_ID,phenotype,broad_category,range_ID,`Case/control`,Control_N,Case_N)
      
      conversion_columns <- x %>% 
        select(!any_of(colnames(fixed_details))) %>% 
        select_if(~ !any(is.na(.)))
      
      all_concepts <- conversion_columns %>% 
        select(starts_with("C_"))
      
      # round about way of getting the Number of concepts in each line
      N_concepts <- lapply(colnames(conversion_columns), function(x) str_count(x,"C_")) %>% 
        reduce(c) %>% 
        sum(.)
      # and number of operators
      N_operators <- lapply(colnames(conversion_columns), function(x) str_count(x,"O_")) %>% 
        reduce(c) %>% 
        sum(.)
      
      # input into function
      concept_seq <- seq(1:N_concepts)
      operator_seq <- seq(1:N_operators)
      
      if(N_operators==0) {
        
        C_N <- paste0("C_",1)
        N_N <- paste0("N_",1)
        G_Na <- paste0("G_",1,"a")
        G_Nb <- paste0("G_",1,"b")
        CBB <- paste0("CBB_",1)
        CBA <- paste0("CBA_",1)   
        
        if((x %>% select(!!G_Na))[[1]]=="minf"){
          gap_a_convert <- -Inf } else {
            gap_a_convert <- as.numeric((x %>% select(!!G_Na))[[1]])
          }
        
        if((x %>% select(!!G_Nb))[[1]]=="minf"){
          gap_b_convert <- -Inf } else {
            gap_b_convert <- as.numeric((x %>% select(!!G_Nb))[[1]])
          }
        
        fully_combined <- list(list(concept_file=list((x %>% select(!!C_N))[[1]]),
                                    n_codes=list((x %>% select(!!N_N))[[1]]),
                                    gap_a=list(gap_a_convert),
                                    gap_b=list(gap_b_convert)))
        
        names(fully_combined) <- "Concept_1"
        
      } else {
        
        # making lists out of predictable inputs that are consecutively numbered
        concept_list_conversion <- function(x,y) {
          
          C_N <- paste0("C_",x)
          N_N <- paste0("N_",x)
          G_Na <- paste0("G_",x,"a")
          G_Nb <- paste0("G_",x,"b")
          CBB <- paste0("CBB_",x)
          CBA <- paste0("CBA_",x)   
          
          if((y %>% select(!!G_Na))[[1]]=="minf"){
            gap_a_convert <- -Inf } else {
              gap_a_convert <- as.numeric((y %>% select(!!G_Na))[[1]])
            }
          
          if((y %>% select(!!G_Nb))[[1]]=="minf"){
            gap_b_convert <- -Inf } else {
              gap_b_convert <- as.numeric((y %>% select(!!G_Nb))[[1]])
            }
          
          concept <- list(concept_file=list((y %>% select(!!C_N))[[1]]),
                          n_codes=list((y %>% select(!!N_N))[[1]]),
                          gap_a=list(gap_a_convert),
                          gap_b=list(gap_b_convert))
          
          if((y %>% select(!!CBB))[[1]]=="NB" & (y %>% select(!!CBA))[[1]]=="NB") {
            return(concept)
          } else if ((y %>% select(!!CBB))[[1]]=="NB" & (y %>% select(!!CBA))[[1]]!="NB") {
            bracket <- c(list(concept_file=list((y %>% select(!!C_N))[[1]]),
                              n_codes=list((y %>% select(!!N_N))[[1]]),
                              gap_a=list(gap_a_convert),
                              gap_b=list(gap_b_convert)),
                         list(bracket=list((y %>% select(!!CBA))[[1]])))
            return(bracket)
          } else if ((y %>% select(!!CBB))[[1]]!="NB" & (y %>% select(!!CBA))[[1]]=="NB") {
            bracket <- list(list(bracket=list((y %>% select(!!CBA))[[1]])),
                            list(concept_file=list((y %>% select(!!C_N))[[1]]),
                                 n_codes=list((y %>% select(!!N_N))[[1]]),
                                 gap_a=list(gap_a_convert),
                                 gap_b=list(gap_b_convert)))
            return(bracket)
          } 
        }
        operator_list_conversion <-  function(x,y) {
          O_N <- paste0("O_",x)
          OBB <- paste0("OBB_",x)
          OBA <- paste0("OBA_",x) 
          operator <- list(Boolean=list((y %>% select(!!O_N))[[1]]))
          
          if((y %>% select(!!OBB))[[1]]=="NB" & (y %>% select(!!OBA))[[1]]=="NB") {
            return(operator)
          } else if ((y %>% select(!!OBB))[[1]]=="NB" & (y %>% select(!!OBA))[[1]]!="NB") {
            bracket <- list(Boolean=list((y %>% select(!!O_N))[[1]]),
                            bracket=list((y %>% select(!!OBA))[[1]]))
            return(bracket)
          } else if ((y %>% select(!!OBB))[[1]]!="NB" & (y %>% select(!!OBA))[[1]]=="NB") {
            bracket <- list(bracket=list((y %>% select(!!OBB))[[1]]),
                            Boolean=list((y %>% select(!!O_N))[[1]]))
            return(bracket)
          } 
        }
        
        concept_conversion <- mapply(concept_list_conversion,concept_seq,MoreArgs = list(y=x),SIMPLIFY = F)
        operator_conversion <- mapply(operator_list_conversion,operator_seq,MoreArgs = list(y=x),SIMPLIFY = F)
        
        # combines the lists again predictable to make it work remove the last concept
        making_case_lists <- function(x,y) {
          part_case <- list(concept_conversion[[x]],operator_conversion[[y]])
          return(part_case)
        }
        
        making_case_names <- function(x,y) {
          part_case <- c(paste0("Concept_",x),paste0("Boolean_",y))
          return(part_case)
        }
        
        concept_seq_edit <- seq(1:(N_concepts-1))
        
        fully_combined <- mapply(making_case_lists,concept_seq_edit,operator_seq, SIMPLIFY = F) %>% 
          reduce(append) %>% 
          append(.,list(concept_conversion[[N_concepts]]))
        
        combined_names <- unlist(mapply(making_case_names,concept_seq_edit,operator_seq, SIMPLIFY = F)) %>% 
          c(paste0("Concept_",N_concepts))
        
        names(fully_combined) <- combined_names 
        
      }
      
      
      # adding the brackets to the correct position in the right list hierarchy
      combined_correct_order <- list()
      
      for(i in 1:length(fully_combined)) {
        if(str_detect(names(fully_combined[i]),"Concept")) {
          list_concept_x <- list(fully_combined[[i]])
          names(list_concept_x) <- names(fully_combined)[i]
          concept_x <- fully_combined[[i]]
          
          if(length(concept_x)>4) {
            bracket_index_concept <- which(names(concept_x)=="bracket")
            concept_list <- list(concept_x[-bracket_index_concept])
            names(concept_list) <- names(fully_combined)[i]
            
            bracket_string <- unlist(concept_x$bracket)
            bracket_names <- rep("bracket",nchar(bracket_string))
            split_brackets <- as.list(strsplit(bracket_string,"")[[1]])
            names(split_brackets) <- bracket_names
            
            if(bracket_index_concept==5) {
              combined_correct_order <- append(combined_correct_order,concept_list)
              for(i in 1:length(split_brackets)) {
                combined_correct_order <- append(combined_correct_order,list(bracket=split_brackets[i]))
              }
            } else {
              for(i in 1:length(split_brackets)) {
                combined_correct_order <- append(combined_correct_order,list(bracket=split_brackets[i]))
              }
              combined_correct_order <- append(combined_correct_order,concept_list)
            } 
          } else {
            combined_correct_order <- append(combined_correct_order,list_concept_x)
          }
          
        } else {
          
          list_operator_x <- list(fully_combined[[i]])
          names(list_operator_x) <- names(fully_combined)[i]
          operator_x <- fully_combined[[i]]
          
          if(length(operator_x)>1) {
            bracket_index_Boolean <- which(names(operator_x)=="bracket")
            operator_list <- list(operator_x[-bracket_index_Boolean])
            names(operator_list) <- names(fully_combined)[i]
            
            bracket_string <- unlist(operator_x$bracket)
            bracket_names <- rep("bracket",nchar(bracket_string))
            split_brackets <- as.list(strsplit(bracket_string,"")[[1]])
            names(split_brackets) <- bracket_names
            
            if(bracket_index_Boolean==2) {
              combined_correct_order <- append(combined_correct_order,operator_list)
              for(i in 1:length(split_brackets)) {
                combined_correct_order <- append(combined_correct_order,list(bracket=split_brackets[i]))
              }
            } else {
              for(i in 1:length(split_brackets)) {
                combined_correct_order <- append(combined_correct_order,list(bracket=split_brackets[i]))
              }
              combined_correct_order <- append(combined_correct_order,operator_list)
            } 
          } else {
            combined_correct_order <- append(combined_correct_order,list_operator_x)
          }
          
        } 
        
      }
      
      # Now the formula is in the correct order including the bracket terms we can evaluate the formula in the correct order using the function below.
      # works by filling and clearing the list items below as the formula is evaluated left to right.
      LHS_list <- list()
      Boolean_list <- list()
      RHS_list <- list()
      start_bracket_list <- list()
      end_bracket_list <- list()
      
      for(i in 1:length(combined_correct_order)) {
        
        if(str_detect(names(combined_correct_order)[[i]], "Concept")) {
          
          concept_x <- combined_correct_order[[i]]
          
          if(length(LHS_list)==0) {
            LHS_list <- append(LHS_list,list(concept=concept_x))
            
            LHS_concept_codes <- unlist(strsplit(LHS_list[[i]]$concept_file[[1]],","))
            LHS_N_codes <- as.numeric(unlist(strsplit(LHS_list[[i]]$n_codes[[1]],",")))
            evaluated_statement <- mapply(combining_codes,LHS_concept_codes,LHS_N_codes,SIMPLIFY = F) %>% 
              reduce(bind_rows) %>% 
              group_by(eid) %>% 
              summarise(any_code=sum(any_code),
                        earliest_date=min(earliest_date))
            
          } else if(length(start_bracket_list)==length(end_bracket_list)) {
            RHS_list <- append(RHS_list,list(concept=concept_x))
          } else if(length(start_bracket_list)>length(end_bracket_list)) {
            
            if(length(LHS_list)>=(length(start_bracket_list)-length(end_bracket_list)+1)) {
              RHS_list <- append(RHS_list,list(concept=concept_x))
            } else {
              LHS_list <- append(LHS_list,list(concept=concept_x))
            }
            
          }
          
        } else if(str_detect(names(combined_correct_order)[[i]], "Boolean")) {
          
          Boolean_x <- combined_correct_order[[i]]
          
          Boolean_list <- append(Boolean_list,list(Boolean_x))
          
          
        }  else if(str_detect(names(combined_correct_order)[[i]], "bracket")) {
          
          bracket_x <- combined_correct_order[[i]]
          
          bracket_extracted <- unlist(bracket_x)
          
          if(bracket_extracted=="(") {
            
            start_bracket_list <- append(start_bracket_list,list(bracket_x))
            
          } else if(bracket_extracted==")") {
            end_bracket_list <- append(end_bracket_list,list(bracket_x))
          }
          
        }
        
        if(length(LHS_list)>=1 & length(RHS_list)>=1 & length(Boolean_list)>=1){
          
          index_LHS <- length(LHS_list)
          index_RHS <- length(RHS_list)
          index_Boolean <- length(Boolean_list)
          
          if(names(LHS_list[index_LHS])=="concept") {
            
            LHS_N_codes <- as.numeric(unlist(strsplit(LHS_list[[index_LHS]]$n_codes[[1]],",")))
            LHS_concept_codes <- unlist(strsplit(LHS_list[[index_LHS]]$concept_file[[1]],","))
            LHS_concept <- mapply(combining_codes,LHS_concept_codes,LHS_N_codes,SIMPLIFY = F) %>% 
              reduce(bind_rows) %>% 
              group_by(eid) %>% 
              summarise(any_code=sum(any_code),
                        earliest_date=min(earliest_date))
            
          } else if(names(LHS_list[index_LHS])=="evaluated_statement") {
            
            LHS_N_codes <- 1
            LHS_concept <- LHS_list[[index_LHS]]$concept_file
          }
          
          if(names(RHS_list[index_RHS])=="concept") {
            
            RHS_N_codes <- as.numeric(unlist(strsplit(RHS_list[[index_RHS]]$n_codes[[1]],",")))
            gap_a <- unlist(RHS_list[[index_RHS]]$gap_a[[1]],",")
            gap_b <- unlist(RHS_list[[index_RHS]]$gap_b[[1]],",")
            RHS_concept_codes <- unlist(strsplit(RHS_list[[index_RHS]]$concept_file[[1]],","))
            RHS_concept <- mapply(combining_codes,RHS_concept_codes,RHS_N_codes,SIMPLIFY = F) %>% 
              reduce(bind_rows) %>% 
              group_by(eid) %>% 
              summarise(any_code=sum(any_code),
                        earliest_date=min(earliest_date))
            
          } else if(names(RHS_list[index_RHS])=="evaluated_statement") {
            
            RHS_N_codes <- 1
            gap_a <- -Inf
            gap_b <- Inf
            RHS_concept <- RHS_list[[index_RHS]]$concept_file
            
          }
          
          Boolean <- unlist(strsplit(Boolean_list[[index_Boolean]]$Boolean[[1]],","))
          
          if(Boolean=="AND") {
            evaluated_statement <- mapply(AND,list(LHS_concept),list(RHS_concept),list(gap_a),list(gap_b),SIMPLIFY = F)[[1]] 
            
          } else if(Boolean=="NOT") {
            evaluated_statement <- mapply(NOT,list(LHS_concept),list(RHS_concept),list(gap_a),list(gap_b),SIMPLIFY = F)[[1]] 
            
          }
          
          evaluated_list <- list(concept_file=evaluated_statement,n_codes=list("1"),gap_a=list("-Inf"),gap_b=list("Inf"))
          
          LHS_list <- LHS_list[-index_LHS]
          Boolean_list <- Boolean_list[-index_Boolean]
          RHS_list <- RHS_list[-index_RHS]
          
          if(length(start_bracket_list)-length(end_bracket_list)==0) {
            LHS_list <- append(LHS_list,list(evaluated_statement=evaluated_list))
            
          } else if(length(start_bracket_list)==length(LHS_list) & str_detect(names(combined_correct_order)[[i+1]], "bracket")) {
            RHS_list <- append(RHS_list,list(evaluated_statement=evaluated_list))
            
          } else {
            LHS_list <- append(LHS_list,list(evaluated_statement=evaluated_list))
          }
        }
        
      }
      return(evaluated_statement)
    }
    
    curated_cases <- mapply(per_line_edit,case_input, SIMPLIFY = F) %>% 
      reduce(bind_rows) %>% 
      group_by(eid) %>% 
      summarise(any_code=sum(any_code),
                earliest_date=min(earliest_date))
    
    if(length(control_input)>0) { 
      
      range_ID <- unique(x$range_ID)
      
      range_ID_full <- phenotypes[[range_ID]] %>% 
        pull(eid)
      
      curated_controls <- mapply(per_line_edit,control_input,SIMPLIFY = F) %>% 
        reduce(bind_rows) %>% 
        distinct(eid, .keep_all = T) %>% 
        mutate(any_code=0) %>% 
        filter(!eid %in% curated_cases$eid,
               !eid %in% range_ID_full)
      
      completed_phenotype <- curated_cases %>% 
        bind_rows(curated_controls)
      
    } else {
      
      completed_phenotype <- curated_cases
      
    }
    return(completed_phenotype)
  }
}
combining_codes <- function(x,y) {
  if(is.character(phenotypes[[x]]$earliest_date)) {
    
    concept_extracted <- phenotypes[[x]] %>% 
      rename(any_code=2) %>% 
      filter(any_code >= y) %>% 
      mutate(earliest_date=ymd(earliest_date))
  } else {
    concept_extracted <- phenotypes[[x]] %>% 
      rename(any_code=2) %>% 
      filter(any_code >= y)
  }
  return(concept_extracted)
}
AND <- function(a,b,c,d) {
  case <- a %>% 
    inner_join(b, by = "eid") %>% 
    mutate(earliest_date.x = ymd(earliest_date.x), earliest_date.y=ymd(earliest_date.y),
           gap = earliest_date.y-earliest_date.x,
           earliest_date=earliest_date.x,
           any_code=any_code.x+any_code.y) %>% 
    filter(between(gap,c,d)) %>% 
    select(eid,any_code,earliest_date)
}
NOT <- function(a,b,c,d)  {
  exclusions <- a %>% 
    inner_join(b, by = "eid") %>% 
    mutate(earliest_date.x = ymd(earliest_date.x), earliest_date.y=ymd(earliest_date.y),
           gap = earliest_date.y-earliest_date.x,
           earliest_date=earliest_date.x,
           any_code=any_code.x+any_code.y) %>% 
    filter(between(gap,c,d)) %>% 
    pull(eid)
  
  case <- a %>% 
    filter(!eid %in% exclusions) 
  
}
# Load in defining variables ----------------------------------------------
# define save location of the data
if(arguments$phenotype_save_file=="data/phenotypes/curated_phenotypes.RDS") {
  save_file_name <- here("data","phenotypes","curated_phenotypes.RDS")
} else {
  save_file_name <- arguments$phenotype_save_file
}
# create control populations to be used in defining curated phenotypes
if(arguments$control_populations=="data/control_populations") {
  control_populations <-  fread(here("data","control_populations"))
} else {
  control_populations <-  fread(arguments$control_populations)
}
# vector used later
list_names <- colnames(control_populations)[-1]
# function to create control populations used in curated_phenotypes
making_control_pop_lists <- function(x) {
  selected_population <- control_populations %>% 
    select(eid,{{x}}) %>% 
    mutate(earliest_date=ymd("2000-01-01")) %>% 
    drop_na() %>% 
    select(eid,{{x}},earliest_date)
  return(selected_population)
}
population_lists <- sapply(list_names,making_control_pop_lists,USE.NAMES = T,simplify = F)
if(arguments$phenotype_folder=="data/phenotypes/"){
  saveRDS(population_lists,here("data","phenotypes","control_populations.RDS")) 
} else {
  saveRDS(population_lists,paste0(arguments$phenotype_folder,"/","control_populations.RDS")) 
}
# number of times the function is iterated over
iteration_N <- as.numeric(arguments$N_iterations)
# read in the data from RDS files from the phenotypes folder and combine
if(arguments$phenotype_folder=="data/phenotypes/") {
  files <- paste0(here("data","phenotypes",list.files(here("data","phenotypes"),pattern = ".RDS")))
} else {
  files <- paste0(arguments$phenotype_folder,list.files(arguments$phenotype_folder,pattern = ".RDS"))
}
phenotypes <- files %>% 
  map(readRDS) %>% 
  reduce(c)
# read in curated_concept maps
if(arguments$curated_phenotype_map=="data/curated_phenotype_map.csv") {
  curated_phenotype_map <- fread(here("data","curated_phenotype_map.csv"), na.strings = "") %>% 
    mutate(across(everything(),as.character))
} else {
  curated_phenotype_map <- fread(arguments$curated_phenotype_map, na.strings = "") %>% 
    mutate(across(everything(),as.character))
}
# option to update existing list
if(arguments$update_list) {
  curated_pheno_list <- readRDS(save_file_name) 
} else{
  curated_pheno_list <- list()
}

# Run functions -----------------------------------------------------------
# some curated phenotypes require the initial creation of curated concepts or phenotypes before they can be created
for(i in 1:iteration_N) {
  curated_phenotype_map_edit <- curated_phenotype_map %>% 
    filter(!PheWAS_ID %in% names(curated_pheno_list),
           !PheWAS_ID %in% names(phenotypes))
  if(nrow(curated_phenotype_map_edit)>0) {
    # input is groups by PheWAS_ID
    PheWAS_ID_groups <- curated_phenotype_map_edit %>% 
      group_split(PheWAS_ID)
    names_PheWAS_ID <- lapply(PheWAS_ID_groups, function(x) input_name <- unique(x$PheWAS_ID)) %>% 
      reduce(c)
    # useful message on log files
    message(names_PheWAS_ID)
    message(length(names_PheWAS_ID))
    # function
    curated_phenotypes <- mapply(curated_phenotype_creation,
                                 names_PheWAS_ID,
                                 PheWAS_ID_groups, 
                                 SIMPLIFY = F,
                                 USE.NAMES = T)
    if(length(which(sapply(curated_phenotypes, is.null)))>0){
      curated_phenotypes = curated_phenotypes[-which(sapply(curated_phenotypes, is.null))]
    }
    # updating lists that will be saved
    curated_pheno_list <- append(curated_pheno_list,curated_phenotypes)
    phenotypes <- append(phenotypes,curated_phenotypes)
  } else {}
}
  # saved as R object
  saveRDS(curated_pheno_list,save_file_name)