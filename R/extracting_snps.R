#! /usr/bin/env Rscript
'
Uses bgenix or plink to extracts SNPs as a pgen file from either bgen, bed or pgen genetic files using the rsid and chromosome via inputted file SNP_list, see user guide for more details of file format of the SNP_list and plink and bgenix software. 

Firstly the user must edit the genetic_file_guide_template showing the location of the genetic files for each of the 23 autosomal and 2 sex chromosomes alongside the location of the corresponding sample/fam/psam file. If chromosome data is missing delete the row of missing chromosome. If genetic data is not stored per chromosome or separation by chromosome does not provide useful then different values can be used in the chromosome coloumn. To ensure the SNPs can be extracted ensure the same values are inputted into SNP_list in the chromosome column. If the values in the chromosome column in SNP_list match the values in the chromosome column in genetic_file_guide_template then the SNPs will be extracted.

The user must select one of plink_input or bgen_input depending on the type of file the genetic data is stored in and a directory that will store the results and temporary files created. The user must also have plink2 and bgenix installed on the system and executable as a single command in bash (normally achieved by adding the executable file directly into usr/bin or something similar) that can be altered from the default values using bgen_exe and plink_exe, see user guide for more information on plink2 and bgenix. The user can select a save name for the pgen file outputted via --variant_save_name=<name> and choose to keep the temporary bgen/plink files using the no_delete_temp flag. 

Usage:
    extracting_snps.R (--genetic_file_guide=<FILE> --SNP_list=<FILE> --analysis_folder=<FOLDER>) (--plink_input | --bgen_input) [--plink_exe=<text> --bgenix_exe=<text> --plink_type=<text> --ref_bgen=<text> --variant_save_name=<name> --no_delete_temp] 
    
Options:
    -h --help  Show this screen.
    -v --version  Show version.
    
    Mandatory inputs
    --genetic_file_guide=<FILE>           Full file path of the completed genetic_file_guide_template.csv
    
    --SNP_list=<FILE>                     Full file path for file containing SNPs to extract. Should be a text file with columns
                                          chromosome, SNP, group, coded_allele, non_coded_allele.
                                          
    --analysis_folder=<FOLDER>            Full path to the folder used to store and create extracted SNPs. Will create sub folder called 
                                          temp_plink within the folder which by default will be deleted upon successful extraction. 
                                          Will be created if it does not already exist.
    Select one of
    --bgen_input                          Input if using bgen files as genetic file input.
    --plink_input                         Input if using plink files as genetic file input.
    
    Options
    --plink_exe=<text>                    Bash command to execute plink2 program [default: plink2] 
    --plink_type=<text>                   if selecting plink the type of plink file, bed or pgen [default: bed]
    --ref_bgen=<text>                     Selecting a flag for bgen files in plink can either be ref-first, ref-last or ref-unknown [default: ref-first]
    --bgenix_exe=<text>                   Command to execute bgenix program [default: bgenix]
    --variant_save_name=<name>            Name of the extracted variants saved pgen files. [default: variants_for_association]
    --no_delete_temp                      Input if not wanting to delete the intermediatory plink and/or bgen files used when extracting SNPs.

' -> doc


suppressMessages(library(docopt))

arguments <- docopt(doc, version = 'v0.1 extracting_snps.R')

suppressMessages(library(data.table))
suppressMessages(library(tidyverse))
library(here)

# Functions ---------------------------------------------------------------
snp_selection <- function(x) {
  
  output_name <- paste0(arguments$analysis_folder,"/temp_plink/",paste0(x,"_temp"))
  
  genetic_file_per_chromo <- snp_guide %>% 
    filter(chromosome==x)
  
  genetic_file <- unique(genetic_file_per_chromo$genetic_file_location)
    
    if (arguments$bgen_input) {
      if(!dir.exists(paste0(arguments$analysis_folder,"/temp_bgen"))) {
        dir.create(paste0(arguments$analysis_folder,"/temp_bgen"),recursive = T)
      }
      temp_bgen <- paste0(arguments$analysis_folder,"/temp_bgen/",paste0(x,"_temp")) 
      bgen_file <- paste0(genetic_file,".bgen")
      ref_bgen <- paste0(arguments$ref_bgen)
      sample_file <-  unique(genetic_file_per_chromo$psam_fam_sample_file_location)
      
      system(paste0(bgenix_exe," -g ",bgen_file," -incl-rsids ",snp_list_file," > ",temp_bgen,"\n",
                    plink_exe," --bgen ",temp_bgen," ",ref_bgen," --sample ",sample_file," --make-pgen --out ",output_name))
      
    } else if (arguments$plink_input) {
      
      if(arguments$plink_type=="bed") {
        bed_file <- paste0(genetic_file,".bed")
        bim_file <- paste0(genetic_file,".bim")
        fam_file <- unique(genetic_file_per_chromo$psam_fam_sample_file_location)
        
        system(paste0(plink_exe," --bed ",bed_file," --bim ",bim_file," --fam ",fam_file," --extract ",snp_list_file," --make-pgen --no-pheno --out ",output_name))
        
      } else if(arguments$plink_type=="pgen") {
        pgen_file <- paste0(genetic_file,".pgen")
        pvar_file <- paste0(genetic_file,".pvar")
        psam_file <- unique(genetic_file_per_chromo$psam_fam_sample_file_location)
        
        system(paste0(plink_exe," --pgen ",pgen_file," --pvar ",pvar_file," --psam ",psam_file," --extract ",snp_list_file," --make-pgen --multiallelics-already-joined --no-psam-pheno --out ",output_name))
        
      }
      
    }
}
recoding_variants <- function(x) {
  
  new_pvar <- fread(x) %>%
    mutate(REF_number=match(str_sub(REF,1,1),LETTERS),
           ALT_number=match(str_sub(ALT,1,1),LETTERS),
           multi=case_when(nchar(REF) > 1 & nchar(ALT)== 1 ~ 1,
                           nchar(ALT) > 1 & nchar(REF)== 1 ~ 2,
                           nchar(REF) > 1 & nchar(ALT) <1 ~ 3),
           allele_order=ifelse(ALT_number>REF_number,paste0(REF,"_",ALT),
                               ifelse(REF_number>ALT_number,paste0(ALT,"_",REF),
                                      ifelse(REF_number==ALT_number,
                                             ifelse(multi==1,paste0(ALT,"_",REF),
                                                    ifelse(multi==2,paste0(REF,"_",ALT),
                                                           ifelse(nchar(REF)>nchar(ALT),paste0(ALT,"_",REF),paste0(REF,"_",ALT)))),NA))),
           ID=paste0(ID,"_",allele_order)) %>% 
    select(`#CHROM`,POS,ID,REF,ALT)
  
  fwrite(new_pvar,x,sep = "\t", na = NA,quote = F)
  
}

# Load in defining variables ----------------------------------------------
# reading in the SNPs to find
snp_search <- fread(arguments$SNP_list) %>% 
  mutate(chromosome=as.character(chromosome))
genetic_file_guide <- fread(arguments$genetic_file_guide) %>% 
  mutate(chromosome=as.character(chromosome))

# PLINK2 + begnix executable commands
plink_exe <- arguments$plink_exe
bgenix_exe <- arguments$bgenix_exe

# folder to store the temporary Plink files
if(!dir.exists(paste0(arguments$analysis_folder,"/temp_plink"))){
  dir.create(paste0(arguments$analysis_folder,"/temp_plink"),recursive = T)
}


# Run functions -----------------------------------------------------------
# selects SNPs from chromosomes
# acceptable chromosome coding
chrom_accept <- genetic_file_guide %>% 
  select(chromosome) %>% 
  pull()

# check for acceptable chromosome input
snp_filter <- snp_search %>% 
  filter(chromosome %in% chrom_accept)

# creates a file of SNPs in a format Plink can use to extract SNPs always goes in the same place in name analysis folder normally
snp_list <- snp_filter %>% 
  select(rsid)
fwrite(snp_list,paste0(arguments$analysis_folder,"/temp_plink/",paste0("snp_list.txt")),col.names = F)

# location of said file is required
snp_list_file <- paste0(arguments$analysis_folder,"/temp_plink/",paste0("snp_list.txt"))

# the required chromosomes to search through requires that genetic data is stored per chromosome
snp_guide <- snp_filter %>% 
  left_join(genetic_file_guide, by="chromosome")

chr_to_search <- unique(snp_guide$chromosome) 

# extract snps per chromosome
sapply(chr_to_search,snp_selection)

# pre_merge_psam edit to create unique ID for SNP input used in later analysis
pvar_names <- list.files(paste0(arguments$analysis_folder,"/temp_plink"),pattern = "temp.pvar")
pvar_file_loc <- lapply(pvar_names, function(x) paste0(arguments$analysis_folder,"/","temp_plink","/",x))

lapply(pvar_file_loc,recoding_variants)

# list all pgen files
plink_merge_list <- list.files(paste0(arguments$analysis_folder,"/temp_plink"),pattern = "_temp.pgen")

if(length(plink_merge_list)==1) {
  
  plink_name <- str_remove(unlist(plink_merge_list),".pgen")
  plink_location <- paste0(arguments$analysis_folder,"/temp_plink/",plink_name)
  plink_output_location <- paste0(arguments$analysis_folder,"/",arguments$variant_save_name)
  system(paste0(plink_exe," --make-pgen --pfile ",plink_location," --out ",plink_output_location))
  
} else {
#convert to df for saving
plink_merge_edit <- data.frame(paste0(paste0(arguments$analysis_folder,"/temp_plink/",(str_remove(plink_merge_list,".pgen")))))
fwrite(plink_merge_edit,paste0(arguments$analysis_folder,"/temp_plink/merge_list"),col.names = F)

merge_file_location <- paste0(arguments$analysis_folder,"/temp_plink/merge_list")
merge_output_location <- paste0(arguments$analysis_folder,"/",arguments$variant_save_name)

#perform merge
system(paste0(plink_exe," --pmerge-list ",merge_file_location," pfile --out ",merge_output_location)) 
}

if(arguments$no_delete_temp) {
  
} else {
unlink(paste0(arguments$analysis_folder,"/temp_plink"),recursive = T)
if(dir.exists(paste0(arguments$analysis_folder,"/temp_bgen"))) {
  unlink(paste0(arguments$analysis_folder,"/temp_bgen"),recursive = T)
}
}
