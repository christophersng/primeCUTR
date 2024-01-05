#' predict.orfs
#' @description Wrapper function to process single VCF input with PrimeCUTR
#'
#' @param single_input Path to a single VEP-annotated .vcf or .vcf.gz file
#' @param output Output path
#'
#' @return
#' @export
#'
#' @examples
predict.orfs <- function(single_input,output){

  #set identifier in the global environment
  identifier <<- "ensembl_transcript_id"

  #create log file
  current_time <- format(Sys.time(), "%Y_%b_%d_%H_%M_%S")
  dir.create(file.path(output),showWarnings = F)
  dir.create(file.path(output,"log"),showWarnings = F)

  logfile <- paste0("log/",study_id,"_",patient,"_log_",current_time,".txt")

  msg <- paste0(
    "My Own Run at ",Sys.time(),"\n\n",
    "#########################\n",
    "Run parameters\n",
    "#########################\n",
    "Patient - ",patient,"\n",
    "Study ID - ",study_id,"\n",
    "Input file - ",single_input,"\n",
    "Output directory - ", output,"\n",
    "Build - ", build,"\n"
  )

  write.table(msg,paste0(output,logfile),row.names = FALSE,col.names = FALSE, quote = FALSE, append = TRUE)

  #read in file
  cat("Reading file ",single_input,"\n")
  vcf <- try(read.delim(gzfile(single_input,"rt"),comment.char = '#',header = FALSE))
  if(class(vcf)=="try-error"){
    stop("Invalid file location provided.")
  }
  colnames(vcf)[1:8] <- c("Chr","bp","rsid","Ref","Alt","Score","Flag","Info")

  if("PASS"%in%vcf$Flag){
    vcf <- vcf[vcf$Flag=="PASS",]
  }

  vcf <-vcf[grep("ENST\\d{11}",vcf$Info),]

  vcf <- distinct(vcf)
  mutation_count <- nrow(vcf)
  transcript_count <- sum(str_count(vcf$Info,"ENST\\d{11}"))

  msg <- paste0("Transcripts identified by ",identifier,"\n",
                "Total mutation count - ",mutation_count,"\n",
                "Total transcripts - ",transcript_count,"\n\n",
                "#########################\n",
                "Mutations\n",
                "#########################\n")

  write.table(msg,paste0(output,logfile),row.names = FALSE,col.names = FALSE, quote = FALSE, append = TRUE)

  #filter out only variants of interest
  vcf <- vcf[grep("5_prime_UTR_premature_start_codon_gain_variant|5_prime_UTR_variant|stop_lost|frameshift_variant|missense_variant|inframe_insertion|inframe_deletion|protein_altering_variant",vcf$Info),]

  if (nrow(vcf)>0){
    extracted <- extract_vep(vcf)
    extracted <- extracted[grep("5_prime_UTR_premature_start_codon_gain_variant|5_prime_UTR_variant|stop_lost|frameshift_variant|missense_variant|inframe_insertion|inframe_deletion|protein_altering_variant",extracted$mutation_effect),]
    row.names(extracted) <- 1:nrow(extracted)
    #load("D:/work/scripts/PrimeCut/extracted_LUAD.rds")
  } else {extracted <- c()}

  cat("\n")

  #create table for the log
  mutation_tbl <- table(extracted$mutation_effect)
  msg <- mutation_tbl
  write.table(msg,paste0(output,logfile),row.names = FALSE,col.names = FALSE, quote = FALSE, append = TRUE, sep = " - ")

  dir.create(file.path(output,"orfs"),showWarnings = F)
  dir.create(file.path(output,"netMHC_inputs"),showWarnings = F)

  run.for <- function(group_name,mutation_effect){
    #set up counter
    group_orfs <- c()
    group_counter <- 0
    group_counter_filtered <- 0

    if(!is.null(extracted)){

      #go through extracted vcf transcript by transcript
      for (row in 1:nrow(extracted)){
        this.transcript <- extracted[row,]

        #reset all
        peptides_msg <<- ""
        this.orf <- ""
        peptide.normal.ends <- ""
        peptide.tryptic.ends <- ""
        this.peptide <- ""


        #Skip any blank rows
        if(any(is.na(c(this.transcript$Chr,this.transcript$Start,this.transcript$ensembl_transcript_id,this.transcript$nuc_change)))){
          warning(paste0("Missing value found in row extracted ",row," for ",patient))
          next
        }

        #Only run on the rows which match the mutation_effect e.g. startgain etc
        if (any(str_detect(this.transcript$mutation_effect,mutation_effect))){

          #Make a unique identifier
          key <- paste0(patient,"--",this.transcript$hgnc_symbol,"--",this.transcript[identifier],"--",this.transcript$Chr,":",this.transcript$Start,":",":",this.transcript$Ref,":",this.transcript$Alt,"--",this.transcript$nuc_change,"--mutant-sequence")

          #Progress bar
          cat(paste0("\r",group_name," - [",
                     paste(rep("=",ceiling(row/nrow(extracted)*70)),collapse=""),
                     paste(rep(" ",floor((1-row/nrow(extracted))*70)),collapse=""),
                     "]"))

          #Run get.peptide
          this.peptide <- get.peptide(this.transcript$ensembl_transcript_id,
                                      this.transcript$nuc_change,
                                      build = build,
                                      check_startgains = ifelse(group_name == "startgain",T,F))

          #extract the peptide sequence and
          if(!is.na(this.peptide$peptide.seq)&&this.peptide$start<=this.peptide$end){
            this.orf <- str_sub(this.peptide$peptide.seq,this.peptide$start,this.peptide$end)

            #hard code in number of missed cleavages
            peptide.tryptic.ends <-  trypsin.trim(this.peptide,missed = 1)
            peptide.normal.ends <- str_sub(this.peptide$peptide.seq,pmax(this.peptide$start-10,0),this.peptide$end+10)

            #if normal ends reaches polyA
            if ((str_detect(peptide.normal.ends,"\\(\\+[012]A\\)"))){
              peptide.normal.ends <- unlist(str_split(peptide.normal.ends,"\\("))[1]
              peptide.normal.ends <- paste0(peptide.normal.ends,paste(rep("K",10),collapse = ""))
            }

          } else {
            this.orf <- ""
            peptide.normal.ends <- ""
            peptide.tryptic.ends <- ""
          }

          group_counter <- group_counter+1
          pass <- FALSE

          if(length(unlist(str_match_all(this.orf, "^[XFLSYCWPHQRIMTNKVADEG]+")))>0){
            group_counter_filtered <- group_counter_filtered+1
            pass <- TRUE

            for (i in 9:11){
              this.peptide.snippet <- str_sub(this.peptide$peptide.seq,ifelse(this.peptide$start-i+1<1,1,this.peptide$start-i+1),this.peptide$end+i-1)

              if (nchar(str_extract(this.peptide.snippet, "^[XFLSYCWPHQRIMTNKVADEG]+"))>=i){
                mers <- get.mers(this.peptide.snippet,i,poly = TRUE)
                #generate as many IDs as length(mers)
                mers_ids <- mapply(paste0,list(">"),as.list(makeRandomID(length(mers)),15),rep(paste0("--",key),length(mers)), SIMPLIFY = FALSE)
                write.table(paste0(c(rbind(mers_ids,mers))),paste0(output,"netMHC_inputs/",study_id,"_",patient,"_",group_name,"_",i,"mer_peptides.",current_time,".txt"),row.names = FALSE,col.names = FALSE, quote = FALSE, append = TRUE)
              }
            }
          }

          if(group_name == "startgain"){

            group_orfs <- rbind(group_orfs,data.frame(
              key = key,
              hgnc_symbol = this.transcript$hgnc_symbol,
              transcript = this.transcript[identifier],
              peptide = this.orf %>% str_remove(.,"\\(\\+[012]A\\)"),
              peptide.tryptic.ends = peptide.tryptic.ends,
              peptide.normal.ends = peptide.normal.ends,
              pass = pass,
              peptides_msg = peptides_msg,
              wt_kz = this.peptide$wt_kz,
              mut_kz = this.peptide$mut_kz,
              relative_strength_native = this.peptide$relative_strength_native,
              atg_count = this.peptide$atg_count,
              overlap_present = this.peptide$overlap_present,
              overlap_pep = this.peptide$overlap_pep,
              perc_overlap = this.peptide$perc_overlap,
              relative_strength_overlap = this.peptide$relative_strength_overlap,
              cds_overlap_expected = this.peptide$cds_overlap_expected,
              overlap_cds_transcripts = this.peptide$overlap_cds_transcripts
            ))
          } else {
            group_orfs <- rbind(group_orfs,data.frame(
              key = key,
              hgnc_symbol = this.transcript$hgnc_symbol,
              transcript = this.transcript[identifier],
              peptide = this.orf %>% str_remove(.,"\\(\\+[012]A\\)"),
              peptide.tryptic.ends = peptide.tryptic.ends,
              peptide.normal.ends = peptide.normal.ends,
              pass = pass,
              peptides_msg))
          }

        }
      }
    }
    cat(paste0("\r",group_name," - [",
               paste(rep("=",70),collapse=""),
               paste(rep(" ",0),collapse=""),
               "]"))
    msg <- paste0("\n\n#########################\n",
                  group_name,"\n",
                  "#########################\n",
                  group_name," filtered from: ",group_counter," transcripts to ",group_counter_filtered," transcripts.\n\n")
    write.table(msg,paste0(output,logfile),row.names = FALSE,col.names = FALSE, quote = FALSE, append = TRUE)

    write.table(group_orfs,paste0(output,"orfs/",paste0(study_id,"_",patient,"_",group_name,"ORFs_",current_time,".tsv")),row.names = FALSE,col.names = TRUE, quote = FALSE, sep = "\t")
  }

  run.for("stoploss","stop_lost")
  cat("\n")
  run.for("startgain",c("5_prime_UTR_premature_start_codon_gain_variant","5_prime_UTR_variant"))
  cat("\n")
  run.for("missense",c("missense_variant","inframe_insertion","inframe_deletion","protein_altering_variant"))
  cat("\n")
  run.for("frameshift",c("frameshift_variant"))
  cat(paste0("\n",patient," - DONE\n"))
}
