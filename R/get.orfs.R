#' get.orfs
#' @description Wrapper function to process single VCF input with PrimeCUTR.
#' Writes three folders to the output directory: ./log/ ./netMHC_inputs/ ./orfs/
#'
#' @param single_input Path to a single VEP-annotated .vcf or .vcf.gz file
#' @param output Output path
#' @param build Reference genome - 37 or 38
#' @param filter Filter vcf FILTER column by a string match e.g. "PASS". Default: NULL
#' @param patient String for patient ID. Default: NULL - will take basename of vcf file
#' @param study_id String for study_id. Default: ""
#'
#' @return returns the log output silently
#' @export
#'
#' @examples
#' input_vcf <- primeCUTR_example("vep_HCC1937.vcf")
#' get.orfs(input_vcf,"./output_dir/",build=38)
get.orfs <- function(single_input,output,build,filter=NULL,patient=NULL,study_id=""){

  identifier <- "ensembl_transcript_id"

  if(is.null(patient)){
    patient <- basename(single_input)
  }

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
  con <- gzfile(single_input,"rt")
  vcf <- try(read.delim(con,comment.char = '#',header = FALSE))
  close(con)

  if(class(vcf)=="try-error"){
    stop("Invalid file location provided.")
  }

  colnames(vcf)[1:8] <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO")

  if(!is.null(filter)){
    orig.rows <- nrow(vcf)
    vcf <- vcf[vcf$FILTER==filter,]
    filt.rows <- nrow(vcf)
    cat("Filtered vcf from,",orig.rows," to ",filt.rows," by FILTER==",filter,".")
  }

  vcf <-vcf[grep("ENST\\d{11}",vcf$INFO),]

  vcf <- distinct(vcf)
  mutation_count <- nrow(vcf)
  transcript_count <- sum(str_count(vcf$INFO,"ENST\\d{11}"))

  msg <- paste0("Transcripts identified by ",identifier,"\n",
                "Total mutation count - ",mutation_count,"\n",
                "Total transcripts - ",transcript_count,"\n\n",
                "#########################\n",
                "Mutations\n",
                "#########################\n")

  write.table(msg,paste0(output,logfile),row.names = FALSE,col.names = FALSE, quote = FALSE, append = TRUE)

  #filter out only variants of interest
  vcf <- vcf[grep("5_prime_UTR_premature_start_codon_gain_variant|5_prime_UTR_variant|stop_lost|frameshift_variant|missense_variant|inframe_insertion|inframe_deletion|protein_altering_variant",vcf$INFO),]

  if (nrow(vcf)>0){
    extracted <- extract_vep(vcf)
    extracted <- extracted[grep("5_prime_UTR_premature_start_codon_gain_variant|5_prime_UTR_variant|stop_lost|frameshift_variant|missense_variant|inframe_insertion|inframe_deletion|protein_altering_variant",extracted$mutation_effect),]
    row.names(extracted) <- 1:nrow(extracted)

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
        this.orf <- ""
        peptide.normal.ends <- ""
        peptide.tryptic.ends <- ""
        this.peptide <- ""


        #Skip any blank rows
        if(any(is.na(c(this.transcript$CHROM,this.transcript$Start,this.transcript$ensembl_transcript_id,this.transcript$nuc_change)))){
          warning(paste0("Missing value found in row extracted ",row," for ",patient))
          next
        }

        #Only run on the rows which match the mutation_effect e.g. startgain etc
        if (any(str_detect(this.transcript$mutation_effect,mutation_effect))){

          #Make a unique identifier
          key <- paste0(patient,"--",this.transcript$hgnc_symbol,"--",this.transcript[identifier],"--",this.transcript$CHROM,":",this.transcript$Start,":",":",this.transcript$REF,":",this.transcript$ALT,"--",this.transcript$nuc_change,"--mutant-sequence")

          print(key)

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
              peptides_msg = this.peptide$peptides_msg,
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
              peptides_msg = this.peptide$peptides_msg))
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

  invisible(logfile)
}


#Helpers --------------
extract_vep <- function(vcf){
  outputtbl <- c()

  check_index <- function(vcf,sample_size=20,loop=0){
    #function to get indexing from the info column

    if(sample_size>nrow(vcf)){
      sample_size = nrow(vcf)
    }

    test <- vcf$INFO[sample(nrow(vcf),sample_size)]
    test <- unlist(str_split(test,";"))
    test <- unlist(str_split(test,","))

    #try to index
    index <- lapply(test,function(line){
      vectors <- str_split(line,"\\|")[[1]]
      variant_loc <- grep("variant|stop|feature|deletion|insertion",vectors)[1]
      enst_loc <- grep("^ENST[0-9\\.]+$",vectors)[1]
      ensg_loc <- grep("ENSG\\d{11}",vectors)[1]
      hgvs_loc <- grep(">|n\\.|c\\.",vectors)[1]
      seq_type_loc <- grep("protein_coding|nonsense_mediated_decay|retained_intron|processed_transcript",vectors)[1]
      return(c(variant_loc,enst_loc,ensg_loc,hgvs_loc,seq_type_loc))
    })

    index <- index[lapply(index,length)==5]

    index.table <- suppressWarnings(data.frame(do.call(rbind,index)))

    #recursively try again if fail
    if(nrow(index.table)==0|ncol(index.table)!=5){

      loop=loop+1

      if(loop>1){
        stop("Unable to extract indexing for INFO column.")
      } else {

        return(check_index(vcf,sample_size = 1000,loop=loop))

      }
    }

    #return the most frequent index
    Mode <- function(vec){
      tab <- table(vec)
      as.numeric(names(tab[tab==max(tab)]))
    }

    return(list(variant_loc = Mode(index.table[1]),
                enst_loc = Mode(index.table[2]),
                ensg_loc = Mode(index.table[3]),
                hgvs_loc = Mode(index.table[4]),
                seq_type_loc = Mode(index.table[5])))
  }

  vcf.index <- check_index(vcf)

  for (i in 1:nrow(vcf)) {
    if (i%%1000 == 0|i == nrow(vcf)){
      cat(paste0("\r extracting row ",i,"/",nrow(vcf)))
    }

    row_data <- vcf[i,]

    #uses the vcf.index where possible
    infoblock <- unlist(str_split(row_data$INFO,";"))
    info <- infoblock[grep("ENSG\\d{11}",infoblock)] %>% str_split(.,"\\|") %>% unlist(.)
    loc <- grep("ENSG\\d{11}",info)

    if(length(loc)==0){
      #some rows lack transcript data in the info cols, skip over
      next
    }

    ensembl_gene_id <- info[loc]
    ensembl_transcript_id <- info[loc+(vcf.index[["enst_loc"]]-vcf.index[["ensg_loc"]])] %>% str_extract(.,"^ENST\\d{11}") %>% unlist(.)
    mutation_effect <- info[loc+(vcf.index[["variant_loc"]]-vcf.index[["ensg_loc"]])]

    hgnc_symbol <- info[loc-1]
    seq_type <- info[loc+(vcf.index[["seq_type_loc"]]-vcf.index[["ensg_loc"]])]
    nuc_change <- info[loc+(vcf.index[["hgvs_loc"]]-vcf.index[["ensg_loc"]])] %>% str_remove(".*:")

    transcript_extract <- cbind(ensembl_transcript_id,
                                mutation_effect,
                                hgnc_symbol,
                                ensembl_gene_id,
                                seq_type,
                                nuc_change)

    transcript_extract <- cbind(row_data,transcript_extract, row.names = NULL)
    outputtbl <- rbind(outputtbl,transcript_extract)
  }
  outputtbl <- dplyr::rename(outputtbl, Start = POS)
  rownames(outputtbl) <- 1:nrow(outputtbl)

  return(outputtbl)
}

trypsin.trim <- function(this.peptide,missed=0){

  if(!is.list(this.peptide)){
    return("")
  }

  n.ends <- str_locate_all(this.peptide$peptide.seq,"R(?!P)|K(?!P)")[[1]][,1]+1

  up.flanks <- n.ends[which(n.ends<=this.peptide$start)]

  if(length(up.flanks)>0){
    this.peptide$start <- up.flanks[ifelse(length(up.flanks)-missed<1,1,length(up.flanks)-missed)]
  } else {this.peptide$start <- 1}

  c.ends <- str_locate_all(this.peptide$peptide.seq,"R(?!P)|K(?!P)")[[1]][,1]

  down.flanks <- c.ends[which(c.ends>=this.peptide$end)]

  if(length(down.flanks)>0){
    this.peptide$end <- down.flanks[ifelse(1+missed>length(down.flanks),length(down.flanks),1+missed)]
  } else {this.peptide$end <- nchar(this.peptide$peptide.seq)}

  return(str_sub(this.peptide$peptide.seq,this.peptide$start,this.peptide$end))
}


get.mers <- function(peptide,length,poly = FALSE){
  peptide <- str_remove_all(peptide,"\\{[\\d]+\\}")

  if (str_detect(peptide,"\\*")){
    peptide <- unlist(str_split(peptide,"\\*"))[1]
  }

  #if the transcript has reached the poly A tail, add a poly-lysine (11) if poly = TRUE
  else if ((str_detect(peptide,"\\(\\+[012]A\\)"))){
    peptide <- unlist(str_split(peptide,"\\("))[1]
    if (poly == TRUE){
      peptide <- paste0(peptide,paste(rep("K",length),collapse = ""))
    }
  }

  #how many peptides will be produced? nchar(peptide)-length+1
  #create empty list holder based on how many mers will be created
  output <- vector(mode = "list",length = nchar(peptide)-length+1)

  for(aa in 1:length(output)){
    output[[aa]] <- str_sub(peptide,aa,aa+length-1)
  }
  return(output)
}

makeRandomID <- function(n=1, len=12){
  randomString <- c(1:n)
  for (i in 1:n)
  {
    randomString[i] <- paste(sample(c(0:9, letters, LETTERS),
                                    len, replace=TRUE),
                             collapse="")
  }
  return(randomString)
}
