#' get.peptide
#' @description Predicts the neopeptide produced for a given mutation and transcript combination
#'
#' @param ensembl_transcript_id The affected ensembl transcript ID
#' @param nuc_change HGVS nomenclature of DNA sequence variant e.g. c.4375C>T
#' @param build Which reference genome to use. Defaults to GRCh38.
#' @param check_startgains Performs checks of kozak strength, checks for overlapping WT uORF or in-frame overlap with coding sequence
#'
#' @return A data frame
#' @export
#'
#' @examples
#' get.peptide("ENST00000539214","c.-61C>T",build = 38,check_startgains = TRUE)
get.peptide <- function(ensembl_transcript_id, nuc_change, build = 38, check_startgains = F){

  peptides_msg <- c()

  changes <- extractNucChange(nuc_change) #parse the hgvs variation

  sequence <- NA

  #identify the location of sequence change
  if(!is.na(changes$start)&&!is.na(changes$stop)&&
     str_detect(changes$start,"^[0-9]+$")&&str_detect(changes$stop,"^[0-9]+$")){
    sequence <- "coding"
  } else if (!is.na(changes$start)&&!is.na(changes$stop)&&
             str_detect(changes$start,"^-[0-9]+$")&&str_detect(changes$stop,"^-[0-9]+$")){
    sequence <- "5utr"
  }

  #NULL outputs
  output <- list(peptide.seq = NA,start = NA,end = NA)

  check_startgain_res <- list(wt_kz = NA,
                              mut_kz = NA,
                              relative_strength_overlap = NA,
                              atg_count = NA,
                              overlap_present = NA,
                              overlap_pep = NA,
                              perc_overlap = NA,
                              relative_strength_native = NA,
                              cds_overlap_expected = NA,
                              overlap_cds_transcripts = NA)

  if(check_startgains == T){
    output <- c(output,check_startgain_res)
  }

  if(is.na(sequence)){
    peptides_msg <<- paste0(peptides_msg,"Incompatible CDS positions (",nuc_change,"). ")
    return(output)
  }

  #read the wt DNA sequence
  wt_sequence <- unlist(getFast(ensembl_transcript_id,sequence,build))

  #check wt_sequence against 3 conditions - if it fails any of the checks, try build 37
  check_sequence <- function(wt_sequence,return_warning = TRUE){
    if(is.null(wt_sequence)){
      if(return_warning==TRUE){
        peptides_msg <<- paste0(peptides_msg,"Transcript ",ensembl_transcript_id," not found")
      }
      return(FALSE)
    } else if (nchar(wt_sequence)<abs(changes$start)||nchar(wt_sequence)<abs(changes$stop)){
      if(return_warning==TRUE){
        peptides_msg <<- paste0(peptides_msg,"Mutation position ",changes$start," not found in ",sequence," sequence. ")
      }
      return(FALSE)
    } else if (str_detect(changes$wt,"[TCGA]+")&&changes$wt!=str_sub(wt_sequence,changes$start,changes$stop)){
      if(return_warning==TRUE){
        peptides_msg <<- paste0(peptides_msg,"CDS supplied is incorrect: ",changes$wt," vs ",str_sub(wt_sequence,changes$start,changes$stop),". ")
      }
      return(FALSE)
    } else {return(TRUE)}
  }

  if(build==37){
    runcheck <- check_sequence(wt_sequence,TRUE)
    if(runcheck==F){
      return(output)
    }
  }

  if(build==38&&check_sequence(wt_sequence,FALSE)==FALSE){
    #try Build 37, if that fails - return blank with warning message

    build = 37
    wt_sequence <- unlist(getFast(ensembl_transcript_id,sequence,37))

    recheck <- check_sequence(wt_sequence,TRUE)
    if(recheck==FALSE){
      return(output)
    } else {
      peptides_msg <<- paste0(peptides_msg,"Build 37 used.")
    }
  }


  #mutate the sequence
  mut_sequence <- wt_sequence

  if (str_detect(changes$wt,"[TCGA]+")&&str_detect(changes$mut,"[TCGA]+")){
    str_sub(mut_sequence,changes$start,changes$stop) <- changes$mut
  }

  if (str_detect(changes$wt,"del")||str_detect(changes$mut,"del")){
    str_sub(mut_sequence,changes$start,changes$stop) <-
      paste(rep("-",changes$stop-changes$start+1),collapse="")
  } #applies "-" as placeholder, to avoid changes in string length at this point

  correction = 0

  if (str_detect(changes$wt,"ins")||str_detect(changes$mut,"ins")){
    str_sub(mut_sequence,changes$start,changes$start) <-
      paste0(str_sub(mut_sequence,changes$start,changes$start),changes$mut)
    correction = 1
    if(str_detect(changes$wt,"delins")){
      correction = 0 #pos in indels refers to the deletion, so the first mutated nucleotide is not 1 forward
    }
  }

  if (changes$mut == "dup"||changes$mut==paste(rep(changes$wt,2),collapse="")){
    str_sub(mut_sequence,changes$start,changes$stop) <- paste0(rep(str_sub(wt_sequence,changes$start,changes$stop),2),collapse = "")
    correction = changes$stop-changes$start+1
  }

  #remove any "-" placeholders
  mut_sequence <- str_remove_all(mut_sequence,"-")

  #now find the neopeptide
  neoStart <- NULL
  neoEnd <- NULL
  frameshift = F

  #coding sequences
  if (sequence == "coding"){
    peptide.seq <- readPeptide(mut_sequence)
    #special case, inframe deletion
    #correction factor is -1 i.e. take the junctional epitope from
    #before and after del
    if(abs(nchar(mut_sequence)-nchar(wt_sequence))%%3==0&&changes$mut=="del"){
      correction = -1
      neoStart <- ceiling((changes$start+correction)/3)
      neoEnd <- ceiling((changes$start)/3)
    }

    #special case, stoploss
    if(nchar(mut_sequence)-changes$stop<=2&&
       !str_detect(peptide.seq,"\\*")){
      peptides_msg <<- paste0(peptides_msg,"Stop loss detected. ")
      neoStart <- ceiling((changes$start+correction)/3)
      frameshift = T #i.e. tell get.peptide to extend the ORF
    }

    #finally everything else
    if(is.null(neoStart)){
      neoStart <- ceiling((changes$start+correction)/3)

      if (abs(nchar(mut_sequence)-nchar(wt_sequence))%%3==0){#if inframe
        neoEnd <- ceiling((changes$start+correction+str_count(changes$mut,"[TCGA]")-1)/3)
        #str_count(changes$mut,"[TCGA]") the length of any insertion
        #this should work for any inframe substitution

        #special case, inframe duplications
        if(changes$mut == "dup"||changes$mut==paste(rep(changes$wt,2),collapse="")){
          neoStart <- ceiling((changes$start+correction-1)/3)
          neoEnd <- ceiling((changes$start+correction)/3) #cover the junction only
        }

      } else {frameshift = T}
    }

  }

  # 5'UTR mutations
  if(sequence == "5utr"){
    frameshift = T #All 5'UTR's referred to as frameshift i.e. not in frame
    length_wt <- nchar(wt_sequence)

    #As 5'utr HGVS uses numbering from -1 i.e. the tail
    #Need to convert to counting from the 5' end
    #This preserves the numbering of the start of the mut with respect to 5' end

    utr_start <- changes$start+length_wt+1 #first modified base
    utr_stop <- utr_start #last modified base

    utr_start <- utr_start+correction
    utr_stop <- utr_start+(nchar(changes$mut))-1

    if(changes$mut=="dup"||changes$mut==paste(rep(changes$wt,2),collapse="")){
      utr_start <- changes$stop+length_wt+1+1
      utr_stop <- utr_start+(changes$stop-changes$start)
    }

    #identify where the neoATG is formed within the mutant peptide
    mut_segment <- str_sub(mut_sequence,utr_start-2,utr_stop+2)
    atg_start <- utr_start+str_locate(mut_segment, "ATG")[1]-3

    if(is.na(atg_start)){
      peptides_msg <<- paste0(peptides_msg,"No ATG gain: ",mut_segment,". ")
      return(output)
    }

    utr_mut_seq <- str_sub(mut_sequence,atg_start)

    peptide.seq <- readPeptide(utr_mut_seq)
    neoStart = 1

    #reaches coding in-frame
    if(str_detect(peptide.seq,"\\(\\+0bp\\)")){ #i.e. reaches coding sequence in frame
      neoEnd = nchar(str_remove(peptide.seq,"\\(\\+0bp\\)")) #stop neoORF before coding sequence
    }
  }

  #Checking the start gains for kozak strength, uORF overlap and CDS overlap
  if(sequence == "5utr"&&check_startgains == T){

    #Subfunction to check kozak_strength
    kozak_strength <- function(input){
      if(!all(c(4,6)%in%unlist(str_locate_all(input,"ATG")))){
        warning("No ATG at pos 4-6")
        return(0)
      }
      if(nchar(input)!=7){
        warning("Input wrong length")
        return(NA)
      }

      score = 1
      if(grepl("A|G",str_sub(input,1,1))){
        score = score + 1
      }

      if(grepl("G",str_sub(input,7,7))){
        score = score + 1
      }

      return(score)

    }

    #Comparing to native WT start position
    wt_coding <- unlist(getFast(ensembl_transcript_id,"coding",build))

    if(!is.null(wt_coding)){
      wt_context <- paste0(str_sub(wt_sequence,-3),str_sub(wt_coding,1,4))

      check_startgain_res$wt_kz <- kozak_strength(wt_context)
    } else {peptides_msg <<- paste0(peptides_msg,"Unable to check native Kozak strength: ",ensembl_transcript_id,". ")}

    #Comparing to other uORFs
    #Other WT ATGs forming potential uORFs
    atgs <- data.frame(str_locate_all(wt_sequence,"ATG")[[1]])

    length <- nchar(wt_sequence)

    if(nrow(atgs)>0){
      atgs <- atgs %>% mutate(tail = start-length-1)
      atgs <- atgs %>% mutate(frame = (-tail)%%3+1)
      atgs <- atgs %>% mutate(pep = sapply(str_sub(wt_sequence,start),readPeptide))
      atgs <- atgs %>% mutate(context = str_sub(wt_sequence,pmax(start-3,0),start+3))
      atgs <- atgs %>% mutate(orflength = 3*nchar(str_remove(pep,"\\(\\+0bp\\)|\\(\\+[012]A\\)")))
      atgs <- atgs %>% mutate(kozak_str = sapply(context,kozak_strength))
      atgs <- atgs %>% mutate(end = start+orflength-1)
    }

    check_startgain_res$atg_count <-  nrow(atgs)

    #Same parameters for the neoORF
    length_mut <- nchar(mut_sequence)
    mut_tail <- atg_start-length_mut-1
    mut_frame <- (-mut_tail)%%3+1
    mut_context <- str_sub(mut_sequence,atg_start-3,atg_start+3)
    mut_end <- atg_start+(3*nchar(str_remove(peptide.seq,"\\(\\+0bp\\)|\\(\\+[012]A\\)")))-1

    check_startgain_res$mut_kz <- kozak_strength(mut_context)
    check_startgain_res$relative_strength_native <- check_startgain_res$mut_kz-check_startgain_res$wt_kz

    #Are there any uORFs overlapping and in the same frame as the mutant?
    check_startgain_res$overlap_present = F

    if(mut_frame%in%atgs$frame){

      atgs_framed <- atgs[atgs$frame==mut_frame,]

      #Calculate the overlap for each row in the data frame
      overlap <- pmin(atgs_framed$end, mut_end) - pmax(atgs_framed$start, atg_start)

      if(any(overlap>0)){
        check_startgain_res$overlap_present = T

        # find the row with the greatest overlap
        max_overlap <- which.max(overlap)

        #find the start of the overlap relative to mut peptide
        #utr_stop+1 accounts for one after the last mutated base
        overlap_start <- pmax((atgs_framed$start[max_overlap]-utr_stop+1)/3+1,1)
        check_startgain_res$overlap_pep <- str_sub(peptide.seq,overlap_start)
        no_overlap_length <- overlap_start-1
        check_startgain_res$relative_strength_overlap = check_startgain_res$mut_kz-atgs_framed$kozak_str[max_overlap]
      }

    }

    #Is there overlap with CDS exon in frame?
    if(str_detect(peptide.seq,"\\(\\+0bp\\)")){
      this.gene <- unlist(getFast(ensembl_transcript_id,"ensembl_gene_id",build))
      this.build <- paste0("build",build,".gff")
      isoforms <- get(this.build)$ensembl_transcript_id[which(get(this.build)$ensembl_gene_id == this.gene)]

      section <- get(this.build)[get(this.build)$ensembl_transcript_id%in%isoforms,]
      section_5utr <- section[section$ensembl_transcript_id==ensembl_transcript_id&section$type=="5utr",]

      if (section_5utr$strand[1] == "+"){
        section_5utr <- section_5utr[order(nrow(section_5utr):1),]
      }

      #recursive function to find the expected genomic position of the ATG start
      utr_cds_to_genomic <- function(position,section_5utr,row=1){
        if(is.list(position)){
          return(unlist(position))
        }

        if(row>nrow(section_5utr)){
          return("")
          warning("Could not locate genomic coordinates.")
        }

        find_pos <- function(cds_pos,start,stop,direction){
          if(direction == "+"){
            downstream = stop
            upstream = start
            mod = 1
          } else {
            downstream = start
            upstream = stop
            mod = -1
          }

          #using mut_tail countr from downstream
          gen_coord <- downstream+mod*(cds_pos+1)

          if(gen_coord-(mod*upstream)<0){
            return(gen_coord-(mod*upstream))
          }

          if(gen_coord-(mod*upstream)>=0){
            return(list(gen_coord = gen_coord))
          }

        }

        res <- find_pos(position,section_5utr$start[row],section_5utr$stop[row],direction=section_5utr$strand[1])

        return(utr_cds_to_genomic(res,section_5utr,row+1))

      }

      atg_genomic_coords <- utr_cds_to_genomic(mut_tail,section_5utr)

      section <- section[section$type=="coding",]

      overlap_cds_sec <- pmin(section$stop, atg_genomic_coords+length_mut-1) -
        pmax(section$start, atg_genomic_coords)
      #this is not perfect as length_mut does not take into account intronic
      #gaps, it will do for now

      if(any(overlap_cds_sec>0)){
        check_startgain_res$cds_overlap_expected <- T
      }

      check_startgain_res$overlap_cds_transcripts <- paste0(section$ensembl_transcript_id[overlap_cds_sec>0],collapse = ",")

    } else {
      check_startgain_res$cds_overlap_expected <- F
    }

  }


  #this if condition below is necessary to prevent get.peptide reading off into
  #the distance if native coding sequence is missing a stop codon
  if(frameshift == T){
    peptide.seq <- next.orf(peptide.seq,mut_sequence,sequence,ensembl_transcript_id,build)
  }

  if(is.null(neoEnd)){
    neoEnd <- nchar(str_remove(peptide.seq,"\\(\\+0bp\\)|\\(\\+[012]A\\)"))
  }

  output <- list(peptide.seq = peptide.seq,start = neoStart,end = neoEnd)


  # if statements check_startgains == T and sequence == "5utr" split
  # therefore if argument check_startgains == T, get.peptide will still return
  # the null outputs of check_startgain_res

  if(check_startgains == T){

    if(sequence == "5utr"){
      peptide_length = nchar(str_remove(peptide.seq,"\\(\\+0bp\\)|\\(\\+[012]A\\)"))

      if(!is.na(check_startgain_res$overlap_present)&&check_startgain_res$overlap_present == T){
        check_startgain_res$perc_overlap = (peptide_length-no_overlap_length)/(peptide_length)
      }
    }

    output <- c(output,check_startgain_res)

  }

  return(output)

}

# Helper --------------------
extractNucChange <- function(hgvs){
  #07/02/2023 added to deal with blanks in some nc transcripts
  if(is.na(hgvs)|hgvs==""){
    return(list(start = "",stop = "",wt = "",mut = ""))
  }

  pos <- unlist(str_extract_all(hgvs,"[0-9-+*]+"))
  change <- str_remove(hgvs,".*\\.") %>% str_extract_all(.,"[A-Z]+|[a-z]+") %>% unlist(.)
  start <- pos[1]
  stop <-  ifelse(!is.na(pos[2]),pos[2],pos[1])

  if(all(str_detect(change,"\\b[TCGA]+\\b")&length(change)>1)){
    wt <- change[1]
    mut <- change[2]
    return(list(start = as.numeric(start),stop = as.numeric(stop),wt = wt,mut = mut))
  } else if (any(change == "ins"|change == "del"|change == "dup"|change == "delins")) {
    if (any(change == "ins")){
      mut <- change[str_which(change,"ins")+1]
      wt <- "ins"
    }
    if (any(change == "del")){
      check <- change[str_which(change,"del")+1]
      wt <- ifelse(is.na(check)|!str_detect(check,"\\b[TCGA]+\\b"),".",check)
      if (!any(change == "ins")){
        mut <- "del"
      }
    }
    if (any(change == "delins")){
      mut <- change[str_which(change,"delins")+1]
      wt <- "delins"
    }
    if (any(change == "dup")){
      check <- change[str_which(change,"dup")+1]
      wt <- ifelse(is.na(check)|!str_detect(check,"\\b[TCGA]+\\b"),".",check)
      mut <- ifelse(is.na(wt)|!str_detect(wt,"\\b[TCGA]+\\b"),"dup",paste0(wt,wt))
    }
    return(list(start = as.numeric(start),stop = as.numeric(stop),wt = wt,mut = mut))
  } else {
    warning(paste0("Unrecognised sequence change: ",hgvs))
    return(list(start = "",stop = "",wt = "",mut = ""))}
}


readPeptide <- function(sequence,report = TRUE){

  if(is.null(sequence)){
    warning("Null sequence input - returning NULL")
    return(NULL)
  }

  #split sequence into triplets
  n <- seq(1, nc <- nchar(sequence), by = 3)
  codons <- str_sub(sequence,n,c(n[-1]-1, nc))

  #find first stop codon
  stop <- grep("TAA|TAG|TGA",codons)[1]
  if(!is.na(stop)){
    codons <- codons[1:stop]
  }

  peptide <- sapply(codons,function(codon){
    aa <- names(triplet.code[codon==triplet.code])

    if(length(aa)==1){
      return(aa)
    } else if(!str_detect(codon,"^[NCTGA]+$")){
      return("ERROR")
    } else if(str_detect(codon,"N")){
      return("X")
    } else if(nchar(codon)<3&nchar(codon)>0&report == T){
      return(paste0("(+",nchar(codon),"bp)"))
    } else if(report == F){
      return("")
    }
  })
  if(any(peptide=="ERROR")){
    warning("Unrecognised codon encountered")
    return("")
  } else {

    if(report == T){
      if(str_detect(peptide[length(peptide)],"^[XFLSYCWPHQRIMTNKVADEG]$")){
        peptide <- paste0(paste(peptide,collapse=""),"(+0bp)")
      }
    }

    peptide <- paste(peptide,collapse="")
    return(peptide)
  }
}

getFast <- function(ensembl_transcript_id,attributes = "all",build){
  outputlist <- list()
  buildname <- paste0("build",build,".gff")

  if (attributes %in% c("5utr","3utr","strand","start","stop","ensembl_gene_id","name","all")){
    section <- get(buildname)[which(get(buildname)$ensembl_transcript_id == ensembl_transcript_id),]
    if (nrow(section) == 0) {
      warning("Ensembl transcript not found")
      return(outputlist)
    }
    if (build == "38"){
      genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    } else if (build == "37") {
      genome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
    }
    if (section$strand[1] == "-"){
      section <- section[order(nrow(section):1),]
    }
  }

  if (any(attributes %in% c("all","5utr"))){
    sec_utr_5 <- section[grep("5utr",section$type,fixed = TRUE),]
    if(nrow(sec_utr_5>0)){
      five_prime_utr <- paste(getSeq(x = genome, names = paste0("chr",sec_utr_5$seqid[1]), start = sec_utr_5$start, end = sec_utr_5$stop, strand = sec_utr_5$strand[1],as.character = TRUE),collapse = "")
    } else {five_prime_utr <- c()}
    outputlist <- c(outputlist,utr5 = five_prime_utr)
  }

  if (any(attributes %in% c("all","coding"))){
    coding <- get(paste0("cds",build,".fasta"))$cds[get(paste0("cds",build,".fasta"))$ensembl_transcript_id==ensembl_transcript_id]
    if(length(coding)==0){coding <- c()}
    outputlist <- c(outputlist,coding = coding)
  }

  if (any(attributes %in% c("all","3utr"))){
    sec_utr_3 <- section[grep("3utr",section$type,fixed = TRUE),]
    if (nrow(sec_utr_3>0)){
      three_prime_utr <- paste(getSeq(x = genome, names = paste0("chr",sec_utr_3$seqid[1]), start = sec_utr_3$start, end = sec_utr_3$stop, strand = sec_utr_3$strand[1],as.character = TRUE),collapse = "")
    } else {three_prime_utr <- c()}
    outputlist <- c(outputlist,utr3 = three_prime_utr)
  }

  if (any(attributes %in% c("all","peptide"))){
    if(!exists("coding")){
      coding <- get(paste0("cds",build,".fasta"))$cds[get(paste0("cds",build,".fasta"))$ensembl_transcript_id==ensembl_transcript_id]
    }
    if (length(coding)>0){
      peptide <- coding %>% readPeptide(.,FALSE)
    } else {peptide <- c()}
    outputlist <- c(outputlist,peptide = peptide)
  }

  if (any(attributes %in% c("all","strand"))){
    outputlist <- c(outputlist,strand = section$strand[1])
  }

  if (any(attributes %in% c("all","start"))){
    outputlist <- c(outputlist,start = as.character(section[section$type == "mRNA",]["start"]))
  }

  if (any(attributes %in% c("all","stop"))){
    outputlist <- c(outputlist,stop = as.character(section[section$type == "mRNA",]["stop"]))
  }

  if (any(attributes %in% c("all","ensembl_gene_id"))){
    outputlist <- c(outputlist,ensembl_gene_id = as.character(section[section$type == "mRNA",]["ensembl_gene_id"]))
  }

  if (any(attributes %in% c("all","name"))){
    outputlist <- c(outputlist,name = as.character(section[section$type == "mRNA",]["name"]))
  }

  return(outputlist)
}

next.orf <- function(peptide,sequence,where,transcript,build){
  steps <- c("5utr","coding","3utr","polyA")
  #if stop codon reached, do nothing
  if (str_detect(peptide,"\\*")){
    return(peptide)
  }

  next_step <- steps[grep(where,steps)+1]
  numberleft_over <- as.numeric(str_extract(peptide,"\\d"))
  peptides_msg <<- paste0(peptides_msg,"ORF reaches ",next_step,". ")

  if (numberleft_over == 0){
    tail_bp = ""
  }

  if (numberleft_over > 0){
    tail_bp <- str_sub(sequence,-numberleft_over)
  }

  if (next_step == "polyA"){
    last_codon <- paste0(tail_bp,strrep("A",(3-numberleft_over)%%3))
    if (nchar(last_codon)==0){
      last_aa <- ""
    } else {
      last_aa <- readPeptide(last_codon,report = FALSE)
    }

    if (last_aa == "*"){
      return(paste0(str_remove(peptide,"\\(\\+[012]bp\\)"),"*"))
    } else {
      return(paste0(str_remove(peptide,"\\(\\+[012]bp\\)"),last_aa,"(+",(3-numberleft_over)%%3,"A)"))
    }
  }

  if (next_step == "coding"&numberleft_over == 0){
    peptides_msg <<- paste0(peptides_msg,"ORF is in-frame. ")
    next_peptide <- unlist(getFast(as.character(transcript),"peptide",build))
    return(paste0(str_remove(peptide,"\\(\\+[012]bp\\)"),next_peptide))
  }

  next_sequence <- unlist(getFast(as.character(transcript),next_step,build))

  #sometimes transcripts are missing UTR data
  if (length(next_sequence)==0){
    peptides_msg <<- paste0(peptides_msg,"No ",next_step," sequence available for ",transcript,". ")
    return("")
  }

  next_peptide <- readPeptide(paste0(tail_bp,next_sequence))

  combined_peptide <- paste0(str_remove(peptide,"\\(\\+[012]bp\\)"),next_peptide)

  return(next.orf(combined_peptide,next_sequence,next_step,transcript,build))

}
