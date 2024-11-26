# 
# BiocManager::install("Biostrings")
# ##install.packages("Peptides", dependencies=TRUE)
#BiocManager::install("dagLogo")
# BiocManager::install("Rdisop")
# install.packages("BiocManager")
# install.packages("protti")
# install.packages("ggridges")
# install.packages("VennDiagram")
# install.packages("hrbrthemes")
# install.packages("hexbin")
# install.packages("plotly")
#install.packages("ggseqlogo")
# install.packages("wesanderson")
#BiocManager::install("rrvgo")
# install.packages("dagLogo")
# install.packages("dbplyr")
# 
# library(devtools)
# devtools::install_github("RobinHankin/Brobdingnag")
# install.packages("Brobdingnag")
save.image()

#libraries
{
  library(tidyverse)
  #library(ggridges)
  library(devtools)
  library(Biostrings)
  library(Peptides)
  library(ggpubr)
  #library(VennDiagram)
  library(RColorBrewer)
  #library(hrbrthemes)
  #library(hexbin)
  library(ggseqlogo)
  library(dagLogo)
  #library(plotly)
  library(wesanderson)
  library(Rdisop)
  library(conflicted)
  library(broom)
  #library(seqRFLP)
  library(OrgMassSpecR) 
  library(scales)
  conflict_scout()
  conflicts_prefer(dplyr::filter())
  conflicts_prefer(dplyr::select())
  conflicts_prefer(dplyr::rename())
  conflicts_prefer(dplyr::first())
  conflicts_prefer(dplyr::desc())
  source("G:/Rcheatsheet/functions/FragmentPeptide2.R")
  library(tidymodels)
  library(ComplexHeatmap)
  library(ggpointdensity)
  library(viridis)
  library(ggVennDiagram)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  Sys.setenv(R_LIBCURL_SSL_REVOKE_BEST_EFFORT=TRUE)
  library(ReactomePA)
  library(rrvgo)
  library(UniProtKeywords)
  library(org.Hs.eg.db)
  library(ggrepel)
}

summaryfolder<-"./Summary/"
peptidefile <- "200821_HELA_MGTG_iNrich_bRP_F_PeptideGroups.txt"
psmfile <- "200821_HELA_MGTG_iNrich_bRP_F_PSMs.txt"
msmsfile <- "200821_HELA_MGTG_iNrich_bRP_F_MSMSSpectrumInfo.txt"

###DB section ##object:: dbtable
{
  ####DB collection
  dbfile <- paste0(summaryfolder,"UP000005640_9606_2302.fasta")
  readAAStringSet(dbfile)->db1 ##Biostrings package
  protein = names(db1)
  sequence = paste(db1)
  data.frame(protein,sequence)->dbtable
  rm(db1)
  dbtable %>% mutate(
    accession = str_extract(protein,"(?<=\\|).+(?=\\|)")
  ) -> dbtable ##db extraction of accession
  dbtable %>% mutate(
    gene = str_extract(protein,"(?<=GN\\=)[a-zA-Z0-9\\-]+")
  ) %>% 
    mutate(
      prot_name = str_extract(protein,"[A-Z0-9]+\\_HUMAN")
    )-> dbtable ##db extraction of gene
  
  
  
  ####DB collection
  dbfile_1 <- paste0("./db/UP000005640_9606_2302.fasta")
  dbfile_2 <- paste0("./db/UP000005640_9606_2302_additional.fasta")
  readAAStringSet(dbfile_1)->db1 ##Biostrings package
  protein = names(db1)
  sequence = paste(db1)
  data.frame(protein,sequence)->dbtable1
  readAAStringSet(dbfile_2)->db2 ##Biostrings package
  protein = names(db2)
  sequence = paste(db2)
  data.frame(protein,sequence)->dbtable2
  bind_rows(dbtable1,dbtable2) -> dbtable_all
  rm(db1,db2)
  
  
  dbtable_all %>% 
    mutate(
      accession = str_extract(protein,"(?<=\\|).+(?=\\|)")
    ) %>% mutate(
      gene = str_extract(protein,"(?<=GN\\=)[a-zA-Z0-9\\-]+")
    ) %>% 
    mutate(
      prot_name = str_extract(protein,"[A-Z0-9]+\\_HUMAN")
    ) -> dbtable_all ##db extraction of gene
  
  rm(dbtable1,dbtable2)
  
  
  
}
##dbtable
##dbtable_all

### read, nmod initiation ##object:: summaryfile
{
  ###Read peptide file
  read_tsv(paste0(summaryfolder,peptidefile)) -> summaryfile
  
  summaryfile %>% distinct(Sequence,Modifications) ## check seq+mod uniqueness 
  
  ##N-term peps:: deconvolution of protein accessions 220407
  summaryfile %>% 
    separate_rows(.,`Positions in Master Proteins`,sep = "; ") -> summaryfile ##deconvolution of accessions
  
  summaryfile %>% 
    arrange(Sequence,Modifications,`Positions in Master Proteins`) -> summaryfile #sort by postion with increasing order
  
  summaryfile %>%
    distinct(Sequence,Modifications, .keep_all = TRUE) -> summaryfile #delete redundant pep+mod
  
  
  ##protein accession, separate first only ## deprecated 220407
  ##summaryfile %>% separate(`Master Protein Accessions`,into="Accessions.Single",sep=";",remove=FALSE)-> summaryfile ###delimiter
  
  ##protein position
  #summaryfile %>% separate(`Positions in Master Proteins`, into="Position", sep=";", remove=FALSE)-> summaryfile ###Position string extraction
  
  ##protein position string
  summaryfile %>% mutate(Accessions.Single = str_extract(`Positions in Master Proteins`, "[a-zA-Z0-9]+(?=[:blank:])")) %>%
    mutate(Position.start = str_extract(`Positions in Master Proteins`,"[0-9]+\\-[0-9]+")) %>% 
    separate(`Position.start`, into=c("Position.start","Position.end"), sep="-", remove=FALSE, convert=TRUE)->summaryfile ###position number extraction
  
  ##nmod start
  summaryfile %>% mutate(nmod = str_extract(Modifications,"1x[a-zA-Z0-9\\>\\:\\-\\(\\)\\+]+ \\[N\\-Term\\]"))->summaryfile ###nmod extraction
  summaryfile %>% mutate(nmod = gsub("1xAcetyl [N-Term]","Acetyl", nmod, fixed=TRUE)) %>%
    mutate(nmod = gsub("1xAcetyl:2H(3) [N-Term]","D3Acetyl", nmod, fixed=TRUE)) %>% 
    mutate(nmod = gsub("1xGlu->pyro-Glu [N-Term]","pyroGlu", nmod, fixed=TRUE)) %>%
    #mutate(nmod = gsub("1xMet-loss+AcetylD3 [N-Term]","MetLossD3Acetyl", nmod, fixed=TRUE)) %>%
    mutate(nmod = gsub("1xArgAcetylD3 [N-Term]","D3AcetylArg", nmod, fixed=TRUE)) %>%
    mutate(nmod = gsub("1xArgAcetylD3NQ [N-Term]","D3AcetylArgDeamid", nmod, fixed=TRUE)) ->summaryfile ##unique(summaryfiletemp$nmod) ## simplifying nmods
  #mutate(nmod = gsub("1xMet-LossArgAcD3NQ [N-Term]","MetLossD3AcetylArgDeamid", nmod, fixed=TRUE)) %>%
  #mutate(nmod = gsub("1xMet-lossArgAcD3 [N-Term]","MetLossD3AcetylArg", nmod, fixed=TRUE)) %>%
  #mutate(nmod = gsub("1xMet-loss [N-Term]",NA, nmod, fixed=TRUE))
  
  
  ##filter contaminant & peptide without reference
  summaryfile %>% 
    dplyr::filter(Contaminant == FALSE) %>%
    dplyr::filter(!is.na(`Master Protein Accessions`)) -> summaryfile
  
  ### db integration 
  left_join(
    summaryfile,dbtable,by=c("Accessions.Single"="accession")
  ) ->summaryfile ##db integration
  
  ### XXX addition
  summaryfile %>% 
    mutate(
      sequence2 = paste0("XXXXX",sequence)
    ) %>%
    mutate(
      p5p1 = str_sub(sequence2,Position.start,Position.start+4)
    ) ->summaryfile
}

## modification setup ##::summaryfile
{
  summaryfile %>%
    mutate(
      arginylation = ifelse(nmod %in% c("D3AcetylArg","D3AcetylArgDeamid"),TRUE,FALSE)
    )-> summaryfile ###arginylated?
  
  summaryfile %>% 
    mutate(
      arginylation_DE = ifelse(str_sub(Sequence,1,1) %in% c("D","E") & nmod == "D3AcetylArg",TRUE,FALSE)
    )-> summaryfile ###arginylated before D,E
  summaryfile %>% 
    mutate(
      arginylation_NQ = ifelse(str_sub(Sequence,1,1) %in% c("N","Q") & nmod == "D3AcetylArgDeamid",TRUE,FALSE)
    )-> summaryfile ###arginylated before N,Q 
  summaryfile %>% 
    mutate(
      arginylation_DENQ = ifelse(arginylation_DE == TRUE | arginylation_NQ == TRUE,TRUE,FALSE)
    )-> summaryfile ###arginylated before D,E,N,Q
  summaryfile %>% 
    mutate(p3p1 = str_sub(p5p1,3,5)) %>%
    mutate(p2p1 = str_sub(p5p1,4,5)) %>%
    mutate(p1p1 = str_sub(p5p1,5,5)) -> summaryfile
  summaryfile  %>%
    mutate(p5p5prime = paste0(
      p5p1, str_sub(Sequence,1,5)
    )
    )-> summaryfile
  summaryfile  %>%
    mutate(p4p4prime = paste0(
      str_sub(p5p1,2,5), str_sub(Sequence,1,4)
    )
    )-> summaryfile
  summaryfile %>% 
    mutate(
      nsite = paste0(Accessions.Single,"_",Position.start,"_",nmod)
    ) -> summaryfile 
  write_tsv(summaryfile,"./Summary/summaryfile_mod_setup.tsv") ##date
}

##for true positive arginylome (GV, VG), (mass ambiguity >> skip)
{
  summaryfile %>%
    mutate(
      arginylation_decision = ifelse(
        p2p1 != "VG" & p2p1 != "GV" & p1p1 != "R" & nmod == "D3AcetylArg",
        "D3Arg",
        ifelse(
          p2p1 != "VG" & p2p1 != "GV" & p1p1 != "R" & nmod == "D3AcetylArgDeamid",
          "D3ArgNQ",
          NA
        )
      )
    ) -> summaryfile
}

##psm file setup ##summaryfile_psm, summaryfile_psm_arg
{
  read_tsv(paste0(summaryfolder,psmfile)) -> summaryfile_psm
  summaryfile_psm %>% 
    mutate(
      index = seq(1:nrow(summaryfile_psm))
    ) %>% 
    separate(`Master Protein Accessions`,into="Accessions_Single",sep=";",remove=FALSE) %>% 
    mutate(
      nmod = str_extract(
        Modifications,
        "N\\-Term[a-zA-Z0-9\\>\\:\\-\\(\\)\\+]+"
      )
    )->summaryfile_psm ###nmod extraction
  
  summaryfile_psm %>% mutate(
    nmod = gsub("N-Term(Acetyl)","Acetyl", nmod, fixed=TRUE)
  ) %>% mutate(
    nmod = gsub("N-Term(Acetyl:2H(3))","D3Acetyl", nmod, fixed=TRUE)
  ) %>% mutate(
    nmod = gsub("N-Term(Glu->pyro-Glu)","pyroGlu", nmod, fixed=TRUE)
  ) %>% mutate(
    nmod = gsub("N-Term(ArgAcetylD3)","D3AcetylArg", nmod, fixed=TRUE)
  ) %>% mutate(
    nmod = gsub("N-Term(ArgAcetylD3NQ)","D3AcetylArgDeamid", nmod, fixed=TRUE)
  ) ->summaryfile_psm ##unique(summaryfiletemp$nmod) ## simplifying nmods
  
  summaryfile_psm %>% 
    select(
      `File ID`,`Spectrum File`
    ) %>% distinct()->file_raw_index
  
  file_raw_index %>% 
    mutate(
      set = str_extract(`File ID`,"^[A-Z0-9]+(?=\\.)")
    ) %>% 
    mutate(
      fraction = str_extract(`File ID`,"(?<=\\.)[A-Z0-9]+") %>% as.integer()
    ) %>% 
    mutate(
      file_name = str_extract(`Spectrum File`,".+(?=.{6}$)")
    ) -> file_raw_index
  
  file_raw_index %>% 
    distinct(set,file_name) %>% 
    bind_cols(
      tibble(
        experiment = c("MGTG","MOCK","MOCK","MG132","MGTG","MG132","MGTG","MOCK","MOCK","MG132","MGTG","MG132")
      )
    ) %>% 
    bind_cols(
      tibble(
        replicate = c("1ST","1ST","2ND","2ND","2ND","1ST","1ST","1ST","2ND","2ND","2ND","1ST")
      )
    ) %>% 
    bind_cols(
      tibble(
        enzyme = c("trypsin","trypsin","trypsin","trypsin","trypsin","trypsin","chymotrypsin","chymotrypsin","chymotrypsin","chymotrypsin","chymotrypsin","chymotrypsin")
      )
    ) -> file_raw_index2
  
  file_raw_index %>% 
    left_join(
      file_raw_index2,
      join_by(
        set==set,
        file_name == file_name
      )
    ) %>% 
    select(
      -c(`Spectrum File`,file_name)
    ) ->file_raw_index
  
  summaryfile_psm %>% 
    left_join(
      file_raw_index,
      join_by(
        `File ID` ==`File ID`
      )
    ) -> summaryfile_psm
  rm(file_raw_index2)
  
  ##make arginylation column
  summaryfile_psm %>% mutate(
    arginylation = ifelse(
      grepl("Arg", nmod),
      TRUE,
      FALSE
    )
  ) -> summaryfile_psm
  
  summaryfile_psm %>% 
    mutate(
      start_with_R = if_else(str_detect(Sequence,"^R"),1,0)
    ) -> summaryfile_psm
  
  summaryfile_psm %>% 
    dplyr::filter(
      arginylation == TRUE
    ) -> summaryfile_psm_arg
}

###peptide new position, from::summaryfile_psm to::reposition_arg_no_dups, reposition_arg
{
  
  
  summaryfile %>% 
    dplyr::filter(arginylation == TRUE) -> summaryfile_arg
  
  summaryfile_psm %>% dplyr::filter(arginylation==TRUE) %>% .$Sequence %>% unique -> arginylated_peptides_all
  
  tibble(
    seq = arginylated_peptides_all
  ) %>% 
    mutate(
      protein = map(seq,\(x) dbtable_all %>% dplyr::filter(str_detect(sequence,x)),.progress=TRUE)
    ) %>% 
    unnest(protein) %>% 
    mutate(
      position = str_locate(sequence,seq)[,"start"]
    ) %>% 
    arrange(
      seq,position
    )  %>% 
    mutate(
      sequence2=paste0("XXXXX",sequence)
    ) %>% 
    mutate(
      p5p1 = str_sub(sequence2,start=position,end=position+4)
    ) %>% 
    mutate(
      reviewed = str_sub(protein,start=1,end=2)
    ) %>% 
    mutate(
      under100 = if_else(
        position <= 100, 1,0
      )
    ) %>% 
    arrange(
      seq,desc(under100),reviewed,position,accession
    ) %>% 
    mutate(
      count = n(),.by=seq
    ) ->reposition_arg
  reposition_arg %>% 
    nest(
      .by=seq
    ) %>% 
    mutate(
      data = map(data,\(x) slice_head(x,n=1))
    ) %>% 
    unnest(data)->reposition_arg_no_dups
  
  
  reposition_arg %>%
    select(-c(sequence,sequence2)) %>% 
    write_tsv("reposition_arg.tsv")
  reposition_arg_no_dups %>%
    select(-c(sequence,sequence2)) %>% 
    write_tsv("reposition_arg_no_dups.tsv")
}

###prepare training file ::summaryfile_psm >> summaryfile_psm_for_training
{
  summaryfile_psm %>% 
    dplyr::filter(start_with_R==1) %>%
    #dplyr::filter(enzyme == "trypsin") %>% for trypsin only
    #dplyr::filter(specfile %in% c("220621_E_1","220621_E_2","220621_E_3")) %>% for single set only
    dplyr::filter(nmod=="D3Acetyl") -> summaryfile_psm_for_training
  
  
  summaryfile_psm_for_training %>%   
    mutate(index=1:nrow(.)) %>% 
    separate_longer_delim(Modifications, delim="; ") %>% 
    mutate(
      mods_aa = case_when(
        str_detect(Modifications,"^N\\-Term")~"Any N-term",
        str_detect(Modifications,"^[A-Z]")~str_extract(Modifications,"^[A-Z]"),
        .default = Modifications
      )
    ) %>% 
    mutate(
      mods_name = case_when(
        str_detect(Modifications,"Acetyl\\)")~"Acetyl",
        str_detect(Modifications,"Acetyl\\:2H\\(3\\)")~"Acetyl:2H(3)",
        str_detect(Modifications,"Oxidation")~"Oxidation",
        str_detect(Modifications,"Carbamidomethyl")~"Carbamidomethyl",
        .default = Modifications
      )
    ) %>% 
    mutate(
      mod_sites = case_when(
        str_detect(Modifications,"^N\\-Term")~"0",
        str_detect(Modifications,"^[A-Z]")~str_extract(Modifications,"(?<=^[A-Z])[0-9]+"),
        .default = Modifications
      )
    ) %>% 
    nest(nestedcol = c(Modifications,mods_name,mods_aa,mod_sites)) %>% 
    mutate(
      Modifications = map_chr(nestedcol,\(x) x$Modifications %>% paste(collapse = ";"))
    ) %>% 
    mutate(
      mods = map_chr(nestedcol,\(x) unite(data = x, col = mods, c(mods_name,mods_aa),sep ="@") %>% .$mods %>% paste(collapse = ";"))
    ) %>% 
    mutate(
      mod_sites = map_chr(nestedcol,\(x) x$mod_sites %>% paste(collapse = ";"))
    )%>% 
    mutate(
      nAA = str_length(Sequence),
      nce = as.integer(27),
      instrument = "QE",
      charge = as.integer(Charge)
    )-> summaryfile_psm_for_training
}

###MS2 read ##don't need to revise after changing training set
{
  folder_pdspectra <- "./annotation_psm/200821_HELA_MGTG_iNrich_bRP_F_startR/"
  ##read xml
  xml2::read_html(
    paste0(folder_pdspectra,"index.html")
  ) -> file_index
  #from xml object read table
  rvest::html_table(
    file_index,
    header = TRUE
  ) ->file_index
  
  ##table trimming
  file_index[[1]] -> file_index #unlist
  file_index[seq(1,nrow(file_index),by=2),seq(2,ncol(file_index),by=2)] -> file_index ## delete unnessary cols and rows
  colnames(file_index)<- as.vector(file_index[1,]) %>% unlist ## colnames
  file_index[c(-1),] -> file_index
  
  #peak list handling
  readLines(paste0(folder_pdspectra,"index.html"))-> file_index_lines
  
  ##add column of peak list
  file_index %>% 
    bind_cols(tibble(peak_file = str_extract(file_index_lines[str_detect(file_index_lines,"Peak List")],"data\\/.+\\.txt"))) -> file_index
  
  file_index %>% 
    mutate(
      peak_table = map(
        peak_file,\(x) read_tsv(paste0(folder_pdspectra,x),show_col_types = FALSE), .progress = TRUE
      )
    )->file_index2
  file_index2->file_index
  rm(file_index2)
  file_index %>% 
    mutate(
      peak_table_dropna = map(
        peak_table,\(x) drop_na(x)
      )
    ) -> file_index
  #### peak_table_dropna_relative intensity!
  file_index %>% 
    mutate(
      peak_table_dropna_wo_modloss = map(peak_table_dropna,\(x) dplyr::filter(x,!str_detect(matches,"NH3|Ox|H2O")))
    ) %>% 
    mutate(
      peak_table_dropna_wo_modloss = map(peak_table_dropna_wo_modloss, \(x) mutate(x,i_relative=i/max(i)))
    ) -> file_index
}

### from file_index>> new  peak_table_alpha_ms2 >>> fragment_inten_df3 ##don't need to revise after changing training set
{
  file_index %>% 
    mutate(
      peak_table_alpha_ms2= map(peak_table_dropna_wo_modloss,\(x) mutate(x, series = str_extract(matches,"b|y")))
    ) %>% 
    mutate(
      peak_table_alpha_ms2= map(
        peak_table_alpha_ms2,
        \(x) mutate(x, ion_num = as.integer(str_extract(matches,"(?<=\\()[0-9]+(?=\\))")))
      )
    ) %>%
    mutate(
      peak_table_alpha_ms2= map(
        peak_table_alpha_ms2,
        \(x) mutate(x, charge = as.integer(str_extract(matches,"(?<=\\()[0-9]+(?=\\+\\))")))
      )
    ) %>% 
    mutate(
      peak_table_alpha_ms2 = map(
        peak_table_alpha_ms2,
        \(x) dplyr::filter(x, charge<=2)
      )
    ) %>% 
    mutate(
      peak_table_alpha_ms2 = map(peak_table_alpha_ms2, \(x) mutate(x,i_relative=i/max(i)))
    ) %>% 
    mutate(
      strip_seq=str_to_upper(str_extract(Sequence,"(?<=\\.)[:alpha:]+(?=\\.)"))
    ) %>% 
    mutate(
      seq_len = str_length(strip_seq)
    ) %>% 
    mutate(
      peak_table_alpha_ms2 = map2(
        peak_table_alpha_ms2,
        seq_len,
        \(x,y) mutate(x,seq_len = y)
      )
    ) %>% 
    mutate(
      peak_table_alpha_ms2 = map(
        peak_table_alpha_ms2,
        \(x) mutate(x,ion_num2 = if_else(series == "y",seq_len - ion_num,ion_num))
      )
    )%>%
    mutate(
      peak_table_alpha_ms2 = map(
        peak_table_alpha_ms2,
        \(x) mutate(x,alpha_col = str_c(series,"_","z",charge))
      )
    ) %>% 
    mutate(
      fragment_inten_df = map(
        peak_table_alpha_ms2,
        \(x) pivot_wider(x,id_cols = ion_num2,names_from = alpha_col,values_from = i_relative)
      )
    ) %>% 
    mutate(
      fragment_inten_df2 = map(seq_len,\(x) tibble(ion_num2=seq(1,x-1)))
    ) %>% 
    mutate(
      fragment_inten_df3 = map2(
        fragment_inten_df,
        fragment_inten_df2,
        \(x,y) full_join(x,y,by=c("ion_num2"="ion_num2"))
      )
    )-> file_index3
  
  colnames(file_index3)
  colnames(file_index3)[16] <-"scan_num"
  colnames(file_index3)
  file_index3 %>% 
    mutate(
      scan_num = as.integer(scan_num)
    ) %>% 
    mutate(
      charge = as.integer(Charge)
    ) %>% 
    mutate(
      XCorr = as.numeric(XCorr)
    ) %>% 
    mutate(
      `Percolator SVMScore` = as.numeric(`Percolator SVMScore`)
    ) ->file_index3
  
  file_index3 %>% 
    mutate(
      Modifications = str_replace_all(Modifications,"[:blank:]","")
    ) -> file_index3
  
  file_index3 %>% 
    distinct(
      strip_seq,scan_num,charge,XCorr,
      .keep_all = TRUE
    ) -> file_index3
  
  file_index3 %>% 
    unnest(fragment_inten_df3) %>% 
    dplyr::rename(percol_svm = `Percolator SVMScore`) %>% 
    select(strip_seq,Modifications,scan_num,XCorr,percol_svm,charge,ion_num2,b_z1,b_z2,y_z1,y_z2) %>%
    mutate(
      across(
        .cols = c(b_z1,b_z2,y_z1,y_z2),
        .fns = \(x) if_else(is.na(x),0,x)
      )
    ) -> file_index4
  ##now strip_seq + scan_num + XCorr + charge
}

##file_index4 + summaryfile_psm_for_training ###230529 Success!
{
  # summaryfile_psm_for_training %>% 
  #   select(Sequence,nAA) %>% print(n=Inf)
  summaryfile_psm_for_training %>% 
    select(-c(nestedcol)) %>% 
    mutate(
      scan_num = as.integer(`First Scan`)
    ) %>% distinct(
      Sequence,scan_num,charge,XCorr,
      .keep_all = TRUE
    ) %>% 
    inner_join(
      file_index4,
      by=join_by(
        Sequence == strip_seq,
        scan_num==scan_num,
        charge==charge,
        XCorr==XCorr
      )
    ) ->summaryfile_psm_for_training2
  # summaryfile_psm_for_training2 %>% 
  #   select(Sequence,nAA) %>% dplyr::filter(nAA==7) %>% print(n=Inf)
  summaryfile_psm_for_training2 %>% 
    arrange(
      index,ion_num2
    ) %>% 
    mutate(
      frag_index = seq(0,nrow(summaryfile_psm_for_training2)-1)
    ) %>%
    nest(frag_df= c(frag_index,ion_num2,b_z1,b_z2,y_z1,y_z2)) -> summaryfile_psm_for_training3 
  
  
  # ###summaryfile_psm_for_training3 is the final table ##dont change it!
  summaryfile_psm_for_training3 %>% 
    mutate(
      index= seq(1,nrow(summaryfile_psm_for_training3))
    ) %>% 
    mutate(
      frag_start_idx = map_dbl(frag_df,\(x) min(x$frag_index))
    ) %>% 
    mutate(
      frag_stop_idx = map_dbl(frag_df,\(x) max(x$frag_index)+1)
    ) %>% 
    dplyr::rename(
      spec_idx = index,
      irt = `RT in min` 
    ) %>% 
    select(
      Sequence,charge,spec_idx,irt,mods,mod_sites,nAA,frag_start_idx,frag_stop_idx,nce,instrument
    ) %>% 
    dplyr::rename(
      sequence = Sequence
    ) -> summaryfile_psm_for_training4
  summaryfile_psm_for_training4 %>% 
    write_tsv(
      "C:/Users/admin/alphapeptdeep/project/folder/alphapeptdeep/precursor_df.tsv"
    )
  ### create frag_intensity_df
  summaryfile_psm_for_training3 %>% 
    select(frag_df) %>% 
    unnest(frag_df) %>% 
    select(
      b_z1,b_z2,y_z1,y_z2
    ) %>% 
    write_tsv(
      "C:/Users/admin/alphapeptdeep/project/folder/alphapeptdeep/frag_intensity_df.tsv"
    )
}

#making a test set! ::summaryfile_psm_for_testing
{
  summaryfile_psm %>% 
    dplyr::filter(
      arginylation == TRUE
    ) -> summaryfile_psm_for_testing
  
  
  summaryfile_psm_for_testing %>%   
    mutate(index=1:nrow(.)) %>% 
    separate_longer_delim(Modifications, delim="; ") %>% 
    mutate(
      mods_aa = case_when(
        str_detect(Modifications,"^N\\-Term")~"Any N-term",
        str_detect(Modifications,"^[A-Z]")~str_extract(Modifications,"^[A-Z]"),
        .default = Modifications
      )
    ) %>% 
    mutate(
      mods_name = case_when(
        str_detect(Modifications,"Acetyl\\)")~"Acetyl",
        str_detect(Modifications,"Acetyl\\:2H\\(3\\)")~"Acetyl:2H(3)",
        str_detect(Modifications,"Oxidation")~"Oxidation",
        str_detect(Modifications,"Carbamidomethyl")~"Carbamidomethyl",
        .default = Modifications
      )
    ) %>% 
    mutate(
      mod_sites = case_when(
        str_detect(Modifications,"^N\\-Term")~"0",
        str_detect(Modifications,"^[A-Z]")~str_extract(Modifications,"(?<=^[A-Z])[0-9]+"),
        .default = Modifications
      )
    ) %>% 
    nest(nestedcol = c(Modifications,mods_name,mods_aa,mod_sites)) %>% 
    mutate(
      Modifications = map_chr(nestedcol,\(x) x$Modifications %>% paste(collapse = ";"))
    ) %>% 
    mutate(
      mods = map_chr(nestedcol,\(x) unite(data = x, col = mods, c(mods_name,mods_aa),sep ="@") %>% .$mods %>% paste(collapse = ";"))
    ) %>% 
    mutate(
      mod_sites = map_chr(nestedcol,\(x) x$mod_sites %>% paste(collapse = ";"))
    )%>% 
    mutate(
      nAA = str_length(Sequence),
      nce = as.integer(27),
      instrument = "QE",
      charge = as.integer(Charge)
    )-> summaryfile_psm_for_testing
  summaryfile_psm_for_testing %>% 
    separate_rows(mod_sites,sep=";") %>% 
    mutate(
      mod_sites = as.integer(mod_sites)
    ) %>%
    mutate(
      mod_sites = case_when(
        mod_sites > 0 ~ mod_sites+1,
        .default = mod_sites
      )
    ) %>%
    mutate(
      nAA = nAA+1
    ) %>%
    nest(
      mod_sites=mod_sites
    ) %>% 
    mutate(
      mod_sites = map_chr(mod_sites,\(x) x$mod_sites %>%  paste(collapse = ";"))
    ) %>% 
    mutate(
      sequence = case_when(
        str_detect(Sequence,"^N") ~ str_replace(Sequence,"^N","RD"),
        str_detect(Sequence,"^Q") ~ str_replace(Sequence,"^Q","RE"),
        str_detect(Sequence,"^E") ~ str_replace(Sequence,"^E","RE"),
        str_detect(Sequence,"^D") ~ str_replace(Sequence,"^D","RD")
      )
    ) %>% 
    mutate(
      irt = `RT in min`
    ) %>% 
    mutate(
      mods = case_when(
        str_detect(mods, "N\\-Term\\(ArgAcetylD3NQ\\)") ~ str_replace(mods,"N\\-Term\\(ArgAcetylD3NQ\\)","Acetyl:2H(3)"),
        str_detect(mods, "N\\-Term\\(ArgAcetylD3\\)") ~ str_replace(mods,"N\\-Term\\(ArgAcetylD3\\)","Acetyl:2H(3)")
      )
    ) -> summaryfile_psm_for_testing
  summaryfile_psm_for_testing %>% 
    dplyr::select(
      sequence,charge,index,mods,mod_sites,nAA,irt
    ) %>% 
    write_tsv(
      "C:/Users/admin/alphapeptdeep/project/folder/alphapeptdeep/input_testing.tsv"
    ) ###this is the testing data
}

###MS2 read for test set! ##>>>file_index_arg4 ##diagnostic_ion_mz + ion_int_relative
###@@@folder input!!!
###b and y ttest!
{
  folder_pdspectra_arg <- "./annotation_psm/200821_HELA_MGTG_iNrich_bRP_F_arg/"
  ##read xml
  xml2::read_html(
    paste0(folder_pdspectra_arg,"index.html")
  ) -> file_index_arg
  #from xml object read table
  rvest::html_table(
    file_index_arg,
    header = TRUE
  ) ->file_index_arg
  
  ##table trimming
  file_index_arg[[1]] -> file_index_arg #unlist
  file_index_arg[seq(1,nrow(file_index_arg),by=2),seq(2,ncol(file_index_arg),by=2)] -> file_index_arg ## delete unnessary cols and rows
  colnames(file_index_arg)<- as.vector(file_index_arg[1,]) %>% unlist ## colnames
  file_index_arg[c(-1),] -> file_index_arg
  
  #peak list handling
  readLines(paste0(folder_pdspectra_arg,"index.html"))-> file_index_lines_arg
  
  ##add column of peak list
  file_index_arg %>% 
    bind_cols(
      tibble(peak_file = str_extract(file_index_lines_arg[str_detect(file_index_lines_arg,"Peak List")],"data\\/.+\\.txt"))
    ) %>% 
    mutate(
      peak_table = map(
        peak_file,\(x) read_tsv(paste0(folder_pdspectra_arg,x),show_col_types = FALSE), .progress = TRUE
      )
    )-> file_index_arg
  
  ###diagnostic ion b (0) (1+)
  # file_index_arg2 %>% 
  #   mutate(
  #     peak_table = map(
  #       peak_table, \(x) mutate(x,matches = if_else(`m/z` < 202.138331+0.005 & `m/z` > 202.138331-0.005, "b (0) (1+)", matches))
  #     )
  #   ) -> file_index_arg3
  
  file_index_arg %>% 
    mutate(
      diagnostic_ion_mz = map(peak_table, \(x) x$`m/z`[str_which(x$matches, "b[:blank:]\\(0")] )
    ) %>% 
    mutate(
      diagnostic_ion_mz = map_dbl(diagnostic_ion_mz,\(x) ifelse(length(x)>0,x,NA))
    ) %>%  
    mutate(
      peak_table_dropna = map(
        peak_table,\(x) drop_na(x)
      )
    ) %>%
    mutate(
      peak_table_dropna_wo_modloss = map(peak_table_dropna,\(x) dplyr::filter(x,!str_detect(matches,"NH3|Ox|H2O")))
    ) %>% 
    mutate(
      peak_table_dropna_wo_modloss = map(peak_table_dropna_wo_modloss, \(x) mutate(x,i_relative=i/max(i)))
    ) %>% 
    mutate(
      peak_table_alpha_ms2= map(peak_table_dropna_wo_modloss,\(x) mutate(x, series = str_extract(matches,"b|y")))
    ) %>% 
    mutate(
      peak_table_alpha_ms2= map(
        peak_table_alpha_ms2,
        \(x) mutate(x, ion_num = as.integer(str_extract(matches,"(?<=\\()[0-9]+(?=\\))")))
      )
    ) %>% ###b ion +1 because of modification! 
    mutate(
      peak_table_alpha_ms2= map(
        peak_table_alpha_ms2,
        \(x) mutate(x, ion_num = if_else(series == "b",ion_num +1,ion_num))
      )
    ) %>% 
    mutate(
      peak_table_alpha_ms2= map(
        peak_table_alpha_ms2,
        \(x) mutate(x, charge = as.integer(str_extract(matches,"(?<=\\()[0-9]+(?=\\+\\))")))
      )
    ) %>% 
    mutate(
      peak_table_alpha_ms2 = map(
        peak_table_alpha_ms2,
        \(x) dplyr::filter(x, charge<=2)
      )
    ) %>% 
    mutate(
      peak_table_alpha_ms2 = map(peak_table_alpha_ms2, \(x) mutate(x,i_relative=i/max(i)))
    ) %>% 
    mutate(
      strip_seq=str_to_upper(str_extract(Sequence,"(?<=\\.)[:alpha:]+(?=\\.)"))
    ) %>% ##strip seq add r
    mutate(
      strip_seq = paste0("R",strip_seq)
    ) %>% 
    mutate(
      seq_len = str_length(strip_seq) ###modification of R! +1!
    ) %>% 
    mutate(
      peak_table_alpha_ms2 = map2(
        peak_table_alpha_ms2,
        seq_len,
        \(x,y) mutate(x,seq_len = y)
      )
    ) %>% 
    mutate(
      peak_table_alpha_ms2 = map(
        peak_table_alpha_ms2,
        \(x) mutate(x,ion_num2 = if_else(series == "y",seq_len - ion_num,ion_num))
      )
    )%>%
    mutate(
      peak_table_alpha_ms2 = map(
        peak_table_alpha_ms2,
        \(x) mutate(x,alpha_col = str_c(series,"_","z",charge))
      )
    ) %>% 
    mutate(
      fragment_inten_df = map(
        peak_table_alpha_ms2,
        \(x) pivot_wider(x,id_cols = ion_num2,names_from = alpha_col,values_from = i_relative)
      )
    ) %>% 
    mutate(
      fragment_inten_df2 = map(seq_len,\(x) tibble(ion_num2=seq(1,x-1)))
    ) %>% 
    mutate(
      fragment_inten_df3 = map2(
        fragment_inten_df,
        fragment_inten_df2,
        \(x,y) full_join(x,y,by=c("ion_num2"="ion_num2"))
      )
    )-> file_index_arg2
  
  
  colnames(file_index_arg2)[16] <-"scan_num"
  
  file_index_arg2 %>% 
    mutate(
      scan_num = as.integer(scan_num)
    ) %>% 
    mutate(
      charge = as.integer(Charge)
    ) %>% 
    mutate(
      XCorr = as.numeric(XCorr)
    ) %>% 
    mutate(
      `Percolator SVMScore` = as.numeric(`Percolator SVMScore`)
    ) %>% 
    mutate(
      Modifications = str_replace_all(Modifications,"[:blank:]","")
    ) -> file_index_arg2
  
  # file_index_arg2 %>% 
  #   distinct(
  #     strip_seq,scan_num,charge,XCorr,
  #     .keep_all = TRUE
  #   ) -> file_index_arg3
  # 
  # file_index_arg3 %>% 
  #   mutate(
  #     diagnostic_ion_int_relative = map(peak_table_alpha_ms2, \(x) x$i_relative[str_which(x$matches, "b[:blank:]\\(0")] )
  #   ) %>% 
  #   mutate(
  #     diagnostic_ion_int_relative = map_dbl(diagnostic_ion_int_relative,\(x) ifelse(length(x)>0,x,NA))
  #   ) -> file_index_arg3
  # 
  #FragmentPeptide()
  ### bion yion ttest
  {
    MonoisotopicMass(
      formula = list(C=8, H=11, N=4, O=2, S=0, x= 3),
      isotopes = list(x = 2.01410177812) ##deuterium
    ) -> mass_ntd3arg
    
    MonoisotopicMass(
      formula = list(C=8, H=11, N=2, O=2, S=0, x= 3),
      isotopes = list(x = 2.01410177812) ##deuterium
    ) -> mass_d3k
    
    MonoisotopicMass(
      formula = list(C=5, H=9, N=1, O=2, S=1)
    ) -> mass_met_ox
    
    file_index_arg2 %>% 
      mutate(
        sequence_for_fragspec = str_extract(Sequence,"(?<=\\.).+(?=\\.)")
      ) %>%
      mutate(
        sequence_for_fragspec = str_replace(sequence_for_fragspec,"^n","D"),
        sequence_for_fragspec = str_replace(sequence_for_fragspec,"^q","E"),
        sequence_for_fragspec = str_replace(sequence_for_fragspec,"^d","D"),
        sequence_for_fragspec = str_replace(sequence_for_fragspec,"^e","E"),
        sequence_for_fragspec = str_replace_all(sequence_for_fragspec,"c","C"),
      ) %>%
      mutate(
        sequence_for_fragspec = paste0("r",sequence_for_fragspec)
      ) %>% 
      mutate(
        frag_theo_spec = map(
          sequence_for_fragspec,
          \(x) OrgMassSpecR::FragmentPeptide(
            x,
            fragments = "by",
            IAA = TRUE,
            custom = list(
              code = c("r","k","m"),
              mass = c(mass_ntd3arg,mass_d3k,mass_met_ox)
            )
          )
        )
      ) -> file_index_arg3
  }
  
  file_index_arg3 %>% 
    mutate(
      frag_theo_spec = map(
        frag_theo_spec,
        \(x) mutate(
          x,
          series = str_extract(ms2type,"[by]"),
          ion_num = as.integer(str_extract(ms2type,"(?<=[by])[0-9]+")),
          charge = as.integer(str_extract(ms2type,"[0-9](?=\\+)"))
        )
      )
    ) -> file_index_arg4
  
  file_index_arg4 %>% 
    mutate(
      peak_table_alpha_ms2_vs_theo_spec = map2(
        peak_table_alpha_ms2,
        frag_theo_spec,
        \(x,y) left_join(x,y, join_by(series == series, ion_num == ion_num, charge == charge))
      )
    ) -> file_index_arg4
  
  file_index_arg4 %>% 
    mutate(
      peak_table_alpha_ms2_vs_theo_spec = map(
        peak_table_alpha_ms2_vs_theo_spec,
        \(x) mutate(x,
                    error = `m/z`-ms2mz
        )
      )
    ) -> file_index_arg4
  
  file_index_arg4 %>% 
    mutate(
      b_ion_count = map_int(
        peak_table_alpha_ms2_vs_theo_spec,
        \(x) filter(x, series == "b") %>% nrow()
      ),
      y_ion_count = map_int(
        peak_table_alpha_ms2_vs_theo_spec,
        \(x) filter(x, series == "y") %>% nrow()
      ),
      b_ion_errors = map(
        peak_table_alpha_ms2_vs_theo_spec,
        \(x) filter(x, series == "b") %>% pull(error)
      ),
      y_ion_errors = map(
        peak_table_alpha_ms2_vs_theo_spec,
        \(x) filter(x, series == "y") %>% pull(error)
      )
    )-> file_index_arg4
  
  file_index_arg4 %>%
    mutate(
      peak_table_alpha_ms2_vs_theo_spec = map(
        peak_table_alpha_ms2_vs_theo_spec,
        \(x) mutate(x,
                    error_ppm = (`m/z`-ms2mz)/`m/z`*1e6
        )
      )
    ) -> file_index_arg4
  
  file_index_arg4 %>% 
    mutate(
      b_ion_errors_ppm = map(
        peak_table_alpha_ms2_vs_theo_spec,
        \(x) filter(x, series == "b") %>% pull(error_ppm)
      ),
      y_ion_errors_ppm = map(
        peak_table_alpha_ms2_vs_theo_spec,
        \(x) filter(x, series == "y") %>% pull(error_ppm)
      ),
      b_ion_mz = map(
        peak_table_alpha_ms2_vs_theo_spec,
        \(x) filter(x, series == "b") %>% pull(`m/z`)
      ),
      y_ion_mz = map(
        peak_table_alpha_ms2_vs_theo_spec,
        \(x) filter(x, series == "y") %>% pull(`m/z`)
      )
    ) -> file_index_arg4
  
  
  file_index_arg4 %>% 
    mutate(
      b_y_ttest = map2_dbl(
        b_ion_errors,
        y_ion_errors,
        \(x,y) if(length(x)>1 & length(y)>1){t.test(unlist(x),unlist(y),var.equal = TRUE)$p.value} else{return(NA)}
      ),
      b_y_ttest_ppm = map2_dbl(
        b_ion_errors_ppm,
        y_ion_errors_ppm,
        \(x,y) if(length(x)>1 & length(y)>1){t.test(unlist(x),unlist(y),var.equal = TRUE)$p.value} else{return(NA)}
      ),
      b_y_ttest_mz = map2_dbl(
        b_ion_mz,
        y_ion_mz,
        \(x,y) if(length(x)>1 & length(y)>1){t.test(unlist(x),unlist(y),var.equal = TRUE)$p.value} else{return(NA)}
      ),
      b_y_ttest_mean_diff = map2_dbl(
        b_ion_mz,
        y_ion_mz,
        \(x,y) if(length(x)>1 & length(y)>1){t.test(unlist(x),unlist(y),var.equal = TRUE) %>% tidy %>% .$estimate} else{return(NA)}
      )
    ) -> file_index_arg5
  
  t.test(c(2,2,2),c(4,4,5))$estimate["mean of x"]
  
  diff(4,5)
  
  file_index_arg5 %>% 
    select(
      b_y_ttest_mean_diff,
      b_y_ttest,
      b_y_ttest_ppm
    ) %>% 
    mutate(
      b_y_ttest_sig = if_else(b_y_ttest<0.05,"smert failed","smert succeed"),
      b_y_ttest_ppm_sig = if_else(b_y_ttest_ppm<0.05,1,0)
    ) %>% 
    drop_na() %>% 
    ggplot(
      aes(x=b_y_ttest_mean_diff,fill=b_y_ttest_sig)
    )+
    geom_density(
      alpha = 0.1
    )-> test1
  
  file_index_arg5 %>% 
    select(
      b_y_ttest_mean_diff,
      b_y_ttest,
      b_y_ttest_ppm
    ) %>% 
    mutate(
      b_y_ttest_sig = if_else(b_y_ttest<0.05,"smert failed","smert succeed"),
      b_y_ttest_ppm_sig = if_else(b_y_ttest_ppm<0.05,"smert failed","smert succeed")
    ) %>% 
    drop_na() %>% 
    ggplot(
      aes(x=b_y_ttest_mean_diff,fill=b_y_ttest_ppm_sig)
    )+
    # geom_histogram(
    #   binwidth = 20,
    #   position = "dodge"
    # )+
    geom_density(
      alpha = 0.1
    ) -> test2
  
  ggarrange(test1,test2, ncol = 1)  
  # lims(y=c(0,70))
  
  file_index_arg5 %>% 
    mutate(
      b_y_ttest_sig = if_else(b_y_ttest<0.05,1,0),
      b_y_ttest_ppm_sig = if_else(b_y_ttest_ppm<0.05,1,0),
      b_y_ttest_mz = if_else(b_y_ttest_mz<0.05,1,0)
    ) %>% 
    summarize(
      .by = c(b_y_ttest_ppm_sig,b_y_ttest_mz),
      count = n()
    ) %>% 
    drop_na()
    
  
  file_index_arg5 %>% 
    mutate(
      b_y_ttest_sig = if_else(b_y_ttest<0.05,"blue",NA),
      b_y_ttest_ppm_sig = if_else(b_y_ttest_ppm<0.05,"red",NA)
    ) %>% 
    summarize(
      .by = c(b_y_ttest_sig,b_y_ttest_ppm_sig),
      count = n()
    ) %>% 
    write_tsv(
      "b_y_ttest.tsv"
    )
  file_index_arg5 %>% 
    mutate(
      b_y_ttest_sig = if_else(b_y_ttest<0.05,"blue",NA),
      b_y_ttest_ppm_sig = if_else(b_y_ttest_ppm<0.05,"red",NA)
    ) %>% 
    ggplot(
      aes(x= b_y_ttest, y= b_y_ttest_ppm)
    )+
    geom_point(
      aes(color = b_y_ttest_sig, fill = b_y_ttest_ppm_sig),
      shape=21,
      size = 1,
      stroke = 1
    )+
    scale_fill_manual(values=c("blue", "black")) + 
    scale_color_manual(values=c("red", "black"))+
    theme_bw()
  ggsave(
    paste0("./export/smert_ppm_mda.png"),
    width = 15, height = 12, dpi=600, units = "cm")
  
  file_index_arg5 %>% 
    mutate(
      data = map(
        peak_table_alpha_ms2_vs_theo_spec,
        \(x) x %>% select(`m/z`,error,error_ppm,series)         
      )
    ) %>% 
    select(data) %>% 
    unnest(data) %>% 
    pivot_longer(
      c(error,error_ppm)
    ) %>% 
    ggplot(
      aes(x=`m/z`,y=value)
    )+
    geom_point(
      alpha = 0.5
    )+
    facet_wrap(
      facets = vars(name,series),
      scales = "free"
    )+
    theme_bw()
  
  
  file_index_arg5 %>%   
    mutate(index=1:nrow(.)) %>% 
    separate_longer_delim(Modifications, delim=";") %>% 
    mutate(
      mods_aa = case_when(
        str_detect(Modifications,"^N\\-Term")~"Any N-term",
        str_detect(Modifications,"^[A-Z]")~str_extract(Modifications,"^[A-Z]"),
        .default = Modifications
      )
    ) %>% 
    mutate(
      mods_name = case_when(
        str_detect(Modifications,"Acetyl\\)")~"Acetyl",
        str_detect(Modifications,"Acetyl\\:2H\\(3\\)")~"Acetyl:2H(3)",
        str_detect(Modifications,"Oxidation")~"Oxidation",
        str_detect(Modifications,"Carbamidomethyl")~"Carbamidomethyl",
        .default = Modifications
      )
    ) %>% 
    mutate(
      mod_sites = case_when(
        str_detect(Modifications,"^N\\-Term")~"0",
        str_detect(Modifications,"^[A-Z]")~as.character(as.integer(str_extract(Modifications,"(?<=^[A-Z])[0-9]+"))+1),
        .default = Modifications
      )
    ) %>% 
    nest(nestedcol = c(Modifications,mods_name,mods_aa,mod_sites)) %>% 
    mutate(
      Modifications = map_chr(nestedcol,\(x) x$Modifications %>% paste(collapse = ";"))
    ) %>% 
    mutate(
      mods = map_chr(nestedcol,\(x) unite(data = x, col = mods, c(mods_name,mods_aa),sep ="@") %>% .$mods %>% paste(collapse = ";"))
    ) %>% 
    mutate(
      mod_sites = map_chr(nestedcol,\(x) x$mod_sites %>% paste(collapse = ";"))
    )%>% 
    mutate(
      mods = case_when(
        str_detect(mods, "N\\-Term\\(ArgAcetylD3NQ\\)") ~ str_replace(mods,"N\\-Term\\(ArgAcetylD3NQ\\)","Acetyl:2H(3)"),
        str_detect(mods, "N\\-Term\\(ArgAcetylD3\\)") ~ str_replace(mods,"N\\-Term\\(ArgAcetylD3\\)","Acetyl:2H(3)")
      )
    ) %>% 
    mutate(
      sequence = str_to_upper(sequence_for_fragspec),
      nAA = seq_len,
      nce = as.integer(27),
      instrument = "QE",
      charge = as.integer(Charge)
    )-> file_index_arg6_ttest
  
  ## need to write down more for exporting tsv from file_index_arg6
}
###file_index_arg6_ttest

####Training setfile: psm_563516.rds >>> psm_311547_trypsin , psm_251969_chymo
####NEW training! 231130!
##psm_311547_trypsin
##psm_251969_chymo
{
  readRDS("./export/psm_563516.rds") -> psm_563516
  
  psm_563516 %>% 
    filter(
      enzyme == "trypsin"
    ) -> psm_311547_trypsin
  
  psm_563516 %>% 
    filter(
      enzyme == "chymotrypsin"
    ) -> psm_251969_chymo
  
  
}
##@psm_311547_trypsin
##@psm_251969_chymo 

###psm_311547_trypsin training
{  
  readRDS("./export/psm_311547_trypsin.rds") -> psm_311547_trypsin
  
  psm_311547_trypsin %>% 
    unnest(fragment_inten_df3) %>% 
    arrange(
      index,ion_num2
    ) %>% 
    mutate(
      frag_index = seq(0,nrow(.)-1)
    ) %>% 
    nest(frag_df= c(frag_index,ion_num2,b_z1,b_z2,y_z1,y_z2)) %>% 
    mutate(
      frag_start_idx = map_dbl(frag_df,\(x) min(x$frag_index))
    ) %>% 
    mutate(
      frag_stop_idx = map_dbl(frag_df,\(x) max(x$frag_index)+1)
    ) %>% 
    dplyr::rename(
      spec_idx = index,
      irt = `RT in min`,
      sequence = strip_seq
    ) -> psm_311547_trypsin
  
  psm_311547_trypsin %>% 
    select(
      sequence,charge,spec_idx,mods,mod_sites,nAA,frag_start_idx,frag_stop_idx,nce,instrument
    ) %>% 
    write_tsv(
      "C:/Users/admin/alphapeptdeep/project/folder/alphapeptdeep/precursor_df_311547.tsv"
    )
  psm_311547_trypsin %>% 
    select(frag_df) %>% 
    unnest(frag_df) %>% 
    select(
      b_z1,b_z2,y_z1,y_z2
    ) %>% 
    mutate(
      across(everything(),\(x) replace_na(x,0))
    ) %>% 
    write_tsv(
      "C:/Users/admin/alphapeptdeep/project/folder/alphapeptdeep/precursor_df_311547_fragdf.tsv"
    )
  
 
}
####
###psm_251969_chymotrypsin training
{  
  readRDS("./export/psm_251969_chymo.rds") -> psm_251969_chymo
  
  psm_251969_chymo %>% 
    unnest(fragment_inten_df3) %>% 
    arrange(
      index,ion_num2
    ) %>% 
    mutate(
      frag_index = seq(0,nrow(.)-1)
    ) %>% 
    nest(frag_df= c(frag_index,ion_num2,b_z1,b_z2,y_z1,y_z2)) %>% 
    mutate(
      frag_start_idx = map_dbl(frag_df,\(x) min(x$frag_index))
    ) %>% 
    mutate(
      frag_stop_idx = map_dbl(frag_df,\(x) max(x$frag_index)+1)
    ) %>% 
    dplyr::rename(
      spec_idx = index,
      irt = `RT in min`,
      sequence = strip_seq
    ) -> psm_251969_chymo
  
  psm_251969_chymo %>% 
    select(
      sequence,charge,spec_idx,mods,mod_sites,nAA,frag_start_idx,frag_stop_idx,nce,instrument
    ) %>% 
    write_tsv(
      "C:/Users/admin/alphapeptdeep/project/folder/alphapeptdeep/psm_251969_chymo.tsv"
    )
  psm_251969_chymo %>% 
    select(frag_df) %>% 
    unnest(frag_df) %>% 
    select(
      b_z1,b_z2,y_z1,y_z2
    ) %>% 
    mutate(
      across(everything(),\(x) replace_na(x,0))
    ) %>% 
    write_tsv(
      "C:/Users/admin/alphapeptdeep/project/folder/alphapeptdeep/psm_251969_chymo_fragdf.tsv"
    )
}

##splitting training set:: trypsin
{
  summaryfile_psm_for_testing3 %>%
    select(-enzyme.y) %>% 
    rename(enzyme = enzyme.x) -> summaryfile_psm_for_testing3
    
  summaryfile_psm_for_testing3 %>% 
    filter(enzyme == "trypsin") %>%
    unnest(frag_df) %>% 
    arrange(
      index,ion_num2
    ) %>% 
    mutate(
      frag_index = seq(0,nrow(.)-1)
    ) %>% 
    nest(frag_df= c(frag_index,ion_num2,b_z1,b_z2,y_z1,y_z2)) %>% 
    mutate(
      frag_start_idx = map_dbl(frag_df,\(x) min(x$frag_index))
    ) %>% 
    mutate(
      frag_stop_idx = map_dbl(frag_df,\(x) max(x$frag_index)+1)
    ) %>% 
    dplyr::rename(
      spec_idx = index
    ) %>% 
    select(
      sequence,charge,spec_idx,mods,mod_sites,nAA,frag_start_idx,frag_stop_idx,nce,instrument
    ) %>% 
    write_tsv(
      "C:/Users/admin/alphapeptdeep/project/folder/alphapeptdeep/precursor_df_arg_1202_trypsin.tsv"
    )
  
  summaryfile_psm_for_testing3 %>% 
    filter(enzyme == "trypsin") %>%
    unnest(frag_df) %>% 
    arrange(
      index,ion_num2
    ) %>% 
    mutate(
      frag_index = seq(0,nrow(.)-1)
    ) %>% 
    nest(frag_df= c(frag_index,ion_num2,b_z1,b_z2,y_z1,y_z2)) %>% 
    select(frag_df) %>% 
    unnest(frag_df) %>% 
    select(
      b_z1,b_z2,y_z1,y_z2
    ) %>% 
    mutate(
      across(everything(),\(x) replace_na(x,0))
    ) %>% 
    write_tsv(
      "C:/Users/admin/alphapeptdeep/project/folder/alphapeptdeep/precursor_df_arg_1202_trypsin_fragdf.tsv"
    )
}
##splitting training set:: chymotrypsin
{
  summaryfile_psm_for_testing3 %>% 
    filter(enzyme == "chymotrypsin") %>%
    unnest(frag_df) %>% 
    arrange(
      index,ion_num2
    ) %>% 
    mutate(
      frag_index = seq(0,nrow(.)-1)
    ) %>% 
    nest(frag_df= c(frag_index,ion_num2,b_z1,b_z2,y_z1,y_z2)) %>% 
    mutate(
      frag_start_idx = map_dbl(frag_df,\(x) min(x$frag_index))
    ) %>% 
    mutate(
      frag_stop_idx = map_dbl(frag_df,\(x) max(x$frag_index)+1)
    ) %>% 
    dplyr::rename(
      spec_idx = index
    ) %>% 
    select(
      sequence,charge,spec_idx,mods,mod_sites,nAA,frag_start_idx,frag_stop_idx,nce,instrument
    ) %>% 
    write_tsv(
      "C:/Users/admin/alphapeptdeep/project/folder/alphapeptdeep/precursor_df_arg_15_chymotrypsin.tsv"
    )
  
  summaryfile_psm_for_testing3 %>% 
    filter(enzyme == "chymotrypsin") %>%
    unnest(frag_df) %>% 
    arrange(
      index,ion_num2
    ) %>% 
    mutate(
      frag_index = seq(0,nrow(.)-1)
    ) %>% 
    nest(frag_df= c(frag_index,ion_num2,b_z1,b_z2,y_z1,y_z2)) %>% 
    select(frag_df) %>% 
    unnest(frag_df) %>% 
    select(
      b_z1,b_z2,y_z1,y_z2
    ) %>% 
    mutate(
      across(everything(),\(x) replace_na(x,0))
    ) %>% 
    write_tsv(
      "C:/Users/admin/alphapeptdeep/project/folder/alphapeptdeep/precursor_df_arg_15_chymotrypsin_fragdf.tsv"
    )
}


####
#### do alphapeptdeep! MS2
####


###output series of 230610!

###MS2 prediction data integration ::
##build data
### summary_prediction_ms2_build, output_prediction_ms2_build
##transfer learning data!
#### summaryfile_psm_for_training3 >> summaryfile_psm_for_training_pcc
#### summaryfile_psm_for_testing3 >> summaryfile_psm_for_testing_pcc
{
  folder_alphapeptdeep = "C:/Users/admin/alphapeptdeep/project/folder/alphapeptdeep/"
  
  
  
  output_prediction_ms2_build <- read_tsv(
    paste0(folder_alphapeptdeep,"output_metrics_describ_230610.tsv")
  )[,-1]
  
  output_prediction_ms2_transfer <- read_tsv(
    paste0(folder_alphapeptdeep,"inrich/output_metrics_describ_for_311547_trypsin_231130.tsv")
  )[,-1] %>% 
    bind_rows(
      read_tsv(
        paste0(folder_alphapeptdeep,"inrich/output_metrics_describ_for_251969_chymo_231130.tsv")
      )[,-1]
    )
  
  output_prediction_ms2_pretrained <- read_tsv(
    paste0(folder_alphapeptdeep,"output_metrics_describ4_230610.tsv")
  )[,-1]
  
  output_prediction_ms2_arg <- read_tsv(
    paste0(folder_alphapeptdeep,"inrich/output/output_metrics_describ_for_arg_trypsin_231201.tsv")
  )[,-1] %>% 
    bind_rows(
      read_tsv(
        paste0(folder_alphapeptdeep,"inrich/output/output_metrics_describ_for_arg_chymotrypsin_231201.tsv")
      )[,-1]
    )
  
  summaryfile_psm_for_training3 %>%
    left_join(
      output_prediction_ms2_transfer,
      join_by(
        index==spec_idx
      )
    ) -> summaryfile_psm_for_training_pcc
  
  summaryfile_psm_for_testing3 %>% 
    left_join(
      output_prediction_ms2_arg %>% select(-c(sequence,charge,irt,mods,mod_sites,nAA,nce,instrument,diagnostic_ion_int_relative)),
      join_by(
        index==spec_idx
      )
    ) -> summaryfile_psm_for_testing3_pcc
}
#summaryfile_psm_for_testing3_pcc #231201!

#231027 fragdf, preddf, ms2_prediction result importing
{
  {
    read_tsv(
      paste0(folder_alphapeptdeep,"inrich/output/output_metrics_describ_for_arg_trypsin_231201.tsv")
    ) %>% 
      select(-...1) %>% 
      mutate(
        false_set = 1
      ) -> fdr_pcc_arg_set
    
    fdr_pcc_arg_set %>% 
      bind_rows(
        read_tsv(
          paste0(folder_alphapeptdeep,"inrich/output/output_metrics_describ_for_arg_chymotrypsin_231201.tsv")
        ) %>% 
          select(-...1) %>% 
          mutate(
            false_set = 1
          )
      ) -> fdr_pcc_arg_set
    
    ###Original SPEC
    read_tsv(
      paste0(folder_alphapeptdeep,"precursor_df_arg_1202_trypsin_fragdf.tsv")
    ) %>% 
      bind_rows(
        read_tsv(
          paste0(folder_alphapeptdeep,"precursor_df_arg_15_chymotrypsin_fragdf.tsv")
        )
      )-> fdr_pcc_arg_set_fragdf
    
    read_tsv(
      paste0(folder_alphapeptdeep,"inrich/output/output_arg_ms2_prediction_for_arg_trypsin_231201.tsv")
    ) %>%
      bind_rows(
        read_tsv(
          paste0(folder_alphapeptdeep,"inrich/output/output_arg_ms2_prediction_for_arg_chymotrypsin_231201.tsv")
        )
      ) %>% 
      select(-c(...1,b_modloss_z1,b_modloss_z2,y_modloss_z1,y_modloss_z2)) %>% 
      rename(
        b_z1_pred = b_z1,
        b_z2_pred = b_z2,
        y_z1_pred = y_z1,
        y_z2_pred = y_z2
      ) -> fdr_pcc_arg_set_preddf
    
    fdr_pcc_arg_set %>% 
      mutate(
        fragindex = map2(frag_start_idx,frag_stop_idx-1,\(x,y) seq(x,y,1))
      ) %>%
      unnest(fragindex) %>% 
      bind_cols(
        fdr_pcc_arg_set_fragdf,
        fdr_pcc_arg_set_preddf
      ) %>% 
      mutate(
        fragindex2=fragindex
      ) %>% 
      nest(
        fragdf = fragindex:y_z2,
        preddf = c(fragindex2,b_z1_pred:y_z2_pred)
      ) -> fdr_pcc_arg_set
    fdr_pcc_arg_set %>% 
      mutate(
        fragdf2 = map(fragdf,\(x) mutate(x,across(.cols=c(b_z1:y_z2),.fns=~na_if(.,y=0)))),
        preddf2 = map(preddf,\(x) mutate(x,across(.cols=c(b_z1_pred:y_z2_pred),.fns=~na_if(.,y=0))))
      ) -> fdr_pcc_arg_set
    
    fdr_pcc_arg_set %>% 
      mutate(
        bion_count = map_dbl(fragdf2,\(x) length(na.omit(append(pull(x,b_z1),pull(x,b_z2)))))
      ) -> fdr_pcc_arg_set
    
    fdr_pcc_arg_set %>%
      mutate(
        bion_pcc = map2_dbl(
          fragdf,
          preddf,
          function(x,y) {
            cor(append(pull(x,b_z1),pull(x,b_z2)),append(pull(y,b_z1_pred),pull(y,b_z2_pred)),method = "pearson")
          },
          .progress = TRUE
        )
      ) %>%
      mutate(
        yion_pcc = map2_dbl(
          fragdf,
          preddf,
          function(x,y) {
            cor(append(pull(x,y_z1),pull(x,y_z2)),append(pull(y,y_z1_pred),pull(y,y_z2_pred)),method = "pearson")
          },
          .progress = TRUE
        )
      ) -> fdr_pcc_arg_set
    
    summaryfile_psm_for_testing3 %>% 
      bind_cols(
        fdr_pcc_arg_set %>% select(spec_idx,PCC:yion_pcc)
      ) -> summaryfile_psm_for_testing3_pcc
  }
}
#summaryfile_psm_for_testing3_pcc



###specFDR specFNR for summaryfile_psm_for_testing3_pcc
{
  summaryfile_psm_for_testing3_pcc %>% 
    mutate(
      fragdf3 = map2(
        fragdf2,
        preddf2,
        \(x,y) left_join(x,y,join_by(fragindex == fragindex2))
      )
    ) -> summaryfile_psm_for_testing3_pcc
  
  summaryfile_psm_for_testing3_pcc %>% 
    mutate(
      fragdf3 = map(
        fragdf3,
        \(x) replace_na(x,list(b_z1=0,b_z2=0,y_z1=0,y_z2=0,b_z1_pred=0,b_z2_pred=0,y_z1_pred=0,y_z2_pred=0))
      )
    ) -> summaryfile_psm_for_testing3_pcc
  
  summaryfile_psm_for_testing3_pcc %>%
    mutate(
      fragdf4 = map(
        fragdf3,
        \(x) mutate(x,
                    b_z1_type = case_when(
                      b_z1 == 0 & b_z1_pred ==0 ~ NA,
                      b_z1 > 0 & b_z1_pred ==0 ~ "FN",
                      b_z1 == 0 & b_z1_pred >0 ~ "FP",
                      b_z1 > 0 & b_z1_pred >0 ~ "TP"
                    ),
                    b_z2_type = case_when(
                      b_z2 == 0 & b_z2_pred ==0 ~ NA,
                      b_z2 > 0 & b_z2_pred ==0 ~ "FN",
                      b_z2 == 0 & b_z2_pred >0 ~ "FP",
                      b_z2 > 0 & b_z2_pred >0 ~ "TP"
                    ),
                    y_z1_type = case_when(
                      y_z1 == 0 & y_z1_pred ==0 ~ NA,
                      y_z1 > 0 & y_z1_pred ==0 ~ "FN",
                      y_z1 == 0 & y_z1_pred >0 ~ "FP",
                      y_z1 > 0 & y_z1_pred >0 ~ "TP"
                    ),
                    y_z2_type = case_when(
                      y_z2 == 0 & y_z2_pred ==0 ~ NA,
                      y_z2 > 0 & y_z2_pred ==0 ~ "FN",
                      y_z2 == 0 & y_z2_pred >0 ~ "FP",
                      y_z2 > 0 & y_z2_pred >0 ~ "TP"
                    )
        ) %>% select(contains("type")) %>% 
          pivot_longer(cols = contains("type"))
      )
    )%>%
    mutate(
      FN = map_dbl(
        fragdf4,
        \(x) sum(pull(x,value) %in% c("FN"))
      )
    ) %>%
    mutate(
      TP = map_dbl(
        fragdf4,
        \(x) sum(pull(x,value) %in% c("TP"))
      )
    ) %>%
    mutate(
      FP = map_dbl(
        fragdf4,
        \(x) sum(pull(x,value) %in% c("FP"))
      )
    ) %>% 
    mutate(
      spec_FDR = FP/(TP+FP),
      spec_FNR = FN/(FN+TP)
    ) -> summaryfile_psm_for_testing3_pcc
}
##out::summaryfile_psm_for_testing3_pcc

###from bidding_pcc_fdr.R, best cutoff is cutoff_pcc or 0.868823 PCC #231212
{
  summaryfile_psm_for_testing3_pcc %>% 
    mutate(
      ms2_filter = if_else(PCC > label_pcc_fdr001, 1, 0)
    ) -> summaryfile_psm_for_testing3_pcc
  
  summaryfile_psm_for_testing3_pcc %>% 
    filter(ms2_filter == 1) ##623 out of 1217 #was 621 at 0.8719604
  
  summaryfile_psm_for_testing3_pcc %>% 
    filter(
      PCC>0.868823 & PCC<0.8719604
    ) %>% pull(`Protein Descriptions`)
  
  summaryfile_psm_for_testing3_pcc %>% 
    filter(gene == "CD55")
}
#summaryfile_psm_for_testing3_pcc@ms2_filter


############
###########

##RT training training REDUX! 231113! + TEST SET! + GVTESTSET! (231214)
##>>@summaryfile_psm_for_testing3
##>@summaryfile_psm_for_training3
##>dummies for training as only 2 peptides were there.
{
  readRDS(file = "./export/psm_563516.rds") -> psm_563516
  
  file_raw_index %>% 
    summarize(
      fraction= first(experiment),
      replicate= first(replicate),
      enzyme = first(enzyme),
      .by = set
    ) %>% 
    bind_cols(
      tibble(
        left_time = c(15,15,10,10,10,10,10,10,10,10,10,10), ###fix:: min for each training data
        right_time = c(140,140,158,158,158,158,158,158,158,158,158,158) ##fix:: 110 duration!,
      )
    ) -> file_raw_index2
  
  file_raw_index2 %>% 
    mutate(
      file_series = 1:12
    ) -> file_raw_index2
  
  psm_563516 %>% 
    mutate(
      nt_arg = if_else(!is.na(nmod) & P1prime == "R",1L,0L)
    ) -> psm_563516
  psm_563516 %>% summarize(
    .by = nt_arg,
    count = n()
  )
  unique(psm_563516$file_id_series)
  for (i in 1:length(unique(psm_563516$file_id_series))){
    psm_563516 %>% 
      filter(file_id_series == unique(psm_563516$file_id_series)[i]) %>% 
      dplyr::rename(
        spec_idx = index 
      ) %>% 
      select(
        strip_seq,spec_idx,mods,mod_sites,nAA,nce,instrument,rt_norm,left_time,right_time,nt_arg
      ) %>% 
      dplyr::rename(
        sequence = strip_seq
      ) %>%
      write_tsv(
        paste0(folder_alphapeptdeep,"input_rt_training_",i,".tsv")
      )
  }
  
  #testset
  for (i in 1:length(unique(summaryfile_psm_for_testing3$set))){
    summaryfile_psm_for_testing3 %>% 
      filter(set == unique(summaryfile_psm_for_testing3$set)[i]) %>% 
      dplyr::rename(
        spec_idx = index 
      ) %>% 
      select(
        sequence,spec_idx,irt,mods,mod_sites,nAA,nce,instrument,rt_norm
      ) %>%
      mutate(
        state = "genuine"
      ) %>% 
      bind_rows(
        summaryfile_psm_for_training3 %>% 
          filter(set == unique(summaryfile_psm_for_training3$set)[i]) %>% 
          dplyr::rename(
            spec_idx = index 
          ) %>% 
          select(
            Sequence,spec_idx,irt,mods,mod_sites,nAA,nce,instrument,rt_norm
          ) %>% 
          dplyr::rename(
            sequence = Sequence
          ) %>%
          mutate(
            state = "dummy"
          )
      ) %>% 
      write_tsv(
        paste0(folder_alphapeptdeep,"input_rt_testing_",i,".tsv")
      )
  }
  
  ##GV testset
  summaryfile_psm_for_testing3_rt_test_gv
  for (i in 1:length(unique(summaryfile_psm_for_testing3_rt_test_gv$set))){
    summaryfile_psm_for_testing3_rt_test_gv %>% 
      filter(set == unique(summaryfile_psm_for_testing3_rt_test_gv$set)[i]) %>% 
      mutate(
        sequence = paste0("GV",Sequence)
      ) %>% 
      mutate(
        nAA = nAA+1
      ) %>% 
      dplyr::rename(
        spec_idx = index 
      ) %>% 
      select(
        sequence,spec_idx,irt,mods,gv_mod_sites,nAA,nce,instrument,rt_norm
      ) %>%
      rename(mod_sites = gv_mod_sites) %>% 
      mutate(
        state = "genuine"
      ) %>% 
      bind_rows(
        summaryfile_psm_for_training3 %>% 
          filter(set == unique(summaryfile_psm_for_training3$set)[i]) %>% 
          dplyr::rename(
            spec_idx = index 
          ) %>% 
          select(
            Sequence,spec_idx,irt,mods,mod_sites,nAA,nce,instrument,rt_norm
          ) %>% 
          dplyr::rename(
            sequence = Sequence
          ) %>%
          mutate(
            state = "dummy"
          )
      ) %>% 
      write_tsv(
        paste0(folder_alphapeptdeep,"input_rt_testing_gvtest_",i,".tsv")
      )
  }
  
  psm_563516 %>% 
    mutate(
      nt_arg = as.character(nt_arg)
    ) %>% 
    summarize(
      count = n(),
      experiment=dplyr::first(fraction),
      enzyme=dplyr::first(enzyme),
      replicate = first(replicate),
      .by = c(file_id_series,nt_arg)
    ) %>% write_tsv("./figures/rt/rt_training_count_renewal.tsv")
  
  
}

###
# do alphapeptdeep rt training
###

####RT Training data summary REDUX!
##>@@summaryfile_psm_for_training3 >> summaryfile_psm_for_training_result
##>
{
  folder_alphapeptdeep
  i=1
  output_prediction_rt_train_result <- tibble()
  output_prediction_rt_test_result <- tibble()
  for (i in 1:12){
    output_prediction_rt_training <- read_tsv(
      paste0(folder_alphapeptdeep,"output_rt_training_",i,".tsv")
    )[,-1] %>% 
      arrange(rt_norm) %>% 
      mutate(file_series= i) %>% 
      filter(nt_arg == 1)
    
    output_prediction_rt_test <- read_tsv(
      paste0(folder_alphapeptdeep,"output_rt_testing_",i,".tsv")
    )[,-1] %>% 
      arrange(rt_norm) %>%
      filter(
        state == "genuine"
      ) %>% 
      mutate(file_series = i)
    
    fit_rt_training <- lm(rt_norm ~ rt_pred, data =output_prediction_rt_training)
    
    output_prediction_rt_training %>%
      bind_cols(
        predict(fit_rt_training, interval = "prediction", level =0.95, newdata = output_prediction_rt_training)
      )->output_prediction_rt_training
    
    
    output_prediction_rt_test %>% 
      bind_cols(
        predict(fit_rt_training, interval = "prediction", level =0.95, newdata = output_prediction_rt_test)
      )->output_prediction_rt_test
    
    
    output_prediction_rt_train_result %>% 
      bind_rows(
        output_prediction_rt_training
      ) -> output_prediction_rt_train_result
    
    output_prediction_rt_test_result %>% 
      bind_rows(
        output_prediction_rt_test
      ) -> output_prediction_rt_test_result
  }
  output_prediction_rt_train_result %>% 
    arrange(
      spec_idx
    ) ->output_prediction_rt_train_result
  
  output_prediction_rt_test_result %>% 
    arrange(
      spec_idx
    ) ->output_prediction_rt_test_result
  
  output_prediction_rt_test_result %>% 
    mutate(
      rt_test = if_else(rt_norm > lwr & rt_norm < upr, 1, 0)
    ) ->output_prediction_rt_test_result
  
  output_prediction_rt_train_result %>% 
    left_join(
      file_raw_index2 %>% select(-c(left_time,right_time)),
      join_by(file_series)
    ) -> output_prediction_rt_train_result
  
  output_prediction_rt_train_result %>% 
    mutate(
      rt_observed = rt_norm*(right_time-left_time)+left_time,
      rt_predicted = rt_pred*(right_time-left_time)+left_time,
      rt_fit = fit*(right_time-left_time)+left_time,
      rt_lwr = lwr*(right_time-left_time)+left_time,
      rt_upr = upr*(right_time-left_time)+left_time
    ) -> output_prediction_rt_train_result
  
  output_prediction_rt_test_result %>% 
    left_join(
      file_raw_index2,
      join_by(file_series)
    ) ->output_prediction_rt_test_result
  
  output_prediction_rt_test_result %>% 
    mutate(
      rt_observed = rt_norm*(right_time-left_time)+left_time,
      rt_predicted = rt_pred*(right_time-left_time)+left_time,
      rt_fit = fit*(right_time-left_time)+left_time,
      rt_lwr = lwr*(right_time-left_time)+left_time,
      rt_upr = upr*(right_time-left_time)+left_time
    ) ->output_prediction_rt_test_result
  
  
  output_prediction_rt_train_result
  output_prediction_rt_test_result
  
  
  ##R, MAE, T90
  {
    
    output_prediction_rt_train_result 
    output_prediction_rt_test_result 
    
    
    output_prediction_rt_train_result %>% 
      filter(file_series==1)
    
    fit_rt_training <- lm(
      rt_observed~ rt_predicted, 
      data =output_prediction_rt_train_result %>% 
        filter(file_series==1))
    
    fit_rt_training %>% tidy()
    
    output_prediction_rt_train_result %>% 
      mutate(t95=rt_upr-rt_lwr) %>% 
      summarise(
        .by = c(set,file_series),
        t95_mean = mean(t95)
      ) -> summary_rtprediction_t95
    summary_rtprediction_t95 %>% write_tsv(
      "./figures/rt/rt_t95.tsv"
    )
    
    summary_rtprediction <- tibble()
    for (i in 1L:12L){
      summary_rtprediction %>% 
        bind_rows(
          output_prediction_rt_train_result %>% 
            filter(file_series==i) %>% 
            metrics(
              truth = rt_observed, estimate = rt_predicted
            ) %>% 
            mutate(
              file_series =i
            )
        )->summary_rtprediction
    }
    summary_rtprediction %>% write_tsv("./figures/rt/summary_rtprediction.tsv")
    
    output_prediction_rt_train_result %>% 
      filter(file_series==i) %>% 
      metrics(
        truth = rt_observed, estimate = rt_predicted
      ) %>% 
      mutate(
        file_series =i
      )
    
  }
  
}
##>>>output_prediction_rt_test_result

##output rt gvtest sum and summarize
{
  output_prediction_rt_gvtest_result <- tibble()
  for (i in 1:12){
    bind_rows(
      output_prediction_rt_gvtest_result,
      read_tsv(
        paste0(folder_alphapeptdeep,"output_rt_testing_gvtest_",i,".tsv")
      )[,-1] %>% 
        arrange(rt_norm) %>%
        filter(
          state == "genuine"
        ) %>% 
        mutate(file_series = i)
    )->output_prediction_rt_gvtest_result
      
  }
  
  output_prediction_rt_gvtest_result %>% 
    left_join(
      output_prediction_rt_test_result %>% select(spec_idx,rt_pred) %>% 
        rename(rt_pred_original = rt_pred),
      join_by(
        spec_idx
      )
    )-> output_prediction_rt_gvtest_result
  
  output_prediction_rt_gvtest_result
}
#output_prediction_rt_gvtest_result

### merge summaryfile_psm_for_testing_pcc+rt_prediction+reposition_arg_no_dups 
### final check! for rt+pcc+diagnostic >> summaryfile_psm_for_testing3_pcc_rt
{
  output_prediction_rt_test_result
  
  summaryfile_psm_for_testing3_pcc %>% colnames()
  summaryfile_psm_for_testing3_pcc %>% 
    left_join(
      output_prediction_rt_test_result %>% select(spec_idx,rt_pred,fit,lwr,upr,rt_test),
      join_by(
        index == spec_idx
      )
    ) %>% 
    left_join(
      dbtable %>% rename(prot_seq_master=sequence) %>% select(accession,prot_seq_master),
      join_by(
        Accessions_Single == accession
      )
    ) %>% 
    rename(accession_master = Accessions_Single) %>% 
    mutate(
      position_master = str_locate(prot_seq_master,Sequence)[,"start"]
    ) %>% mutate(
      prot_pos_master = paste0(accession_master,"|",position_master)
    ) %>% 
    left_join(
      reposition_arg_no_dups %>%
        mutate(
          prot_name = str_extract(protein,"(?<=\\|)[A-Z0-9]+\\_HUMAN")
        ) %>% 
        select(seq,gene,prot_name,accession,position,p5p1,reviewed,under100),
      join_by(
        Sequence == seq
      )
    ) %>% 
    mutate(
      p1 = str_sub(p5p1,start=5,end=5)
    ) %>% 
    mutate(
      p2 = str_sub(p5p1,start=4,end=5)
    ) %>% 
    mutate(
      p3 = str_sub(p5p1,start=3,end=5)
    ) %>%
    mutate(
      prot_pos_reposition = paste0(gene,"|",accession,"|",prot_name,"|",position)
    ) %>% 
    select(-prot_seq_master)-> summaryfile_psm_for_testing3_pcc_rt
  
  summaryfile_psm_for_testing3_pcc_rt %>% 
    relocate(
      rt_norm, .before = rt_pred
    ) -> summaryfile_psm_for_testing3_pcc_rt
  
  summaryfile_psm_for_testing3_pcc_rt %>% write_tsv("./export/summaryfile_psm_for_testing3_pcc_rt_231201.tsv")
}
##>>>summaryfile_psm_for_testing3_pcc_rt


####ms2 t-test or homogeneity of fragment species mass  SMERT
#>>>summaryfile_psm_for_testing3_pcc_rt, file_index_arg6_ttest
{
  
  summaryfile_psm_for_testing3_pcc_rt %>% 
    bind_cols(
      file_index_arg6_ttest %>% 
        select(sequence_for_fragspec:b_y_ttest)
    ) -> summaryfile_psm_for_testing3_pcc_rt_hfsm
  summaryfile_psm_for_testing3_pcc_rt_hfsm %>% 
    mutate(b_y_ttest = replace_na(b_y_ttest,0)) -> summaryfile_psm_for_testing3_pcc_rt_hfsm
}
#>>>summaryfile_psm_for_testing3_pcc_rt_hfsm

### validation purpose! arginylome_validation_matrix
## arginylome_validation_matrix2 << w/o gv
{
  ##finalcheck
  summaryfile_psm_for_testing3_pcc_rt_hfsm %>%
    mutate(
      pcc_check = if_else(PCC>cutoff_pcc,1,0),###PCC cutoff
      rt_check = rt_test,
      b_y_ttest_check = if_else(b_y_ttest>=0.05,1,0),
      #diagnostic_ion_check = if_else(!is.na(diagnostic_ion_int_relative),1,0),
      p2gv_check = if_else(p2 == "GV" | p2 == "VG",1,0),
      p1r_check = if_else(p1 == "R",1,0),
    ) -> summaryfile_psm_for_testing3_pcc_rt_hfsm
  
  summaryfile_psm_for_testing3_pcc_rt_hfsm %>% 
    write_tsv("./export/summaryfile_psm_for_testing3_pcc_rt_hfsm_231201.tsv")
  
  summaryfile_psm_for_testing3_pcc_rt_hfsm %>% 
    select(prot_pos_master,prot_pos_reposition,sequence,pcc_check,rt_check,b_y_ttest_check,p2gv_check,p1r_check) -> arginylome_validation_matrix0
  
  summaryfile_psm_for_testing3_pcc_rt_hfsm %>%
    mutate(
      prot_pos_new = paste0(accession,"|",position)
    ) %>%
    left_join(
      output_cleavage %>% 
        rename(
          cleavage_name = name,
          cleavage_score = score,
          cleavage_type = type
        ) %>% 
        distinct(value,.keep_all=TRUE), ### 
      join_by(
        prot_pos_master == value
      )
    ) -> summaryfile_psm_for_testing3_pcc_rt_hfsm2
  
  summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
    mutate(
      dl_score = pcc_check+rt_check+b_y_ttest_check
    )-> summaryfile_psm_for_testing3_pcc_rt_hfsm2
  
  summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
    select(prot_pos_master,prot_pos_reposition,sequence,pcc_check,rt_check,b_y_ttest_check,p2gv_check,p1r_check) -> arginylome_validation_matrix0
  
  arginylome_validation_matrix0 %>% 
    mutate(
      score_sum = pcc_check+rt_check+b_y_ttest_check#+diagnostic_ion_check
    ) %>% 
    arrange(desc(score_sum)) %>% 
    distinct(prot_pos_master,.keep_all = TRUE) -> arginylome_validation_matrix0
  
  arginylome_validation_matrix0 %>% 
    dplyr::filter(p1r_check == 0) -> arginylome_validation_matrix1
  
  arginylome_validation_matrix1 %>% 
    dplyr::filter(p2gv_check !=1) -> arginylome_validation_matrix2
  
}
##arginylome_validation_matrix0::without filter
##arginylome_validation_matrix1::p1r removed
##arginylome_validation_matrix2::p2gv removed
#summaryfile_psm_for_testing3_pcc_rt_hfsm2
#arginylome_validation_matrix0
#arginylome_validation_matrix1
#arginylome_validation_matrix2

#####============================================
#####============================================
#####============================================

####fetching alphafold :: @summaryfile_psm_for_testing_pcc_rt >>alphafold_fetched_tidied
{
  library(protti)
  fetch_alphafold_prediction(
    uniprot_ids = distinct(summaryfile_psm_for_testing_pcc_rt,accession_master)$accession_master,
    version = "v4",
    timeout = 3600,
    return_data_frame = TRUE,
    show_progress = TRUE
  ) -> alphafold_fetched ##raw data
  
  readRDS("./export/alphafold_fetched.rds") -> alphafold_fetched
  
  alphafold_fetched %>% colnames
  alphafold_fetched %>% 
    distinct(
      uniprot_id,auth_seq_id,auth_comp_id,prediction_score,score_quality
    ) %>% 
    mutate(
      auth_seq_id = as.integer(auth_seq_id)
    ) ->alphafold_fetched_tidied
  
  
  
  alphafold_fetched_tidied %>% 
    nest(.by=uniprot_id) %>% 
    left_join(
      distinct(summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered,accession_master,position_master,Sequence,prot_pos_master),
      join_by(
        uniprot_id == accession_master
      )
    ) %>% drop_na() %>% 
    mutate(
      alphafold_subset = map(position_master,\(x) tibble(position = seq(x-10,x+10),mod_position = seq(-10,10)))
    ) %>% 
    mutate(
      alphafold_subset = map(alphafold_subset, \(x) dplyr::filter(x,position>0))
    ) %>%
    mutate(
      alphafold_subset = map2(alphafold_subset, data, \(x,y) left_join(x,y,join_by(position == auth_seq_id)))
    ) -> alphafold_fetched_tidied
  
  alphafold_fetched_tidied %>% 
    distinct(prot_pos_master,.keep_all = TRUE) %>% 
    inner_join(
      arginylome_validation_matrix0,join_by(prot_pos_master == prot_pos_master)
    ) %>% 
    select(-data) %>% 
    unnest(alphafold_subset) -> alphafold_fetched_tidied2
  
  remove(alphafold_fetched)
  
  alphafold_fetched_tidied2 %>% 
    mutate(
      gene = str_split_i(prot_pos_reposition,"\\|",1)
    ) %>% 
    left_join(
      arginylome_validation_matrix3_sequence_distinct %>% 
        select(gene,loc_first) %>% 
        distinct(gene,.keep_all = TRUE)
    ) -> alphafold_fetched_tidied2
}
###alphafold_fetched_tidied
###alphafold_fetched_tidied2 is for figures

###signalp + targetp+procleave + deeploc
###prot_list_all
###df.fasta
### >>> output_cleavage
### >>> output_deeploc
{
  base::union(
    summaryfile_psm_for_testing_pcc_rt$accession_master,
    summaryfile_psm_for_testing_pcc_rt$accession
  ) -> prot_list_all
  
  tibble(
    accession = prot_list_all
  ) %>% 
    left_join(
      dbtable_all,
      join_by(
        accession == accession
      )
    ) %>% 
    mutate(
      seq100 = str_sub(sequence,1,100)
    ) -> input_signalp
  
  ##input for signalp or targetp
  input_signalp %>% 
    select(accession,seq100) %>% 
    as.data.frame() %>% 
    dataframe2fas(file="./export/df.fasta")
  
  ###read signalp
  read_tsv(
    paste0(
      "./signalp/output.gff3"
    ),
    skip = 1,
    col_names = c("accession","program","type","start","end","score")
  )->output_signalp
  
  output_signalp %>% 
    mutate(
      type = "SP"
    ) %>% 
    select(
      c("accession","program","type","start","end","score")
    ) -> output_signalp
  
  #read targetp
  read_tsv(
    paste0(
      "./targetp/output_protein_type.txt"
    ),
    skip = 1
  ) %>% 
    filter(
      Prediction == "mTP"
    ) %>% 
    mutate(
      program = "TargetP-2.0",
      type = "mTP",
      score = mTP,
      start = 1L,
      end = str_extract(`CS Position`, "[0-9]+(?=\\-)") %>% as.integer(),
      accession = ID
    ) %>% 
    select(
      c("accession","program","type","start","end","score")
    ) -> output_targetp
  
  bind_rows(
    output_signalp,
    output_targetp
  ) -> output_cleavage
  
  ##procleave
  
  summaryfile_psm_for_testing_pcc_rt %>% 
    mutate(
      p4p4prime = paste0(str_sub(p5p1,2,5),str_sub(Sequence,1,4))
    ) %>% 
    select(
      prot_pos_master,p4p4prime
    ) %>% 
    distinct() %>% 
    as.data.frame() %>%
    dataframe2fas(file="./export/df_procleave.fasta")
  
  #read procleave output
  output_procleave <- tibble()
  for(i in 1:length(dir("./procleave/",pattern = "csv"))){
    read_csv(
      paste0(
        "./procleave/",dir("./procleave/",pattern = "csv")[i]
      )
    ) -> output_procleave_temp
    
    bind_rows(
      output_procleave,
      output_procleave_temp
    )->output_procleave
    
    rm(output_procleave_temp)
  }
  
  
  output_procleave %>% 
    mutate(
      p4p4prime = str_replace(`P4-P4' Site`,"\\†","")
    ) %>% 
    filter(
      Score >= 0.9 & Score <= 1
    ) %>% 
    arrange(
      desc(Score)
    ) %>% 
    distinct(
      p4p4prime,.keep_all = TRUE
    ) -> output_procleave
  
  output_procleave %>% 
    left_join(
      summaryfile_psm_for_testing_pcc_rt %>% 
        mutate(
          p4p4prime = paste0(str_sub(p5p1,2,5),str_sub(Sequence,1,4))
        ) %>% 
        select(
          prot_pos_master,p4p4prime
        ) %>% distinct(),
      join_by(p4p4prime == p4p4prime)
    ) %>% 
    select(
      prot_pos_master,p4p4prime,Score,Family
    ) -> output_procleave
  output_procleave %>% 
    rename(
      type=Family,
      score=Score,
      value = prot_pos_master
    )
  
  output_cleavage %>% 
    mutate(
      prot_pos_0 = paste0(accession,"|",end+1),
      prot_pos_neg1 = paste0(accession,"|",end),
      prot_pos_neg2 = paste0(accession,"|",end-1),
      prot_pos_pos1 = paste0(accession,"|",end+2),
      prot_pos_pos2 = paste0(accession,"|",end+3),
    ) %>% 
    select(type,score,prot_pos_0:prot_pos_pos2) %>% 
    pivot_longer(
      prot_pos_0:prot_pos_pos2
    ) %>% 
    bind_rows(
      output_procleave %>% 
        rename(
          type=Family,
          score=Score,
          value = prot_pos_master
        ) %>% 
        select(-p4p4prime) %>% 
        mutate(
          name = "procleave"
        )
    ) -> output_cleavage
  
  ##deeploc
  
  tibble(
    accession = prot_list_all
  ) %>% 
    left_join(
      dbtable_all,
      join_by(
        accession == accession
      )
    ) %>% 
    mutate(
      seq6000 = str_sub(sequence,1,6000)
    ) -> input_deeploc
  
  input_deeploc %>% 
    select(accession,seq6000) %>% 
    as.data.frame() %>% 
    dataframe2fas(file="./export/df_deeploc.fasta")
  
  read_csv(
    "./deeploc/results_6485A661000022229ACDEA12.csv"
  ) -> output_deeploc
}

###summaryfile_psm_for_testing_pcc_rt+output_cleavage
##arginylome_validation_matrix2
###SP>>MTP>>Cleavage
{
  
  
  summaryfile_psm_for_testing3_pcc_rt_hfsm %>% colnames
  
  
  
  summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
    left_join(
      output_deeploc %>% 
        select(
          Protein_ID,Localizations
        ),
      join_by(
        accession == Protein_ID
      )
    ) -> summaryfile_psm_for_testing3_pcc_rt_hfsm2
  
  arginylome_validation_matrix2 %>% 
    mutate(
      prot = str_extract(prot_pos_master,".+(?=\\|)")
    ) %>% 
    left_join(
      output_deeploc %>% 
        select(
          Protein_ID,Localizations
        ),
      join_by(
        prot == Protein_ID
      )
    ) %>% 
    left_join(
      output_cleavage %>% 
        rename(
          cleavage_name = name,
          cleavage_score = score,
          cleavage_type = type
        ) %>% 
        distinct(value,.keep_all=TRUE), ### 
      join_by(
        prot_pos_master == value
      )
    ) %>% 
    mutate(
      site_type = case_when(
        cleavage_type == "SP" ~ "SP",
        cleavage_type == "mTP" ~ "mTP",
        cleavage_name == "procleave" ~ "Protease"
      )
    )->arginylome_validation_matrix3
}
###summaryfile_psm_for_testing3_pcc_rt_hfsm2
###arginylome_validation_matrix3

###uniprot signal transit
{
  folder_uniprot <- "C:/Users/admin/Desktop/Data/DB/UniprotRefHuman/202302/"
  tibble(
    origin =  read_file(
      paste0(folder_uniprot,"UP000005640_9606.dat")
    )
  )%>% 
    separate_longer_delim(
      origin,
      delim = regex("\\/\\/\\\n")
    ) %>% 
    mutate(
      id = str_extract(origin,"ID.+(?=\\\n)"),
      ac = str_extract(origin,"(?<=AC[:blank:]{3}).+(?=\\;\\\n)"),
      prem_ac = str_split_i(ac,";",1),
      gn = str_extract(origin,"(?<=GN[:blank:]{3}Name\\=).+(?=[:blank:])"),
      signal= str_extract(origin,"(?<=FT[:blank:]{3}SIGNAL[:blank:]{10}).+(?=\\\n)"),
      transit= str_extract(origin,"(?<=FT[:blank:]{3}TRANSIT[:blank:]{9}).+(?=\\\n)"),
      propep= str_extract_all(origin,"(?<=FT[:blank:]{3}PROPEP[:blank:]{10}).+(?=\\\n)"),
    ) %>% 
    unnest(
      propep, keep_empty = TRUE
    ) %>% 
    filter(
      !is.na(signal) | !is.na(transit) | !is.na(propep)
    )%>% 
    mutate(
      signal_start = as.integer(str_extract(signal,"^[0-9]+")),
      signal_end = as.integer(str_extract(signal,"[0-9]+$")),
      transit_start = as.integer(str_extract(transit,"^[0-9]+")),
      transit_end = as.integer(str_extract(transit,"[0-9]+$")),
      propep_start = as.integer(str_extract(propep,"^[0-9]+")),
      propep_end = as.integer(str_extract(propep,"[0-9]+$"))
    ) %>% 
    select(-origin) -> uniprot_cleavage
  output_cleavage %>% distinct(type)
  
  uniprot_cleavage %>% mutate(
    type = case_when(
      !is.na(signal) ~ "SP",
      !is.na(transit) ~ "mTP",
      !is.na(propep) ~ "ProP",
      .default = NA
    )
  ) %>% 
    mutate(
      end = case_when(
        type == "SP" ~ signal_end,
        type == "mTP" ~ transit_end,
        type == "ProP" ~ propep_end,
        .default=NA_integer_
      )
    ) %>% 
    select(type,ac,end) %>% 
    mutate(
      accession = str_split_i(ac,regex("\\;[:blank:]"),1)
    ) %>% 
    mutate(
      score = 1.00,
      source_cleavage = "uniprot"
    ) %>% 
    select(type,score,accession,end) %>% 
    mutate(
      prot_pos_0 = paste0(accession,"|",end+1),
      prot_pos_neg1 = paste0(accession,"|",end),
      prot_pos_neg2 = paste0(accession,"|",end-1),
      prot_pos_pos1 = paste0(accession,"|",end+2),
      prot_pos_pos2 = paste0(accession,"|",end+3),
    ) %>% 
    select(type,score,prot_pos_0:prot_pos_pos2) %>% 
    pivot_longer(
      prot_pos_0:prot_pos_pos2
    ) %>% 
    drop_na()->uniprot_cleavage2
  
  uniprot_cleavage2 %>% 
    mutate(
      source_cleavage = "uniprot"
    ) -> uniprot_cleavage2
  
  #output_cleavage %>% write_rds("output_cleavage_240322.rds")
  
  #read_rds("output_cleavage_240322.rds") -> output_cleavage
  
  uniprot_cleavage2 %>% 
    bind_rows(
      output_cleavage
    ) %>% 
    arrange(desc(score)) -> output_cleavage
  
  #output_cleavage %>% write_rds("output_cleavage_240401.rds")
}
#@uniprot_cleavage
#@output_cleavage
#@output_cleavage_240401.rds

### value generation:: gene, prot
###@summaryfile_psm_for_testing3_pcc_rt_hfsm2 >> summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered
##>arginylome_genes
##>arginylome_uniprot_master
##>arginylome_uniprot_reposition
##>arginylome_prot_pos_master
{
  summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
    filter(
      p1r_check ==0,
      p2gv_check == 0
    ) -> summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered
  
  summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered %>% 
    select(
      gene
    ) %>% 
    distinct %>% .$gene -> arginylome_genes
  
  summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered %>% 
    select(
      `Master Protein Accessions`
    ) %>% 
    distinct %>% .$`Master Protein Accessions` -> arginylome_uniprot_master
  
  summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered %>% 
    select(
      accession
    ) %>% 
    distinct %>% .$accession -> arginylome_uniprot_reposition
  
  summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered %>% 
    select(
      prot_pos_master
    ) %>% 
    distinct %>% .$prot_pos_master -> arginylome_prot_pos_master
}
#summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered

####arg ##finding argpairs
#>>>summaryfile_psm
#>>>arginylome_pairs_prot_pos_master
#>>>summaryfile_psm_for_testing_pcc_rt2
{
  summaryfile_psm %>% colnames
  summaryfile_psm %>% pull(accession_master)
  
  
  ###DONE!
  # summaryfile_psm %>% 
  #   left_join(
  #     dbtable %>% rename(prot_seq_master=sequence) %>% select(accession,gene,prot_seq_master),
  #     join_by(
  #       accession_master == accession
  #     )
  #   ) %>% 
  #   #rename(accession_master = Accessions_Single) %>% 
  #   mutate(
  #     position_master = str_locate(prot_seq_master,Sequence)[,"start"]
  #   ) %>% mutate(
  #     prot_pos_master = paste0(accession_master,"|",position_master)
  #   ) %>% 
  #   mutate(
  #     p5p1 = str_sub(paste0("XXXXX",prot_seq_master),start=position_master,end=position_master+4)
  #   ) %>% 
  #   mutate(
  #     p1 = str_sub(p5p1,start=5,end=5)
  #   ) %>% 
  #   mutate(
  #     p2 = str_sub(p5p1,start=4,end=5)
  #   ) %>% 
  #   mutate(
  #     p3 = str_sub(p5p1,start=3,end=5)
  #   ) -> summaryfile_psm
  
  summaryfile_psm %>% 
    filter(prot_pos_master %in% arginylome_prot_pos_master) -> arginylome_psm_pairs
  arginylome_psm_pairs %>%
    filter(
      !str_detect(Modifications,"Arg")
    ) -> arginylome_psm_pairs
  
  arginylome_psm_pairs %>% 
    distinct(prot_pos_master) %>% .$prot_pos_master -> arginylome_pairs_prot_pos_master
  
  # summaryfile_psm_for_testing_pcc_rt2 %>% 
  #   mutate(
  #     arg_pair = if_else(prot_pos_master %in% arginylome_pairs_prot_pos_master,1,0)
  #   ) -> summaryfile_psm_for_testing_pcc_rt2
  
  summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered %>% 
    mutate(
      arg_pair = if_else(prot_pos_master %in% arginylome_pairs_prot_pos_master,1,0)
    ) -> summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered
  
  summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered %>% 
    filter(dl_score == 3) -> summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered_trinity
}

##cleavage
{
  table_protease_all### cleavage name source
  arginylome_psm_pcc_rt_di_basic_filtered %>% 
    filter(pcc_check == 1 & rt_check == 1) %>% 
    filter(cleavage_name == "procleave") %>% 
    distinct(prot_pos_reposition,.keep_all = TRUE) %>% 
    group_by(cleavage_type) %>% 
    summarize(
      site_count = n_distinct(prot_pos_reposition),
      site = paste0(prot_pos_reposition,collapse = ",")
    ) %>% 
    left_join(
      table_protease_all %>% select(cleavage_type,protease)
    ) ->table_protease
  
  table_protease %>% write_tsv("./export/protease.tsv")
}
#table_protease

##@arginylome_validation_matrix2 >> arginylome_validation_matrix3 (w/o DiagIon)
##summaryfile_psm_for_testing_pcc_rt2 >> table_nsite_score2_experiment_psmcount 
##arginylome_validation_matrix3 + table_nsite_score2_experiment_psmcount >> arginylome_validation_matrix_sequence (for quantitation matching!)
{
  arginylome_validation_matrix3 %>% 
    filter(
      score_sum == 3
    )
  
  arginylome_validation_matrix3 %>% 
    left_join(
      summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered %>% select(prot_pos_master,gene,Sequence)
      # ) %>% 
      # filter(
      #   pcc_check == 1 & rt_check == 1
    ) -> arginylome_validation_matrix3_sequence
  
  arginylome_validation_matrix3 %>% 
    left_join(
      arginylome_validation_matrix3_sequence %>% 
        group_by(prot_pos_master,Sequence) %>% 
        summarise(
          count = n()
        ) 
    ) -> arginylome_validation_matrix3_sequence
  
  summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered %>% select(
    prot_pos_master,prot_pos_reposition
  ) %>% distinct() -> table_reposition
  
  arginylome_validation_matrix3_sequence %>% 
    left_join(
      table_reposition
    ) %>% 
    mutate(
      gene = str_extract(prot_pos_reposition,"^[^\\|]+")
    )-> arginylome_validation_matrix3_sequence
  
  arginylome_validation_matrix3_sequence %>% 
    arrange(
      prot_pos_master, Sequence
    ) %>% 
    distinct(
      prot_pos_master,.keep_all = TRUE
    ) -> arginylome_validation_matrix3_sequence_distinct
  
  summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered %>% 
    group_by(prot_pos_master,experiment) %>% 
    # filter(
    #   pcc_check ==1 & rt_check == 1 & b_y_ttest_check == 1 & p2gv_check == 0 & p1r_check == 0 
    # ) %>% 
    summarize(
      count = n()
    ) %>% 
    pivot_wider(
      id_cols = prot_pos_master,
      names_from = experiment,
      values_from = count
    ) %>% 
    select(
      prot_pos_master,MOCK,MG132,MGTG
    ) %>% ungroup() %>% 
    mutate(
      across(c(MOCK,MG132,MGTG), \(x) replace_na(x,0))  
    )-> table_nsite_score3_experiment_psmcount
  
  arginylome_psm_pairs %>% 
    group_by(prot_pos_master,experiment) %>% 
    summarize(
      count = n()
    ) %>% pivot_wider(
      id_cols = prot_pos_master,
      names_from = experiment,
      values_from = count
    ) %>% 
    select(
      prot_pos_master,MOCK,MG132,MGTG
    ) %>% ungroup() %>% 
    mutate(
      across(c(MOCK,MG132,MGTG), \(x) replace_na(x,0))  
    ) %>% 
    rename(
      argpair_MOCK = MOCK,
      argpair_MG132 = MG132,
      argpair_MGTG = MGTG,
    )-> table_nsite_arg_pair_experiment_psmcount
  
  arginylome_validation_matrix3_sequence_distinct %>% 
    # filter(
    #   pcc_check ==1 & rt_check == 1 & b_y_ttest_check == 1
    # ) %>% 
    left_join(
      table_nsite_score3_experiment_psmcount
    ) %>% 
    left_join(
      table_nsite_arg_pair_experiment_psmcount
    ) -> arginylome_validation_matrix3_sequence_distinct
  
  arginylome_validation_matrix3_sequence_distinct %>% 
    mutate(
      across(c(argpair_MOCK,argpair_MG132,argpair_MGTG), \(x) replace_na(x,0))  
    ) -> arginylome_validation_matrix3_sequence_distinct
  
  arginylome_validation_matrix3_sequence_distinct %>% write_tsv("./export/arginylome_validation_matrix3_sequence_distinct_240118.tsv")
  
}

##HPA
#HPA
{
  read_tsv("C:/Users/admin/Desktop/Data/HumanProteinAtlas/proteinatlas.tsv") -> proteinatlas
  proteinatlas$Gene
  
  arginylome_validation_matrix3_sequence_distinct %>% 
    left_join(
      proteinatlas,
      join_by(gene == Gene)
    ) -> table_hpa_arg_all
  
  arginylome_validation_matrix3_sequence_distinct %>% 
    filter(
      pcc_check == 1 & rt_check == 1
    ) %>% 
    left_join(
      proteinatlas,
      join_by(gene == Gene)
    ) -> table_hpa_arg_trinity
  
  table_hpa_arg_trinity %>% write_tsv("./export/table_hpa_arg_trinity.tsv")
  
  proteinatlas %>% colnames()
  
  ##Protein Class
  {
    proteinatlas %>% 
      select(
        Gene,Ensembl,`Protein class`,`CCD Protein`
      ) %>% 
      separate_rows(`Protein class`,sep = ", ") %>% 
      pivot_wider(
        id_cols = c(Ensembl,Gene,`CCD Protein`),
        names_from = `Protein class`,
        values_from = `Protein class`
      ) %>% 
      mutate(
        across(!c(Ensembl,Gene,`CCD Protein`),\(x) !is.na(x))
      ) -> proteinatlas_proteinclass
    
    proteinatlas_proteinclass %>% 
      mutate(
        `CCD Protein` = if_else(`CCD Protein`=="yes",TRUE,FALSE)
      ) -> proteinatlas_proteinclass
    
    proteinatlas_proteinclass %>% 
      pivot_longer(
        cols = colnames(proteinatlas_proteinclass)[3:ncol(proteinatlas_proteinclass)],
      ) %>%
      filter(
        value == TRUE
      ) %>% 
      group_by(name) %>% 
      summarise(
        count = n_distinct(Ensembl)
      ) %>% 
      mutate(
        perc = count/sum(count)*100
      ) -> table_hpa_proteinclass_count
    
    arginylome_validation_matrix3_sequence_distinct %>% 
      left_join(
        proteinatlas_proteinclass,
        join_by(gene == Gene)
      ) -> arginylome_validation_matrix3_sequence_distinct_hpa_class
    
    
    arginylome_validation_matrix3_sequence_distinct_hpa_class %>% 
      filter(pcc_check == 1 & rt_check == 1 & b_y_ttest_check == 1) %>% 
      distinct(gene,.keep_all = TRUE) %>% 
      pivot_longer(
        cols = colnames(proteinatlas_proteinclass)[3:ncol(proteinatlas_proteinclass)],
      ) %>%
      filter(
        value == TRUE
      ) %>% 
      group_by(name) %>% 
      summarise(
        count = n_distinct(prot_pos_master),
        genes = paste0(gene,collapse=";")
      ) %>% 
      mutate(
        perc = count/sum(count)*100
      ) -> table_hpa_proteinclass_count_arg_trinity
    
    table_hpa_proteinclass_count %>% 
      left_join(
        table_hpa_proteinclass_count_arg_trinity,
        join_by(name == name)
      ) -> table_hpa_proteinclass_count_total_arg
  }
  
  ##Secretome location
  {
    proteinatlas %>% 
      select(
        Gene,Ensembl,`Secretome location`
      ) %>% 
      separate_rows(`Secretome location`,sep = ", ") %>% 
      pivot_wider(
        id_cols = c(Ensembl,Gene,),
        names_from = `Secretome location`,
        values_from = `Secretome location`
      ) %>% 
      mutate(
        across(!c(Ensembl,Gene),\(x) !is.na(x))
      ) -> proteinatlas_secretomelocation
    
    proteinatlas_secretomelocation %>% 
      pivot_longer(
        cols = colnames(proteinatlas_secretomelocation)[3:ncol(proteinatlas_secretomelocation)],
      ) %>%
      filter(
        value == TRUE
      ) %>% 
      group_by(name) %>% 
      summarise(
        count = n_distinct(Ensembl)
      ) %>% 
      mutate(
        perc = count/sum(count)*100
      ) -> table_hpa_secretomelocation_count
    
    arginylome_validation_matrix3_sequence_distinct %>% 
      left_join(
        proteinatlas_secretomelocation,
        join_by(gene == Gene)
      ) -> arginylome_validation_matrix3_sequence_distinct_hpa_secretomelocation
    
    arginylome_validation_matrix3_sequence_distinct_hpa_secretomelocation %>% 
      filter(pcc_check == 1 & rt_check == 1 & b_y_ttest_check == 1) %>% 
      distinct(gene,.keep_all = TRUE) %>% 
      pivot_longer(
        cols = colnames(proteinatlas_secretomelocation)[3:ncol(proteinatlas_secretomelocation)],
      ) %>%
      filter(
        value == TRUE
      ) %>% 
      group_by(name) %>% 
      summarise(
        count = n_distinct(prot_pos_master),
        genes = paste0(gene,collapse=";")
      ) %>% 
      mutate(
        perc = count/sum(count)*100
      ) -> table_hpa_secretomelocation_count_arg_trinity
    
    table_hpa_secretomelocation_count %>% 
      left_join(
        table_hpa_secretomelocation_count_arg_trinity,
        join_by(name == name)
      ) -> table_hpa_secretomelocation_count_total_arg
  }
  
  ##Secretome Function
  {
    proteinatlas %>% 
      select(
        Gene,Ensembl,`Secretome function`
      ) %>% 
      separate_rows(`Secretome function`,sep = ", ") %>% 
      pivot_wider(
        id_cols = c(Ensembl,Gene,),
        names_from = `Secretome function`,
        values_from = `Secretome function`
      ) %>% 
      mutate(
        across(!c(Ensembl,Gene),\(x) !is.na(x))
      ) -> proteinatlas_secretomefunction
    
    proteinatlas_secretomefunction %>% 
      pivot_longer(
        cols = colnames(proteinatlas_secretomefunction)[3:ncol(proteinatlas_secretomefunction)],
      ) %>%
      filter(
        value == TRUE
      ) %>% 
      group_by(name) %>% 
      summarise(
        count = n_distinct(Ensembl)
      ) %>% 
      mutate(
        perc = count/sum(count)*100
      ) -> table_hpa_secretomefunction_count
    
    arginylome_validation_matrix3_sequence_distinct_hpa_secretomelocation %>% 
      left_join(
        proteinatlas_secretomefunction %>% select(-`NA`),
        join_by(gene == Gene)
      ) -> arginylome_validation_matrix3_sequence_distinct_hpa_secretome
    
    arginylome_validation_matrix3_sequence_distinct_hpa_secretome %>% 
      filter(score_sum==3) %>% 
      distinct(gene,.keep_all = TRUE) %>% 
      pivot_longer(
        cols = colnames(proteinatlas_secretomefunction)[3:ncol(proteinatlas_secretomefunction)],
      ) %>%
      filter(
        value == TRUE
      ) %>% 
      group_by(name) %>% 
      summarise(
        count = n_distinct(prot_pos_master),
        genes = paste0(gene,collapse=";")
      ) %>% 
      mutate(
        perc = count/sum(count)*100
      ) -> table_hpa_secretomefunction_count_arg_trinity
    
    table_hpa_secretomefunction_count %>% 
      left_join(
        table_hpa_secretomefunction_count_arg_trinity,
        join_by(name == name)
      ) -> table_hpa_secretomefunction_count_total_arg
  }
  table_hpa_secretomefunction_count_arg_ms2_rt
  table_hpa_proteinclass_count_arg_ms2_rt
  table_hpa_secretomelocation_count_arg_ms2_rt
  
  table_hpa_proteinclass_count_total_arg
  table_hpa_secretomelocation_count_total_arg
  table_hpa_secretomefunction_count_total_arg
  
  
  table_hpa_arg_trinity %>% 
    distinct(gene,.keep_all = TRUE) %>% 
    separate_rows(`Biological process`,sep=", ") %>% 
    group_by(
      `Biological process`
    ) %>% 
    summarize(
      count = n_distinct(gene)
    ) %>% 
    arrange(desc(count)) %>%
    drop_na() %>% print(n=Inf)
  
  table_hpa_arg_trinity %>% 
    distinct(gene,.keep_all = TRUE) %>% 
    separate_rows(`Subcellular main location`,sep=", ") %>% 
    group_by(
      `Subcellular main location`
    ) %>% 
    summarize(
      count = n_distinct(gene)
    ) %>% 
    arrange(desc(count)) %>% 
    drop_na() %>% print(n=Inf)
  
}
###TRINITY
#arginylome_validation_matrix3_sequence_distinct_hpa_secretome
{
  arginylome_validation_matrix3_sequence_distinct_hpa_secretome %>% 
    mutate(
      trinity = if_else(score_sum == 3, "Passed","Failed")
    ) -> arginylome_validation_matrix3_sequence_distinct_hpa_secretome
  
  arginylome_validation_matrix3_sequence_distinct_hpa_secretome %>% filter(
    trinity == "Passed"
  ) -> arginylome_validation_matrix3_sequence_distinct_trinity
  
  arginylome_validation_matrix3_sequence_distinct_trinity %>% 
    write_tsv(
      "./export/arginylome_validation_matrix3_sequence_distinct_trinity_134_240118.tsv"
    )
}
#@arginylome_validation_matrix3_sequence_distinct_trinity

###peptide new position, from::summaryfile_psm to::reposition_arg_no_dups, reposition_arg
{
  file_index_arg6_ttest_pcc_rt1 %>% 
    distinct(strip_seq) %>% 
    pull(strip_seq) -> arginylated_peptides_all
  
  tibble(
    seq = arginylated_peptides_all
  ) %>% 
    mutate(
      protein = map(seq,\(x) dbtable_all %>% dplyr::filter(str_detect(sequence,x)),.progress=TRUE)
    ) %>% 
    unnest(protein) %>% 
    mutate(
      position = str_locate(sequence,seq)[,"start"]
    ) %>% 
    arrange(
      seq,position
    )  %>% 
    mutate(
      sequence2=paste0("XXXXX",sequence)
    ) %>% 
    mutate(
      p5p1 = str_sub(sequence2,start=position,end=position+4)
    ) %>% 
    mutate(
      reviewed = str_sub(protein,start=1,end=2)
    ) %>% 
    mutate(
      under100 = if_else(
        position <= 100, 1,0
      )
    ) %>% 
    arrange(
      seq,desc(under100),reviewed,position,accession
    ) %>% 
    mutate(
      count = n(),.by=seq
    ) %>% 
    mutate(
      p2p1 = str_sub(p5p1,start=4),
      p1p1 = str_sub(p5p1,start=5),
      p2gv = if_else(p2p1 == "GV"|p2p1 == "VG", TRUE,FALSE),
      p1r = if_else(p1p1 == "R", TRUE,FALSE)
    )->reposition_arg
  
  reposition_arg %>% 
    arrange(
      seq,reviewed,position,gene
    ) %>% 
    mutate(
      prot_pos_reposition = paste0(gene,"|",accession,"|",position)
    ) ->  reposition_arg
  
  reposition_arg %>% 
    mutate(
      trinity = if_else(seq %in% (arginylome_validation_matrix3_sequence_distinct_trinity %>% pull(Sequence)),TRUE,FALSE)
    ) -> reposition_arg
  
  reposition_arg %>% 
    filter(trinity == TRUE) -> reposition_arg_trinity
  
  reposition_arg_trinity %>%
    select(-sequence, -sequence2) %>% 
    write_tsv("./export/reposition_arg_trinity.tsv")
  
  reposition_arg %>% 
    filter(position<=2)
  
  reposition_arg %>%
    distinct(seq) ##392
  
  dbtable %>% pull(accession)
  
  reposition_arg %>% 
    mutate(
      dbtable = if_else(
        accession %in% (dbtable %>% pull(accession)), "ref","additional"
      )
    )->reposition_arg
  
  reposition_arg %>% 
    mutate(
      prot_nterm = if_else(position %in% c(1,2),TRUE,FALSE)
    ) %>% 
    filter(
      dbtable == "ref" | prot_nterm == TRUE
    ) %>% 
    arrange(seq,desc(prot_nterm),desc(dbtable),gene) -> reposition_arg
  
  reposition_arg %>% 
    distinct(
      seq,.keep_all = TRUE
    ) -> reposition_arg_no_dups
  
  reposition_arg %>%
    select(-c(sequence,sequence2)) %>% 
    write_tsv("reposition_arg.tsv")
  reposition_arg_no_dups %>%
    select(-c(sequence,sequence2)) %>% 
    write_tsv("reposition_arg_no_dups.tsv")
}
#@reposition_arg
#@reposition_arg_no_dups


##GENE POSITIon!!!!!!
#summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered_trinity gene position
{
  
  arginylome_validation_matrix3_sequence_distinct_trinity %>% 
    mutate(
      gene_position = paste0(gene,"|",str_split_i(prot_pos_master,"\\|",2))
    ) -> arginylome_validation_matrix3_sequence_distinct_trinity
  
  arginylome_validation_matrix3_sequence_distinct_trinity %>% 
    left_join(
      proteinatlas,
      join_by(
        prot == Uniprot
      )
    ) -> arginylome_validation_matrix3_sequence_distinct_trinity
  
  arginylome_validation_matrix3_sequence_distinct_trinity %>% 
    left_join(
      peptide_lfq_trinity_distinct_arg %>% 
        select(-c(mean,stdev)) %>% 
        rename_with(
          \(x) paste0(x,"_arg"), .cols = c(MG132:MOCK)
        )
    ) %>% 
    left_join(
      peptide_lfq_trinity_distinct_free %>% 
        select(-c(mean,stdev)) %>% 
        rename_with(
          \(x) paste0(x,"_free"), .cols = c(MG132:MOCK)
        )
    ) %>% 
    write_tsv(
      "./export/arginylome_validation_matrix3_sequence_distinct_trinity.tsv"
    )
  
  
  arginylome_validation_matrix3_sequence_distinct_trinity %>% 
    select(Sequence,gene_position) -> table_gene_position
  
  proteinatlas %>% write_tsv(
    "proteinatlas.tsv"
  )
}
#@summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered
#@summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered_trinity
#@arginylome_validation_matrix3_sequence_distinct_trinity

###updating arginylome table if cleavage info has been updated
#output_cleavage
{
  output_cleavage %>% 
    rename(
      cleavage_type = type,
      cleavage_score = score,
      cleavage_pos = name
    ) -> output_cleavage
  
  output_cleavage %>% 
    arrange(desc(cleavage_score)) %>% 
    distinct(value,.keep_all = TRUE) -> output_cleavage
  
  arginylome_validation_matrix3_sequence_distinct_trinity %>% 
    select(
      -c(cleavage_type,cleavage_score,cleavage_pos,source_cleavage)
    ) %>% 
    left_join(
      output_cleavage,
      join_by(
        prot_pos_master == value
      )
    ) -> arginylome_validation_matrix3_sequence_distinct_trinity
  
  
  arginylome_validation_matrix3_sequence_distinct_trinity %>% 
    select(-c(protease)) %>% 
    left_join(
      table_protease_all %>% select(cleavage_type,protease)
    ) -> arginylome_validation_matrix3_sequence_distinct_trinity
    
    
  arginylome_validation_matrix3_sequence_distinct_trinity %>% 
    write_tsv(
      "./export/arginylome_validation_matrix3_sequence_distinct_trinity_240401.tsv"
    )
 
  
  # read_tsv(
  #   "./export/arginylome_validation_matrix3_sequence_distinct_trinity.tsv"
  # ) -> arginylome_validation_matrix3_sequence_distinct_trinity
}
#@output_cleavage
#@arginylome_validation_matrix3_sequence_distinct_trinity


#########################ADDITIONAL



##diagnostic ions in training set
{
  read_rds(
    "./export/psm_563516.rds"
  ) -> psm_563516
  
  psm_563516 %>% 
    mutate(
      peak_table = map(
        peak_table, \(x) mutate(x,matches = if_else(`m/z` < 202.138331+0.005 & `m/z` > 202.138331-0.005, "diagnostic", matches)),
        .progress = TRUE
      )
    ) -> psm_563516
  
  psm_563516 %>% 
    mutate(
      ion_int_relative_max = map_dbl(
        peak_table_alpha_ms2,
        \(x) max(pull(x,i))
      )
    ) %>% 
    mutate(
      diagnostic_ion_mz = map(peak_table, \(x) x$`m/z`[str_which(x$matches, "diagnostic")] )
    ) %>% 
    mutate(
      diagnostic_ion_mz = map_dbl(diagnostic_ion_mz,\(x) ifelse(length(x)>0,x,NA))
    ) %>% 
    mutate(
      diagnostic_ion_int = map(peak_table, \(x) x$i[str_which(x$matches, "diagnostic")] )
    ) %>%
    mutate(
      diagnostic_ion_int = map_dbl(diagnostic_ion_int,\(x) ifelse(length(x)>0,x,NA))
    ) %>% 
    mutate(
      diagnostic_ion_int_relative = diagnostic_ion_int/ion_int_relative_max
    ) -> psm_563516
  
  psm_563516 %>% 
    mutate(
      artefact_arg = if_else(
        is.na(artefact_arg),FALSE,artefact_arg
      )
    ) -> psm_563516
  
  psm_563516 %>% 
    write_rds(
      "./export/psm_563516.rds"
    )
}
#psm_563516@ion_int_relative_max
#psm_563516@diagnostic_ion_mz
#psm_563516@diagnostic_ion_int
#psm_563516@diagnostic_ion_int_relative




###frag_theo_spec >>>> peak_table_forgraph
##+peak_table from file_index_arg6_ttest_pcc_rt1
#summaryfile_psm_for_testing3_pcc_rt_hfsm2$peak_table_alpha_ms2_vs_theo_spec
#summaryfile_psm_for_testing3_pcc_rt_hfsm2
{
  head(file_index_arg6_ttest_pcc_rt1)$peak_table
  
  ##peak_table object
  summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>%
    bind_cols(
      select(file_index_arg6_ttest_pcc_rt1,peak_table)
    ) -> summaryfile_psm_for_testing3_pcc_rt_hfsm2
  
  ##to check peak_table
  summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
    head() %>% 
    pull(peak_table)
    
  
  ##to create peak_table_for_graph
  summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>%  ##wo modloss!
    mutate(
      peak_table_for_graph = map(peak_table,\(x) mutate(x,matches = if_else(!str_detect(matches,"NH3|Ox|H2O"),matches,NA)))
    ) %>% 
    mutate(
      peak_table_for_graph= map(peak_table_for_graph,\(x) mutate(x, series = str_extract(matches,"b|y")))
    ) %>% 
    mutate(
      peak_table_for_graph= map(
        peak_table_for_graph,
        \(x) mutate(x, ion_num = as.integer(str_extract(matches,"(?<=\\()[0-9]+(?=\\))")))
      )
    ) %>% ###b ion +1 because of modification! 
    mutate(
      peak_table_for_graph= map(
        peak_table_for_graph,
        \(x) mutate(x, ion_num = if_else(series == "b",ion_num +1,ion_num))
      )
    ) %>%
    mutate(
      peak_table_for_graph= map(
        peak_table_for_graph,
        \(x) mutate(x, charge = as.integer(str_extract(matches,"(?<=\\()[0-9]+(?=\\+\\))")))
      )
    ) -> summaryfile_psm_for_testing3_pcc_rt_hfsm2
  
  
  ##combination #peak_table_for_graph+peak_table_alpha_ms2_vs_theo_spec+fragdf3(fragdf+preddf)
  summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
    head() %>% 
    pull(peak_table_for_graph)
  
  summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
    head() %>% 
    pull(peak_table_alpha_ms2_vs_theo_spec)
  
  summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
    mutate(
      alpha_ms2_frag_pred = map(
        fragdf3,
        \(x) mutate(x, bion_num = 1:nrow(x)) %>% mutate(yion_num = max(bion_num)-bion_num+1)
      )
    ) %>%
    mutate(
      alpha_ms2_frag_pred = map(
        alpha_ms2_frag_pred,
        \(x) pivot_longer(x, cols = c(b_z1:y_z2_pred)) %>% 
          mutate(prediction = if_else(str_detect(name,"pred"),"pred","frag")) %>% 
          mutate(name = str_replace(name,"\\_pred","")) %>% 
          rename(type = "name")
      )
    ) %>% 
    mutate(
      alpha_ms2_frag_pred = map(
        alpha_ms2_frag_pred,
        \(x) pivot_wider(x, id_cols = c(fragindex,bion_num,yion_num,type),names_from = prediction) %>% 
          mutate(
            series = str_extract(type,"^."),
            ion_num = as.integer(if_else(series == "b",bion_num,yion_num)),
            charge = as.integer(str_extract(type,"(?<=\\_z)[0-9]"))
          ) 
      )
    ) %>%
    mutate(
      peak_table_for_graph = map2(
        peak_table_for_graph,
        alpha_ms2_frag_pred,
        \(x,y) left_join(x,y %>% select(series,ion_num,charge,frag,pred), join_by(series,ion_num,charge))
      )
    ) ->summaryfile_psm_for_testing3_pcc_rt_hfsm2
  
  summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
    mutate(
      peak_table_for_graph = map2(
        peak_table_for_graph,
        frag_theo_spec,
        \(x,y) left_join(x,y %>% select(series,ion_num,charge,ms2mz), join_by(series,ion_num,charge))
      )
    ) %>% 
    mutate(
      peak_table_for_graph=map(
        peak_table_for_graph,
        \(x) mutate(x,ms2error = `m/z`-ms2mz) %>% mutate(
          ms2error_ppm = ms2error/`m/z`*1e6
        )
      )
    ) -> summaryfile_psm_for_testing3_pcc_rt_hfsm2
  
  summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
    mutate(
      theo_spec_preddf = map2(
        frag_theo_spec,
        alpha_ms2_frag_pred,
        \(x,y) left_join(x,y %>% select(series,ion_num,charge,pred),join_by(series,ion_num,charge)) %>% mutate(pred = -pred)
      )
    )%>%
    mutate(
      peak_table_for_graph = map2(
        peak_table_for_graph,
        theo_spec_preddf,
        \(x,y) bind_rows(x %>% mutate(frag = i/max(i)),y %>% filter(pred !=0) %>% rename(`m/z`=ms2mz,frag = pred) %>% select(`m/z`,series,ion_num,charge,frag))
      )
    ) -> summaryfile_psm_for_testing3_pcc_rt_hfsm2
  
  summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
    mutate(
      peak_table_for_graph = map(
        peak_table_for_graph,
        \(x) mutate(x, ion_type = paste0(series,ion_num," ",charge,"+"))
      )
    ) -> summaryfile_psm_for_testing3_pcc_rt_hfsm2
  
  summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
    head(n=1) %>% 
    pull(frag_theo_spec)
  
  summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
    head(n=1) %>% 
    pull(alpha_ms2_frag_pred)
  
  summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
    head(n=1) %>% 
    pull(peak_table_for_graph)
  
  write_rds(
    summaryfile_psm_for_testing3_pcc_rt_hfsm2,
    "./export/summaryfile_psm_for_testing3_pcc_rt_hfsm2.rds"
  )
}
#summaryfile_psm_for_testing3_pcc_rt_hfsm2$peak_table_forgraph
#summaryfile_psm_for_testing3_pcc_rt_hfsm2$theo_spec_preddf
#summaryfile_psm_for_testing3_pcc_rt_hfsm2$alpha_ms2_frag_pred


####231212 UPDATE if PCC cut off have changed
{
  ##new prot_pos_reposition
  {
    reposition_arg
    reposition_arg_no_dups
    reposition_arg_trinity
  }
  
  #new cutoff: label_pcc_fdr001 0.868823+ chymotrypsin PCC >0.9
  summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
    filter(enzyme == "chymotrypsin") %>%
    select(prot_pos_reposition,PCC)
  
  #to change 1: PSM level. 231212: trypsin=newfdr, chymo>=0.9
  {
    summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
      mutate(
        pcc_check = case_when(
          enzyme == "trypsin" & PCC > label_pcc_fdr001 ~ 1,
          enzyme == "chymotrypsin" & PCC >= 0.9 ~ 1,
          .default = 0
        )
      ) %>% 
      mutate(
        dl_score = pcc_check + rt_check + b_y_ttest_check,
        score_sum = pcc_check + rt_check + b_y_ttest_check
      ) -> summaryfile_psm_for_testing3_pcc_rt_hfsm2
    summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered %>% 
      mutate(
        pcc_check = case_when(
          enzyme == "trypsin" & PCC > label_pcc_fdr001 ~ 1,
          enzyme == "chymotrypsin" & PCC >= 0.9 ~ 1,
          .default = 0
        )
      ) %>% 
      mutate(
        dl_score = pcc_check + rt_check + b_y_ttest_check,
        score_sum = pcc_check + rt_check + b_y_ttest_check
      ) -> summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered
    
    summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered %>% 
      filter(
        score_sum == 3
      ) -> summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered_trinity
  }
  #summaryfile_psm_for_testing3_pcc_rt_hfsm2
  #summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered
  #summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered_trinity
  
  #to change 2: site level.
  {
    summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
      mutate(
        pcc_check = case_when(
          enzyme == "trypsin" & PCC > label_pcc_fdr001 ~ 1,
          enzyme == "chymotrypsin" & PCC >= 0.9 ~ 1,
          .default = 0
        )
      ) %>% 
      mutate(
        dl_score = pcc_check + rt_check + b_y_ttest_check,
        score_sum = pcc_check + rt_check + b_y_ttest_check
      ) -> summaryfile_psm_for_testing3_pcc_rt_hfsm2
    summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered %>% 
      mutate(
        pcc_check = case_when(
          enzyme == "trypsin" & PCC > label_pcc_fdr001 ~ 1,
          enzyme == "chymotrypsin" & PCC >= 0.9 ~ 1,
          .default = 0
        )
      ) %>% 
      mutate(
        dl_score = pcc_check + rt_check + b_y_ttest_check,
        score_sum = pcc_check + rt_check + b_y_ttest_check
      ) -> summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered
    
    summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered %>% 
      filter(
        score_sum == 3
      ) -> summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered_trinity
  }
  #summaryfile_psm_for_testing3_pcc_rt_hfsm2
  #summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered
  #summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered_trinity
}

####231219 UPDATE for Th error to Da error for SMERT
{
  summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
    mutate(
      peak_table_alpha_ms2_vs_theo_spec = map(
        peak_table_alpha_ms2_vs_theo_spec,
        \(x) mutate(x, error = error * charge)
      )
    ) -> summaryfile_psm_for_testing3_pcc_rt_hfsm2
    
  summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
    mutate(
      b_ion_count = map_int(
        peak_table_alpha_ms2_vs_theo_spec,
        \(x) filter(x, series == "b") %>% nrow()
      ),
      y_ion_count = map_int(
        peak_table_alpha_ms2_vs_theo_spec,
        \(x) filter(x, series == "y") %>% nrow()
      ),
      b_ion_errors = map(
        peak_table_alpha_ms2_vs_theo_spec,
        \(x) filter(x, series == "b") %>% pull(error)
      ),
      y_ion_errors = map(
        peak_table_alpha_ms2_vs_theo_spec,
        \(x) filter(x, series == "y") %>% pull(error)
      )
    )-> summaryfile_psm_for_testing3_pcc_rt_hfsm2
  
  
  summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
    mutate(
      b_y_ttest = map2_dbl(
        b_ion_errors,
        y_ion_errors,
        \(x,y) if(length(x)>1 & length(y)>1){t.test(unlist(x),unlist(y),var.equal = TRUE)$p.value} else{return(NA)}
      )
    ) -> summaryfile_psm_for_testing3_pcc_rt_hfsm2
  
  summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
    mutate(
      b_y_ttest_check = if_else(b_y_ttest>=0.05,1,0),
      b_y_ttest_check = if_else(is.na(b_y_ttest_check),0,b_y_ttest_check)
    ) -> summaryfile_psm_for_testing3_pcc_rt_hfsm2
  
  summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
    mutate(
      dl_score = pcc_check+rt_check+b_y_ttest_check
    ) -> summaryfile_psm_for_testing3_pcc_rt_hfsm2
  
  summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
    filter(
      p1r_check != 1 & p2gv_check !=1
    ) -> summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered
  
  summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered %>% 
    filter(
      dl_score == 3
    ) -> summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered_trinity
  
  summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
    select(prot_pos_master,prot_pos_reposition,sequence,pcc_check,rt_check,b_y_ttest_check,p2gv_check,p1r_check) -> arginylome_validation_matrix0
  
  arginylome_validation_matrix0 %>% 
    mutate(
      score_sum = pcc_check+rt_check+b_y_ttest_check#+diagnostic_ion_check
    ) %>% 
    arrange(desc(score_sum)) %>% 
    distinct(prot_pos_master,.keep_all = TRUE) -> arginylome_validation_matrix0
  
  arginylome_validation_matrix0 %>% 
    dplyr::filter(p1r_check == 0) -> arginylome_validation_matrix1
  
  arginylome_validation_matrix1 %>% 
    dplyr::filter(p2gv_check !=1) -> arginylome_validation_matrix2
  
  arginylome_validation_matrix2 %>% 
    select(b_y_ttest_check)
  
  ####
  arginylome_validation_matrix3 %>% 
    select(
      -b_y_ttest_check
    ) %>% 
    bind_cols(
      arginylome_validation_matrix2 %>% 
        select(b_y_ttest_check)
    ) -> arginylome_validation_matrix3
  
  arginylome_validation_matrix3 %>% 
    left_join(
      arginylome_validation_matrix3_sequence %>% 
        group_by(prot_pos_master,Sequence) %>% 
        summarise(
          count = n()
        ) 
    ) %>% 
    left_join(
      table_reposition
    ) %>% 
    mutate(
      gene = str_extract(prot_pos_reposition,"^[^\\|]+")
    )-> arginylome_validation_matrix3_sequence
  
  arginylome_validation_matrix3_sequence %>% 
    mutate(
      score_sum = pcc_check + rt_check + b_y_ttest_check 
    ) -> arginylome_validation_matrix3_sequence
  #####
  
  ####
  arginylome_validation_matrix3_sequence_distinct %>% 
    select(
      -b_y_ttest_check
    ) %>% 
    bind_cols(
      arginylome_validation_matrix3_sequence %>% 
        arrange(
          prot_pos_master, Sequence
        ) %>% 
        distinct(
          prot_pos_master,.keep_all = TRUE
        ) %>% 
        select(
          b_y_ttest_check
        )
    ) -> arginylome_validation_matrix3_sequence_distinct
  
  arginylome_validation_matrix3_sequence_distinct %>% 
    mutate(
     score_sum = pcc_check + rt_check + b_y_ttest_check 
    ) -> arginylome_validation_matrix3_sequence_distinct
  ####
  
  arginylome_validation_matrix3_sequence_distinct %>% 
    mutate(
      trinity = if_else(score_sum == 3,"Passed","Failed")
    ) -> arginylome_validation_matrix3_sequence_distinct
  
  ###
  arginylome_validation_matrix3_sequence_distinct_hpa_secretome %>% 
    select(
      -b_y_ttest_check
    ) %>% 
    bind_cols(
      arginylome_validation_matrix3_sequence %>% 
        arrange(
          prot_pos_master, Sequence
        ) %>% 
        distinct(
          prot_pos_master,.keep_all = TRUE
        ) %>% 
        select(
          b_y_ttest_check
        )
    ) %>% 
    mutate(
      score_sum = pcc_check + rt_check + b_y_ttest_check 
    ) -> arginylome_validation_matrix3_sequence_distinct_hpa_secretome
  ####
  
  ####
  arginylome_validation_matrix3_sequence_distinct_hpa_secretome %>% 
    filter(
      score_sum == 3
    ) -> arginylome_validation_matrix3_sequence_distinct_trinity
  ####
  
  ##p5p5prime
  summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
    select(
      Sequence,p5p5prime
    ) %>% 
    distinct(Sequence,.keep_all = TRUE)-> table_p5p5prime
  
  arginylome_validation_matrix3_sequence_distinct_trinity %>% 
    left_join(
      table_p5p5prime,
      join_by(Sequence)
    ) -> arginylome_validation_matrix3_sequence_distinct_trinity
  
  arginylome_validation_matrix3_sequence_distinct %>% 
    left_join(
      table_p5p5prime,
      join_by(Sequence)
    ) -> arginylome_validation_matrix3_sequence_distinct
  
  
  arginylome_validation_matrix3_sequence_distinct %>% 
    mutate(
      loc_first = str_split_i(Localizations,"\\|",1)
    ) -> arginylome_validation_matrix3_sequence_distinct
  
  arginylome_validation_matrix3_sequence_distinct_trinity %>% 
    mutate(
      loc_first = str_split_i(Localizations,"\\|",1)
    ) -> arginylome_validation_matrix3_sequence_distinct_trinity
  
  arginylome_validation_matrix3_sequence_distinct_trinity %>% 
    mutate(
      protease_site = case_when(
        cleavage_type == "C14.001" ~ "Caspase-1",
        cleavage_type == "C14.005" ~ "Caspase-6",
        cleavage_type == "C14.003" ~ "Caspase-3",
        cleavage_type == "M12.033" ~ "LAST_MAM peptidase",
        cleavage_type == "SP" ~ "Signal Peptide",
        cleavage_type == "mTP" ~ "Transit Peptide",
        cleavage_type == "M12.004" ~ "meprin beta subunit",
        cleavage_type == "A01.009" ~ "cathepsin D",
        cleavage_type == "M12.002" ~ "matrix metallopeptidase-8",
        cleavage_type == "S01.010" ~ "granzyme B",
        .default = cleavage_type
      )
    ) -> arginylome_validation_matrix3_sequence_distinct_trinity
}

##IUPRED3 analysis
{
  #input file
  dbtable %>% 
    filter(
      accession %in% unique(arginylome_validation_matrix3_sequence_distinct_trinity$prot)
    ) %>% 
    mutate(merged = paste0(protein,"$",sequence)) %>% 
    select(merged) %>% 
    mutate(merged = paste0(">",merged)) %>% 
    tidyr::separate_rows(merged,sep = "\\$") %>% 
    write_tsv(paste0("./export/target_fasta.fasta"),col_names = FALSE)
}

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
###END------------------------------------------------------------------------------------









####230728 added
{
  #output_metrics_describ3_230610.tsv
  file_index_arg6_ttest %>% 
    bind_cols(
      read_tsv(
        paste0(folder_alphapeptdeep,"output_metrics_describ3_230610.tsv")
      ) %>% select(
        PCC:SPC
      )
    ) -> file_index_arg6_ttest_pcc
  file_index_arg6_ttest_pcc %>% 
    bind_cols(
      output_prediction_rt_test_result %>% select(rt_norm:rt_test) 
    ) -> file_index_arg6_ttest_pcc_rt
  
  file_index_arg6_ttest_pcc_rt %>% 
    mutate(
      pcc_passed = if_else(PCC>= label_pcc_fdr001, 1L,0L)
    ) %>% 
    mutate(
      rt_passed = rt_test
    ) %>% 
    mutate(
      ttest_passed = if_else(!is.na(b_y_ttest) & b_y_ttest >= 0.05, 1L, 0L)
    ) %>% 
    mutate(
      pass_count = fct(as.character(pcc_passed+rt_passed+ttest_passed),levels=c("0","1","2","3"))
    )-> file_index_arg6_ttest_pcc_rt1
  
  
  file_index_arg6_ttest_pcc_rt1 %>% 
    mutate(
      strip_seq=str_sub(strip_seq,start=2L)
    )-> file_index_arg6_ttest_pcc_rt1
  
  
}



##from file_index_arg6_ttest_pcc_rt1
{
  file_index_arg6_ttest_pcc_rt1 %>% 
    left_join(
      reposition_arg_no_dups %>% select(
        seq,accession,gene,prot_name,position,p5p1,reviewed,p2p1,p2p1_list,prot_pos_reposition
      ),
      join_by(strip_seq == seq)
    ) -> file_index_arg6_ttest_pcc_rt2
  
  file_index_arg6_ttest_pcc_rt2 %>% 
    mutate(
      p1r_passed = if_else(!str_detect(p2p1_list,"[^\\;]R"),1L,0L)
    ) %>% 
    mutate(
      p2gv_passed = if_else(!str_detect(p2p1_list,"[VG]{2}"),1L,0L)
    ) %>% 
    mutate(
      db_pass_count = p1r_passed+p2gv_passed
    )-> file_index_arg6_ttest_pcc_rt2
  
  file_index_arg6_ttest_pcc_rt2 %>% write_tsv("./export/inrich_file_index_arg6_ttest_pcc_rt2_230728.tsv") 
  
}


####
#563516-17669 MS2 prediction model
####
{
  remove(psm_536660)
  remove(psm_17670)
  ##536660
  psm_563516 %>% 
    filter(
      artefact_arg == FALSE
    ) %>% 
    unnest(fragment_inten_df3) %>% 
    arrange(
      index,ion_num2
    ) %>% 
    mutate(
      frag_index = seq(0,nrow(.)-1)
    ) %>% 
    nest(frag_df= c(frag_index,ion_num2,b_z1,b_z2,y_z1,y_z2)) %>% 
    mutate(
      frag_start_idx = map_dbl(frag_df,\(x) min(x$frag_index))
    ) %>% 
    mutate(
      frag_stop_idx = map_dbl(frag_df,\(x) max(x$frag_index)+1)
    ) %>% 
    dplyr::rename(
      spec_idx = index,
      irt = `RT in min`,
      sequence = strip_seq
    ) -> psm_536653
  
  ###17670>> only artefact R
  psm_563516 %>% 
    filter(
      artefact_arg == TRUE
    ) %>% 
    unnest(fragment_inten_df3) %>% 
    arrange(
      index,ion_num2
    ) %>% 
    mutate(
      frag_index = seq(0,nrow(.)-1)
    ) %>% 
    nest(frag_df= c(frag_index,ion_num2,b_z1,b_z2,y_z1,y_z2)) %>% 
    mutate(
      frag_start_idx = map_dbl(frag_df,\(x) min(x$frag_index))
    ) %>% 
    mutate(
      frag_stop_idx = map_dbl(frag_df,\(x) max(x$frag_index)+1)
    ) %>% 
    dplyr::rename(
      spec_idx = index,
      irt = `RT in min`,
      sequence = strip_seq
    ) -> psm_17669
  
  psm_536653 %>% 
    select(
      sequence,charge,spec_idx,mods,mod_sites,nAA,frag_start_idx,frag_stop_idx,nce,instrument
    ) %>% 
    write_tsv(
      "C:/Users/admin/alphapeptdeep/project/folder/alphapeptdeep/precursor_df_536653.tsv"
    )
  
  psm_536653 %>% 
    select(frag_df) %>% 
    unnest(frag_df) %>% 
    select(
      b_z1,b_z2,y_z1,y_z2
    ) %>% 
    mutate(
      across(everything(),\(x) replace_na(x,0))
    ) %>% 
    write_tsv(
      "C:/Users/admin/alphapeptdeep/project/folder/alphapeptdeep/precursor_df_536653_fragdf.tsv"
    )
  
  read_rds(
    "./export/psm_17669.rds"
  ) -> psm_17669
  
  psm_17669 %>% 
    select(
      sequence,charge,spec_idx,mods,mod_sites,nAA,frag_start_idx,frag_stop_idx,nce,instrument
    ) %>% 
    write_tsv(
      "C:/Users/admin/alphapeptdeep/project/folder/alphapeptdeep/precursor_df_17669.tsv"
    )
  
  psm_17669 %>% 
    select(frag_df) %>% 
    unnest(frag_df) %>% 
    select(
      b_z1,b_z2,y_z1,y_z2
    ) %>% 
    mutate(
      across(everything(),\(x) replace_na(x,0))
    ) %>% 
    write_tsv(
      "C:/Users/admin/alphapeptdeep/project/folder/alphapeptdeep/precursor_df_17669_fragdf.tsv"
    )
  
  
  psm_17669 %>% 
    filter(
      enzyme == "trypsin"
    ) %>% 
    select(
      sequence,charge,spec_idx,mods,mod_sites,nAA,frag_start_idx,frag_stop_idx,nce,instrument
    ) %>% 
    write_tsv(
      "C:/Users/admin/alphapeptdeep/project/folder/alphapeptdeep/precursor_df_9569_trypsin.tsv"
    )
  
  psm_17669 %>% 
    filter(
      enzyme == "trypsin"
    ) %>% 
    select(frag_df) %>% 
    unnest(frag_df) %>% 
    select(
      b_z1,b_z2,y_z1,y_z2
    ) %>% 
    mutate(
      across(everything(),\(x) replace_na(x,0))
    ) %>% 
    write_tsv(
      "C:/Users/admin/alphapeptdeep/project/folder/alphapeptdeep/precursor_df_9569_trypsin_fragdf.tsv"
    )
  
  psm_17669 %>% 
    filter(
      enzyme == "chymotrypsin"
    ) %>% 
    unnest(frag_df) %>% 
    arrange(
      frag_index,ion_num2
    ) %>% 
    mutate(
      frag_index = seq(0,nrow(.)-1)
    ) %>% 
    nest(frag_df= c(frag_index,ion_num2,b_z1,b_z2,y_z1,y_z2)) %>% 
    mutate(
      frag_start_idx = map_dbl(frag_df,\(x) min(x$frag_index))
    ) %>% 
    mutate(
      frag_stop_idx = map_dbl(frag_df,\(x) max(x$frag_index)+1)
    )-> psm_17669_chymo
  
  
  psm_17669_chymo %>% 
    select(
      sequence,charge,spec_idx,mods,mod_sites,nAA,frag_start_idx,frag_stop_idx,nce,instrument
    ) %>% 
    write_tsv(
      "C:/Users/admin/alphapeptdeep/project/folder/alphapeptdeep/precursor_df_8100_chymotrypsin.tsv"
    )
  
  psm_17669_chymo %>%
    select(frag_df) %>% 
    unnest(frag_df) %>% 
    select(
      b_z1,b_z2,y_z1,y_z2
    ) %>% 
    mutate(
      across(everything(),\(x) replace_na(x,0))
    ) %>% 
    write_tsv(
      "C:/Users/admin/alphapeptdeep/project/folder/alphapeptdeep/precursor_df_8100_chymotrypsin_fragdf.tsv"
    )
  
  
}
#precursor_df_9569_trypsin.tsv
#precursor_df_9569_trypsin_fragdf.tsv
#precursor_df_8100_chymotrypsin.tsv
#precursor_df_8100_chymotrypsin_fragdf.tsv

##############################

##############################


###Chpater1
###1mg hct116 arg search non tandem 
{
  readxl::read_xlsx(
    "./chapter1data/160929_iNrich_160926_HPH_non_tandem.xlsx",
    sheet="Sheet1"
  ) -> data1
  
  data1 %>% 
    summarize(
      .by = `p1`,
      count = n()
    ) %>% 
    mutate(
      count = count/sum(count)
    ) %>% 
    arrange(count)
  
  data1 %>% 
    summarize(
      .by = `p1'`,
      count = n()
    ) %>% 
    mutate(
      count = count/sum(count)
    ) %>% 
    arrange(count)
  
  data1 %>% 
    summarize(
      .by = `p1'`,
      count = n()
    ) %>% 
    mutate(
      `p1'` = fct(`p1'`,levels = LETTERS)
    ) %>% 
    ggplot(
      aes(x= `p1'`, y = count)
    )+
    geom_bar(stat="identity")+
    labs(
      x="P1 residue",
      y="Count"
    )+
    theme_classic2()
  ggsave("./figures/1_01_p1prime.png",width = 8,height = 8, units = "cm", dpi = 600)
  data1 %>% 
    summarize(
      .by = p1,
      count = n()
    ) %>% 
    mutate(
      p1 = fct(p1,levels = c(LETTERS,"-"))
    ) %>% 
    ggplot(
      aes(x= p1, y = count)
    )+
    geom_bar(stat="identity")+
    labs(
      x="P1 residue",
      y="Count"
    )+
    theme_classic2()
  ggsave("./figures/1_01_p1.png",width = 8,height = 8, units = "cm", dpi = 600)
  
  ggseqlogo(
    data=readxl::read_xlsx(
      "./chapter1data/160929_iNrich_160926_HPH_non_tandem.xlsx",
      sheet="Sheet1"
    )$logo
  )+
    scale_x_discrete(
      labels = c("test1","test2")
    )+
    theme(
      legend.position = "none"
    )
  ggsave("./figures/1_01_logo.png",width = 8,height = 8, units = "cm", dpi = 600)
  
  data1 %>% 
    filter(
      `p1'` %in% c("E","D")
    ) %>% 
    summarize(
      .by = p1,
      count = n()
    ) %>% 
    mutate(
      perc = count/sum(count)
    )
  
  data1 %>% 
    filter(
      `p1'` %in% c("E","D")
    ) %>% 
    summarize(
      .by = p1,
      count = n()
    ) %>% 
    mutate(
      p1 = fct(p1,levels = c(LETTERS,"-"))
    ) %>% 
    ggplot(
      aes(x= p1, y = count)
    )+
    geom_bar(stat="identity")+
    labs(
      x="P1 residue",
      y="Count"
    )+
    theme_classic2()
  ggsave("./figures/1_01_p1_ed.png",width = 8,height = 8, units = "cm", dpi = 600)
}
#@ data1

###1mg hct116 arg search tandem
{
  ###data2 = arg peptide sheet
  readxl::read_xlsx(
    "./chapter1data/160929_iNrich_160926_HPH_F.xlsx",
    sheet="PeptideGroups"
  ) -> data2
  
  data2 %>% 
    mutate(
      prot_position = str_split_i(`Positions in Master Proteins`,";",1)
    ) %>%
    mutate(
      prot_pos_master = paste0(str_split_i(prot_position," ",1),"_",str_extract(prot_position,"(?<=\\[)[0-9]+"))
    ) %>%
    mutate(
      prot_master = str_split_i(prot_pos_master,"_",1)
    ) %>% 
    mutate(
      prot_position = str_split_i(prot_pos_master,"_",2)
    ) -> data2
  
  data2 %>% 
    left_join(
      dbtable, 
      join_by(prot_master==accession)
    ) %>% 
    mutate(
      sequence = paste0("XXXXX",sequence)
    ) %>% 
    rename(
      prot_seq = sequence
    ) %>% 
    mutate(
      prot_position = as.integer(prot_position)
    ) %>% 
    mutate(
      p5p1 = str_sub(prot_seq,prot_position,prot_position+4)
    ) %>% 
    mutate(
      p5p5prime = paste0(
        p5p1,
        str_sub(Sequence,1,5)
      )
    ) %>% 
    mutate(
      p2p1 = str_sub(prot_seq,prot_position+3,prot_position+4)
    ) %>% 
    mutate(
      p1p1 = str_sub(prot_seq,prot_position+4,prot_position+4)
    ) %>% 
    mutate(
      mass_ambiguity_passed = if_else(
        p2p1 == "GV" | p2p1 == "VG" | p1p1 =="R", FALSE,TRUE
      )
    )-> data2
  
  data2 %>% 
    mutate(
      p3p3 = str_sub(prot_seq,prot_position+2,prot_position+2)
    ) %>% 
    mutate(
      p2p2 = str_sub(prot_seq,prot_position+3,prot_position+3)
    ) %>% 
    mutate(
      index = 1:nrow(.)
    )-> data2
  
  data2 %>% 
    filter(
      p3p3 == "R"
    )
  data2 %>% 
    filter(
      p2p2 == "R"
    )
  
  data2 %>% 
    filter(
      mass_ambiguity_passed == TRUE
    ) %>% 
    pull(p5p5prime) %>% 
    ggseqlogo(
    )+
    theme(
      legend.position = "none"
    )
  ggsave("./figures/1_01_logo_ednq_p5p5prime_wo_massambig.png",width = 8,height = 8, units = "cm", dpi = 600)
  
  data2 %>% 
    summarize(
      .by = `p1`,
      count = n()
    ) %>% 
    mutate(
      count = count/sum(count)
    ) %>% 
    arrange(count)
  
  data2 %>% 
    filter(
      mass_ambiguity_passed == TRUE
    ) %>%
    select(
      gene,prot_pos_master
    ) %>% view()
  
  ##daqlogo analysis
  {
    logo_proteome <- prepareProteomeByFTP(source = NULL, species = "Homo sapiens", 
                                          fastaFile="UP000005640_9606_2301.fasta")
    formatSequence(seq = data2 %>% 
                     filter(
                       mass_ambiguity_passed == TRUE
                     ) %>% 
                     pull(p5p5prime),
                   proteome = logo_proteome, upstreamOffset = 5,
                   downstreamOffset = 5) -> dag_seq
    
    
    bg_fisher <- buildBackgroundModel(dag_seq, background = "wholeProteome", 
                                      proteome = logo_proteome, testType = "fisher")
    bg_ztest <- buildBackgroundModel(dag_seq, background = "wholeProteome", 
                                     proteome = logo_proteome, testType = "ztest")
    
    png("./figures/1_01_dag_heatmap.png", width = 10, height = 10, units = "cm", res = 600)
    dagHeatmap(
      testDAU(dag_seq, dagBackground = bg_ztest),
      type="diff",
      labels_col = c("P5","P4","P3","P2","P1","P1'","P2'","P3'","P4'","P5'"),
      display_numbers = TRUE
    )
    dev.off()
    }
}

#####Tandem Search Result of MOCK, MG132, MGTG
{
  read_tsv(
    "F:/PDdata/newargnrich/200821_HELA_MGTG_iNrich_bRP_F_PSMs.txt"
  ) -> psm_1015276
  
  psm_1015276 %>% 
    mutate(
      nmod = str_extract(Modifications,"N\\-Term\\([^\\;]+"),
      P1 = str_extract(`Annotated Sequence`,"(?<=^\\[)[A-Z]+"),
      P1prime = str_extract(Sequence,"^[A-Z]"),
      file_id_series = str_extract(`File ID`,"^F[0-9]+"),
      enzyme = case_when(
        file_id_series %in% c("F1","F2","F5","F6","F7","F12") ~ "trypsin",
        file_id_series %in% c("F3","F4","F8","F9","F10","F11") ~ "chymotrypsin"
      )
    ) -> psm_1015276
  
  psm_1015276 %>% distinct(nmod)
  
  psm_1015276 %>% 
    summarize(
      .by = c(nmod,P1prime),
      count = n()
    ) %>% 
    ggplot(
      aes(y=P1prime,x=nmod)
    )+
    geom_text(aes(label = count))+
    theme_bw()
  
  psm_1015276 %>% 
    summarize(
      .by = c(nmod,enzyme,P1prime),
      count = n()
    ) %>% print(n=Inf)
  
  psm_1015276 %>%
    filter(P1prime == "R" & nmod == "N-Term(Acetyl:2H(3))")
}

### peptide analysis of THREE
{
  summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
    select(
      Sequence,p2gv_check,p1r_check
    ) %>% 
    pivot_longer(
      p2gv_check:p1r_check
    ) %>% 
    filter(value == 1) %>% 
    distinct(Sequence)
  
  summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
    select(
      Sequence,p2gv_check,p1r_check
    ) %>% 
    pivot_longer(
      p2gv_check:p1r_check
    ) %>% 
    filter(value == 1) %>% 
    distinct(Sequence)
  
  summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
    mutate(
      p5p5prime = paste0(p5p1,str_sub(Sequence,start=1,end=5))
    ) -> summaryfile_psm_for_testing3_pcc_rt_hfsm2
  
  summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered %>% 
    mutate(
      p5p5prime = paste0(p5p1,str_sub(Sequence,start=1,end=5))
    ) -> summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered
  
  
  summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
    distinct(Sequence,.keep_all = TRUE) %>% #339
    pull(p5p5prime) %>% 
    ggseqlogo(
    )+
    theme(
      legend.position = "none"
    )
  ggsave("./figures/1_02_logo1_ednq_p5p5prime_all.png",width = 8,height = 8, units = "cm", dpi = 600)
  
  summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered %>% 
    distinct(Sequence,.keep_all = TRUE) %>% #339
    pull(p5p5prime) %>% 
    ggseqlogo(
    )+
    theme(
      legend.position = "none"
    )
  ggsave("./figures/1_02_logo2_ednq_p5p5prime_p1r_p2gv_removed.png",width = 8,height = 8, units = "cm", dpi = 600)
  
  ##daqlogo analysis
  {
    logo_proteome <- prepareProteomeByFTP(source = NULL, species = "Homo sapiens", 
                                          fastaFile="UP000005640_9606_2301.fasta")
    formatSequence(seq = summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
                     distinct(Sequence,.keep_all = TRUE) %>% pull(p5p1),
                   proteome = logo_proteome, upstreamOffset = 5,
                   downstreamOffset = 1) -> dag_seq2
    
    
    bg_fisher2 <- buildBackgroundModel(dag_seq2, background = "wholeProteome", 
                                       proteome = logo_proteome, testType = "fisher")
    #bg_ztest2 <- buildBackgroundModel(dag_seq2, background = "wholeProteome", 
                                      #proteome = logo_proteome, testType = "ztest")
    
    png("./figures/1_02_dag_heatmap_all.png", width = 10, height = 10, units = "cm", res = 600)
    dagHeatmap(
      testDAU(dag_seq2, dagBackground = bg_fisher2),
      type="diff",
      labels_col = c("P5","P4","P3","P2","P1"),
      display_numbers = TRUE
    )
    dev.off()
  }
  
  ##daqlogo analysis
  {
    logo_proteome <- prepareProteomeByFTP(source = NULL, species = "Homo sapiens", 
                                          fastaFile="UP000005640_9606_2301.fasta")
    formatSequence(seq = summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered %>% 
                     distinct(Sequence,.keep_all = TRUE) %>% pull(p5p1),
                   proteome = logo_proteome, upstreamOffset = 5,
                   downstreamOffset = 1) -> dag_seq3
    
    
    bg_fisher3 <- buildBackgroundModel(dag_seq3, background = "wholeProteome", 
                                       proteome = logo_proteome, testType = "fisher")
    #bg_ztest2 <- buildBackgroundModel(dag_seq2, background = "wholeProteome", 
    #proteome = logo_proteome, testType = "ztest")
    
    png("./figures/1_02_dag_heatmap_p2gv_removed.png", width = 10, height = 10, units = "cm", res = 600)
    dagHeatmap(
      testDAU(dag_seq3, dagBackground = bg_fisher3),
      type="diff",
      labels_col = c("P5","P4","P3","P2","P1"),
      display_numbers = TRUE
    )
    dev.off()
  }
}
#1_02_logo1_ednq_p5p5prime_all.png
#1_02_logo2_ednq_p5p5prime_p1r_p2gv_removed.png
#1_02_dag_heatmap_all.png

#logo analysis of P1R
{
  summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
    distinct(Sequence,.keep_all = TRUE) %>% #339
    filter(p1r_check==1) %>% 
    filter(enzyme =="trypsin") %>% 
    pull(p5p1) %>% 
    ggseqlogo(
    )+
    theme(
      legend.position = "none"
    )->fig_logo
  
  fig_logo$scales$scales[[1]]$labels <- c("P5","P4","P3","P2","P1")
  
  fig_logo
  
  ggsave("./figures/1_02_logo2_p1r.png",width = 8,height = 8, units = "cm", dpi = 600)
  
  summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
    distinct(Sequence,.keep_all = TRUE) %>% #339
    filter(p2gv_check==1) %>% 
    filter(enzyme =="trypsin") %>% 
    pull(p5p1) -> seqs_logo
    ggplot(
    )+
    geom_logo(
      data = seqs_logo
    )+
    theme_logo(
    )+
    theme(
      legend.position = "none"
    ) ->fig_logo
    
    fig_logo$scales$scales[[1]]$labels <- c("P5","P4","P3","P2","P1")
    
    fig_logo
  ggsave("./figures/1_02_logo2_p2gv.png",width = 8,height = 8, units = "cm", dpi = 600)
}
#1_02_logo2_p1r.png

file_index_arg %>% distinct(Sequence)

###
####Chapter 2 TrinityNmod Method Development
{
  ##arginylome >> tryp vs. chymo
  {
    summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
      summarize(
        .by = enzyme,
        count =n()
      )
    
    summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
      filter(PCC>=0.9) %>% 
      summarize(
        .by = enzyme,
        count =n()
      )
    
    summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
      mutate(
        enzyme = fct(enzyme,levels = c("trypsin","chymotrypsin"))
      ) %>% 
      ggplot(
        aes(x=enzyme,y=PCC,fill = enzyme)
      )+
      geom_boxplot()+
      scale_x_discrete(
        labels = c("Trypsin","Chymotrypsin")
      )+
      labs(
        x="Protease"
      )+
      theme_classic2()+
      theme(
        legend.position = "none",
        axis.text = element_text(color="black")
      )
    ggsave("./figures/1_pcc_arginylome_protease_result.png",width=10,height=10,dpi = 600,units="cm")
  }
  #1_pcc_arginylome_protease_result.png
  
  ##missed cleavage peptides count
  {
    file_index %>% 
      distinct(Sequence)
  }
  
  ##chapter missed cleavage peptides
  ###psm_563516
  {
    readRDS("230814_file_index_trainingset2_563516.rds")->psm_563524
    
    read_tsv(
      "F:/PDdata/newargnrich/200821_trainingset-(1)_PSMs.txt"
    ) -> psm_563524_table
    
    psm_563524_table %>% 
      mutate(
        nmod = str_extract(Modifications,"N\\-Term\\([^\\;]+"),
        P1 = str_extract(`Annotated Sequence`,"(?<=^\\[)[A-Z]+"),
        P1prime = str_extract(Sequence,"^[A-Z]"),
        file_id_series = str_extract(`File ID`,"^F[0-9]+"),
        enzyme = case_when(
          file_id_series %in% c("F1","F2","F5","F6","F7","F12") ~ "trypsin",
          file_id_series %in% c("F3","F4","F8","F9","F10","F11") ~ "chymotrypsin"
        )
      ) -> psm_563524_table
    
    psm_563524_table %>% distinct(
      `Annotated Sequence`,XCorr,`First Scan`
    )###check!! 563524 -> 563516
    
    psm_563524 %>% distinct(
      Sequence,XCorr, scan_num
    )
    
    psm_563524 %>%
      left_join(
        psm_563524_table,
        by = join_by(Sequence == `Annotated Sequence`, XCorr==XCorr,scan_num == `First Scan`,rt==`RT in min`),
        suffix = c("",".y")
      ) %>%
      dplyr::select(-contains(".y")) -> psm_563516
    
    remove(psm_563524_table)
    remove(psm_563524)
    
    
    
    
    psm_563516 %>% mutate(
      bion_int = map_dbl(peak_table_alpha_ms2, \(x) filter(x,series == "b") %>% .$i %>% sum()),
      yion_int = map_dbl(peak_table_alpha_ms2, \(x) filter(x,series == "y") %>% .$i %>% sum()),
      bion_int_perc = bion_int/(bion_int+yion_int),
      yion_int_perc = yion_int/(bion_int+yion_int),
    )-> psm_563516
    
    
    
    psm_563516 %>% 
      filter(
        P1prime != "R" & nmod =="Acetyl:2H(3" & enzyme == "trypsin"
      ) %>% 
      select(bion_int_perc,yion_int_perc) %>% 
      pivot_longer(everything()) %>% 
      drop_na() %>% 
      summarise(
        .by = name,
        grp_mean = mean(value)
      )
    
    readRDS(
      file = "./export/psm_563516.rds"
    )->psm_563516
    
    ##9570
    psm_563516 %>% filter(
      P1prime == "R" & nmod =="Acetyl:2H(3" & enzyme == "trypsin"
    ) %>% 
      select(bion_int_perc,yion_int_perc) %>% 
      pivot_longer(everything()) %>% 
      ggplot(
      )+ 
      geom_histogram(
        aes(x=value,color = name),
        fill = "white",alpha=0.5,binwidth = 0.01
      )+
      geom_vline(
        data = psm_563516 %>% 
          filter(
            P1prime == "R" & nmod =="Acetyl:2H(3" & enzyme == "trypsin"
          ) %>% 
          select(bion_int_perc,yion_int_perc) %>% 
          pivot_longer(everything()) %>% 
          drop_na() %>% 
          summarise(
            .by = name,
            grp_mean = mean(value)
          ),
        aes(xintercept = grp_mean),color="#4DAF4A", linetype="dashed"
      )+
      scale_color_manual(values = c("#377EB8","#E41A1C"))+
      scale_x_continuous(labels = scales::label_percent())+
      labs(x="Area Percentage",y="PSM count")+
      theme_classic2()+
      theme(legend.position = "none")+
      facet_grid(name ~.)
    ggsave("./figures/2_02_ion_series_9570.png",width =4, height = 4,dpi = 600)
    
    psm_563516 %>%
      filter(
        P1prime == "R" & nmod =="Acetyl:2H(3" & enzyme == "trypsin"
      ) %>% 
      select(bion_int_perc,yion_int_perc) %>% 
      pivot_longer(everything()) %>% 
      summarize(
        .by = name,
        mean = mean(value)
      )
    
    RColorBrewer::brewer.pal(3,"Set1")
    
    ##156125
    psm_563516 %>% filter(
      P1prime != "R" & nmod =="Acetyl:2H(3" & enzyme == "trypsin"
    ) %>% 
      select(bion_int_perc,yion_int_perc) %>% 
      pivot_longer(everything()) %>% 
      ggplot(
      )+ 
      geom_histogram(
        aes(x=value,color = name),
        fill = "white",alpha=0.5,binwidth = 0.01
      )+
      geom_vline(
        data = psm_563516 %>% 
          filter(
            P1prime != "R" & nmod =="Acetyl:2H(3" & enzyme == "trypsin"
          ) %>% 
          select(bion_int_perc,yion_int_perc) %>% 
          pivot_longer(everything()) %>% 
          drop_na() %>% 
          summarise(
            .by = name,
            grp_mean = mean(value)
          ),
        aes(xintercept = grp_mean),color="#4DAF4A", linetype="dashed"
      )+
      scale_color_manual(values = c("#377EB8","#E41A1C"))+
      scale_x_continuous(labels = scales::label_percent())+
      labs(x="Area Percentage",y="PSM count")+
      theme_classic2()+
      theme(legend.position = "none")+
      facet_grid(name ~.)
    ggsave("./figures/2_02_ion_series_156125.png",width =4, height = 4,dpi = 600)
    
    psm_563516 %>%
      filter(
        P1prime != "R" & nmod =="Acetyl:2H(3" & enzyme == "trypsin"
      ) %>% 
      select(bion_int_perc,yion_int_perc) %>% 
      pivot_longer(everything()) %>% 
      drop_na() %>% 
      summarize(
        .by = name,
        mean = mean(value),
        count = n(),
        sd = sd()
      )
    
    psm_563516 %>%
      filter(
        P1prime != "R" & nmod =="Acetyl:2H(3" & enzyme == "trypsin"
      ) %>% 
      select(bion_int_perc,yion_int_perc) %>% 
      drop_na()
    t.test(
      psm_563516 %>%
        filter(
          P1prime != "R" & nmod =="Acetyl:2H(3" & enzyme == "trypsin"
        ) %>% 
        select(bion_int_perc,yion_int_perc) %>% 
        drop_na() %>% .$bion_int_perc,
      psm_563516 %>%
        filter(
          P1prime != "R" & nmod =="Acetyl:2H(3" & enzyme == "trypsin"
        ) %>% 
        select(bion_int_perc,yion_int_perc) %>% 
        drop_na() %>% .$yion_int_perc,
      var.equal = TRUE
    ) ##2.2e-16
    
    t.test(
      psm_563516 %>%
        filter(
          P1prime == "R" & nmod =="Acetyl:2H(3" & enzyme == "trypsin"
        ) %>% 
        select(bion_int_perc,yion_int_perc) %>% 
        drop_na() %>% .$bion_int_perc,
      psm_563516 %>%
        filter(
          P1prime == "R" & nmod =="Acetyl:2H(3" & enzyme == "trypsin"
        ) %>% 
        select(bion_int_perc,yion_int_perc) %>% 
        drop_na() %>% .$yion_int_perc,
      var.equal = TRUE
    ) ##2.2e-16
    
    
    
    ###low RT distribution
    psm_563516 %>% filter(
      nmod =="Acetyl:2H(3" & enzyme == "trypsin"
    ) %>% 
      select(P1prime,rt_norm) %>% 
      mutate(
        p1prime = if_else(P1prime =="R","R","Others")
      ) %>% 
      select(p1prime,rt_norm) %>% 
      ggplot()+
      geom_density(aes(x=rt_norm,fill = p1prime),alpha = 0.5)+
      labs(x="Normalized RT",y="Density")+
      theme_classic2()+
      theme(legend.position = "none")
    ggsave("./figures/2_02_rt.png",width =4, height = 4,dpi = 600)
    
    
    t.test(
      psm_563516 %>% filter(
        nmod =="Acetyl:2H(3" & enzyme == "trypsin"
      ) %>% 
        select(P1prime,rt_norm) %>% 
        mutate(
          p1prime = if_else(P1prime =="R","R","Others")
        ) %>% 
        select(p1prime,rt_norm) %>% filter(p1prime == "Others") %>% .$rt_norm,
      psm_563516 %>% filter(
        nmod =="Acetyl:2H(3" & enzyme == "trypsin"
      ) %>% 
        select(P1prime,rt_norm) %>% 
        mutate(
          p1prime = if_else(P1prime =="R","R","Others")
        ) %>% 
        select(p1prime,rt_norm) %>% filter(p1prime == "R") %>% .$rt_norm,
      var.equal = TRUE
    )
    
    psm_563516 %>% filter(
      nmod =="Acetyl:2H(3" & enzyme == "trypsin"
    ) %>% 
      select(P1prime,rt_norm) %>% 
      mutate(
        p1prime = if_else(P1prime =="R","R","Others")
      ) %>% 
      select(p1prime,rt_norm) %>% filter(p1prime == "Others") %>% .$rt_norm
  }
  #2_02_ion_series_9570.png
  #2_02_ion_series_156125.png
  #2_02_rt.png
  #psm_563516
  #psm_1015276
  
  ##chapter ms2 prediction model
  ###output_prediction_comparison ##9569
  ##output_prediction_comparison_limit
  {
    ##read alphapeptdeep prediction result of trypsin 9569
    read_tsv(
      "C:/Users/admin/alphapeptdeep/project/folder/alphapeptdeep/inrich/output/prediction_17669_output_metrics_pretrained_model_trypsin_231206.tsv"
    ) %>% 
      mutate(
        model = "pretrained"
      ) %>% 
      bind_rows(
        read_tsv(
          "C:/Users/admin/alphapeptdeep/project/folder/alphapeptdeep/inrich/output/prediction_17669_output_metrics_built_model_trypsin_231206.tsv"
        ) %>% 
          mutate(
            model = "built"
          )
      ) %>% 
      bind_rows(
        read_tsv(
          "C:/Users/admin/alphapeptdeep/project/folder/alphapeptdeep/inrich/output/prediction_17669_output_metrics_finetuned_model_trypsin_231206.tsv"
        ) %>% 
          mutate(
            model = "finetuned"
          )
      ) %>% 
      bind_rows(
        read_tsv(
          "C:/Users/admin/alphapeptdeep/project/folder/alphapeptdeep/inrich/output/prediction_9569_output_metrics_311547woarg_model_trypsin_231213.tsv"
        ) %>% 
          mutate(
            model = "finetunedwoarg"
          )
      )-> output_prediction_comparison_trypsin
    
    ##read alphapeptdeep prediction result of chymotrypsin 8100
    read_tsv(
      "C:/Users/admin/alphapeptdeep/project/folder/alphapeptdeep/inrich/output/prediction_17669_output_metrics_pretrained_model_chymotrypsin_231206.tsv"
    ) %>% 
      mutate(
        model = "pretrained"
      ) %>% 
      bind_rows(
        read_tsv(
          "C:/Users/admin/alphapeptdeep/project/folder/alphapeptdeep/inrich/output/prediction_17669_output_metrics_built_model_chymotrypsin_231206.tsv"
        ) %>% 
          mutate(
            model = "built"
          )
      ) %>% 
      bind_rows(
        read_tsv(
          "C:/Users/admin/alphapeptdeep/project/folder/alphapeptdeep/inrich/output/prediction_17669_output_metrics_finetuned_model_chymotrypsin_231206.tsv"
        ) %>% 
          mutate(
            model = "finetuned"
          )
      ) -> output_prediction_comparison_chymotrypsin
      
    
    RColorBrewer::display.brewer.all()
    RColorBrewer::brewer.pal(3,"Set1")
    
    output_prediction_comparison_trypsin %>% write_tsv(
      "fig2_1.tsv"
    )
    output_prediction_comparison_trypsin %>%
      mutate(
        model = fct(model,levels=c("","built","finetuned","finetunedwoarg"))
      ) %>% 
      ggplot(aes(x=model,y=PCC))+
      geom_violin(aes(fill = model),linewidth = 0.1,scale = "width")+
      labs(y="PCC")+
      geom_boxplot(width = 0.1,outlier.shape = NA,linewidth = 0.1)+
      scale_x_discrete(labels =c("Pretrained","Built","Fine-tuned","Fine-tuned w/o\nR-starting peptides"))+
      scale_fill_manual(
        values = c("#4DAF4A","#377EB8","#E41A1C","#984EA3")
      )+
      scale_y_continuous(
        limits=c(0,1)
      )+
      theme_classic2()+
      theme(
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = "black",size = 7),
        axis.text.x = element_blank(),
        axis.text.y = element_text(color = "black",size = 7)
      )
    ggsave("./figures/2_03_finetuning1.png",width =4.5, height = 4,dpi = 600,units = "cm")
    
    output_prediction_comparison_trypsin %>% 
      mutate(
        model = fct(model,levels=c("pretrained","built","finetuned","finetunedwoarg"))
      ) %>% 
      filter(PCC >= 0.9) %>% 
      summarize(
        .by = model,
        count = n()
      ) %>% 
      mutate(
        perc = count/9569
      ) %>% 
      ggplot(aes(x=model,y=perc,fill=model))+
      geom_bar(stat = "identity")+
      geom_text(
        aes(label = scales::percent(perc,accuracy = 0.1)),vjust =-0.3,
        size =2
      )+
      labs(y="PCC \u2265 0.9")+
      scale_y_continuous(labels = scales::percent,breaks = seq(0,1,0.25),limits = c(0,1))+
      #scale_x_discrete(labels =c("Built","Pretrained","Fine-tuned","Fine-tuned w/o\nR-starting peptides"))+
      scale_fill_manual(
        values = c("#4DAF4A","#377EB8","#E41A1C","#984EA3")
      )+
      theme_classic2()+
      theme(
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = "black",size = 7),
        axis.text.x = element_blank(),
        axis.text.y = element_text(color = "black",size = 7)
      )
    ggsave("./figures/2_03_finetuning2.png",width =4.5, height = 4,dpi = 600,units = "cm")
  }
  #2_03_finetuning1.png
  #2_03_finetuning2.png
  
  ##chapter ms2 prediction model chymotrypsin
  ###output_prediction_comparison ##8100
  ##output_prediction_comparison_limit
  {
    ##read alphapeptdeep prediction result of chymotrypsin 8100
    read_tsv(
      "C:/Users/admin/alphapeptdeep/project/folder/alphapeptdeep/inrich/output/prediction_17669_output_metrics_pretrained_model_chymotrypsin_231206.tsv"
    ) %>% 
      mutate(
        model = "pretrained"
      ) %>% 
      bind_rows(
        read_tsv(
          "C:/Users/admin/alphapeptdeep/project/folder/alphapeptdeep/inrich/output/prediction_17669_output_metrics_built_model_chymotrypsin_231206.tsv"
        ) %>% 
          mutate(
            model = "built"
          )
      ) %>% 
      bind_rows(
        read_tsv(
          "C:/Users/admin/alphapeptdeep/project/folder/alphapeptdeep/inrich/output/prediction_17669_output_metrics_finetuned_model_chymotrypsin_231206.tsv"
        ) %>% 
          mutate(
            model = "finetuned"
          )
      ) -> output_prediction_comparison_chymotrypsin
    
    
    RColorBrewer::display.brewer.all()
    RColorBrewer::brewer.pal(3,"Set1")
    
    output_prediction_comparison_chymotrypsin %>%
      mutate(
        model = fct(model,levels=c("built","pretrained","finetuned"))
      ) %>% 
      ggplot(aes(x=model,y=PCC))+
      geom_violin(aes(fill = model),linewidth = 0.1,scale = "width")+
      labs(y="PCC")+
      geom_boxplot(width = 0.1,outlier.shape = NA,linewidth = 0.1)+
      scale_x_discrete(labels =c("Built","Pretrained","Fine-tuned"))+
      scale_fill_manual(
        values = c("#4DAF4A","#377EB8","#E41A1C")
      )+
      theme_classic2()+
      theme(
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black",size = 9,angle = 45,hjust=1),
        axis.text.y = element_text(color = "black",size = 9)
      )
    ggsave("./figures/2_03_finetuning3.png",width =5, height = 8,dpi = 600,units = "cm")
    
    output_prediction_comparison_chymotrypsin %>% 
      mutate(
        model = fct(model,levels=c("built","pretrained","finetuned"))
      ) %>% 
      filter(PCC >= 0.9) %>% 
      summarize(
        .by = model,
        count = n()
      ) %>% 
      mutate(
        perc = count/9569
      ) %>% 
      ggplot(aes(x=model,y=perc,fill=model))+
      geom_bar(stat = "identity")+
      geom_text(
        aes(label = scales::percent(perc,accuracy = 0.1)),vjust =-0.3,
        size =3
      )+
      labs(y="PCC \u2265 0.9")+
      scale_y_continuous(labels = scales::percent,breaks = seq(0,1,0.1),limits = c(0,1))+
      scale_x_discrete(labels =c("Built","Pretrained","Fine-tuned"))+
      scale_fill_manual(
        values = c("#4DAF4A","#377EB8","#E41A1C")
      )+
      theme_classic2()+
      theme(
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black",size = 9,angle = 45,hjust=1),
        axis.text.y = element_text(color = "black",size = 9)
      )
    ggsave("./figures/2_03_finetuning4.png",width =5, height = 8,dpi = 600,units = "cm")
  }
  #2_03_finetuning3.png
  #2_03_finetuning4.png
    
  ##deprecated!
  ####training limit
  {
    read_tsv(
      "C:/Users/admin/alphapeptdeep/project/folder/alphapeptdeep/inrich/output_metrics_describ_230917_250.tsv"
    ) %>% 
      select(PCC) %>% 
      mutate(training_model = "250") %>% 
      bind_rows(
        read_tsv(
          "C:/Users/admin/alphapeptdeep/project/folder/alphapeptdeep/inrich/output_metrics_describ_230917_500.tsv"
        ) %>% 
          select(PCC) %>% 
          mutate(training_model = "500"),
        read_tsv(
          "C:/Users/admin/alphapeptdeep/project/folder/alphapeptdeep/inrich/output_metrics_describ_230917_1000.tsv"
        ) %>% 
          select(PCC) %>% 
          mutate(training_model = "1000"),
        read_tsv(
          "C:/Users/admin/alphapeptdeep/project/folder/alphapeptdeep/inrich/output_metrics_describ_230917_5000.tsv"
        ) %>% 
          select(PCC) %>% 
          mutate(training_model = "5000"),
        read_tsv(
          "C:/Users/admin/alphapeptdeep/project/folder/alphapeptdeep/inrich/output_metrics_describ_230917_10000.tsv"
        ) %>% 
          select(PCC) %>% 
          mutate(training_model = "10000"),
        read_tsv(
          "C:/Users/admin/alphapeptdeep/project/folder/alphapeptdeep/inrich/output_metrics_describ_230917_15000.tsv"
        ) %>% 
          select(PCC) %>% 
          mutate(training_model = "15000"),
        read_tsv(
          "C:/Users/admin/alphapeptdeep/project/folder/alphapeptdeep/inrich/output_metrics_describ_230917_17669.tsv"
        ) %>% 
          select(PCC) %>% 
          mutate(training_model = "17669"),
        read_tsv(
          "C:/Users/admin/alphapeptdeep/project/folder/alphapeptdeep/output_metrics_describ_training_230728.tsv"
        ) %>% 
          filter(str_detect(sequence,"^R")&str_detect(mod_sites,"^0")&str_detect(mods,"^Acetyl:2")) %>% 
          select(PCC) %>% 
          mutate(training_model = "563516")
      ) -> output_prediction_comparison_limit
    
    output_prediction_comparison_limit %>% 
      mutate(
        training_model = fct(training_model,levels=c("250","500","1000","5000","10000","15000","17669","563516"))
      ) %>% 
      ggplot(
        aes(x=training_model,y=PCC)
      )+
      geom_violin(aes(fill = training_model),linewidth = 0.1,scale = "area")+
      labs(y="PCC",x="Number of MS/MS scans into fine-tuning")+
      geom_boxplot(width = 0.05,outlier.shape = NA)+
      scale_y_continuous(limits = c(0.5,1))+
      theme_classic2()+
      theme(
        legend.position = "none",
        #axis.title.x = element_blank(),
        axis.text = element_text(color = "black")
      )
    ggsave("./figures/2_03_finetuning_limit_pcc.png",width =8, height = 4,dpi = 600)
    
    output_prediction_comparison_limit %>%
      mutate(
        training_model = fct(training_model,levels=c("250","500","1000","5000","10000","15000","17669","563516"))
      ) %>% 
      filter(PCC >= 0.9) %>% 
      summarize(
        .by = training_model,
        count = n(),
        count = count/17669
      ) %>% 
      ggplot(aes(x=training_model,y=count))+
      geom_bar(aes(fill = training_model),stat = "identity")+
      geom_text(aes(label = scales::percent(count)),vjust =-0.3)+
      labs(y="PCC",x="Number of MS/MS scans into fine-tuning")+
      scale_y_continuous(labels = scales::percent,breaks = seq(0,1,0.1))+
      theme_classic2()+
      theme(
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text = element_text(color = "black")
      )
    ggsave("./figures/2_03_finetuning_limit_pcc90.png",width =8, height = 4,dpi = 600)
  }
  #2_03_finetuning_limit_pcc90.png
  
  #### bion yion pcc
  {
    summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
      select(
        prot_pos_reposition,PCC,bion_pcc,yion_pcc
      ) %>% 
      mutate(
        group = case_when(
          PCC>=0.9 ~ fct("1",levels = c("1","2","3","4")),
          PCC<0.9 & PCC >=0.6 ~ fct("2",levels = c("1","2","3","4")),
          PCC<0.6 ~ fct("3",levels = c("1","2","3","4"))
        )
      ) %>% 
      summarise(
        .by = group,
        count = n()
      ) ## count!
    
    summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
      select(
        prot_pos_reposition,PCC,bion_pcc,yion_pcc
      ) %>% 
      pivot_longer(cols=c(bion_pcc,yion_pcc)) %>% 
      mutate(
        group = case_when(
          PCC>=0.9 ~ fct("1",levels = c("1","2","3","4")),
          PCC<0.9 & PCC >=0.6 ~ fct("2",levels = c("1","2","3","4")),
          PCC<0.6 ~ fct("3",levels = c("1","2","3","4"))
        )
      ) %>% 
      ggplot(
        aes(y=value,x=group,fill=name)
      )+
      geom_boxplot(
        size = 0.2,
        outlier.size = 0.2,
      )+
      annotate(
        "text",
        x=1:3,
        y = 1.08,
        label = paste("N = ",c(194,446,577)),
        size = 2
      )+
      scale_x_discrete(
        limits = c("3","2","1"),
        labels = c(expression("PCC"<"0.6", paste("0.6"<="","PCC",""<"0.9"),"0.9"<="PCC"))
      )+
      scale_fill_manual(labels = c("b-ion","y-ion"),values = c(wes_palette("Zissou1")[1],wes_palette("Zissou1")[5]))+
      labs(y="PCC")+
      theme_classic()+
      theme(
        legend.position = "none",
        legend.title = element_blank(),
        text= element_text(color = "black", size = 7),
        axis.text.y = element_text(color = "black", size = 7, hjust=1),
        axis.text.x = element_text(color = "black", size = 7,angle = 20, hjust =1),
        axis.title.x = element_blank(),
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm")
      )
    ggsave("./figures/2_pcc_ions.png",width=4.5,height=4,dpi = 600,units="cm")
  }
  #2_pcc_ions.png
  
  ###summary
  {
    summaryfile_psm_for_testing3_pcc_rt_hfsm %>% 
      summarize(
        .by = pcc_check,
        count = n()
      )
    
    summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
      mutate(
        histogram_color = case_when(
          PCC>=0.9 ~ fct("1",levels = c("1","2","3","4")),
          PCC<0.9 & PCC >=0.6 ~ fct("2",levels = c("1","2","3","4")),
          PCC<0.6 ~ fct("3",levels = c("1","2","3","4"))
        )
      ) %>% 
      gghistogram(
        x = "PCC", 
        add = "median", rug = TRUE,
        fill = "histogram_color", palette = c("#66C2A5", "#FC8D62","#8DA0CB"),
        binwidth = 0.025,
        #add_density = TRUE,
        position="stack",
        size = 0.1
      )+
      labs(x="PCC",y="Count")+
      scale_x_continuous(limits = c(-0.2,1))+
      theme(
        legend.position = "none",
        #axis.title.x = element_blank(),
        axis.title = element_text(color = "black",size = 7),
        #axis.text.x = element_blank(),
        axis.text = element_text(color = "black",size = 6),
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm")
      )
    ggsave("./figures/2_03_pcc_summary.png",width=4.5,height=4.5, units = "cm",dpi = 600)
    
    
    summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
      mutate(
        pcc_check = fct(as.character(pcc_check))
      ) %>% 
      gghistogram(
        x = "PCC", 
        add = "median", rug = TRUE,
        fill = "pcc_check", palette = c("#0868AC", "#A8DDB5"),
        binwidth = 0.025,
        position="stack"
      )+
      labs(x="PCC",y="Count")+
      scale_x_continuous(limits = c(-0.2,1))+
      theme(
        legend.position = "none"
      )
    ggsave("./figures/2_03_pcc_summary_fdr.png",width=8,height=8, units = "cm",dpi = 600)
  }
  #2_03_pcc_summary.png
  display.brewer.all()
  brewer.pal(5,"Set2")
  brewer.pal(9,"GnBu")
  
  ##2A >> 311547-9569, psm_251969_chymo - 8100 prediction data
  {
    readRDS(
      "./export/psm_311547_trypsin.rds"
    )->psm_311547_trypsin
    
    readRDS(
      "./export/psm_251969_chymo.rds"
    )->psm_251969_chymo
    
    #311547-9569
    {
      psm_311547_trypsin
      psm_311547_trypsin %>% 
        filter(artefact_arg != TRUE) %>% 
        unnest(fragment_inten_df3) %>% 
        arrange(
          index,ion_num2
        ) %>% 
        mutate(
          frag_index = seq(0,nrow(.)-1)
        ) %>% 
        nest(frag_df= c(frag_index,ion_num2,b_z1,b_z2,y_z1,y_z2)) %>% 
        mutate(
          frag_start_idx = map_dbl(frag_df,\(x) min(x$frag_index))
        ) %>% 
        mutate(
          frag_stop_idx = map_dbl(frag_df,\(x) max(x$frag_index)+1)
        ) %>% 
        dplyr::rename(
          spec_idx = index,
          irt = `RT in min`,
          sequence = strip_seq
        ) -> psm_299486_trypsin
      
      psm_299486_trypsin %>% 
        select(
          sequence,charge,spec_idx,mods,mod_sites,nAA,frag_start_idx,frag_stop_idx,nce,instrument
        ) %>% 
        write_tsv(
          "C:/Users/admin/alphapeptdeep/project/folder/alphapeptdeep/inrich/input/precursor_df_311547_wo_arg.tsv"
        )
      psm_299486_trypsin %>% 
        select(frag_df) %>% 
        unnest(frag_df) %>% 
        select(
          b_z1,b_z2,y_z1,y_z2
        ) %>% 
        mutate(
          across(everything(),\(x) replace_na(x,0))
        ) %>% 
        write_tsv(
          "C:/Users/admin/alphapeptdeep/project/folder/alphapeptdeep/inrich/input/precursor_df_311547_wo_arg_fragdf.tsv"
        )
        
    }
    
    #251969-8100
    {
      psm_251969_chymo %>% 
        filter(artefact_arg != TRUE) %>% 
        unnest(fragment_inten_df3) %>% 
        arrange(
          index,ion_num2
        ) %>% 
        mutate(
          frag_index = seq(0,nrow(.)-1)
        ) %>% 
        nest(frag_df= c(frag_index,ion_num2,b_z1,b_z2,y_z1,y_z2)) %>% 
        mutate(
          frag_start_idx = map_dbl(frag_df,\(x) min(x$frag_index))
        ) %>% 
        mutate(
          frag_stop_idx = map_dbl(frag_df,\(x) max(x$frag_index)+1)
        ) %>% 
        dplyr::rename(
          spec_idx = index,
          irt = `RT in min`,
          sequence = strip_seq
        ) -> psm_237167_chymotrypsin
      
      psm_237167_chymotrypsin %>% 
        select(
          sequence,charge,spec_idx,mods,mod_sites,nAA,frag_start_idx,frag_stop_idx,nce,instrument
        ) %>% 
        write_tsv(
          "C:/Users/admin/alphapeptdeep/project/folder/alphapeptdeep/inrich/input/precursor_df_251969_wo_arg.tsv"
        )
      psm_237167_chymotrypsin %>% 
        select(frag_df) %>% 
        unnest(frag_df) %>% 
        select(
          b_z1,b_z2,y_z1,y_z2
        ) %>% 
        mutate(
          across(everything(),\(x) replace_na(x,0))
        ) %>% 
        write_tsv(
          "C:/Users/admin/alphapeptdeep/project/folder/alphapeptdeep/inrich/input/precursor_df_251969_wo_arg_fragdf.tsv"
        )
      
    }
    ##save rds
    {
      write_rds(
        psm_299486_trypsin,
        "./export/psm_299486_trypsin.rds"
      )
      
      write_rds(
        psm_237167_chymotrypsin,
        "./export/psm_237167_chymotrypsin.rds"
      )
    }
    
    #alphapeptdeep!
    
    ##figure
    {
      
    }
  }
  
  ###Prediction vs. measured spectra
  ##ms2 comparison figure! @summaryfile_psm_for_testing3_pcc_rt_hfsm2
  {
    ##figure testbed
    {
      summaryfile_psm_for_testing3_pcc_rt_hfsm2$index[1] -> target_scan #1 >> j
      summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
        dplyr::filter(
          index == target_scan
        ) -> spectrum_file_plot
      
      spectrum_file_plot$sequence_for_fragspec -> title1
      spectrum_file_plot$index -> subtitle1 ##index!
      round(spectrum_file_plot$PCC,3) -> title2 ##pcc
      spectrum_file_plot$`MHplus in Da` -> spec_plot_xmax
      spectrum_file_plot$prot_pos_reposition -> title3
      spectrum_file_plot %>% select(peak_table_for_graph) %>%
        unnest(cols=peak_table_for_graph) %>% 
        ggplot( 
          aes(x=`m/z`,y=frag, label = ion_type)
        )+
        geom_bar(stat="identity",width = 0.2)+
        geom_col(aes(fill = series),width = spec_plot_xmax/1000)+
        geom_text(aes(y=if_else(frag<0,frag-0.08,frag),label=ion_type,color= series),angle=90,hjust=-0.05,vjust=0,size=2, na.rm = TRUE)+
        geom_hline(yintercept = 0, color = "black")+
        scale_fill_manual(
          values = c("b"=RColorBrewer::brewer.pal(n=9,"Set1")[2],"y"=RColorBrewer::brewer.pal(n=9,"Set1")[1]),
          na.value = NA
        )+
        scale_color_manual(
          values = c("b"=RColorBrewer::brewer.pal(n=9,"Set1")[2],"y"=RColorBrewer::brewer.pal(n=9,"Set1")[1]),
          na.value = NA
        )+
        scale_x_continuous(limits=c(0, spec_plot_xmax),expand=c(0,0))+ #limits=c(250,800),
        scale_y_continuous(breaks = c(1,0.5,0,-0.5,-1),labels = scales::percent(c(1,0.5,0,0.5,1)), limits=c(-1,1),expand=c(0,0))+
        theme_classic()+
        labs(x="m/z",y="Normalized Intensity", fill = "Ion Type", color = "Ion Type", title = paste0(title1,"  ",title3), subtitle = paste0("PCC = ",title2,"  index: ",subtitle1))+
        theme(
          legend.position = "none",
          plot.subtitle = element_text(size=10, color="black")
        )->mirror_plot1
      
      
      spectrum_file_plot %>% select(peak_table_for_graph) %>%
        unnest(cols=peak_table_for_graph) %>% 
        filter(!is.na(ms2error)) %>% 
        ggplot(
          aes(x=`m/z`,y=ms2error,color = series)
        )+
        geom_point(size=1)+
        geom_hline(yintercept = 0, color = "black")+
        scale_color_manual(
          values = c("b"=RColorBrewer::brewer.pal(n=9,"Set1")[2],"y"=RColorBrewer::brewer.pal(n=9,"Set1")[1])
        )+
        scale_x_continuous(
          limits=c(0, spec_plot_xmax),expand=c(0,0)
        )+
        scale_y_continuous(
          limits = c(-0.005,0.005),
          labels = \(x) x*1000
        )+
        labs(
          y="Mass error (mDa)"
        )+
        theme_bw()+
        theme(
          legend.position = "none"
        )->mirror_plot2
      ggarrange(
        mirror_plot1,mirror_plot2,
        ncol=1,
        heights = c(2,1),
        align = "v"
      )
      ggsave(paste0("./ms2figures/",j,".png"), width=20, height = 12, units = "cm", dpi = 600)
    }
    
    
    ###figure generation
    for (j in 1:nrow(summaryfile_psm_for_testing3_pcc_rt_hfsm2)){
      summaryfile_psm_for_testing3_pcc_rt_hfsm2$index[j] -> target_scan #1 >> j
      summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
        dplyr::filter(
          index == target_scan
        ) -> spectrum_file_plot
      
      spectrum_file_plot$sequence_for_fragspec -> title1
      spectrum_file_plot$index -> subtitle1 ##index!
      round(spectrum_file_plot$PCC,3) -> title2 ##pcc
      spectrum_file_plot$`MHplus in Da` -> spec_plot_xmax
      spectrum_file_plot$prot_pos_reposition -> title3
      #summary_filtered_PRM_table_final_DL_library_merged_spec_DLspec_crossed$pearson_correl_coeff[j] |> round(digits = 4) -> subtitle2
      #summary_filtered_PRM_table_final_DL_library_merged_spec_DLspec_crossed$cosine_similarity [j] |> round(digits = 4) -> subtitle3
      spectrum_file_plot %>% select(peak_table_for_graph) %>%
        unnest(cols=peak_table_for_graph) %>% 
        ggplot( 
          aes(x=`m/z`,y=frag, label = ion_type)
        )+
        geom_bar(stat="identity",width = 0.2)+
        geom_col(aes(fill = series),width = spec_plot_xmax/1000)+
        geom_text(aes(y=if_else(frag<0,frag-0.08,frag),label=ion_type,color= series),angle=90,hjust=-0.05,vjust=0,size=2, na.rm = TRUE)+
        # ggrepel::geom_text_repel(
        #   size=2,
        #   alpha = 0.2,
        #   color = "black",
        #   angle= 90,
        #   nudge_x = 10,
        #   #nudge_y = 1000,
        #   box.padding = 0.1,
        #   na.rm = TRUE,
        #   direction = "x",
        #   min.segment.length = 0.5,
      #   ylim = c(0.4,0.8)
      # )+
      geom_hline(yintercept = 0, color = "black")+
        scale_fill_manual(
          values = c("b"=RColorBrewer::brewer.pal(n=9,"Set1")[2],"y"=RColorBrewer::brewer.pal(n=9,"Set1")[1]),
          na.value = NA
        )+
        scale_color_manual(
          values = c("b"=RColorBrewer::brewer.pal(n=9,"Set1")[2],"y"=RColorBrewer::brewer.pal(n=9,"Set1")[1]),
          na.value = NA
        )+
        scale_x_continuous(limits=c(0, spec_plot_xmax),expand=c(0,0))+ #limits=c(250,800),
        scale_y_continuous(breaks = c(1,0.5,0,-0.5,-1),labels = scales::percent(c(1,0.5,0,0.5,1)), limits=c(-1.2,1.2),expand=c(0,0))+
        theme_classic()+
        labs(x="m/z",y="Normalized Intensity", fill = "Ion Type", color = "Ion Type", title = paste0(title1,"  ",title3), subtitle = paste0("PCC = ",title2,"  index: ",subtitle1))+
        theme(
          legend.position = "none",
          plot.subtitle = element_text(size=10, color="black")
        )->mirror_plot1
      
      
      spectrum_file_plot %>% select(peak_table_for_graph) %>%
        unnest(cols=peak_table_for_graph) %>% 
        filter(!is.na(ms2error)) %>% 
        ggplot(
          aes(x=`m/z`,y=ms2error,color = series)
        )+
        geom_point(size=1)+
        geom_hline(yintercept = 0, color = "black")+
        scale_color_manual(
          values = c("b"=RColorBrewer::brewer.pal(n=9,"Set1")[2],"y"=RColorBrewer::brewer.pal(n=9,"Set1")[1])
        )+
        scale_x_continuous(
          limits=c(0, spec_plot_xmax),expand=c(0,0)
        )+
        scale_y_continuous(
          limits = c(-0.05,0.05),
          labels = \(x) x*1000
        )+
        labs(
          y="Mass error (mDa)"
        )+
        theme_bw()+
        theme(
          legend.position = "none"
        )->mirror_plot2
      ggarrange(
        mirror_plot1,mirror_plot2,
        ncol=1,
        heights = c(2,1),
        align = "v"
      )
      ggsave(paste0("./ms2figures/",j,".png"), width=20, height = 12, units = "cm", dpi = 600)
    }
  }
  
  ##Prediction vs. measured spectra best case
  {
    summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
      filter(gene %in% c("FBLN1","CALR","P4HB")) %>% 
      arrange(
        prot_pos_reposition,desc(PCC)
      ) %>% 
      select(
        prot_pos_reposition,PCC,index
      ) %>% print(n=Inf)
    
    #CALR: 1028
    #FBLN1: 267
    #P4HB: 366
  }
  
  
  ###FDR: false TTEST
  {
    fdr_pcc_negative_set7 %>%
      filter(enzyme == "trypsin") %>% 
      summarize(
        .by=fdr,
        count =n(),
        percol_score= mean(percolator_score),
        #.by = moddb_result
      )
    fdr_pcc_negative_set7 %>% 
      filter(enzyme == "trypsin") %>% 
      ggplot(
        aes(x=as.character(fdr),y=percolator_score,fill=as.character(fdr))
      )+
      geom_boxplot()+
      stat_compare_means(
        aes(label = paste0("p= ",after_stat(p.format))),
        comparisons = list(c("1","0")),
        method="t.test",
        p.adjust.methods="BH"
      )+
      labs(y="Percolator Score")+
      scale_x_discrete(
        labels= c("TRUE","FALSE")
      )+
      ylim(c(-2.5,5.5))+
      scale_fill_manual(
        values = c("#4DAF4A","#E41A1C")
      )+
      theme_classic2()+
      theme(
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black",size = 9),
        axis.text.y = element_text(color = "black",size = 9)
      )
    ggsave(
      "./figures/2_03_false_ttest.png",width=8,height=8, units = "cm",dpi = 600
    )
    
  }
  #2_03_false_ttest.png
  
  ###FDR: result after fdr
  {
    summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
      summarize(
        .by = pcc_check,
        count = n(),
        bion_pcc = median(bion_pcc,na.rm= TRUE),
        yion_pcc = median(yion_pcc,na.rm= TRUE)
      )
    summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
      select(
        prot_pos_reposition,PCC,bion_pcc,yion_pcc,pcc_check
      ) %>% 
      pivot_longer(cols=c(bion_pcc,yion_pcc)) %>% 
      mutate(
        group = case_when(
          pcc_check==1 ~ fct("1",levels = c("1","0")),
          pcc_check==0 ~ fct("0",levels = c("1","0"))
        )
      ) %>% 
      ggplot(
        aes(y=value,x=group,fill=name)
      )+
      geom_boxplot(
        size = 0.2,
        outlier.size = 0.2,
      )+
      annotate(
        "text",
        x=1:2,
        y = 1.08,
        label = paste("N = ",c(594,623)),
        size = 2
      )+
      scale_x_discrete(
        limits = c("0","1"),
        labels = c(expression("FDR">"1%", "FDR"<="1%"))
      )+
      scale_fill_manual(labels = c("b-ion","y-ion"),values = c(wes_palette("Zissou1")[1],wes_palette("Zissou1")[5]))+
      labs(y="PCC")+
      theme_classic()+
      theme(
        legend.position = "none",
        legend.title = element_blank(),
        text= element_text(color = "black", size = 7),
        axis.text.y = element_text(color = "black", size = 7, hjust=1),
        axis.text.x = element_text(color = "black", size = 7),
        axis.title.x = element_blank(),
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm")
      )
    ggsave("./figures/2_07_fdr_ions.png",width=4.5,height=4,dpi = 600,units="cm")
  }
  
  ###chapter RT prediction model
  ##RT prediction figure ::@summaryfile_psm_for_testing_pcc_rt
  {
    output_prediction_rt_train_result
    summary_rtprediction_t95
    summary_rtprediction
    
    output_prediction_rt_test_result %>% filter(file_series==1) 
    
    
    ##figure
    i=1
    
    paste0(
      "Deltat[95] = ",str_sub(paste0(summary_rtprediction_t95[i,3]),1,5)
    )
    
    i=1
    for (i in 1L:12L){
      output_prediction_rt_train_result %>% 
        filter(file_series == i) %>% 
        ggplot(aes(x=rt_observed, y=rt_predicted))+
        geom_point(size=0.5,alpha=0.5) +
        annotate(
          "text",x=-Inf,y=Inf,hjust=0,vjust=1.2,
          label=bquote(
            ~Delta*italic(t)[95*'%'] == .(str_sub(paste0(summary_rtprediction_t95[i,3]),1,5))*","~R^2==.(summary_rtprediction %>% filter(file_series==i,.metric =="rsq") %>% pull(.estimate) %>% round(3))
          ),size =4
        )+
        annotate(
          "text",x=-Inf,y=Inf,hjust=0,vjust=3.2,
          label=bquote(
            ~MAE==.(summary_rtprediction %>% filter(file_series==i,.metric =="mae") %>% pull(.estimate) %>% round(3))*","~N==.(output_prediction_rt_train_result %>% filter(file_series==i) %>% nrow())
          ),size =4
        )+
        geom_line(aes(x=rt_lwr), color = "red", linetype = "dashed")+
        geom_line(aes(x=rt_upr), color = "red", linetype = "dashed")+
        geom_smooth(method=lm, se=TRUE)+
        #lims(x=c(0,1),y=c(0,1))+
        labs(
          x="Observed RT (min)",
          y="Predicted RT (min)"
        )+
        theme_classic2()+
        theme(
          axis.text = element_text(color = "black",size=10),
          axis.title = element_text(color = "black",size=11)
        )
      ggsave(
        paste0("./figures/rt/rt_prediction_training_set",i,".png"),width=10,height=10,dpi = 600,units="cm"
      )
    }
    i=1
    
    #rt mainfigure
    {
      output_prediction_rt_train_result %>% 
        filter(file_series == i) %>% 
        ggplot(aes(x=rt_observed, y=rt_predicted))+
        geom_point(size=0.1,alpha=0.25) +
        annotate(
          "text",x=-Inf,y=Inf,hjust=0,vjust=1.2,
          label=bquote(
            ~Delta*italic(t)[95*'%'] == .(str_sub(paste0(summary_rtprediction_t95[i,3]),1,5))*","~R^2==.(summary_rtprediction %>% filter(file_series==i,.metric =="rsq") %>% pull(.estimate) %>% round(3))
          ),size =2
        )+
        annotate(
          "text",x=-Inf,y=Inf,hjust=0,vjust=3.2,
          label=bquote(
            ~MAE==.(summary_rtprediction %>% filter(file_series==i,.metric =="mae") %>% pull(.estimate) %>% round(3))*","~N==.(output_prediction_rt_train_result %>% filter(file_series==i) %>% nrow())
          ),size =2
        )+
        geom_line(aes(x=rt_lwr), color = "red", linetype = "dashed",size =0.5)+
        geom_line(aes(x=rt_upr), color = "red", linetype = "dashed",size =0.5)+
        geom_smooth(method=lm, se=TRUE, size = 0.5)+
        #lims(x=c(0,1),y=c(0,1))+
        labs(
          x="Observed RT (min)",
          y="Predicted RT (min)"
        )+
        theme_classic2()+
        theme(
          legend.position = "none",
          #axis.title.x = element_blank(),
          axis.title = element_text(color = "black",size = 7),
          #axis.text.x = element_blank(),
          axis.text = element_text(color = "black",size = 6),
          plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm")
        )
      ggsave(
        paste0("./figures/2_08_rt_prediction_training_set_mainfigure",i,".png"),width=4.5,height=4.5,dpi = 600,units="cm"
      )
      
      paste0(
        output_prediction_rt_test_result %>% 
          filter(file_series == i) %>% 
          summarise(
            .by=rt_test,
            count = n()
          ) %>% 
          mutate(sum = sum(count)) %>% pull(count) %>% .[2],"/",
        output_prediction_rt_test_result %>% 
          filter(file_series == i) %>% 
          summarise(
            .by=rt_test,
            count = n()
          ) %>% 
          mutate(sum = sum(count)) %>% pull(sum) %>% .[1]
      ) -> count_for_rt
      
      
      output_prediction_rt_test_result %>% 
        filter(file_series == i) %>% 
        ggplot(aes(x=rt_observed, y=rt_predicted))+
        geom_point(
          aes(color = fct(as.character(rt_test))),
          size=0.1
        )+
        annotate(
          "text",x=-Inf,y=Inf,hjust=-0.2,vjust=2,
          label=bquote(
            Passed ==.(count_for_rt)
          ),size =2
        )+
        geom_line(aes(x=rt_lwr), color = "red", linetype = "dashed",size =0.5)+
        geom_line(aes(x=rt_upr), color = "red", linetype = "dashed",size =0.5)+
        #geom_smooth(method=lm, se=TRUE)+
        #lims(x=c(0,1),y=c(0,1))+
        scale_color_manual(
          values = c("red","black")
        )+
        labs(
          x="Observed RT (min)",
          y="Predicted RT (min)"
        )+
        theme_classic2()+
        theme(
          legend.position = "none",
          #axis.title.x = element_blank(),
          axis.title = element_text(color = "black",size = 7),
          #axis.text.x = element_blank(),
          axis.text = element_text(color = "black",size = 7),
          plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm")
        )
      ggsave(
        paste0("./figures/2_09_rt2.png"),width=4.5,height=4.5,dpi = 600,units="cm"
      )
    }
    
    i=1
    for (i in 1:12){
      paste0(
        output_prediction_rt_test_result %>% 
          filter(file_series == i) %>% 
          summarise(
            .by=rt_test,
            count = n()
          ) %>% 
          mutate(sum = sum(count)) %>% pull(count) %>% .[2],"/",
        output_prediction_rt_test_result %>% 
          filter(file_series == i) %>% 
          summarise(
            .by=rt_test,
            count = n()
          ) %>% 
          mutate(sum = sum(count)) %>% pull(sum) %>% .[1]
      ) -> count_for_rt
      
      
      output_prediction_rt_test_result %>% 
        filter(file_series == i) %>% 
        ggplot(aes(x=rt_observed, y=rt_predicted))+
        geom_point(
          aes(color = fct(as.character(rt_test))),
          size=0.5,alpha=0.5
        )+
        annotate(
          "text",x=-Inf,y=Inf,hjust=-0.2,vjust=2,
          label=bquote(
            Passed ==.(count_for_rt)
          ),size =4
        )+
        geom_line(aes(y=rt_lwr), color = "red", linetype = "dashed")+
        geom_line(aes(y=rt_upr), color = "red", linetype = "dashed")+
        #geom_smooth(method=lm, se=TRUE)+
        #lims(x=c(0,1),y=c(0,1))+
        scale_color_manual(
          values = c("black","red")
        )+
        labs(
          x="Observed RT (min)",
          y="Predicted RT (min)"
        )+
        theme_classic2()+
        theme(
          legend.position = "none",
          axis.text = element_text(color = "black",size=10),
          axis.title = element_text(color = "black",size=11)
        )
      ggsave(
        paste0("./figures/rt/rt_prediction_test_set",i,".png"),width=10,height=10,dpi = 600,units="cm"
      )
    }
    
    ##summary
    {
      summaryfile_psm_for_testing3_pcc_rt_hfsm %>% 
        summarize(
          .by = rt_check,
          count = n()
        )
      summaryfile_psm_for_testing3_pcc_rt_hfsm %>% 
        mutate(
          rt_check = fct(as.character(rt_check)),
          rt_dev = rt_pred-rt_norm
        ) %>% 
        gghistogram(
          x = "rt_dev", 
          add = "mean", rug = TRUE,
          fill = "rt_check", palette = c("#D53E4F", "#3288BD"),
          binwidth = 0.025
        )+
        labs(x="Normalized RT Deviation",y="Count")+
        scale_x_continuous(limits = c(-0.7,0.7))+
        theme(
          legend.position = "none"
        )
      ggsave("./figures/2_04_rt_summary.png",width=8,height=8, units = "cm",dpi = 600)
      
    }
    display.brewer.all()
    brewer.pal(11,"Spectral")
  }
  
  
  ##RT Additional GV VG training
  {
    readRDS(file = "./export/psm_563516.rds")->psm_563516
    
    file_raw_index2
    
    psm_563516 %>% 
      filter(str_detect(strip_seq,"^GV")|str_detect(strip_seq,"^VG")) %>% 
      filter(str_detect(mods,"3\\)\\@Any")) -> training_set_for_rt_gv ##2586
    
    training_set_for_rt_gv %>% 
      summarize(
        .by = file_id_series,
        count = n()
      ) ###F1,F2 target
    
    training_set_for_rt_gv %>% 
      summarize(
        .by = file_id_series,
        count = n()
      )
    
    ##training set
    {
      training_set_for_rt_gv %>% 
        filter(file_id_series == "F1") %>% 
        select(
          strip_seq,index,mods,mod_sites,nAA,nce,instrument,rt_norm,file_id_series
        ) %>% 
        rename(
          sequence = strip_seq
        ) %>% 
        write_tsv(
          paste0(folder_alphapeptdeep,"inrich/input/input_rt_gv_model_f1_231109.tsv")
        )
      
      training_set_for_rt_gv %>% 
        filter(file_id_series == "F2") %>% 
        select(
          strip_seq,index,mods,mod_sites,nAA,nce,instrument,rt_norm,file_id_series
        ) %>% 
        rename(
          sequence = strip_seq
        ) %>% 
        write_tsv(
          paste0(folder_alphapeptdeep,"inrich/input/input_rt_gv_model_f2_231109.tsv")
        )
      
      training_set_for_rt_arg %>% 
        filter(file_id_series == "F1") %>% 
        select(
          strip_seq,index,mods,mod_sites,nAA,nce,instrument,rt_norm,file_id_series
        ) %>% 
        rename(
          sequence = strip_seq
        ) %>% 
        write_tsv(
          paste0(folder_alphapeptdeep,"inrich/input/input_rt_gv_model_f1_231109.tsv")
        )
      
      training_set_for_rt_gv %>% 
        filter(file_id_series == "F2") %>% 
        select(
          strip_seq,index,mods,mod_sites,nAA,nce,instrument,rt_norm,file_id_series
        ) %>% 
        rename(
          sequence = strip_seq
        ) %>% 
        write_tsv(
          paste0(folder_alphapeptdeep,"inrich/input/input_rt_gv_model_f2_231109.tsv")
        )
      }
    
    ###test set generation
    {
      summaryfile_psm_for_testing3 %>% view()
      
      summaryfile_psm_for_testing3 %>% 
        left_join(
          summaryfile_psm_for_testing3 %>% 
            separate_rows(
              mod_sites,sep = "\\;"
            ) %>% 
            mutate(
              mod_sites = as.integer(mod_sites)
            ) %>% 
            mutate(
              mod_sites = if_else(mod_sites == 0,mod_sites,mod_sites+1)
            ) %>% 
            nest(
              gv_mod_sites=mod_sites,
              .by = index
            ) %>% 
            mutate(
              gv_mod_sites = map_chr(
                gv_mod_sites,
                \(x) pull(x,mod_sites) %>% paste0(collapse = ";")
              )
            ),
          join_by(index)
        ) -> summaryfile_psm_for_testing3_rt_test_gv
      
      
      summaryfile_psm_for_testing3_rt_test_gv %>% 
        filter(
          str_detect(Sequence,"^[DE]")
        ) %>% 
        filter(
          set == "F1"
        ) %>% 
        mutate(
          nAA = nAA+1
        ) %>% 
        mutate(
          seq_gv_introduced = paste0("GV",Sequence)
        ) %>% 
        select(
          seq_gv_introduced,index,mods,gv_mod_sites,nAA,nce,instrument,rt_norm,set
        ) %>% 
        rename(
          sequence = seq_gv_introduced,
          mod_sites = gv_mod_sites
        ) %>% 
        write_tsv(
          paste0(folder_alphapeptdeep,"inrich/input/input_rt_gv_test_f1_231109.tsv")
        )
      
      ### GV_F2
      summaryfile_psm_for_testing3_rt_test_gv %>% 
        filter(
          str_detect(Sequence,"^[DE]")
        ) %>% 
        filter(
          set == "F2"
        ) %>% 
        mutate(
          nAA = nAA+1
        ) %>% 
        mutate(
          seq_gv_introduced = paste0("GV",Sequence)
        ) %>% 
        select(
          seq_gv_introduced,index,mods,gv_mod_sites,nAA,nce,instrument,rt_norm,set
        ) %>% 
        rename(
          sequence = seq_gv_introduced,
          mod_sites = gv_mod_sites
        ) %>% 
        write_tsv(
          paste0(folder_alphapeptdeep,"inrich/input/input_rt_gv_test_f2_231109.tsv")
        )
      
      ###normal F1
      summaryfile_psm_for_testing3 %>% 
        filter(
          str_detect(Sequence,"^[DE]")
        ) %>% 
        filter(
          set == "F1"
        ) %>%
        select(
          sequence,index,mods,mod_sites,nAA,nce,instrument,rt_norm,set
        ) %>% 
        write_tsv(
          paste0(folder_alphapeptdeep,"inrich/input/input_rt_test_f1_231109.tsv")
        )
      
      ###normal F2
      summaryfile_psm_for_testing3 %>% 
        filter(
          str_detect(Sequence,"^[DE]")
        ) %>% 
        filter(
          set == "F2"
        ) %>%
        select(
          sequence,index,mods,mod_sites,nAA,nce,instrument,rt_norm,set
        ) %>% 
        write_tsv(
          paste0(folder_alphapeptdeep,"inrich/input/input_rt_test_f2_231109.tsv")
        )
      
    }
    
    ###alphapeptdeep::RTpredict_Transfer.ipynb
    
    ###training arg model
    {
      folder_alphapeptdeep
      i=1
      output_prediction_rt_train_result <- tibble()
      output_prediction_rt_test_result <- tibble()
      for (i in 1:2){
        output_prediction_rt_training <- read_tsv(
          paste0(folder_alphapeptdeep,"inrich/output/output_rt_training_f",i,".tsv")
        )[,-1] %>% 
          arrange(rt_norm)
        
        output_prediction_rt_test <- read_tsv(
          paste0(folder_alphapeptdeep,"inrich/output/output_rt_normal_arg_f",i,"_231109.tsv")
        )[,-1] %>% 
          arrange(rt_norm) 
        
        fit_rt_training <- lm(rt_pred~ rt_norm, data =output_prediction_rt_training)
        
        output_prediction_rt_training %>%
          bind_cols(
            predict(fit_rt_training, interval = "prediction", level =0.90, newdata = output_prediction_rt_training)
          )->output_prediction_rt_training
        
        
        output_prediction_rt_test %>% 
          bind_cols(
            predict(fit_rt_training, interval = "prediction", level =0.90, newdata = output_prediction_rt_test)
          )->output_prediction_rt_test
        
        output_prediction_rt_training %>% 
          write_tsv(
            paste0("./export/output_rt_training_prediction_f",i,".tsv")
          )
        
        output_prediction_rt_test %>% 
          write_tsv(
            paste0("./export/output_rt_test_prediction_f",i,".tsv")
          )
      }
    }
    ###gv vg model
    {
      folder_alphapeptdeep
      i=1
      output_prediction_rt_train_result <- tibble()
      output_prediction_rt_test_result <- tibble()
      for (i in 1:2){
        output_prediction_rt_training <- read_tsv(
          paste0(folder_alphapeptdeep,"inrich/output/output_rt_gv_model_f",i,"_231109.tsv")
        )[,-1] %>% 
          arrange(rt_norm)
        
        output_prediction_rt_test <- read_tsv(
          paste0(folder_alphapeptdeep,"inrich/output/output_rt_gv_test_f",i,"_231109.tsv")
        )[,-1] %>% 
          arrange(rt_norm) 
        
        fit_rt_training <- lm(rt_pred~ rt_norm, data =output_prediction_rt_training)
        
        output_prediction_rt_training %>%
          bind_cols(
            predict(fit_rt_training, interval = "prediction", level =0.90, newdata = output_prediction_rt_training)
          )->output_prediction_rt_training
        
        
        output_prediction_rt_test %>% 
          bind_cols(
            predict(fit_rt_training, interval = "prediction", level =0.90, newdata = output_prediction_rt_test)
          )->output_prediction_rt_test
        
        output_prediction_rt_training %>% 
          write_tsv(
            paste0("./export/output_rt_training_gv_prediction_f",i,".tsv")
          )
        
        output_prediction_rt_test %>% 
          write_tsv(
            paste0("./export/output_rt_test_gv_prediction_f",i,".tsv")
          )
      }
    }
    
    
    
  }
  
  ##RT data interpretation
  {
    read_tsv(
      "./export/output_rt_training_prediction_f1.tsv"
    ) %>% mutate(
      file_id = "F1"
    ) %>% 
      bind_rows(
        read_tsv(
          "./export/output_rt_training_prediction_f2.tsv"
        ) %>% mutate(
          file_id = "F2"
        )
      ) -> output_prediction_rt_training_arg_231109
    
    read_tsv(
      "./export/output_rt_test_prediction_f1.tsv"
    ) %>% mutate(
      file_id = "F1"
    ) %>% 
      bind_rows(
        read_tsv(
          "./export/output_rt_test_prediction_f2.tsv"
        ) %>% mutate(
          file_id = "F2"
        )
      ) -> output_prediction_rt_test_arg_231109
    
    read_tsv(
      "./export/output_rt_training_gv_prediction_f1.tsv"
    ) %>% mutate(
      file_id = "F1"
    ) %>% 
      bind_rows(
        read_tsv(
          "./export/output_rt_training_gv_prediction_f2.tsv"
        ) %>% mutate(
          file_id = "F2"
        )
      ) -> output_prediction_rt_training_gv_231109
    
    read_tsv(
      "./export/output_rt_test_gv_prediction_f1.tsv"
    ) %>% mutate(
      file_id = "F1"
    ) %>% 
      bind_rows(
        read_tsv(
          "./export/output_rt_test_gv_prediction_f2.tsv"
        ) %>% mutate(
          file_id = "F2"
        )
      ) -> output_prediction_rt_test_gv_231109
    
    
    
  }
  #output_prediction_rt_training_arg_231109
  #output_prediction_rt_test_arg_231109
  #output_prediction_rt_training_gv_231109
  #output_prediction_rt_test_gv_231109
  
  ###GV TEST RESULT ###COUNT!
  #@output_prediction_rt_test_gv_231109
  {
    output_prediction_rt_test_gv_231109 %>% 
      mutate(
        gv_test = if_else(
          rt_pred >= lwr & rt_pred <= upr,1,0
        )
      ) %>% 
      select(index,gv_test) %>% 
      inner_join(
        summaryfile_psm_for_testing3_pcc_rt_hfsm %>% select(index,rt_check),
        join_by(index)
      ) -> table_rt_stage_gv_test
    
    table_rt_stage_gv_test %>% 
      mutate(
        test_result = case_when(
          gv_test == 0 & rt_check == 1 ~ "Passed",
          gv_test == 1 & rt_check == 0 ~ "Passed",
          gv_test == 1 & rt_check == 1 ~ "Failed",
          gv_test == 0 & rt_check == 0 ~ "Passed",
          .default = NA
        )
      ) %>% 
      summarize(
        .by = test_result,
        count = n()
      ) %>% write_tsv("./figures/rt/gv_vg_test_f1tof2.result")
  }
  
  ##Figs
  {
    output_prediction_rt_training_arg_231109 %>% 
      filter(file_id == "F1") %>% 
      ggplot(aes(x=rt_norm, y=rt_pred))+
      geom_point(
        size=1,alpha=0.5
      ) +
      geom_line(aes(y=lwr), color = "red", linetype = "dashed")+
      geom_line(aes(y=upr), color = "red", linetype = "dashed")+
      #geom_smooth(method=lm, se=TRUE)+
      scale_color_manual(labels = c("true","false"),values = c("red","black"))+
      lims(x=c(0,1),y=c(0,1))+
      labs(
        x="Observed Normalized RT",
        y="Predicted Normalized RT",
        color = "FDR >= 0.01"
      )+
      theme_classic2()+
      theme(
        legend.position = "top"
      )
    ggsave("./figures/2_04_rt_training_arg_f1.png",width=8,height=8, units = "cm",dpi = 600)
    
    output_prediction_rt_training_arg_231109 %>% 
      filter(file_id == "F2") %>% 
      ggplot(aes(x=rt_norm, y=rt_pred))+
      geom_point(
        size=1,alpha=0.5
      ) +
      geom_line(aes(y=lwr), color = "red", linetype = "dashed")+
      geom_line(aes(y=upr), color = "red", linetype = "dashed")+
      #geom_smooth(method=lm, se=TRUE)+
      scale_color_manual(labels = c("true","false"),values = c("red","black"))+
      lims(x=c(0,1),y=c(0,1))+
      labs(
        x="Observed Normalized RT",
        y="Predicted Normalized RT",
        color = "FDR >= 0.01"
      )+
      theme_classic2()+
      theme(
        legend.position = "top"
      )
    ggsave("./figures/2_04_rt_training_arg_f2.png",width=8,height=8, units = "cm",dpi = 600)
    
    output_prediction_rt_test_arg_231109 %>% 
      filter(file_id == "F1") %>% 
      ggplot(aes(x=rt_norm, y=rt_pred))+
      geom_point(
        size=1,alpha=0.5
      ) +
      geom_line(aes(y=lwr), color = "red", linetype = "dashed")+
      geom_line(aes(y=upr), color = "red", linetype = "dashed")+
      #geom_smooth(method=lm, se=TRUE)+
      scale_color_manual(labels = c("true","false"),values = c("red","black"))+
      lims(x=c(0,1),y=c(0,1))+
      labs(
        x="Observed Normalized RT",
        y="Predicted Normalized RT",
        color = "FDR >= 0.01"
      )+
      theme_classic2()+
      theme(
        legend.position = "top"
      )
    ggsave("./figures/2_04_rt_test_arg_f1.png",width=8,height=8, units = "cm",dpi = 600)
    
    output_prediction_rt_test_arg_231109 %>% 
      filter(file_id == "F2") %>% 
      ggplot(aes(x=rt_norm, y=rt_pred))+
      geom_point(
        size=1,alpha=0.5
      ) +
      geom_line(aes(y=lwr), color = "red", linetype = "dashed")+
      geom_line(aes(y=upr), color = "red", linetype = "dashed")+
      #geom_smooth(method=lm, se=TRUE)+
      scale_color_manual(labels = c("true","false"),values = c("red","black"))+
      lims(x=c(0,1),y=c(0,1))+
      labs(
        x="Observed Normalized RT",
        y="Predicted Normalized RT",
        color = "FDR >= 0.01"
      )+
      theme_classic2()+
      theme(
        legend.position = "top"
      )
    ggsave("./figures/2_04_rt_test_arg_f2.png",width=8,height=8, units = "cm",dpi = 600)
    
    ###TRAINING GV F1 264
    output_prediction_rt_training_gv_231109 %>% 
      filter(file_id == "F1") %>% 
      ggplot(aes(x=rt_norm, y=rt_pred))+
      geom_point(
        size=1,alpha=0.5
      ) +
      geom_line(aes(y=lwr), color = "red", linetype = "dashed")+
      geom_line(aes(y=upr), color = "red", linetype = "dashed")+
      #geom_smooth(method=lm, se=TRUE)+
      scale_color_manual(labels = c("true","false"),values = c("red","black"))+
      lims(x=c(0,1),y=c(0,1))+
      labs(
        x="Observed Normalized RT",
        y="Predicted Normalized RT",
        color = "FDR >= 0.01"
      )+
      theme_classic2()+
      theme(
        legend.position = "top"
      )
    ggsave("./figures/2_04_rt_training_gv_f1.png",width=8,height=8, units = "cm",dpi = 600)
    
    
    ##TRAINING GV 240
    output_prediction_rt_training_gv_231109 %>% 
      filter(file_id == "F2") %>% 
      ggplot(aes(x=rt_norm, y=rt_pred))+
      geom_point(
        size=1,alpha=0.5
      ) +
      geom_line(aes(y=lwr), color = "red", linetype = "dashed")+
      geom_line(aes(y=upr), color = "red", linetype = "dashed")+
      #geom_smooth(method=lm, se=TRUE)+
      scale_color_manual(labels = c("true","false"),values = c("red","black"))+
      lims(x=c(0,1),y=c(0,1))+
      labs(
        x="Observed Normalized RT",
        y="Predicted Normalized RT",
        color = "FDR >= 0.01"
      )+
      theme_classic2()+
      theme(
        legend.position = "top"
      )
    ggsave("./figures/2_04_rt_training_gv_f2.png",width=8,height=8, units = "cm",dpi = 600)
    
    output_prediction_rt_test_gv_231109 %>% 
      filter(file_id == "F1") %>% 
      ggplot(aes(x=rt_norm, y=rt_pred))+
      geom_point(
        size=1,alpha=0.5
      ) +
      geom_line(aes(y=lwr), color = "red", linetype = "dashed")+
      geom_line(aes(y=upr), color = "red", linetype = "dashed")+
      #geom_smooth(method=lm, se=TRUE)+
      scale_color_manual(labels = c("true","false"),values = c("red","black"))+
      lims(x=c(0,1),y=c(0,1))+
      labs(
        x="Observed Normalized RT",
        y="Predicted Normalized RT",
        color = "FDR >= 0.01"
      )+
      theme_classic2()+
      theme(
        legend.position = "top"
      )
    ggsave("./figures/2_04_rt_test_gv_f1.png",width=8,height=8, units = "cm",dpi = 600)
    
    output_prediction_rt_test_gv_231109 %>% 
      filter(file_id == "F2") %>% 
      ggplot(aes(x=rt_norm, y=rt_pred))+
      geom_point(
        size=1,alpha=0.5
      ) +
      geom_line(aes(y=lwr), color = "red", linetype = "dashed")+
      geom_line(aes(y=upr), color = "red", linetype = "dashed")+
      #geom_smooth(method=lm, se=TRUE)+
      scale_color_manual(labels = c("true","false"),values = c("red","black"))+
      lims(x=c(0,1),y=c(0,1))+
      labs(
        x="Observed Normalized RT",
        y="Predicted Normalized RT",
        color = "FDR >= 0.01"
      )+
      theme_classic2()+
      theme(
        legend.position = "top"
      )
    ggsave("./figures/2_04_rt_test_gv_f2.png",width=8,height=8, units = "cm",dpi = 600)
    
  }
  #2_04_rt_training_arg_f1.png
  
  ###RT sample figure
  {
    summaryfile_psm_for_testing3_pcc_rt_hfsm %>% 
      filter(set == "F1") %>% 
      ggplot(aes(x=rt_norm, y=rt_pred))+
      geom_point(
        aes(color = as.character(rt_check)),
        size=1,alpha=0.5
      ) +
      geom_line(aes(y=lwr), color = "red", linetype = "dashed")+
      geom_line(aes(y=upr), color = "red", linetype = "dashed")+
      #geom_smooth(method=lm, se=TRUE)+
      scale_color_manual(labels = c("1","0"),values = c("black","red"))+
      lims(x=c(0.2,1),y=c(0.2,1))+
      labs(
        x="Observed Normalized RT",
        y="Predicted Normalized RT"
      )+
      theme_classic2()+
      theme(
        legend.position = "none"
      )
    ggsave("./figures/rt/rt_result_f1.png",width =10,height = 10, dpi = 600, units = "cm")
    
    summaryfile_psm_for_testing3_pcc_rt_hfsm %>% 
      filter(set == "F1") %>%  #333
      filter(rt_check == 1) ##156
    
    summaryfile_psm_for_testing3_pcc_rt_hfsm %>% 
      filter(set == "F1") %>%  #333
      filter(pcc_check == 1 & rt_check == 1) ##137
    
    summaryfile_psm_for_testing3_pcc_rt_hfsm %>% 
      filter(set == "F1") %>%  #333
      filter(pcc_check == 1) ##217
    
    summaryfile_psm_for_testing3_pcc_rt_hfsm %>% 
      filter(set == "F1") %>% 
      ggplot(aes(x=rt_norm, y=rt_pred))+
      geom_point(
        aes(color = as.character(pcc_check)),
        size=1,alpha=0.5
      ) +
      geom_line(aes(y=lwr), color = "red", linetype = "dashed")+
      geom_line(aes(y=upr), color = "red", linetype = "dashed")+
      #geom_smooth(method=lm, se=TRUE)+
      scale_color_manual(labels = c("1","0"),values = c("black","blue"))+
      lims(x=c(0.2,1),y=c(0.2,1))+
      labs(
        x="Observed Normalized RT",
        y="Predicted Normalized RT"
      )+
      theme_classic2()+
      theme(
        legend.position = "none"
      )
    ggsave("./figures/rt/rt_result_f1_pcc.png",width =10,height = 10, dpi = 600, units = "cm")
  }
  #rt_result_f1_pcc.png
  
  ##2K ditribution of rt
  {
    summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
      mutate(
        rt_deviation = rt_norm-rt_pred,
        rt_test = case_when(
          rt_deviation > 0 & rt_norm > upr ~ "increase",
          rt_deviation < 0 & rt_norm < lwr ~ "decrease",
          .default = "none"
        )
      ) %>% 
      summarize(
        .by = rt_test,
        count = n()
      )
    
    summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
      mutate(
        rt_deviation = rt_norm-rt_pred,
        rt_test = case_when(
          rt_deviation > 0 & rt_norm > upr ~ "increase",
          rt_deviation < 0 & rt_norm < lwr ~ "decrease",
          .default = "none"
        )
      ) %>% 
      gghistogram(
        x = "rt_deviation",
        add="median",
        fill = "rt_test",
        position = "stack",
        size = 0.1
      )+
      labs(
        x="RT Deviation",
        y="PSM count",
        fill = "RT prediction",
        color = "RT prediction"
      )+
      scale_fill_manual(
        labels = c("Decrease","Increase",expression(" "*Under~Delta*t[95*"%"])),
        values = c("#66C2A5","#FC8D62","#8DA0CB")
      )+
      scale_color_manual(
        labels = c("Decrease","Increase",expression(" "*Under~Delta*t[95*"%"])),
        values = c("#66C2A5","#FC8D62","#8DA0CB")
      )+
      theme(
        #legend.position = "none",
        # legend.margin = margin(0.1,0.1,0.1,0.1, "cm"),
        # legend.spacing = unit(c(0,0,0,0), "cm"),
        # legend.title = element_blank(),
        # legend.key.height = unit(0.1,"cm"),
        # legend.key.width = unit(0.1,"cm"),
        # legend.text = element_text(color = "black",size = 7),
        #axis.title.x = element_blank(),
        axis.title = element_text(color = "black",size = 7),
        #axis.text.x = element_blank(),
        axis.text = element_text(color = "black",size = 7),
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm")
      )
    ggsave("./figures/2_04_rt_result_all.png",width =4.5,height = 4.5, dpi = 600, units = "cm")
    RColorBrewer::display.brewer.all()
    RColorBrewer::brewer.pal(8,"Set2")
  }
  #2_04_rt_result_all.png
  
  ###GV figure
  {
    summaryfile_psm_for_testing3_pcc_rt_hfsm %>% 
      inner_join(
        output_prediction_rt_test_gv_231109 %>% 
          mutate(
            gv_test = if_else(
              rt_pred >= lwr & rt_pred <= upr,1,0
            )
          ) %>% 
          select(index,gv_test),
        join_by(index)
      ) %>% 
      filter(set == "F1") %>% 
      ggplot(aes(x=rt_norm, y=rt_pred))+
      geom_point(
        aes(color = as.character(gv_test)),
        size=1,alpha=0.5
      ) +
      geom_line(aes(y=lwr), color = "red", linetype = "dashed")+
      geom_line(aes(y=upr), color = "red", linetype = "dashed")+
      #geom_smooth(method=lm, se=TRUE)+
      scale_color_manual(labels = c("1","0"),values = c("black","red"))+
      lims(x=c(0.2,1),y=c(0.2,1))+
      labs(
        x="Observed Normalized RT",
        y="Predicted Normalized RT"
      )+
      theme_classic2()+
      theme(
        legend.position = "none"
      )
    ggsave("./figures/rt/rt_gv_test_result_f1.png",width =10,height = 10, dpi = 600, units = "cm")
    
    summaryfile_psm_for_testing3_pcc_rt_hfsm %>% 
      inner_join(
        output_prediction_rt_test_gv_231109 %>% 
          mutate(
            gv_test = if_else(
              rt_pred >= lwr & rt_pred <= upr,1,0
            )
          ) %>% 
          select(index,gv_test),
        join_by(index)
      ) %>% 
      filter(set == "F1") %>%  ##245
      #filter(rt_check == 1) ##105
      filter(gv_test == 1 & rt_check == 1) ##5
    
    summaryfile_psm_for_testing3_pcc_rt_hfsm %>% 
      inner_join(
        output_prediction_rt_test_gv_231109 %>% 
          mutate(
            gv_test = if_else(
              rt_pred >= lwr & rt_pred <= upr,1,0
            )
          ) %>% 
          select(index,gv_test),
        join_by(index)
      ) %>% 
      filter(set == "F2") %>% 
      ggplot(aes(x=rt_norm, y=rt_pred))+
      geom_point(
        aes(color = as.character(gv_test)),
        size=1,alpha=0.5
      ) +
      geom_line(aes(y=lwr), color = "red", linetype = "dashed")+
      geom_line(aes(y=upr), color = "red", linetype = "dashed")+
      #geom_smooth(method=lm, se=TRUE)+
      scale_color_manual(labels = c("1","0"),values = c("black","red"))+
      lims(x=c(0.2,1),y=c(0.2,1))+
      labs(
        x="Observed Normalized RT",
        y="Predicted Normalized RT"
      )+
      theme_classic2()+
      theme(
        legend.position = "none"
      )
    ggsave("./figures/rt/rt_gv_test_result_f2.png",width =10,height = 10, dpi = 600, units = "cm")
    
    output_prediction_rt_test_gv_231109 %>% 
      mutate(
        gv_test = if_else(
          rt_pred >= lwr & rt_pred <= upr,1,0
        )
      ) %>% 
      filter(file_id == "F1") %>% 
      ggplot(aes(x=rt_norm, y=rt_pred))+
      geom_point(
        aes(color = as.character(gv_test)),
        size=1,alpha=0.5
      ) +
      geom_line(aes(y=lwr), color = "blue", linetype = "dashed")+
      geom_line(aes(y=upr), color = "blue", linetype = "dashed")+
      #geom_smooth(method=lm, se=TRUE)+
      scale_color_manual(labels = c("1","0"),values = c("black","blue"))+
      lims(x=c(0.2,1),y=c(0.2,1))+
      labs(
        x="Observed Normalized RT",
        y="Predicted Normalized RT"
      )+
      theme_classic2()+
      theme(
        legend.position = "none"
      )
    ggsave("./figures/rt/rt_gv_test_result_f1_gv_data.png",width =10,height = 10, dpi = 600, units = "cm")
    
    
    output_prediction_rt_test_gv_231109 %>% 
      mutate(
        gv_test = if_else(
          rt_pred >= lwr & rt_pred <= upr,1,0
        )
      ) %>% 
      filter(file_id == "F1") %>% 
      ggplot(aes(x=rt_norm, y=rt_pred))+
      geom_point(
        aes(color = as.character(gv_test)),
        size=1,alpha=0.5
      ) +
      geom_line(aes(y=lwr), color = "blue", linetype = "dashed")+
      geom_line(aes(y=upr), color = "blue", linetype = "dashed")+
      #geom_smooth(method=lm, se=TRUE)+
      scale_color_manual(labels = c("1","0"),values = c("black","blue"))+
      lims(x=c(0.2,1),y=c(0.2,1))+
      labs(
        x="Observed Normalized RT",
        y="Predicted Normalized RT"
      )+
      theme_classic2()+
      theme(
        legend.position = "none"
      )
    ggsave("./figures/rt/rt_gv_test_result_f1_gv_data.png",width =10,height = 10, dpi = 600, units = "cm")
    
    
    
    summaryfile_psm_for_testing3_pcc_rt_hfsm %>% 
      inner_join(
        output_prediction_rt_test_gv_231109 %>% 
          mutate(
            gv_test = if_else(
              rt_pred >= lwr & rt_pred <= upr,1,0
            )
          ) %>% 
          select(index,gv_test,rt_pred),
        join_by(index)
      ) %>% 
      filter(set == "F1") %>% 
      ggplot(aes(x=rt_norm))+
      geom_point(
        aes(y= rt_pred.x),
        size=1,color = "red",alpha=0.5
      )+
      geom_point(
        aes(y= rt_pred.y),
        size=1, color = "blue", alpha =0.5
      ) +
      geom_line(aes(y=lwr), color = "red", linetype = "dashed")+
      geom_line(aes(y=upr), color = "red", linetype = "dashed")+
      #geom_smooth(method=lm, se=TRUE)+
      #scale_color_manual(labels = c("1","0"),values = c("black","red"))+
      lims(x=c(0.2,1),y=c(0.2,1))+
      labs(
        x="Observed Normalized RT",
        y="Predicted Normalized RT"
      )+
      theme_classic2()+
      theme(
        legend.position = "none"
      )
    ggsave("./figures/rt/rt_gv_test_result_f1_arg_plus_gv.png",width =5,height = 5, dpi = 600, units = "cm")
    
    
    output_prediction_rt_training_gv_231109
    output_prediction_rt_test_gv_231109
    output_prediction_rt_training_arg_231109
    output_prediction_rt_test_arg_231109
  }
  
  ###GV figure histogram
  {
    output_prediction_rt_gvtest_result %>% 
      arrange(spec_idx) %>% 
      mutate(
        rt_deviation = rt_pred-rt_pred_original
      ) -> output_prediction_rt_gvtest_result
    
    summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
      bind_cols(
        output_prediction_rt_gvtest_result %>% select(rt_pred) %>% 
          rename(rt_pred_gv = rt_pred)
      ) -> summaryfile_psm_for_testing3_pcc_rt_hfsm2_rtgvtest
    
    summaryfile_psm_for_testing3_pcc_rt_hfsm2_rtgvtest %>% 
      mutate(
        avg_lwr = fit-lwr,
        avg_upr = upr-fit
      ) %>% 
      summarize(
        avg_lwr = mean(avg_lwr),
        avg_upr = mean(avg_upr)
      ) #0.0476
    
    summaryfile_psm_for_testing3_pcc_rt_hfsm2_rtgvtest %>% filter(p2gv_check == 1) %>% print(width =Inf)
    
    summaryfile_psm_for_testing3_pcc_rt_hfsm2_rtgvtest %>% 
      mutate(
        rt_deviation_gv = rt_pred_gv-rt_pred,
        rt_test_gv = case_when(
          rt_deviation_gv > 0 & rt_pred_gv > upr ~ "increase",
          rt_deviation_gv < 0 & rt_pred_gv < lwr ~ "decrease",
          .default = "none"
        )
      ) %>% 
      summarize(
        .by = rt_test_gv,
        count = n()
      )#63
    
    summaryfile_psm_for_testing3_pcc_rt_hfsm2_rtgvtest %>% 
      mutate(
        rt_deviation_gv = rt_pred_gv-rt_pred,
        rt_test_gv = case_when(
          rt_deviation_gv > 0 & rt_pred_gv > upr ~ "increase",
          rt_deviation_gv < 0 & rt_pred_gv < lwr ~ "decrease",
          .default = "none"
        )
      ) %>% 
      summarize(
        avg = mean(rt_deviation_gv)
      )#0.105
    
    summaryfile_psm_for_testing3_pcc_rt_hfsm2_rtgvtest %>% 
      mutate(
        rt_deviation_gv = rt_pred_gv-rt_pred,
        rt_test_gv = case_when(
          rt_deviation_gv > 0 & rt_pred_gv > upr ~ "increase",
          rt_deviation_gv < 0 & rt_pred_gv < lwr ~ "decrease",
          .default = "none"
        )
      ) %>% 
      ggplot(aes(x=rt_pred,y=rt_pred_gv,))+
      geom_point(
        aes(color = rt_test_gv),alpha =0.5,size = 0.25
      )+
      geom_abline(
        linewidth = 0.25
      )+
      # geom_abline(intercept = 0.0476,linetype = "dashed", color="red")+
      # geom_abline(intercept = -0.0476,linetype = "dashed", color="red")+
      geom_smooth(method = "lm",linewidth = 0.25)+
      labs(
        x="Predicted RT as Nt-arginylation",
        y="Predicted RT as Nt-Gly-Val"
      )+
      scale_color_manual(
        labels = c(expression(Increased),expression(" "*Under~Delta*t[95*"%"])),
        values = c("black","red")
      )+
      theme_classic2()+
      theme(
        legend.position = "none",
        #axis.title.x = element_blank(),
        axis.title = element_text(color = "black",size = 7),
        #axis.text.x = element_blank(),
        axis.text = element_text(color = "black",size = 7),
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm")
      )
    ggsave("./figures/2_09_rt_gv_prediction.png",width =4.5,height = 4.5, dpi = 600, units = "cm")
    
    summaryfile_psm_for_testing3_pcc_rt_hfsm2_rtgvtest %>% 
      mutate(
        rt_deviation_gv = rt_pred_gv-rt_pred,
        rt_test_gv = case_when(
          rt_deviation_gv > 0 & rt_pred_gv > upr ~ "increase",
          rt_deviation_gv < 0 & rt_pred_gv < lwr ~ "decrease",
          .default = "none"
        )
      ) %>% 
      filter(rt_test_gv=="none") %>% 
      ggplot(
        aes(x=rt_norm)
      )+
      geom_histogram(
        aes(y=after_stat(density)),
        fill="white",
        color="black"
      )+
      geom_density()+
      scale_x_continuous(
        breaks=seq(0,1,0.1)
      )+
      labs(
        x="Normalized Observed RT",
        y="PSM Count"
      )+
      theme_classic2()+
      theme(
        axis.text = element_text(color = "black")
      )
    ggsave("./figures/2_04_rt_gv_prediction3.png",width =10,height = 10, dpi = 600, units = "cm")
    
    
    summaryfile_psm_for_testing3_pcc_rt_hfsm2_rtgvtest %>% 
      filter(p2gv_check == 1) %>% 
      mutate(
        rt_deviation_gv = rt_pred_gv-rt_pred,
        rt_test_gv = case_when(
          rt_deviation_gv > 0 & rt_pred_gv > upr ~ "increase",
          rt_deviation_gv < 0 & rt_pred_gv < lwr ~ "decrease",
          .default = "none"
        )
      ) %>% 
      ggplot(aes(x=rt_norm,y=rt_pred,))+
      geom_point(
        aes(color = rt_test_gv),alpha =0.5
      )+
      geom_abline()+
      # geom_abline(intercept = 0.0476,linetype = "dashed", color="red")+
      # geom_abline(intercept = -0.0476,linetype = "dashed", color="red")+
      geom_smooth(method = "lm")+
      labs(
        x="Observed Normalized RT",
        y="Predicted Normalized RT"
      )+
      scale_color_manual(
        labels = c(expression(Increased),expression(" "*Under~Delta*t[95*"%"])),
        values = c("black","red")
      )+
      theme_classic2()+
      theme(
        legend.position = "none",
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_blank()
      )
    ggsave("./figures/2_04_rt_gv_prediction.png",width =10,height = 10, dpi = 600, units = "cm")
    
    summaryfile_psm_for_testing3_pcc_rt_hfsm2_rtgvtest %>% 
      filter(p2gv_check == 1 & rt_check == 1) %>% 
      select(p5p1,rt_norm)
      # ggplot(aes(x=rt_pred,y=rt_pred_gv,))+
      # geom_point(
      #   aes(color = rt_test_gv),alpha =0.5
      # )+
      # geom_abline()+
      # # geom_abline(intercept = 0.0476,linetype = "dashed", color="red")+
      # # geom_abline(intercept = -0.0476,linetype = "dashed", color="red")+
      # geom_smooth(method = "lm")+
      # labs(
      #   x="Predicted Normalized RT as Nt-arginylation",
      #   y="Predicted Normalized RT as Mass Ambiguity"
      # )+
      # scale_color_manual(
      #   labels = c(expression(Increased),expression(" "*Under~Delta*t[95*"%"])),
      #   values = c("black","red")
      # )+
      # theme_classic2()+
      # theme(
      #   legend.position = "top",
      #   axis.title = element_text(size = 10),
      #   legend.text = element_text(size = 10),
      #   legend.title = element_blank()
      # )
    
    RColorBrewer::display.brewer.all()
    RColorBrewer::display.brewer.all()
    RColorBrewer::brewer.pal(8,"Set2")
  }
  #2_04_rt_gv_prediction.png
  
  ###chapter SMERT
  ##SMERT cases
  {
    summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
      filter(
        p1r_check == 0 & p2gv_check == 0 & rt_check == 1 & pcc_check == 1 & b_y_ttest_check ==0
      ) %>%
      write_tsv(
        "./export/failed_only_smert.tsv"
      )
  }
  #failed_only_smert.tsv
  
  ##SMERT MASS functions.. deprecated
  {
    summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
      filter(str_detect(p3,"^R")) %>%
      mutate(
        p2 = str_extract(p3,"..$")
      ) %>% 
      mutate(
        mod = Peptides::mw(p2,monoisotopic = TRUE)-MonoisotopicMass(list(H=2,O=1)),
        mod_acetyl = mod + MonoisotopicMass(list(H=2,C=2,O=1)),
        mod_d3acetyl = mod + MonoisotopicMass(list(H=-1,x=3,C=2,O=1),isotopes = list(x = 2.01410177812)),
        mod_minus_d3arg = mod- MonoisotopicMass(list(H=11,x=3,C=8,O=2,N=4),isotopes = list(x = 2.01410177812)),
        mod_acetyl_minus_d3arg = mod_acetyl- MonoisotopicMass(list(H=11,x=3,C=8,O=2,N=4),isotopes = list(x = 2.01410177812)),
        mod_d3acetyl_minus_d3arg = mod_d3acetyl- MonoisotopicMass(list(H=11,x=3,C=8,O=2,N=4),isotopes = list(x = 2.01410177812))
      ) %>% 
      ggplot(
        aes(x=mod)
      )+
      geom_histogram(
        binwidth = 1
      )+
      scale_y_continuous(expand = c(0, 0))+
      labs(
        x="Mass (Da)",
        y="Count"
      )+
      theme_classic2()
    ggsave("./figures/2_05_smert_mass_error.png",width =10,height = 10, dpi = 600, units = "cm")
      
      
    Peptides::
    
    Peptides::mw("VT")-MonoisotopicMass(list(H=2,O=1))
    Peptides::mw("GL")-MonoisotopicMass(list(H=2,O=1))+MonoisotopicMass(list(H=2,C=2,O=1))
    Peptides::mw("LS")-MonoisotopicMass(list(H=2,O=1))
    Peptides::mw("VG")-MonoisotopicMass(list(H=2,O=1))
    Peptides::mw("AG")-MonoisotopicMass(list(H=2,O=1))
    Peptides::mw("IA")-MonoisotopicMass(list(H=2,O=1))
    Peptides::mw("SL")-MonoisotopicMass(list(H=2,O=1))
    Peptides::mw("GV",monoisotopic = TRUE)-MonoisotopicMass(list(H=2,O=1))
    Peptides::mw("R",monoisotopic = TRUE)-MonoisotopicMass(list(H=2,O=1))
    mw("R",monoisotopic = TRUE)
    MonoisotopicMass()
  }
  
  ##bion yion error distribution
  {
    summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% write_tsv("./export/summaryfile_psm_for_testing3_pcc_rt_hfsm2.tsv")
    summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>%colnames
      filter(str_detect(p3,"^R")) %>% 
      gghistogram(
        x="Delta M in ppm"
      )
    
    summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>%
      filter(str_detect(p3,"^R")) %>% 
      summarize(
        .by = c(pcc_check,rt_check,b_y_ttest_check),
        count = n()
      )
      
    summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>%
      mutate(
        p3r=if_else(str_detect(p3,"^R"),TRUE,FALSE)
      ) %>% 
      ggplot(
        aes(x=p3r,y=`Delta M in ppm`,fill=p3r)
      )+
      geom_boxplot()+
      stat_compare_means(method = "t.test")
    
    summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>%
      mutate(
        p3r=if_else(str_detect(p3,"^R"),TRUE,FALSE)
      ) %>% 
      unnest(b_ion_errors) %>% 
      ggplot(
        aes(x=p3r,y=b_ion_errors,fill=p3r)
      )+
      geom_boxplot()+
      stat_compare_means(method = "t.test")
      
    
    summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>%
      mutate(
        p3r=if_else(str_detect(p3,"^R"),TRUE,FALSE)
      ) %>% 
      unnest(y_ion_errors) %>% 
      ggplot(
        aes(x=p3r,y=y_ion_errors,fill=p3r)
      )+
      geom_boxplot()+
      stat_compare_means(method = "t.test")
  }
  
  ###frag t-test result @summaryfile_psm_for_testing3_pcc_rt_hfsm 
  {
    ##target PLIN3 #index 259 ##0.000321
    {
      summaryfile_psm_for_testing3_pcc_rt_hfsm %>% 
        filter(gene=="PLIN3") %>% 
        print(width = Inf)
      
      ##PLIN3 + legends
      tibble(
        error = summaryfile_psm_for_testing3_pcc_rt_hfsm %>% 
          filter(gene=="PLIN3") %>% .$b_ion_errors %>% unlist(),name ="b-ion"
      ) %>% 
        bind_rows(
          tibble(
            error = summaryfile_psm_for_testing3_pcc_rt_hfsm %>% 
              filter(gene=="PLIN3") %>% .$y_ion_errors %>% unlist(),name ="y-ion"
          )
        ) %>% 
        gghistogram(
          x = "error", 
          add = "mean", rug = TRUE,
          fill = "name", palette = c("#377EB8", "#E41A1C")
        )
      ggsave("./figures/2_05_masserror_plin3_legend.png",width=7,height=7, units = "cm",dpi = 600)
      
      ##PLIN3
      tibble(
        error = summaryfile_psm_for_testing3_pcc_rt_hfsm %>% 
          filter(gene=="PLIN3") %>% .$b_ion_errors %>% unlist(),name ="b-ion"
      ) %>% 
        bind_rows(
          tibble(
            error = summaryfile_psm_for_testing3_pcc_rt_hfsm %>% 
              filter(gene=="PLIN3") %>% .$y_ion_errors %>% unlist(),name ="y-ion"
          )
        ) %>% 
        gghistogram(
          x = "error", 
          add = "mean", rug = TRUE,
          fill = "name", palette = c("#377EB8", "#E41A1C")
        )+
        labs(x="Mass error (Da)",y="Count")+
        scale_x_continuous(limits = c(-0.05,0.05))+
        theme(
          legend.position = "none"
        )
      ggsave("./figures/2_05_masserror_plin3.png",width=8,height=8, units = "cm",dpi = 600)
    }
    
    ##CALR #0.888
    {
      summaryfile_psm_for_testing3_pcc_rt_hfsm %>% 
        filter(gene=="CALR") %>% 
        print(width = Inf)
      
      ##CALR index 233
      tibble(
        error = summaryfile_psm_for_testing3_pcc_rt_hfsm %>% 
          filter(index==233) %>% .$b_ion_errors %>% unlist(),name ="b-ion"
      ) %>% 
        bind_rows(
          tibble(
            error = summaryfile_psm_for_testing3_pcc_rt_hfsm %>% 
              filter(index==233) %>% .$y_ion_errors %>% unlist(),name ="y-ion"
          )
        ) %>% 
        gghistogram(
          x = "error", 
          add = "mean", rug = TRUE,
          fill = "name", palette = c("#377EB8", "#E41A1C")
        )+
        labs(x="Mass error (Da)",y="Count")+
        scale_x_continuous(limits = c(-0.05,0.05))+
        theme(
          legend.position = "none"
        )
      ggsave("./figures/2_05_masserror_calr.png",width=8,height=8, units = "cm",dpi = 600)
    }
    
    #distribution 2d density
    {
      summaryfile_psm_for_testing3_pcc_rt_hfsm %>% 
        mutate(
          b_ion_errors = map_dbl(b_ion_errors,\(x) mean(x)),
          y_ion_errors = map_dbl(y_ion_errors,\(x) mean(x))
        ) %>% 
        ggplot(
          aes(x=b_ion_errors*1000,y=y_ion_errors*1000)
        )+
        geom_pointdensity(adjust=5,size = 0.25)+
        scale_color_viridis()+
        scale_y_continuous(
          limits = c(-50,50)
        )+
        scale_x_continuous(
          limits = c(-50,50)
        )+
        labs(
          x="b-ion Errors (mDa)",
          y="y-ion Errors (mDa)",
        )+
        theme_bw()+
        theme(
          legend.position = "none",
          #axis.title.x = element_blank(),
          axis.title = element_text(color = "black",size = 7),
          #axis.text.x = element_blank(),
          axis.text = element_text(color = "black",size = 7),
          plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm")
        )
      ggsave("./figures/2_smert_error_dist.png",width=4.5,height=4.5,dpi = 600,units="cm")
      
      summaryfile_psm_for_testing3_pcc_rt_hfsm %>% 
        mutate(
          b_ion_errors = map_dbl(b_ion_errors,\(x) mean(x)),
          y_ion_errors = map_dbl(y_ion_errors,\(x) mean(x))
        ) %>% 
        filter(
          b_y_ttest >=0.05
        ) %>% 
        select(b_ion_errors,y_ion_errors) %>% 
        ggplot(
          aes(x=b_ion_errors,y=y_ion_errors)
        )+
        geom_pointdensity(adjust=0.005)+
        scale_color_viridis()+
        scale_y_continuous(
          limits = c(-0.05,0.05)
        )+
        scale_x_continuous(
          limits = c(-0.05,0.05)
        )+
        labs(
          x="b-ion Errors (Da)",
          y="y-ion Errors (Da)",
        )+
        theme_bw()+
        theme(
          legend.position = "none",
          axis.text = element_text(color="black",size=10),
          axis.title = element_text(size=12)
        )
      ggsave("./figures/2_smert_error_dist_p005_over.png",width=10,height=10,dpi = 600,units="cm")
      
      summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
        mutate(
          b_ion_errors = map_dbl(b_ion_errors,\(x) mean(x)),
          y_ion_errors = map_dbl(y_ion_errors,\(x) mean(x))
        ) %>% 
        filter(
          b_y_ttest <0.05
        ) %>% 
        select(b_ion_errors,y_ion_errors) %>% 
        ggplot(
          aes(x=b_ion_errors,y=y_ion_errors)
        )+
        geom_pointdensity(adjust=0.005)+
        scale_color_viridis()+
        scale_y_continuous(
          limits = c(-0.05,0.05)
        )+
        scale_x_continuous(
          limits = c(-0.05,0.05)
        )+
        labs(
          x="b-ion Errors (Da)",
          y="y-ion Errors (Da)",
        )+
        theme_bw()+
        theme(
          legend.position = "none",
          axis.text = element_text(color="black",size=10),
          axis.title = element_text(size=12)
        )
      ggsave("./figures/2_smert_error_dist_p005_under.png",width=10,height=10,dpi = 600,units="cm")
      
      summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
        mutate(
          b_ion_errors = map_dbl(b_ion_errors,\(x) mean(x)),
          y_ion_errors = map_dbl(y_ion_errors,\(x) mean(x))
        ) %>% 
        filter(
          pcc_check ==1
        ) %>% 
        select(b_ion_errors,y_ion_errors) %>% 
        ggplot(
          aes(x=b_ion_errors,y=y_ion_errors)
        )+
        geom_pointdensity(adjust=0.005)+
        scale_color_viridis()+
        scale_y_continuous(
          limits = c(-0.05,0.05)
        )+
        scale_x_continuous(
          limits = c(-0.05,0.05)
        )+
        labs(
          x="b-ion Errors (Da)",
          y="y-ion Errors (Da)",
        )+
        theme_bw()+
        theme(
          legend.position = "none",
          axis.text = element_text(color="black",size=10),
          axis.title = element_text(size=12)
        )
      ggsave("./figures/2_smert_density_pcc.png",width=10,height=10,dpi = 600,units="cm")
      
      summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
        mutate(
          b_ion_errors = map_dbl(b_ion_errors,\(x) mean(x)),
          y_ion_errors = map_dbl(y_ion_errors,\(x) mean(x))
        ) %>% 
        filter(
          rt_check ==1
        ) %>% 
        select(b_ion_errors,y_ion_errors) %>% 
        ggplot(
          aes(x=b_ion_errors,y=y_ion_errors)
        )+
        geom_pointdensity(adjust=0.005)+
        scale_color_viridis()+
        scale_y_continuous(
          limits = c(-0.05,0.05)
        )+
        scale_x_continuous(
          limits = c(-0.05,0.05)
        )+
        labs(
          x="b-ion Errors (Da)",
          y="y-ion Errors (Da)",
        )+
        theme_bw()+
        theme(
          legend.position = "none",
          axis.text = element_text(color="black",size=10),
          axis.title = element_text(size=12)
        )
      ggsave("./figures/2_smert_density_rt.png",width=10,height=10,dpi = 600,units="cm")
      
      summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
        mutate(
          b_ion_errors = map_dbl(b_ion_errors,\(x) mean(x)),
          y_ion_errors = map_dbl(y_ion_errors,\(x) mean(x))
        ) %>% 
        filter(
          dl_score ==3
        ) %>% 
        select(b_ion_errors,y_ion_errors) %>% 
        ggplot(
          aes(x=b_ion_errors,y=y_ion_errors)
        )+
        geom_pointdensity(adjust=0.005)+
        scale_color_viridis()+
        scale_y_continuous(
          limits = c(-0.05,0.05)
        )+
        scale_x_continuous(
          limits = c(-0.05,0.05)
        )+
        labs(
          x="b-ion Errors (Da)",
          y="y-ion Errors (Da)",
        )+
        theme_bw()+
        theme(
          legend.position = "none",
          axis.text = element_text(color="black",size=10),
          axis.title = element_text(size=12)
        )
      ggsave("./figures/2_smert_density_trinity.png",width=10,height=10,dpi = 600,units="cm")
    }
    #2_smert_error_dist.png
    
    ##histograms
    {
      summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
        select(
          b_ion_errors
        ) %>% 
        unnest(b_ion_errors) %>% #pull(b_ion_errors) %>% sd() #-0.001370094 sd 0.01418807
        gghistogram(
          x= "b_ion_errors",
          add = "mean",
          binwidth = 0.001,
          fill = c("#377EB8"),
          size= 0.2
        )+
        labs(
          x="b-ion errors (Da)",
          y="Count"
        )+
        theme(
          axis.text = element_text(size = 7),
          axis.title = element_text(size = 9)
        )
      ggsave("./figures/2_smert_b_error_dist.png",width=10,height=8,dpi = 600,units="cm")
      RColorBrewer::brewer.pal(3,"Set1")
      
      summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
        select(
          y_ion_errors
        ) %>% 
        unnest(y_ion_errors) %>% #pull(y_ion_errors) %>% sd() #3.831874e-05 #0.008929471
        gghistogram(
          x= "y_ion_errors",
          add = "mean",
          binwidth = 0.001,
          fill = c("#E41A1C"),
          size= 0.2
        )+
        labs(
          x="y-ion errors (Da)",
          y="Count"
        )+
        theme(
          axis.text = element_text(size = 7),
          axis.title = element_text(size = 9)
        )
      ggsave("./figures/2_smert_y_error_dist.png",width=10,height=8,dpi = 600,units="cm")
      
      summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
        filter(pcc_check==0) %>% 
        select(
          b_ion_errors
        ) %>% 
        unnest(b_ion_errors) %>% 
        gghistogram(
          x= "b_ion_errors",
          add = "mean",
          binwidth = 0.001,
          fill = c("#377EB8"),
          size= 0.2
        )
      ggsave("./figures/2_smert_b_error_dist_pcc0.png",width=10,height=10,dpi = 600,units="cm")
      
      summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
        filter(pcc_check==1) %>% 
        select(
          b_ion_errors
        ) %>% 
        unnest(b_ion_errors) %>% 
        gghistogram(
          x= "b_ion_errors",
          add = "mean",
          binwidth = 0.001,
          fill = c("#377EB8"),
          size= 0.2
        )
      ggsave("./figures/2_smert_b_error_dist_pcc1.png",width=10,height=10,dpi = 600,units="cm")
      
      summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
        filter(rt_check==0) %>% 
        select(
          b_ion_errors
        ) %>% 
        unnest(b_ion_errors) %>% 
        gghistogram(
          x= "b_ion_errors",
          add = "mean",
          binwidth = 0.001,
          fill = c("#377EB8"),
          size= 0.2
        )
      ggsave("./figures/2_smert_b_error_dist_rt0.png",width=10,height=10,dpi = 600,units="cm")
      
      summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
        filter(rt_check==1) %>% 
        select(
          b_ion_errors
        ) %>% 
        unnest(b_ion_errors) %>% 
        gghistogram(
          x= "b_ion_errors",
          add = "mean",
          binwidth = 0.001,
          fill = c("#377EB8"),
          size= 0.2
        )
      ggsave("./figures/2_smert_b_error_dist_rt1.png",width=10,height=10,dpi = 600,units="cm")
      
      summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
        filter(pcc_check==0) %>% 
        select(
          y_ion_errors
        ) %>% 
        unnest(y_ion_errors) %>% 
        gghistogram(
          x= "y_ion_errors",
          add = "mean",
          binwidth = 0.001,
          fill = c("#E41A1C"),
          size= 0.2
        )
      ggsave("./figures/2_smert_y_error_dist_pcc0.png",width=10,height=10,dpi = 600,units="cm")
      
      summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
        filter(pcc_check==1) %>% 
        select(
          y_ion_errors
        ) %>% 
        unnest(y_ion_errors) %>% 
        gghistogram(
          x= "y_ion_errors",
          add = "mean",
          binwidth = 0.001,
          fill = c("#E41A1C"),
          size= 0.2
        )
      ggsave("./figures/2_smert_y_error_dist_pcc1.png",width=10,height=10,dpi = 600,units="cm")
      
      summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
        filter(rt_check==0) %>% 
        select(
          y_ion_errors
        ) %>% 
        unnest(y_ion_errors) %>% 
        gghistogram(
          x= "y_ion_errors",
          add = "mean",
          binwidth = 0.001,
          fill = c("#E41A1C"),
          size= 0.2
        )
      ggsave("./figures/2_smert_y_error_dist_rt0.png",width=10,height=10,dpi = 600,units="cm")
      summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
        filter(rt_check==1) %>% 
        select(
          y_ion_errors
        ) %>% 
        unnest(y_ion_errors) %>% 
        gghistogram(
          x= "y_ion_errors",
          add = "mean",
          binwidth = 0.001,
          fill = c("#E41A1C"),
          size= 0.2
        )
      ggsave("./figures/2_smert_y_error_dist_rt1.png",width=10,height=10,dpi = 600,units="cm")
    }
    #2_smert_error_dist.png
    
    ##summary
    {
      summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
        summarize(
          .by = b_y_ttest_check,
          count = n()
        )
      summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
        mutate(
          b_y_ttest_check = fct(as.character(b_y_ttest_check)),
          #b_y_ttest = -log(b_y_ttest)
        ) %>% 
        gghistogram(
          x = "b_y_ttest", 
          add = "median", rug = TRUE,
          fill = "b_y_ttest_check", palette = c("#9E0142", "#5E4FA2"),
          binwidth = 0.025,
          position = "stack",
          size = 0.25
        )+
        labs(x=expression(italic(P)-value),y="Count")+
        #scale_x_continuous(limits = c(-0.05,0.05))+
        theme(
          legend.position = "none",
          #axis.title.x = element_blank(),
          axis.title = element_text(color = "black",size = 7),
          #axis.text.x = element_blank(),
          axis.text = element_text(color = "black",size = 7),
          plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm")
        )
      ggsave("./figures/2_05_masserror_all.png",width=4.5,height=4.5, units = "cm",dpi = 600)
      
      summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
        unnest(b_ion_errors) %>% 
        summarize(
          .by = pcc_check,
          b_mean = mean(b_ion_errors),
          b_sd = sd(b_ion_errors)
        ) %>% 
        bind_cols(
          summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
            unnest(y_ion_errors) %>% 
            summarize(
              .by = pcc_check,
              y_mean = mean(y_ion_errors),
              y_sd = sd(y_ion_errors)
            ) %>% select(y_mean,y_sd)
        )->table1
      
      summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
        unnest(b_ion_errors) %>% 
        summarize(
          .by = rt_check,
          b_mean = mean(b_ion_errors),
          b_sd = sd(b_ion_errors)
        ) %>% 
        bind_cols(
          summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
            unnest(y_ion_errors) %>% 
            summarize(
              .by = rt_check,
              y_mean = mean(y_ion_errors),
              y_sd = sd(y_ion_errors)
            ) %>% select(y_mean,y_sd)
        )->table2
      
      summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
        unnest(b_ion_errors) %>% 
        summarize(
          .by = b_y_ttest_check,
          b_mean = mean(b_ion_errors),
          b_sd = sd(b_ion_errors)
        ) %>% 
        bind_cols(
          summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
            unnest(y_ion_errors) %>% 
            summarize(
              .by = b_y_ttest_check,
              y_mean = mean(y_ion_errors),
              y_sd = sd(y_ion_errors)
            ) %>% select(y_mean,y_sd)
        )->table3
      
      bind_rows(
        table1,
        table2,
        table3
      ) %>% 
        mutate(
          result = case_when(
            pcc_check == 1 ~ "pcc_passed",
            pcc_check == 0 ~ "pcc_failed",
            rt_check == 1 ~ "rt_passed",
            rt_check == 0 ~ "rt_failed",
            b_y_ttest_check == 1 ~ "smert_passed",
            b_y_ttest_check == 0 ~ "smert_failed"
          )
        ) %>% 
        select(result, b_mean,b_sd,y_mean,y_sd) -> table_smert_result
      remove(table1,table2, table3)
      table_smert_result %>% 
        pivot_longer(
          cols = c(b_mean,y_mean)
        ) %>% 
        filter(
          str_detect(result, "passed")
        ) %>% 
        mutate(
          sd = if_else(str_detect(name,"^b"),b_sd,y_sd)
        ) %>% 
        ggplot(
          aes(x=result, y=value,fill=name)
        )+
        geom_bar(
          stat="identity",
          position = "dodge"
        )+
        geom_hline(
          aes(yintercept = 0)
        )+
        scale_y_continuous(
          labels = \(x) x*1000
        )+
        scale_x_discrete(
          labels = c("MS2 model\n passed", "RT model\n passed", "SMERT\n passed")
        )+
        scale_fill_discrete(
          labels = c("b-ion","y-ion")
        )+
        labs(
          y="Average Mass Error (mDa)"
        )+
        theme_bw()+
        theme(
          axis.title.x = element_blank(),
          legend.title = element_blank(),
          legend.position = "top",
          axis.text = element_text(color = "black")
        )
      ggsave("./figures/2_05_masserror_mean.png",width=10,height=10, units = "cm",dpi = 600)
      
      summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
        unnest(
          b_ion_errors
        ) %>% 
        filter(
          pcc_check == 1
        ) %>% 
        select(b_ion_errors) %>% 
        mutate(
          ion = "b"
        ) %>% 
        rename(error=b_ion_errors) %>% 
        bind_rows(
          summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
            unnest(
              y_ion_errors
            ) %>% 
            filter(
              pcc_check == 1
            ) %>% 
            select(y_ion_errors) %>% 
            mutate(
              ion = "y"
            ) %>% 
            rename(error=y_ion_errors)
        ) %>% 
        ggplot(
          aes(x=ion, y=error,fill=ion)
        )+
        geom_boxplot()+
        stat_compare_means(
          method = "t.test",
          hjust = -0.5,
          label.y = 0.06
        )+
        scale_y_continuous(limits = c(-0.05,0.06),breaks = seq(-0.05,0.05,0.025))+
        labs(
          y="Average Mass Error (mDa)"
        )+
        theme_classic2()+
        theme(
          legend.position = "none",
          axis.text = element_text(color="black",size = 12),
          axis.title.x = element_blank()
        )
      ggsave("./figures/2_05_pcc_error_ttest.png",width=10,height=10, units = "cm",dpi = 600)
      
      summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
        unnest(
          b_ion_errors
        ) %>% 
        filter(
          rt_check == 1
        ) %>% 
        select(b_ion_errors) %>% 
        mutate(
          ion = "b"
        ) %>% 
        rename(error=b_ion_errors) %>% 
        bind_rows(
          summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
            unnest(
              y_ion_errors
            ) %>% 
            filter(
              rt_check == 1
            ) %>% 
            select(y_ion_errors) %>% 
            mutate(
              ion = "y"
            ) %>% 
            rename(error=y_ion_errors)
        ) %>% 
        ggplot(
          aes(x=ion, y=error,fill=ion)
        )+
        geom_boxplot()+
        stat_compare_means(
          method = "t.test",
          hjust = -0.5,
          label.y = 0.06
        )+
        scale_y_continuous(limits = c(-0.05,0.06),breaks = seq(-0.05,0.05,0.025))+
        labs(
          y="Average Mass Error (mDa)"
        )+
        theme_classic2()+
        theme(
          legend.position = "none",
          axis.text = element_text(color="black",size = 12),
          axis.title.x = element_blank()
        )
      ggsave("./figures/2_05_rt_error_ttest.png",width=10,height=10, units = "cm",dpi = 600)
      
      summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
        unnest(
          b_ion_errors
        ) %>% 
        filter(
          b_y_ttest_check == 1
        ) %>% 
        select(b_ion_errors) %>% 
        mutate(
          ion = "b"
        ) %>% 
        rename(error=b_ion_errors) %>% 
        bind_rows(
          summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
            unnest(
              y_ion_errors
            ) %>% 
            filter(
              b_y_ttest_check == 1
            ) %>% 
            select(y_ion_errors) %>% 
            mutate(
              ion = "y"
            ) %>% 
            rename(error=y_ion_errors)
        ) %>% 
        ggplot(
          aes(x=ion, y=error,fill=ion)
        )+
        geom_boxplot()+
        stat_compare_means(
          method = "t.test",
          hjust = -0.5,
          label.y = 0.06
        )+
        scale_y_continuous(limits = c(-0.05,0.06),breaks = seq(-0.05,0.05,0.025))+
        labs(
          y="Average Mass Error (mDa)"
        )+
        theme_classic2()+
        theme(
          legend.position = "none",
          axis.text = element_text(color="black",size = 12),
          axis.title.x = element_blank()
        )
      ggsave("./figures/2_05_smert_error_ttest.png",width=10,height=10, units = "cm",dpi = 600)
      
    }
    #2_05_masserror_all.png
    #2_05_smert_error_ttest.png
    
    #SMERT LIST
    {
      arginylome_validation_matrix3_sequence_distinct %>% 
        filter(
          pcc_check == 1 & b_y_ttest_check == 0
        ) %>%
        left_join(
          reposition_arg_no_dups %>% 
            select(seq,p5p1),
          join_by(Sequence == seq)
        ) %>% 
        write_tsv(
          "./export/arginylome_validation_matrix3_failed_SMERT.tsv"
        )
      
      arginylome_validation_matrix3_sequence_distinct %>% 
        filter(
          pcc_check == 1 & b_y_ttest_check == 0
        ) %>%
        left_join(
          reposition_arg_no_dups %>% 
            select(seq,p5p1),
          join_by(Sequence == seq)
        ) -> arginylome_validation_matrix3_smert_only_failed
      
      arginylome_validation_matrix3_sequence_distinct_trinity %>% 
        left_join(
          reposition_arg_no_dups %>% 
            select(seq,p5p1),
          join_by(Sequence == seq)
        ) -> arginylome_validation_matrix3_sequence_distinct_trinity
        write_tsv(
          "./export/arginylome_validation_matrix3_failed_SMERT.tsv"
        )
      
      ##daqlogo analysis
      {
        logo_proteome <- prepareProteomeByFTP(source = NULL, species = "Homo sapiens", 
                                              fastaFile="UP000005640_9606_2301.fasta")
        formatSequence(seq = arginylome_validation_matrix3_smert_only_failed %>% 
                         pull(p5p1),
                       proteome = logo_proteome, upstreamOffset = 5,
                       downstreamOffset = 1) -> dag_seq
        
        
        bg_fisher <- buildBackgroundModel(dag_seq, background = "wholeProteome", 
                                          proteome = logo_proteome, testType = "fisher")
        bg_ztest <- buildBackgroundModel(dag_seq, background = "wholeProteome", 
                                         proteome = logo_proteome, testType = "ztest")
        
        png("./figures/2_smert_dag_heatmap.png", width = 10, height = 10, units = "cm", res = 600)
        dagHeatmap(
          testDAU(dag_seq, dagBackground = bg_fisher),
          type="diff",
          labels_col = c("P5","P4","P3","P2","P1"),
          display_numbers = TRUE
        )
        dev.off()
      }
      
      ##daqlogo analysis2 trinity
      {
        logo_proteome <- prepareProteomeByFTP(source = NULL, species = "Homo sapiens", 
                                              fastaFile="UP000005640_9606_2301.fasta")
        formatSequence(seq = arginylome_validation_matrix3_sequence_distinct_trinity %>% 
                         pull(p5p1),
                       proteome = logo_proteome, upstreamOffset = 5,
                       downstreamOffset = 1) -> dag_seq
        
        
        bg_fisher <- buildBackgroundModel(dag_seq, background = "wholeProteome", 
                                          proteome = logo_proteome, testType = "fisher")
        bg_ztest <- buildBackgroundModel(dag_seq, background = "wholeProteome", 
                                         proteome = logo_proteome, testType = "ztest")
        
        png("./figures/2_smert_trinity_dag_heatmap.png", width = 10, height = 10, units = "cm", res = 600)
        dagHeatmap(
          testDAU(dag_seq, dagBackground = bg_fisher),
          type="diff",
          labels_col = c("P5","P4","P3","P2","P1"),
          display_numbers = TRUE
        )
        dev.off()
      }
      
      
      reposition_arg_no_dups %>% 
        select(p5p1)
    }
    
    display.brewer.all()
    brewer.pal(11,"Spectral")
  }
  
  ###chapter 2 finish >> venn diagram
  {
    ##spectrum 1217
    {
      writeClipboard(summaryfile_psm_for_testing3_pcc_rt_hfsm %>% 
                       select(index,pcc_check,rt_check,b_y_ttest_check) %>% 
                       mutate(index = as.character(index)) %>% 
                       filter(pcc_check ==1) %>% .$index)
      
      writeClipboard(summaryfile_psm_for_testing3_pcc_rt_hfsm %>% 
                       select(index,pcc_check,rt_check,b_y_ttest_check) %>% 
                       mutate(index = as.character(index)) %>% 
                       filter(rt_check ==1) %>% .$index)
      
      writeClipboard(summaryfile_psm_for_testing3_pcc_rt_hfsm %>% 
                       select(index,pcc_check,rt_check,b_y_ttest_check) %>% 
                       mutate(index = as.character(index)) %>% 
                       filter(b_y_ttest_check ==1) %>% .$index)
    }
    
    ##spectrum 1217
    {
      writeClipboard(arginylome_validation_matrix2 %>% 
                       select(prot_pos_reposition,pcc_check,rt_check,b_y_ttest_check) %>% 
                       filter(pcc_check ==1) %>% .$prot_pos_reposition)
      
      writeClipboard(arginylome_validation_matrix2 %>% 
                       select(prot_pos_reposition,pcc_check,rt_check,b_y_ttest_check) %>% 
                       filter(rt_check ==1) %>% .$prot_pos_reposition)
      
      writeClipboard(arginylome_validation_matrix2 %>% 
                       select(prot_pos_reposition,pcc_check,rt_check,b_y_ttest_check) %>% 
                       filter(b_y_ttest_check ==1) %>% .$prot_pos_reposition)
    }
    
    summaryfile_psm_for_testing3_pcc_rt_hfsm %>% 
      select(index,pcc_check,rt_check,b_y_ttest_check) %>% 
      mutate(index = as.character(index)) %>% 
      filter(
        pcc_check == 1|
          rt_check == 1|
          b_y_ttest_check ==1
      ) %>% .$index %>% length()
    
    
  }
  
  ###Diagnostic ions
  #3_05_diagnostic_tests.png
  ##6. diagnostic ions in training set
  {
    install.packages("ggmosaic")
    library(ggmosaic)
    read_rds(
      "./export/psm_563516.rds"
    ) -> psm_563516
    
    psm_563516 %>% 
      select(
        index,enzyme,artefact_arg,diagnostic_ion_int_relative
      ) %>% 
      mutate(
        diagnostic_ion_presence = if_else(is.na(diagnostic_ion_int_relative),FALSE,TRUE)
      ) %>% 
      summarise(
        .by = c(artefact_arg,enzyme,diagnostic_ion_presence),
        count = n()
      ) %>% 
      write_tsv(
        "./export/3_06_diagnostic_ion_summary.tsv"
      )
    filter(
      enzyme == "trypsin"
    ) %>% 
      pull(count) %>% 
      matrix(nrow =2, ncol =2, dimnames = list(c("none","artifact_arg"),c("none","diag"))) %>% 
      fisher.test()
    
    psm_563516 %>% 
      select(
        index,enzyme,artefact_arg,diagnostic_ion_int_relative
      ) %>% 
      mutate(
        diagnostic_ion_presence = if_else(is.na(diagnostic_ion_int_relative),FALSE,TRUE)
      ) %>% 
      filter(
        enzyme == "trypsin"
      ) %>% 
      ggplot()+
      geom_mosaic(
        aes(x=product(artefact_arg),fill = diagnostic_ion_presence)
      )+
      labs(
        x="Artifact Nt-arginyl peptide",
        y="Presence of diagnostic ion",
        fill = "Presence of\ndiagnostic ion"
      )+
      theme_classic2()+
      theme(
        #axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text=element_text(color ="black"),
        legend.position = "none"
      )
    ggsave("./figures/3_06_diagnostic_ion_trypsin.png",width = 8, height = 8, units = "cm",dpi=600)
    
    psm_563516 %>% 
      select(
        index,enzyme,artefact_arg,diagnostic_ion_int_relative
      ) %>% 
      mutate(
        diagnostic_ion_presence = if_else(is.na(diagnostic_ion_int_relative),FALSE,TRUE)
      ) %>% 
      filter(
        enzyme == "chymotrypsin"
      ) %>% 
      ggplot()+
      geom_mosaic(
        aes(x=product(artefact_arg),fill = diagnostic_ion_presence)
      )+
      labs(
        x="Artifact Nt-arginyl peptide",
        y="Presence of diagnostic ion",
        fill = "Presence of\ndiagnostic ion"
      )+
      theme_classic2()+
      theme(
        #axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text=element_text(color ="black"),
        legend.position = "none"
      )
    ggsave("./figures/3_06_diagnostic_ion_chymotrypsin.png",width = 8, height = 8, units = "cm",dpi=600)
    
    psm_563516 %>% 
      select(
        index,enzyme,artefact_arg,diagnostic_ion_int_relative,peak_table
      ) %>% 
      filter(artefact_arg == TRUE) %>% 
      mutate(
        peak_table_diag = map(
          peak_table,
          \(x) mutate(x,rel_i = i/max(i)) %>% filter(`m/z`>202.128331 & `m/z`<202.148331)
        )
      ) %>% 
      select(-peak_table) %>% 
      unnest(peak_table_diag) ->psm_563516_diag
    
    psm_563516_diag %>% 
      filter(enzyme == "trypsin") %>% 
      ggplot(
        aes(x=`m/z`)
      )+
      geom_histogram(
        aes(y=after_stat(density)),
        binwidth = 0.0001,
        fill="black",
        color="white",
        alpha = 0.5,
        linewidth = 0.1
      )+
      geom_density(
        color="red"
      )+
      geom_vline(
        aes(xintercept = 202.143331),linetype="dashed"
      )+
      geom_vline(
        aes(xintercept = 202.133331),linetype="dashed"
      )+
      scale_x_continuous(
        limits = c(202.128331,202.148331),expand=c(0,0),
        breaks= seq(202.123331,202.148331,0.005),
        labels = scales::number_format(accuracy = 0.000001)
        )+
      scale_y_continuous(expand=c(0,0))+
      labs(y="Density")+
      theme_bw()+
      theme(
        axis.text = element_text(color = "black")
      )
    ggsave(
      "./figures/2_06_diagnostic_ions_mz_trypsin.png",width =10, height =10,dpi = 600, units ="cm"
    )
    
    psm_563516_diag %>% 
      filter(enzyme == "chymotrypsin") %>% 
      ggplot(
        aes(x=`m/z`)
      )+
      geom_histogram(
        aes(y=after_stat(density)),
        binwidth = 0.0001,
        fill="black",
        color="white",
        alpha = 0.5,
        linewidth = 0.1
      )+
      geom_density(
        color="red"
      )+
      geom_vline(
        aes(xintercept = 202.143331),linetype="dashed"
      )+
      geom_vline(
        aes(xintercept = 202.133331),linetype="dashed"
      )+
      scale_x_continuous(
        limits = c(202.128331,202.148331),expand=c(0,0),
        breaks= seq(202.123331,202.148331,0.005),
        labels = scales::number_format(accuracy = 0.000001)
      )+
      scale_y_continuous(expand=c(0,0))+
      labs(y="Density")+
      theme_bw()+
      theme(
        axis.text = element_text(color = "black")
      )
    ggsave(
      "./figures/2_06_diagnostic_ions_mz_chymotrypsin.png",width =10, height =10,dpi = 600, units ="cm"
    )
  }
  #2_06_diagnostic_ions_mz_chymotrypsin.png
  #3_06_diagnostic_ion_trypsin.png
  #3_06_diagnostic_ion_chymotrypsin.png
  
  ##5. diagnostic ions_trinity
  #@summaryfile_psm_for_testing3_pcc_rt_hfsm2
  {
    summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
      pull(diagnostic_ion_mz)
    summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
      pull(diagnostic_ion_int_relative)
    
    summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
      mutate(
        trinity = if_else(dl_score == 3,1,0)
      ) %>% 
      mutate(
        trinity_failed = if_else(dl_score != 3,1,0)
      ) %>%
      select(
        diagnostic_ion_int_relative,trinity_failed,pcc_check,rt_check,b_y_ttest_check,trinity
      ) %>% 
      mutate(
        diagnostic_ion_int_relative = replace_na(diagnostic_ion_int_relative,0)
      ) %>% 
      pivot_longer(
        trinity_failed:trinity
      ) %>% 
      filter(value == 1) %>% 
      summarize(
        .by = name,
        median = median(diagnostic_ion_int_relative)
      )
    
    summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
      filter(pcc_check == 0 & rt_check == 1) %>% 
      arrange(
        desc(diagnostic_ion_int_relative)
      ) %>% view()
    
    summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>%
      filter(
        gene == "PRSS21"
      ) %>% 
      arrange(
        desc(PCC)
      ) %>%
      pull(index)
    
    summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
      mutate(
        trinity = if_else(dl_score == 3,1,0)
      ) %>% 
      mutate(
        trinity_failed = if_else(dl_score != 3,1,0)
      ) %>%
      select(
        diagnostic_ion_int_relative,trinity_failed,pcc_check,rt_check,b_y_ttest_check,trinity
      ) %>% 
      mutate(
        diagnostic_ion_int_relative = replace_na(diagnostic_ion_int_relative,0)
      ) %>% 
      pivot_longer(
        trinity_failed:trinity
      ) %>% 
      filter(value == 1) %>% 
      ggplot(
        aes(x=fct(name,levels = c("trinity_failed","pcc_check","rt_check","b_y_ttest_check","trinity")),y=diagnostic_ion_int_relative)
      )+
      geom_violin(
        scale = "width",
        aes(fill = fct(name,levels = c("trinity_failed","pcc_check","rt_check","b_y_ttest_check","trinity"))),
        linewidth = 0.25
      )+
      geom_boxplot(width = 0.2,outlier.size = 0.25,linewidth = 0.25)+
      scale_y_continuous(labels=label_percent())+
      scale_x_discrete(
        labels = c("None","PCC","RT","SMERT","All")
      )+
      scale_fill_manual(
        values = wes_palette("Darjeeling1",5,type="discrete")
      )+
      stat_compare_means(
        comparisons = list(c("pcc_check","trinity"),c("rt_check","trinity"),c("b_y_ttest_check","trinity")), label.y = c(0.9,0.8,0.7),
        p.adjust.methods="fdr",
        size = 2
      )+
      labs(
        x="Trinity Tests",y="Diagnostic Ion Intensity %"
      )+
      scale_y_continuous(
        breaks = c(0,0.25,0.5,0.75,1),
        limits = c(0,1)
      )+
      theme_classic2()+
      theme(
        legend.position = "none",
        #axis.title.x = element_blank(),
        axis.title = element_text(color = "black",size = 7),
        #axis.text.x = element_blank(),
        axis.text = element_text(color = "black",size = 7),
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm")
      )
    ggsave("./figures/3_05_diagnostic_tests.png",width = 4.5, height = 4.5, units = "cm",dpi=600)
    
    
  }
  #3_05_diagnostic_tests.png
  
  ##6. diagnostic ions in training set
  {
    install.packages("ggmosaic")
    library(ggmosaic)
    read_rds(
      "./export/psm_563516.rds"
    ) -> psm_563516
    
    psm_563516 %>% 
      select(
        index,enzyme,artefact_arg,diagnostic_ion_int_relative
      ) %>% 
      mutate(
        diagnostic_ion_presence = if_else(is.na(diagnostic_ion_int_relative),FALSE,TRUE)
      ) %>% 
      summarise(
        .by = c(artefact_arg,enzyme,diagnostic_ion_presence),
        count = n(),
        median = median(diagnostic_ion_int_relative,na.rm = TRUE)
      ) %>% 
      write_tsv(
        "./export/2_06_diagnostic_ions_summary.tsv"
      )
    
    psm_563516 %>% 
      select(
        index,enzyme,artefact_arg,diagnostic_ion_int_relative
      ) %>% 
      mutate(
        diagnostic_ion_presence = if_else(is.na(diagnostic_ion_int_relative),FALSE,TRUE)
      ) %>% 
      summarise(
        .by = c(artefact_arg,enzyme,diagnostic_ion_presence),
        count = n()
      ) %>% 
      filter(
        enzyme == "chymotrypsin"
      ) %>% 
      pull(count) %>% 
      matrix(nrow =2, ncol =2, dimnames = list(c("none","artifact_arg"),c("none","diag"))) %>% 
      fisher.test()
    
    psm_563516 %>% 
      select(
        index,enzyme,artefact_arg,diagnostic_ion_int_relative
      ) %>% 
      mutate(
        diagnostic_ion_presence = if_else(is.na(diagnostic_ion_int_relative),FALSE,TRUE)
      ) %>% 
      filter(
        enzyme == "trypsin"
      ) %>% 
      ggplot()+
      geom_mosaic(
        aes(x=product(artefact_arg),fill = diagnostic_ion_presence)
      )+
      labs(
        x="Artifact Nt-arginyl peptide",
        y="Presence of diagnostic ion",
        fill = "Presence of\ndiagnostic ion"
      )+
      theme_classic2()+
      theme(
        #axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text=element_text(color ="black"),
        legend.position = "none"
      )
    ggsave("./figures/3_06_diagnostic_ion_trypsin.png",width = 8, height = 8, units = "cm",dpi=600)
    
    psm_563516 %>% 
      select(
        index,enzyme,artefact_arg,diagnostic_ion_int_relative
      ) %>% 
      mutate(
        diagnostic_ion_presence = if_else(is.na(diagnostic_ion_int_relative),FALSE,TRUE)
      ) %>% 
      filter(
        enzyme == "chymotrypsin"
      ) %>% 
      ggplot()+
      geom_mosaic(
        aes(x=product(artefact_arg),fill = diagnostic_ion_presence)
      )+
      labs(
        x="Artifact Nt-arginyl peptide",
        y="Presence of diagnostic ion",
        fill = "Presence of\ndiagnostic ion"
      )+
      theme_classic2()+
      theme(
        #axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text=element_text(color ="black"),
        legend.position = "none"
      )
    ggsave("./figures/3_06_diagnostic_ion_chymotrypsin.png",width = 8, height = 8, units = "cm",dpi=600)
  }
  #3_06_diagnostic_ion_trypsin.png
  #3_06_diagnostic_ion_chymotrypsin.png
  
  
  ##7. diagnostic ions in training set: trypsin, chymotrypsin
  {
    psm_563516 %>% 
      mutate(
        diagnostic_ion_int_relative = if_else(is.na(diagnostic_ion_int_relative),0,diagnostic_ion_int_relative)
      ) %>% 
      filter(
        artefact_arg==TRUE
      ) %>% 
      select(
        enzyme, diagnostic_ion_int_relative
      ) %>% 
      mutate(
        set = "training"
      ) %>%
      ggplot(
        aes(y=diagnostic_ion_int_relative,x=enzyme,fill = enzyme)
      )+
      geom_boxplot()+
      scale_x_discrete(
        limits = c("trypsin","chymotrypsin"),
        labels = c("Trypsin","Chymotrypsin")
      )+
      scale_y_continuous(
        labels = \(x) scales::percent(x)
      )+
      labs(
        x="Protease",
        y="Diagnostic Ion Intensity (%)",
        fill = "Protease"
      )+
      theme_classic2()+
      theme(
        axis.text = element_text(color = "black"),
        legend.position = "none",
        legend.just = "left"
      )
    ggsave("./figures/3_07_diagnostic_ion_int_protease.png",width = 10, height = 10, units = "cm",dpi=600)
  }
  #3_07_diagnostic_ion_int_protease.png
  
  ##1. VG GV
  {
    summaryfile_psm_for_testing3_pcc_rt_hfsm %>% 
      filter(pcc_check == 1 | p2gv_check ==1) %>% summarise(.by = c(pcc_check, p2gv_check),count = n())
    summaryfile_psm_for_testing3_pcc_rt_hfsm %>% 
      filter(pcc_check == 1 | p2gv_check ==1) %>%
      mutate(
        p2gv_check = fct(as.character(p2gv_check),levels = c("1","0"))
      ) %>% 
      ggplot(
        aes(y=PCC,x = p2gv_check)
      )+
      geom_boxplot()+
      geom_hline(
        yintercept = cutoff_pcc,
        linetype = "dashed",
        color = "red"
      )+
      scale_x_discrete(
        labels = c("P2-P1 = GV or VG","MS2 model passed")
      )+theme_classic2()+
      theme(
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text = element_text(color = "black")
      )
    ggsave("./figures/3_01_gv_vg_pcc.png",width =10, height = 10,units = "cm",dpi = 600)
    
    summaryfile_psm_for_testing3_pcc_rt_hfsm
    summaryfile_psm_for_testing3_pcc_rt_hfsm %>% 
      filter(rt_check == 1 | p2gv_check ==1) %>% summarise(.by = c(rt_check, p2gv_check),count = n())
    
    summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
      filter(rt_check == 1 | p2gv_check ==1) %>% 
      mutate(
        p2gv_check = fct(as.character(p2gv_check),levels = c("1","0")),
        rt_dev = rt_pred-rt_norm
      ) %>% 
      ggplot(
        aes(y=rt_dev,x = p2gv_check)
      )+
      geom_boxplot()+
      geom_hline(
        yintercept = 0,
        linetype = "dashed",
        color = "red"
      )+
      scale_x_discrete(
        labels = c("P2-P1 = GV or VG","RT model passed")
      )+
      labs(y="Normalized RT Deviation (Predicted RT - Observed RT)")+
      theme_classic2()+
      theme(
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10),
        axis.text = element_text(color = "black")
      )
    ggsave("./figures/3_01_gv_vg_rt.png",width =10, height = 10,units = "cm",dpi = 600)
  }
  
}
#######

####
####Chapter 3 Nt-arginylome ####
{
  ##1. Venn diagram
  {
    summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
      filter(p1r_check ==1 | p2gv_check ==1) %>% 
      pull(index)
    
    ggVennDiagram(
      list(
        # A=summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
        #   filter(p1r_check ==1 | p2gv_check ==1) %>% 
        #   pull(index),
        B=summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered %>% 
          filter(pcc_check == 1) %>% 
          pull(index),
        C=summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered %>% 
          filter(rt_check == 1) %>% 
          pull(index),
        D=summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered %>% 
          filter(b_y_ttest_check == 1) %>% 
          pull(index)
      ),
      category.names = c("MS2","RT","SMERT"),
      set_size = 4
    )+
      theme_void() +
      scale_fill_distiller(palette = "Reds",direction=1)+
      scale_color_manual(
        values = c("black","black","black")
      )+
      coord_sf(clip="off")+
      theme(
        legend.position = "none",
        plot.margin = margin(rep(0.5,4),unit="cm")
      )
    ggsave("./figures/3_01_venn.png",width = 10, height = 10, units = "cm",dpi=600)
    
    arginylome_validation_matrix3_sequence_distinct %>% 
      filter(pcc_check == 1) %>% pull(prot_pos_master)
    ggVennDiagram(
      list(
        # A=summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
        #   filter(p1r_check ==1 | p2gv_check ==1) %>% 
        #   pull(index),
        B=arginylome_validation_matrix3_sequence_distinct %>% 
          filter(pcc_check == 1) %>% pull(prot_pos_master),
        C=arginylome_validation_matrix3_sequence_distinct %>% 
          filter(rt_check == 1) %>% pull(prot_pos_master),
        D=arginylome_validation_matrix3_sequence_distinct %>% 
          filter(b_y_ttest_check == 1) %>% pull(prot_pos_master)
      ),
      category.names = c("MS2","RT","SMERT"),
      set_size = 4
    )+
      theme_void() +
      scale_fill_distiller(palette = "Reds",direction=1)+
      scale_color_manual(
        values = c("black","black","black")
      )+
      coord_sf(clip="off")+
      theme(
        legend.position = "none",
        plot.margin = margin(rep(0.5,4),unit="cm")
      )
    ggsave("./figures/3_01_venn_sites.png",width = 10, height = 10, units = "cm",dpi=600)
  }
  
  #removed arg in drug model
  #@summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered
  #@arginylome_validation_matrix3_sequence_distinct
  {
    sum(
      sum(arginylome_validation_matrix3_sequence_distinct$MOCK), ##205
      sum(arginylome_validation_matrix3_sequence_distinct$MG132), ##370
      sum(arginylome_validation_matrix3_sequence_distinct$MGTG) ##533
    )
    sum(
      sum(arginylome_validation_matrix3_sequence_distinct_trinity$MOCK), ##27
      sum(arginylome_validation_matrix3_sequence_distinct_trinity$MG132), ##184
      sum(arginylome_validation_matrix3_sequence_distinct_trinity$MGTG) ##340
    )
    
    arginylome_validation_matrix3_sequence_distinct %>% 
      dplyr::select(trinity,MOCK,MG132,MGTG) %>% 
      pivot_longer(c(MOCK,MG132,MGTG)) %>% 
      summarise(
        .by = c(name,trinity),
        count = sum(value)
      ) %>% 
      ggplot(
        aes(x=fct(name,level = c("MOCK","MG132","MGTG")),y=count,fill=trinity)
      )+
      geom_bar(stat="identity",color="black")+
      geom_text(aes(label=count),position = position_stack(vjust=0.5),size=3)+
      scale_x_discrete(
        label = c("MOCK","MG132","MG132 + TG")
      )+
      scale_fill_manual(
        values = c(wesanderson::wes_palette("Royal1")[2],wesanderson::wes_palette("Royal1")[3])
      )+
      labs(x="Experiment",y="PSM Count",fill="TrinityNmod")+
      theme_classic2()+
      theme(
        legend.position = "none",
        axis.text = element_text(color="black")
      )
    ggsave("./figures/3_02_removed_arg_drug_model.png",width=8,height=8, units = "cm",dpi = 600)
    
    arginylome_validation_matrix3_sequence_distinct %>% 
      dplyr::select(trinity,MOCK,MG132,MGTG) %>% 
      pivot_longer(c(MOCK,MG132,MGTG)) %>% 
      summarise(
        .by = c(name,trinity),
        count = sum(value)
      ) %>% 
      ggplot(
        aes(x=fct(name,level = c("MOCK","MG132","MGTG")),y=count,fill=trinity)
      )+
      geom_bar(stat="identity",color="black")+
      geom_text(aes(label=count),position = position_stack(vjust=0.5),size=3)+
      scale_x_discrete(
        label = c("MOCK","MG132","MG132 + TG")
      )+
      scale_fill_manual(
        values = c(wesanderson::wes_palette("Royal1")[2],wesanderson::wes_palette("Royal1")[3])
      )+
      labs(x="Experiment",y="Count",fill="TrinityNmod")+
      theme_classic2()+
      theme(
        axis.text = element_text(color="black")
      )
    ggsave("./figures/3_02_removed_arg_drug_model_legend.png",width=8,height=8, units = "cm",dpi = 600)
    
  }
  #3_02_removed_arg_drug_model.png
  
  ##4. logo
  {
    library(ggseqlogo)
    
    {
      arginylome_validation_matrix3_sequence_distinct_trinity %>% select(p5p5prime) %>%
        ggplot()+geom_logo(data = arginylome_validation_matrix3_sequence_distinct_trinity$p5p5prime)+
        theme_logo()+
        theme(
          legend.position = "none"
        )->fig_logo
      
      fig_logo$scales$scales[[1]]$labels <- c("P5","P4","P3","P2","P1","P1'","P2'","P3'","P4'","P5'")
      
      fig_logo
      ggsave("./figures/3_04_logo_trinity.png",width=10,height=10,dpi = 600,units="cm")
    }
    #3_04_logo_trinity.png
    
    #trinity DAU
    {
      logo_proteome <- prepareProteomeByFTP(source = NULL, species = "Homo sapiens", 
                                            fastaFile="UP000005640_9606_2301.fasta")
      formatSequence(seq = arginylome_validation_matrix3_sequence_distinct_trinity %>% select(p5p5prime) %>% mutate(p5p1 = str_sub(p5p5prime,1,5)) %>% pull(p5p1),
                     proteome = logo_proteome, upstreamOffset = 5,
                     downstreamOffset = 1) -> dag_seq
      
      
      bg_fisher <- buildBackgroundModel(dag_seq, background = "wholeProteome", 
                                         proteome = logo_proteome, testType = "fisher")
      #bg_ztest2 <- buildBackgroundModel(dag_seq2, background = "wholeProteome", 
      #proteome = logo_proteome, testType = "ztest")
      
      png("./figures/3_02_dag_heatmap_trinity.png", width = 10, height = 10, units = "cm", res = 600)
      dagHeatmap(
        testDAU(dag_seq, dagBackground = bg_fisher),
        type="diff",
        labels_col = c("P5","P4","P3","P2","P1"),
        display_numbers = TRUE
      )
      dev.off()
    }
    #3_02_dag_heatmap_trinity.png
    
    {
      arginylome_validation_matrix3_sequence_distinct %>% filter(trinity == "Failed")%>% select(p5p5prime) -> logo_p5p1_nontrinity
      arginylome_validation_matrix3_sequence_distinct %>% filter(trinity == "Failed")%>% select(p5p5prime) %>% 
        ggplot()+geom_logo(data = logo_p5p1_nontrinity$p5p5prime)+
        theme_logo()+
        theme(
          legend.position = "none"
        )->fig_logo
      
      fig_logo$scales$scales[[1]]$labels <- c("P5","P4","P3","P2","P1","P1'","P2'","P3'","P4'","P5'")
      
      fig_logo
      ggsave("./figures/3_04_logo_nontrinity.png",width=10,height=10,dpi = 600,units="cm")
    }
    #3_04_logo_nontrinity.png
    
    #trinity DAU
    {
      logo_proteome <- prepareProteomeByFTP(source = NULL, species = "Homo sapiens", 
                                            fastaFile="UP000005640_9606_2301.fasta")
      formatSequence(seq = arginylome_validation_matrix3_sequence_distinct %>% filter(trinity == "Failed") %>% select(p5p5prime) %>% mutate(p5p1 = str_sub(p5p5prime,1,5)) %>% pull(p5p1),
                     proteome = logo_proteome, upstreamOffset = 5,
                     downstreamOffset = 1) -> dag_seq
      
      
      bg_fisher <- buildBackgroundModel(dag_seq, background = "wholeProteome", 
                                        proteome = logo_proteome, testType = "fisher")
      #bg_ztest2 <- buildBackgroundModel(dag_seq2, background = "wholeProteome", 
      #proteome = logo_proteome, testType = "ztest")
      
      png("./figures/3_02_dag_heatmap_nontrinity.png", width = 10, height = 10, units = "cm", res = 600)
      dagHeatmap(
        testDAU(dag_seq, dagBackground = bg_fisher),
        type="diff",
        labels_col = c("P5","P4","P3","P2","P1"),
        display_numbers = TRUE
      )
      dev.off()
    }
    #3_02_dag_heatmap_trinity.png
  }
  
  ##2. subcellular localization
  {
    
    #subcellular localization all
    {
      arginylome_validation_matrix3_sequence_distinct %>% 
        mutate(
          subloc=str_split_i(Localizations,"\\|",1)
        ) %>%
        summarise(
          .by=c(subloc,trinity),
          count =n()
        ) %>% 
        pivot_wider(names_from = trinity,id_cols = subloc,values_from = count) %>% 
        mutate(
          subloc = fct(subloc,level = c("Extracellular","Cell membrane","Lysosome/Vacuole","Golgi apparatus","Endoplasmic reticulum","Cytoplasm","Nucleus","Mitochondrion")),
          Failed = replace_na(Failed,0),
          sum = Failed+Passed,
          perc = percent(Passed/sum,accuracy=1),
          text = paste0(perc)
        ) %>% arrange(subloc) %>% .$text -> fig_3_02_text
      
      arginylome_validation_matrix3_sequence_distinct %>% 
        mutate(
          subloc=str_split_i(Localizations,"\\|",1)
        ) %>%
        summarise(
          .by=c(subloc,trinity),
          count =n()
        ) %>% 
        ggplot(
          aes(x=fct(subloc,level = c("Extracellular","Cell membrane","Lysosome/Vacuole","Golgi apparatus","Endoplasmic reticulum","Cytoplasm","Nucleus","Mitochondrion")),y=count,fill=trinity)
        )+
        geom_bar(stat="identity",color="black",linewidth = 0.1)+
        annotate(
          geom = "text",
          x=1:8,
          y=220,
          label = fig_3_02_text,
          size = 3
        )+
        #geom_text(aes(label=count),position = position_stack(vjust=0.5),size=3)+
        scale_fill_manual(
          values = c(wesanderson::wes_palette("Royal1")[2],wesanderson::wes_palette("Royal1")[3])
        )+
        labs(x="Subcellular Localization",y="Site Count",fill="TrinityNmod")+
        theme_classic2()+
        theme(
          legend.position = "none",
          axis.text = element_text(color="black",size=8),
          axis.text.x = element_text(angle=45,hjust=1),
          axis.title = element_text(color="black",size = 10)
        )
      ggsave("./figures/3_02_removed_arg_subloc.png",width=8,height=8, units = "cm",dpi = 600)
    }
    #3_02_removed_arg_subloc.png
    
    ####MOCK
    { 
      arginylome_validation_matrix3_sequence_distinct %>% 
        filter(MOCK>0) %>% 
        mutate(
          subloc=str_split_i(Localizations,"\\|",1)
        ) %>%
        summarise(
          .by=c(subloc,trinity),
          count =n()
        ) %>% 
        pivot_wider(names_from = trinity,id_cols = subloc,values_from = count) %>% 
        mutate(
          subloc = fct(subloc,level = c("Extracellular","Cell membrane","Lysosome/Vacuole","Golgi apparatus","Endoplasmic reticulum","Cytoplasm","Nucleus","Mitochondrion")),
          Failed = replace_na(Failed,0),
          Passed = replace_na(Passed,0),
          sum = Failed+Passed,
          perc = percent(Passed/sum,accuracy=0.1),
          text = paste0(perc)
        ) %>% arrange(subloc) %>% .$text -> fig_3_02_text_mock
      
      arginylome_validation_matrix3_sequence_distinct %>% 
        filter(MOCK>0) %>% 
        mutate(
          subloc=str_split_i(Localizations,"\\|",1)
        ) %>%
        summarise(
          .by=c(subloc,trinity),
          count =n()
        ) %>% 
        ggplot(
          aes(x=fct(subloc,level = c("Extracellular","Cell membrane","Lysosome/Vacuole","Golgi apparatus","Endoplasmic reticulum","Cytoplasm","Nucleus","Mitochondrion")),y=count,fill=trinity)
        )+
        geom_bar(stat="identity",color="black",linewidth = 0.1)+
        annotate(
          geom = "text",
          x=1:5,
          y=65,
          label = fig_3_02_text_mock,
          size = 2.3
        )+
        #geom_text(aes(label=count),position = position_stack(vjust=0.5),size=3)+
        scale_fill_manual(
          values = c(wesanderson::wes_palette("Royal1")[2],wesanderson::wes_palette("Royal1")[3])
        )+
        labs(x="Subcellular Localization",y="Count",fill="TrinityNmod")+
        theme_classic2()+
        theme(
          legend.position = "none",
          axis.text = element_text(color="black",size=8),
          axis.text.x = element_text(angle=45,hjust=1),
          axis.title = element_text(color="black",size = 10)
        )
      ggsave("./figures/3_02_removed_arg_subloc_mock.png",width=8,height=8, units = "cm",dpi = 600)
    }
    #3_02_removed_arg_subloc_mock.png
    
    ####MG132
    {
      arginylome_validation_matrix3_sequence_distinct %>% 
        filter(MG132>0) %>% 
        mutate(
          subloc=str_split_i(Localizations,"\\|",1)
        ) %>%
        summarise(
          .by=c(subloc,trinity),
          count =n()
        ) %>% 
        pivot_wider(names_from = trinity,id_cols = subloc,values_from = count) %>% 
        mutate(
          subloc = fct(subloc,level = c("Extracellular","Cell membrane","Lysosome/Vacuole","Golgi apparatus","Endoplasmic reticulum","Cytoplasm","Nucleus","Mitochondrion")),
          Failed = replace_na(Failed,0),
          Passed = replace_na(Passed,0),
          sum = Failed+Passed,
          perc = percent(Passed/sum,accuracy=0.1),
          text = paste0(perc)
        ) %>% arrange(subloc) %>% .$text -> fig_3_02_text_mg132
      
      arginylome_validation_matrix3_sequence_distinct %>% 
        filter(MG132>0) %>% 
        mutate(
          subloc=str_split_i(Localizations,"\\|",1)
        ) %>%
        summarise(
          .by=c(subloc,trinity),
          count =n()
        ) %>% 
        ggplot(
          aes(x=fct(subloc,level = c("Extracellular","Cell membrane","Lysosome/Vacuole","Golgi apparatus","Endoplasmic reticulum","Cytoplasm","Nucleus","Mitochondrion")),y=count,fill=trinity)
        )+
        geom_bar(stat="identity",color="black",linewidth = 0.1)+
        annotate(
          geom = "text",
          x=1:8,
          y=90,
          label = fig_3_02_text_mg132,
          size = 2.3
        )+
        #geom_text(aes(label=count),position = position_stack(vjust=0.5),size=3)+
        scale_fill_manual(
          values = c(wesanderson::wes_palette("Royal1")[2],wesanderson::wes_palette("Royal1")[3])
        )+
        labs(x="Subcellular Localization",y="Count",fill="TrinityNmod")+
        theme_classic2()+
        theme(
          legend.position = "none",
          axis.text = element_text(color="black",size=8),
          axis.text.x = element_text(angle=45,hjust=1),
          axis.title = element_text(color="black",size = 10)
        )
      ggsave("./figures/3_02_removed_arg_subloc_mg132.png",width=8,height=8, units = "cm",dpi = 600)
    }
    #3_02_removed_arg_subloc_mg132.png
    
    ####MGTG
    {
      arginylome_validation_matrix3_sequence_distinct %>% 
        filter(MGTG>0) %>% 
        mutate(
          subloc=str_split_i(Localizations,"\\|",1)
        ) %>%
        summarise(
          .by=c(subloc,trinity),
          count =n()
        ) %>% 
        pivot_wider(names_from = trinity,id_cols = subloc,values_from = count) %>% 
        mutate(
          subloc = fct(subloc,level = c("Extracellular","Cell membrane","Lysosome/Vacuole","Golgi apparatus","Endoplasmic reticulum","Cytoplasm","Nucleus","Mitochondrion")),
          Failed = replace_na(Failed,0),
          Passed = replace_na(Passed,0),
          sum = Failed+Passed,
          perc = percent(Passed/sum,accuracy=0.1),
          text = paste0(perc)
        ) %>% arrange(subloc) %>% .$text -> fig_3_02_text_mgtg
      
      arginylome_validation_matrix3_sequence_distinct %>% 
        filter(MGTG>0) %>% 
        mutate(
          subloc=str_split_i(Localizations,"\\|",1)
        ) %>%
        summarise(
          .by=c(subloc,trinity),
          count =n()
        ) %>% 
        ggplot(
          aes(x=fct(subloc,level = c("Extracellular","Cell membrane","Lysosome/Vacuole","Golgi apparatus","Endoplasmic reticulum","Cytoplasm","Nucleus","Mitochondrion")),y=count,fill=trinity)
        )+
        geom_bar(stat="identity",color="black",linewidth = 0.1)+
        annotate(
          geom = "text",
          x=1:7,
          y=125,
          label = fig_3_02_text_mgtg,
          size = 2.3
        )+
        #geom_text(aes(label=count),position = position_stack(vjust=0.5),size=3)+
        scale_fill_manual(
          values = c(wesanderson::wes_palette("Royal1")[2],wesanderson::wes_palette("Royal1")[3])
        )+
        labs(x="Subcellular Localization",y="Count",fill="TrinityNmod")+
        theme_classic2()+
        theme(
          legend.position = "none",
          axis.text = element_text(color="black",size=8),
          axis.text.x = element_text(angle=45,hjust=1),
          axis.title = element_text(color="black",size = 10)
        )
      ggsave("./figures/3_02_removed_arg_subloc_mgtg.png",width=8,height=8, units = "cm",dpi = 600)
    }
    #3_02_removed_arg_subloc_mgtg.png
    
    
    
    wesanderson::wes_palette("Royal1")[3]
    wesanderson::wes_palette()
    
    arginylome_validation_matrix3_sequence_distinct
    arginylome_validation_matrix3_sequence_distinct %>% 
      filter(score_sum ==3)
  }
  
  ##3. Alphafold
  ###Alphafold
  {
    alphafold_fetched_tidied2 %>% 
      dplyr::filter(score_sum == 3) %>% 
      ggplot(
        aes(
          x=fct(as.character(mod_position)),
          y=prediction_score,
          fill = fct(as.character(mod_position))
        )
      )+
      geom_boxplot(
        linewidth = 0.25
      )+
      lims(y=c(20,100))+
      labs(
        x = "Arginylation Position", y = "AlphaFold pLDDT score"
      )+
      scale_x_discrete(
        labels = c(
          "-10","",
          "-8","",
          "-6","",
          "-4","",
          "-2","",
          "0","",
          "2","",
          "4","",
          "6","",
          "8","",
          "10"
        )
      )+
      theme_classic()+
      theme( 
        legend.position = "none",
        #axis.title.x = element_blank(),
        axis.title = element_text(color = "black",size = 7),
        #axis.text.x = element_blank(),
        axis.text = element_text(color = "black",size = 7),
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm")
      )
    ggsave("./figures/3_02_alphafold_trinity.png",width=4.5,height=4.5,dpi = 600,units="cm")
    
    alphafold_fetched_tidied2 %>% 
      dplyr::filter(score_sum != 3) %>% 
      ggplot(
        aes(
          x=fct(as.character(mod_position)),
          y=prediction_score,
          fill = fct(as.character(mod_position))
        )
      )+
      geom_boxplot(
        linewidth = 0.25,
        outlier.size = 0.25
      )+
      lims(y=c(20,100))+
      labs(
        x = "Arginylation Position", y = "AlphaFold pLDDT score"
      )+
      scale_x_discrete(
        labels = c(
          "-10","",
          "-8","",
          "-6","",
          "-4","",
          "-2","",
          "0","",
          "2","",
          "4","",
          "6","",
          "8","",
          "10"
        )
      )+
      theme_classic()+
      theme( 
        legend.position = "none",
        #axis.title.x = element_blank(),
        axis.title = element_text(color = "black",size = 7),
        #axis.text.x = element_blank(),
        axis.text = element_text(color = "black",size = 7),
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm")
      )
    ggsave("./figures/3_02_alphafold_non_trinity.png",width=4.5,height=4.5,dpi = 600,units="cm")
    
    alphafold_fetched_tidied2 %>% 
      dplyr::filter(score_sum == 3) %>% 
      filter(
        loc_first == "Cytoplasm"
      ) %>% 
      ggplot(
        aes(
          x=fct(as.character(mod_position)),
          y=prediction_score,
          fill = fct(as.character(mod_position))
        )
      )+
      geom_boxplot()+
      lims(y=c(20,100))+
      labs(
        x = "Arginylation Position", y = "AlphaFold pLDDT score"
      )+
      theme_classic()+
      theme(
        legend.position = "none",
        axis.text = element_text(color="black")
      )
    ggsave("./figures/3_02_alphafold_cytoplasm.png",width=10,height=10,dpi = 600,units="cm")
    
    alphafold_fetched_tidied2 %>% 
      dplyr::filter(score_sum == 3) %>% 
      filter(
        loc_first == "Nucleus"
      ) %>% 
      ggplot(
        aes(
          x=fct(as.character(mod_position)),
          y=prediction_score,
          fill = fct(as.character(mod_position))
        )
      )+
      geom_boxplot()+
      lims(y=c(20,100))+
      labs(
        x = "Arginylation Position", y = "AlphaFold pLDDT score"
      )+
      theme_classic()+
      theme(
        legend.position = "none",
        axis.text = element_text(color="black")
      )
    ggsave("./figures/3_02_alphafold_nucleus.png",width=10,height=10,dpi = 600,units="cm")
    
    alphafold_fetched_tidied2 %>% 
      dplyr::filter(score_sum == 3) %>% 
      filter(
        loc_first == "Mitochondrion"
      ) %>% 
      ggplot(
        aes(
          x=fct(as.character(mod_position)),
          y=prediction_score,
          fill = fct(as.character(mod_position))
        )
      )+
      geom_boxplot()+
      lims(y=c(20,100))+
      labs(
        x = "Arginylation Position", y = "AlphaFold pLDDT score"
      )+
      theme_classic()+
      theme(
        legend.position = "none",
        axis.text = element_text(color="black")
      )
    ggsave("./figures/3_02_alphafold_mito.png",width=10,height=10,dpi = 600,units="cm")
    
    alphafold_fetched_tidied2 %>% 
      dplyr::filter(score_sum == 3) %>% 
      filter(
        loc_first == "Endoplasmic reticulum"
      ) %>% 
      ggplot(
        aes(
          x=fct(as.character(mod_position)),
          y=prediction_score,
          fill = fct(as.character(mod_position))
        )
      )+
      geom_boxplot()+
      lims(y=c(20,100))+
      labs(
        x = "Arginylation Position", y = "AlphaFold pLDDT score"
      )+
      theme_classic()+
      theme(
        legend.position = "none",
        axis.text = element_text(color="black")
      )
    ggsave("./figures/3_02_alphafold_er.png",width=10,height=10,dpi = 600,units="cm")
    
    alphafold_fetched_tidied2 %>% 
      dplyr::filter(score_sum == 3) %>% 
      filter(
        loc_first == "Golgi apparatus"
      ) %>% 
      ggplot(
        aes(
          x=fct(as.character(mod_position)),
          y=prediction_score,
          fill = fct(as.character(mod_position))
        )
      )+
      geom_boxplot()+
      lims(y=c(20,100))+
      labs(
        x = "Arginylation Position", y = "AlphaFold pLDDT score"
      )+
      theme_classic()+
      theme(
        legend.position = "none",
        axis.text = element_text(color="black")
      )
    ggsave("./figures/3_02_alphafold_golgi.png",width=10,height=10,dpi = 600,units="cm")
    
    alphafold_fetched_tidied2 %>% 
      dplyr::filter(score_sum == 3) %>% 
      filter(
        loc_first == "Lysosome/Vacuole"
      ) %>% 
      ggplot(
        aes(
          x=fct(as.character(mod_position)),
          y=prediction_score,
          fill = fct(as.character(mod_position))
        )
      )+
      geom_boxplot()+
      lims(y=c(20,100))+
      labs(
        x = "Arginylation Position", y = "AlphaFold pLDDT score"
      )+
      theme_classic()+
      theme(
        legend.position = "none",
        axis.text = element_text(color="black")
      )
    ggsave("./figures/3_02_alphafold_lysosome.png",width=10,height=10,dpi = 600,units="cm")
    
    alphafold_fetched_tidied2 %>% 
      dplyr::filter(score_sum == 3) %>% 
      filter(
        loc_first == "Cell membrane"
      ) %>% 
      ggplot(
        aes(
          x=fct(as.character(mod_position)),
          y=prediction_score,
          fill = fct(as.character(mod_position))
        )
      )+
      geom_boxplot()+
      lims(y=c(20,100))+
      labs(
        x = "Arginylation Position", y = "AlphaFold pLDDT score"
      )+
      theme_classic()+
      theme(
        legend.position = "none",
        axis.text = element_text(color="black")
      )
    ggsave("./figures/3_02_alphafold_membrane.png",width=10,height=10,dpi = 600,units="cm")
    
    alphafold_fetched_tidied2 %>% 
      dplyr::filter(score_sum == 3) %>% 
      filter(
        loc_first == "Extracellular"
      ) %>% 
      ggplot(
        aes(
          x=fct(as.character(mod_position)),
          y=prediction_score,
          fill = fct(as.character(mod_position))
        )
      )+
      geom_boxplot()+
      lims(y=c(20,100))+
      labs(
        x = "Arginylation Position", y = "AlphaFold pLDDT score"
      )+
      theme_classic()+
      theme(
        legend.position = "none",
        axis.text = element_text(color="black")
      )
    ggsave("./figures/3_02_alphafold_extracellular.png",width=10,height=10,dpi = 600,units="cm")
    
    ##t.test modposition 0 
    t.test(
      alphafold_fetched_tidied2 %>% 
        dplyr::filter(score_sum == 3) %>% 
        filter(mod_position == 0) %>% 
        pull(prediction_score),
      alphafold_fetched_tidied2 %>% 
        dplyr::filter(score_sum != 3) %>% 
        filter(mod_position == 0) %>% 
        pull(prediction_score),
      var.equal = TRUE
    )
    
    sd(alphafold_fetched_tidied2 %>% 
         dplyr::filter(score_sum == 3) %>% 
         filter(mod_position == 0) %>% 
         pull(prediction_score))
    
    sd(alphafold_fetched_tidied2 %>% 
         dplyr::filter(score_sum != 3) %>% 
         filter(mod_position == 0) %>% 
         pull(prediction_score))
    
    ### Summary p-value
    alphafold_fetched_tidied2 %>% 
      dplyr::filter(score_sum == 3) %>% head() %>% print(width = Inf)
    
    ##average #deprecated
    {
      alphafold_fetched_tidied2 %>% 
        mutate(
          trinity = if_else(score_sum == 3,TRUE,FALSE)
        ) %>% 
        filter(mod_position<=0) %>% 
        ggplot(
          aes(x=trinity,y=prediction_score)
        )+
        geom_boxplot()+
        scale_x_discrete(
          limits = c(TRUE,FALSE)
        )+
        labs(
          x="TrinityNmod Passed",
          y="Average pLDDT Score"
        )+
        theme_classic2()+
        theme(
          axis.text = element_text(color = "black")
        )
      ggsave(
        "./figures/3_02_alphafold_comparison_average2.png",width=5,height=10,dpi = 600,units="cm"
      )
        
      
      t.test(
        alphafold_fetched_tidied2 %>% 
          dplyr::filter(score_sum == 3) %>% 
          filter(mod_position<=0) %>% 
          .$prediction_score,
        alphafold_fetched_tidied2 %>% 
          dplyr::filter(score_sum != 3) %>% 
          filter(mod_position<=0) %>% 
          .$prediction_score,
        alternative = "two.sided",
        var.equal = TRUE
      ) %>% tidy() %>% mutate(test = "trinity") %>% 
        bind_rows(
          t.test(
            alphafold_fetched_tidied2 %>% 
              dplyr::filter(pcc_check == 1) %>% 
              filter(mod_position<=0) %>% 
              .$prediction_score,
            alphafold_fetched_tidied2 %>% 
              dplyr::filter(pcc_check != 1) %>% 
              filter(mod_position<=0) %>% 
              .$prediction_score,
            alternative = "two.sided",
            var.equal = TRUE
          ) %>% tidy() %>% mutate(test = "ms2")
        ) %>% 
        bind_rows(
          t.test(
            alphafold_fetched_tidied2 %>% 
              dplyr::filter(rt_check == 1) %>% 
              filter(mod_position<=0) %>% 
              .$prediction_score,
            alphafold_fetched_tidied2 %>% 
              dplyr::filter(rt_check != 1) %>% 
              filter(mod_position<=0) %>% 
              .$prediction_score,
            alternative = "two.sided",
            var.equal = TRUE
          ) %>% tidy() %>% mutate(test = "rt")
        )%>% 
        bind_rows(
          t.test(
            alphafold_fetched_tidied2 %>% 
              dplyr::filter(b_y_ttest_check == 1) %>% 
              filter(mod_position<=0) %>% 
              .$prediction_score,
            alphafold_fetched_tidied2 %>% 
              dplyr::filter(b_y_ttest_check != 1) %>% 
              filter(mod_position<=0) %>% 
              .$prediction_score,
            alternative = "two.sided",
            var.equal = TRUE
          ) %>% tidy() %>% mutate(test = "met")
        ) %>% 
        select(
          test, estimate
        )
        filter(test =="trinity")
        ggplot(
          aes(x=fct(test,levels = c("ms2","rt","met","trinity")),y=estimate,
              fill=fct(test,levels = c("ms2","rt","met","trinity")))
        )+
        geom_bar(stat="identity",color = "black",width=0.5)+ 
        geom_line(aes(y=0), color='black')+
        scale_x_discrete(
          #labels = c(expression(MS^2),"RT","MET","ALL"),
          limits=c("trinity")
        )+
        labs(
          y="Average pDLLT Score Change",
          x="TrinityNmod Stage"
        )+
        scale_fill_brewer(type="div",palette = "Set3")+
        theme_classic2()+
        theme(
          legend.position = "none",
          legend.title = element_blank(),
          axis.text = element_text(color = "black", size = 8),
          axis.title = element_text(color = "black", size = 10)
        )
      ggsave(
        "./figures/3_02_alphafold_comparison_average.png",width=5,height=10,dpi = 600,units="cm"
      )
      }
    #3_02_alphafold_comparison_average2.png
    
    ##pvalue #deprecated
    {
      t.test(
        alphafold_fetched_tidied2 %>% 
          dplyr::filter(score_sum == 3) %>% 
          filter(mod_position<=0) %>% 
          .$prediction_score,
        alphafold_fetched_tidied2 %>% 
          dplyr::filter(score_sum != 3) %>% 
          filter(mod_position<=0) %>% 
          .$prediction_score,
        alternative = "two.sided",
        var.equal = TRUE
      ) %>% tidy() %>% mutate(test = "trinity") %>% 
        bind_rows(
          t.test(
            alphafold_fetched_tidied2 %>% 
              dplyr::filter(pcc_check == 1) %>% 
              filter(mod_position<=0) %>% 
              .$prediction_score,
            alphafold_fetched_tidied2 %>% 
              dplyr::filter(pcc_check != 1) %>% 
              filter(mod_position<=0) %>% 
              .$prediction_score,
            alternative = "two.sided",
            var.equal = TRUE
          ) %>% tidy() %>% mutate(test = "ms2")
        ) %>% 
        bind_rows(
          t.test(
            alphafold_fetched_tidied2 %>% 
              dplyr::filter(rt_check == 1) %>% 
              filter(mod_position<=0) %>% 
              .$prediction_score,
            alphafold_fetched_tidied2 %>% 
              dplyr::filter(rt_check != 1) %>% 
              filter(mod_position<=0) %>% 
              .$prediction_score,
            alternative = "two.sided",
            var.equal = TRUE
          ) %>% tidy() %>% mutate(test = "rt")
        )%>% 
        bind_rows(
          t.test(
            alphafold_fetched_tidied2 %>% 
              dplyr::filter(b_y_ttest_check == 1) %>% 
              filter(mod_position<=0) %>% 
              .$prediction_score,
            alphafold_fetched_tidied2 %>% 
              dplyr::filter(b_y_ttest_check != 1) %>% 
              filter(mod_position<=0) %>% 
              .$prediction_score,
            alternative = "two.sided",
            var.equal = TRUE
          ) %>% tidy() %>% mutate(test = "met")
        ) %>% 
        select(
          test, p.value
        ) %>% 
        ggplot(
          aes(x=fct(test,levels = c("ms2","rt","met","trinity")),y=-log(p.value),
              fill=fct(test,levels = c("ms2","rt","met","trinity")))
        )+
        geom_bar(stat="identity",color = "black",width=0.5)+ 
        scale_x_discrete(
          labels = c(expression(MS^2),"RT","MET","ALL")
        )+
        labs(
          y=expression(-Log[10]~italic(P)-Value),
          x="TrinityNmod Stage"
        )+
        scale_fill_brewer(type="div",palette = "Set3")+
        theme_classic2()+
        theme(
          legend.position = "none",
          legend.title = element_blank(),
          axis.text = element_text(color = "black", size = 8),
          axis.title = element_text(color = "black", size = 10)
        )
      ggsave(
        "./figures/3_02_alphafold_comparison_pvalue.png",width=5,height=10,dpi = 600,units="cm"
      )
    }
  }
  
  ##4. Cleavage
  {
    #trinity vs. Cleavage
    {
      arginylome_validation_matrix3_sequence_distinct %>% 
        summarize(
          .by = c(site_type,trinity),
          count =n()
        ) %>% 
        mutate(site_type=replace_na(site_type,"None")) %>% 
        ggplot(
          aes(x=fct(site_type,level = c("SP","mTP","Protease","None")),y=count,fill=trinity)
        )+
        geom_bar(position="fill",stat="identity",color="black",linewidth = 0.1)+
        geom_text(aes(label=count),position = position_fill(vjust=0.5),size=4)+
        scale_x_discrete(label = c("Signal","Transit","Protease","None"))+
        scale_fill_manual(
          values = c(wesanderson::wes_palette("Royal1")[2],wesanderson::wes_palette("Royal1")[3])
        )+
        scale_y_continuous(label=label_percent())+
        labs(x="Known Cleavage Site",y="Percent",fill="TrinityNmod")+
        theme_classic2()+
        theme(
          legend.position = "none",
          axis.text = element_text(color="black",size=10),
          #axis.text.x = element_text(angle=45,hjust=1),
          axis.title = element_text(color="black",size = 10),
          axis.title.y = element_blank()
        )
      ggsave("./figures/3_04_cleavage_site_trinity.png",width=8,height=8, units = "cm",dpi = 600)
    }
    #3_04_cleavage_site_trinity.png
    
  }

  ##6. Cleavage LOGO #geomlogo
  #@arginylome_validation_matrix3_sequence_distinct_trinity
  #@summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered_trinity
  {
    library(ggseqlogo)
    
    
    
    
    arginylome_validation_matrix3_sequence_distinct_trinity %>% 
      summarise(
        .by = protease_site,
        count = n()
      ) %>% 
      arrange(-count) %>% 
      filter(!is.na(protease_site)) %>% write_tsv("./export/4_06_protease_site.tsv")
    
    arginylome_validation_matrix3_sequence_distinct_trinity %>% 
      summarize(
        .by = protease_site,
        count = n()
      )
  
    #signal  
    {
      arginylome_validation_matrix3_sequence_distinct_trinity %>% 
        filter(protease_site=="Signal Peptide") %>% select(p5p5prime) -> logo_data
      ggplot()+geom_logo(data = logo_data)+
        theme_logo()+
        theme(
          legend.position = "none"
        )->fig_logo
      
      fig_logo$scales$scales[[1]]$labels <- c("P5","P4","P3","P2","P1","P1'","P2'","P3'","P4'","P5'")
      
      fig_logo
      ggsave("./figures/3_05_logo_trinity_signal.png",width=10,height=10,dpi = 600,units="cm")
      }
    
    ##transit
    
    {
      arginylome_validation_matrix3_sequence_distinct_trinity %>% 
        filter(protease_site=="Transit Peptide") %>% select(p5p5prime) -> logo_data
      ggplot()+geom_logo(data = logo_data)+
        theme_logo()+
        theme(
          legend.position = "none"
        )->fig_logo
      
      fig_logo$scales$scales[[1]]$labels <- c("P5","P4","P3","P2","P1","P1'","P2'","P3'","P4'","P5'")
      
      fig_logo
      ggsave("./figures/3_05_logo_trinity_transit.png",width=10,height=10,dpi = 600,units="cm")
    }
    
    ##caspase3
    
    {
      arginylome_validation_matrix3_sequence_distinct_trinity %>% 
        filter(protease_site=="Caspase-3") %>% select(p5p5prime) -> logo_data
      ggplot()+geom_logo(data = logo_data)+
        theme_logo()+
        theme(
          legend.position = "none"
        )->fig_logo
      
      fig_logo$scales$scales[[1]]$labels <- c("P5","P4","P3","P2","P1","P1'","P2'","P3'","P4'","P5'")
      
      fig_logo
      ggsave("./figures/3_05_logo_trinity_caspase3.png",width=10,height=10,dpi = 600,units="cm")
    }
    
    ##NA
    
    {
      arginylome_validation_matrix3_sequence_distinct_trinity %>% 
        filter(is.na(protease_site)) %>% select(p5p5prime) -> logo_data
      ggplot()+geom_logo(data = logo_data)+
        theme_logo()+
        theme(
          legend.position = "none"
        )->fig_logo
      
      fig_logo$scales$scales[[1]]$labels <- c("P5","P4","P3","P2","P1","P1'","P2'","P3'","P4'","P5'")
      
      fig_logo
      ggsave("./figures/3_05_logo_trinity_NA.png",width=10,height=10,dpi = 600,units="cm")
    }
    
    
    
    arginylome_validation_matrix3_sequence_distinct_trinity %>% 
      filter(protease_site=="Signal Peptide") %>% 
      select(p5p5prime) -> logo_site_signal
    ggplot()+geom_logo(data = logo_site_signal$p5p5prime)+
      theme_logo()+
      theme(
        legend.position = "none"
      )
    ggsave("./figures/4_06_logo_signal.png",width=10,height=10,dpi = 600,units="cm")
    
    arginylome_validation_matrix3_sequence_distinct_trinity %>% 
      filter(protease_site=="Caspase-3") %>% 
      select(p5p5prime) -> logo_site_casp3
    ggplot()+geom_logo(data = logo_site_casp3$p5p5prime)+
      theme_logo()+
      theme(
        legend.position = "none"
      )
    ggsave("./figures/4_06_logo_casp3.png",width=10,height=10,dpi = 600,units="cm")
  }
}
##cleavage violin
{
  
  arginylome_validation_matrix3_sequence_distinct_trinity %>% 
    filter(protease_site=="Caspase-3") %>% pull(gene_position) -> gene_position_casp3
  
  arginylome_validation_matrix3_sequence_distinct_trinity %>% 
    filter(protease_site=="Signal Peptide") %>% pull(gene_position) -> gene_position_signal
  
  arginylome_validation_matrix3_sequence_distinct_trinity %>% 
    filter(protease_site=="Transit Peptide") %>% pull(gene_position) -> gene_position_transit
  
  arginylome_validation_matrix3_sequence_distinct_trinity %>% 
    filter(protease_site=="Caspase-6") %>% pull(gene_position) -> gene_position_casp6
  
  ##ARG
  {
    peptide_lfq_trinity_distinct_arg %>% 
      filter(
        gene_position %in% gene_position_casp3
      ) %>% 
      pivot_longer(MG132:MOCK) %>% 
      ggplot()+
      aes(x= factor(name, levels = c("MOCK","MG132","MGTG")), y = value, fill = factor(name, levels = c("MOCK","MG132","MGTG")))+
      geom_violin(linewidth = 0.2)+
      geom_boxplot(width = 0.1,linewidth = 0.2,outlier.size = 0)+
      stat_compare_means(
        comparisons = list(c("MOCK","MG132"),c("MG132","MGTG"),c("MOCK","MGTG")),
        label = "p.signif",
        method="t.test",
        #label.y = c(2,2.7,3.4),
        size = 2
      )+
      labs(x="Experiments", y="Z-score")+
      scale_fill_manual(values =c(wes_palette("Zissou1",5)[c(1,3,5)]))+
      lims(y = c(-2,3))+
      theme_classic2()+
      theme(
        legend.position = "none",
        axis.title = element_text(size=7),
        axis.text = element_text(size=6,color = "black"),
        strip.text = element_text(size=6)
      )
    ggsave("./figures/4_06_cleavage_violin_casp3.png",width=4,height=4,dpi = 600,units="cm")
    
    peptide_lfq_trinity_distinct_arg %>% 
      filter(
        gene_position %in% gene_position_signal
      ) %>% 
      pivot_longer(MG132:MOCK) %>% 
      ggplot()+
      aes(x= factor(name, levels = c("MOCK","MG132","MGTG")), y = value, fill = factor(name, levels = c("MOCK","MG132","MGTG")))+
      geom_violin(linewidth = 0.2)+
      geom_boxplot(width = 0.1,linewidth = 0.2,outlier.size = 0)+
      stat_compare_means(
        comparisons = list(c("MOCK","MG132"),c("MG132","MGTG"),c("MOCK","MGTG")),
        label = "p.signif",
        method="t.test",
        #label.y = c(2,2.7,3.4),
        size = 2
      )+
      labs(x="Experiments", y="Z-score")+
      scale_fill_manual(values =c(wes_palette("Zissou1",5)[c(1,3,5)]))+
      lims(y = c(-2,3))+
      theme_classic2()+
      theme(
        legend.position = "none",
        axis.title = element_text(size=7),
        axis.text = element_text(size=6,color = "black"),
        strip.text = element_text(size=6)
      )
    ggsave("./figures/4_06_cleavage_violin_signal.png",width=4,height=4,dpi = 600,units="cm")
    
    peptide_lfq_trinity_distinct_arg %>% 
      filter(
        gene_position %in% gene_position_transit
      ) %>% 
      pivot_longer(MG132:MOCK) %>% 
      ggplot()+
      aes(x= factor(name, levels = c("MOCK","MG132","MGTG")), y = value, fill = factor(name, levels = c("MOCK","MG132","MGTG")))+
      geom_violin(linewidth = 0.2)+
      geom_boxplot(width = 0.1,linewidth = 0.2,outlier.size = 0)+
      stat_compare_means(
        comparisons = list(c("MOCK","MG132"),c("MG132","MGTG"),c("MOCK","MGTG")),
        label = "p.signif",
        method="t.test",
        #label.y = c(2,2.7,3.4),
        size = 2
      )+
      labs(x="Experiments", y="Z-score")+
      scale_fill_manual(values =c(wes_palette("Zissou1",5)[c(1,3,5)]))+
      lims(y = c(-2,3))+
      theme_classic2()+
      theme(
        legend.position = "none",
        axis.title = element_text(size=7),
        axis.text = element_text(size=6,color = "black"),
        strip.text = element_text(size=6)
      )
    ggsave("./figures/4_06_cleavage_violin_transit.png",width=4,height=4,dpi = 600,units="cm")
    
    peptide_lfq_trinity_distinct_arg %>% 
      filter(
        gene_position %in% gene_position_casp6
      ) %>% 
      pivot_longer(MG132:MOCK) %>% 
      ggplot()+
      aes(x= factor(name, levels = c("MOCK","MG132","MGTG")), y = value, fill = factor(name, levels = c("MOCK","MG132","MGTG")))+
      geom_violin(linewidth = 0.2)+
      geom_boxplot(width = 0.1,linewidth = 0.2,outlier.size = 0)+
      stat_compare_means(
        comparisons = list(c("MOCK","MG132"),c("MG132","MGTG"),c("MOCK","MGTG")),
        label = "p.signif",
        method="t.test",
        #label.y = c(2,2.7,3.4),
        size = 2
      )+
      labs(x="Experiments", y="Z-score")+
      scale_fill_manual(values =c(wes_palette("Zissou1",5)[c(1,3,5)]))+
      lims(y = c(-2,3))+
      theme_classic2()+
      theme(
        legend.position = "none",
        axis.title = element_text(size=7),
        axis.text = element_text(size=6,color = "black"),
        strip.text = element_text(size=6)
      )
    ggsave("./figures/4_06_cleavage_violin_casp6.png",width=4,height=4,dpi = 600,units="cm")
  }
  
  ##FREE
  {
    peptide_lfq_trinity_distinct_free %>% 
      filter(
        gene_position %in% gene_position_signal
      ) %>% 
      pivot_longer(MG132:MOCK) %>% 
      ggplot()+
      aes(x= factor(name, levels = c("MOCK","MG132","MGTG")), y = value, fill = factor(name, levels = c("MOCK","MG132","MGTG")))+
      geom_violin(linewidth = 0.2)+
      geom_boxplot(width = 0.1,linewidth = 0.2,outlier.size = 0)+
      stat_compare_means(
        comparisons = list(c("MOCK","MG132"),c("MG132","MGTG"),c("MOCK","MGTG")),
        label = "p.signif",
        method="t.test",
        #label.y = c(2,2.7,3.4),
        size = 2
      )+
      labs(x="Experiments", y="Z-score")+
      scale_fill_manual(values =c(wes_palette("Zissou1",5)[c(1,3,5)]))+
      lims(y = c(-2,3))+
      theme_classic2()+
      theme(
        legend.position = "none",
        axis.title = element_text(size=7),
        axis.text = element_text(size=6,color = "black"),
        strip.text = element_text(size=6)
      )
    ggsave("./figures/4_06_cleavage_violin_signal_free.png",width=4,height=4,dpi = 600,units="cm")
    
    peptide_lfq_trinity_distinct_free %>% 
      filter(
        gene_position %in% gene_position_transit
      ) %>% 
      pivot_longer(MG132:MOCK) %>% 
      ggplot()+
      aes(x= factor(name, levels = c("MOCK","MG132","MGTG")), y = value, fill = factor(name, levels = c("MOCK","MG132","MGTG")))+
      geom_violin(linewidth = 0.2)+
      geom_boxplot(width = 0.1,linewidth = 0.2,outlier.size = 0)+
      stat_compare_means(
        comparisons = list(c("MOCK","MG132"),c("MG132","MGTG"),c("MOCK","MGTG")),
        label = "p.signif",
        method="t.test",
        #label.y = c(2,2.7,3.4),
        size = 2
      )+
      labs(x="Experiments", y="Z-score")+
      scale_fill_manual(values =c(wes_palette("Zissou1",5)[c(1,3,5)]))+
      lims(y = c(-2,3))+
      theme_classic2()+
      theme(
        legend.position = "none",
        axis.title = element_text(size=7),
        axis.text = element_text(size=6,color = "black"),
        strip.text = element_text(size=6)
      )
    ggsave("./figures/4_06_cleavage_violin_transit_free.png",width=4,height=4,dpi = 600,units="cm")
  }
}



####
#PreChapter 4 >> peptide level quantifications
#####
##summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered -> summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered_trinity
{
  summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered %>% filter(dl_score ==3) %>% colnames ###323
  
  summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered %>% 
    mutate(
      gene_position = paste0(gene,"_",position)
    ) -> summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered
  
  summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered %>% filter(dl_score ==3) -> summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered_trinity
  
  summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered_trinity %>% 
    distinct(Sequence) %>% 
    pull(Sequence)-> temp_seq
  
  arginylome_validation_matrix3_sequence_distinct%>% 
    filter(
      Sequence %in% temp_seq
    )
  
  
  
  
  summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered %>% filter(dl_score ==3) %>% 
    filter(enzyme == "trypsin") %>% 
    distinct(Sequence) %>% 
    pull(Sequence) -> stripseq_arg_trinity_trypsin
  
  
  summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered %>% filter(dl_score ==3) %>% 
    filter(enzyme == "chymotrypsin") %>% 
    distinct(Sequence) %>% 
    pull(Sequence) -> stripseq_arg_trinity_chymotrypsin
  
  summaryfile %>% colnames()
  ###filter peptide file with arg seq (not constrainted by mod)
  
  ###update::240711
  summaryfile %>% # peptide file
    rename(
      "MG132_1"="Abundances Normalized F12 Sample MG132 1ST trypsin",
      "MG132_2"="Abundances Normalized F6 Sample MG132 2ND trypsin",
      "MGTG_1"="Abundances Normalized F1 Sample MGTG 1ST trypsin",
      "MGTG_2"="Abundances Normalized F7 Sample MGTG 2ND trypsin",
      "MOCK_1"="Abundances Normalized F2 Sample MOCK 1ST trypsin",
      "MOCK_2"="Abundances Normalized F5 Sample MOCK 2ND trypsin"
    ) %>% 
    mutate(
      across(
        .cols = c("MG132_1","MG132_2","MGTG_1","MGTG_2","MOCK_1","MOCK_2"),
        .fns = \(x) log(x,10)
      )
    ) %>% 
    select(
      Sequence,nmod,Modifications,MG132_1,MG132_2,MGTG_1,MGTG_2,MOCK_1,MOCK_2
    ) %>% 
    mutate(
      across(c("MG132_1","MG132_2","MGTG_1","MGTG_2","MOCK_1","MOCK_2"),~na_if(.x,0))
    ) %>% 
    mutate(
      across(
        c("MG132_1","MG132_2","MGTG_1","MGTG_2","MOCK_1","MOCK_2"),
        \(x) (x-mean(x,na.rm=TRUE))/sd(x,na.rm = TRUE)
      )
    ) %>%
    filter(
      Sequence %in% stripseq_arg_trinity_trypsin
    ) %>% 
    bind_rows(
      summaryfile %>% # peptide file
        rename(
          "MG132_1"="Abundances Normalized F11 Sample MG132 1ST chymotrypsin",
          "MG132_2"="Abundances Normalized F9 Sample MG132 2ND chymotrypsin",
          "MGTG_1"="Abundances Normalized F3 Sample MGTG 1ST chymotrypsin",
          "MGTG_2"="Abundances Normalized F10 Sample MGTG 2ND chymotrypsin",
          "MOCK_1"="Abundances Normalized F4 Sample MOCK 1ST chymotrypsin",
          "MOCK_2"="Abundances Normalized F8 Sample MOCK 2ND chymotrypsin"
        ) %>% 
        mutate(
          across(
            .cols = c("MG132_1","MG132_2","MGTG_1","MGTG_2","MOCK_1","MOCK_2"),
            .fns = \(x) log(x,10)
          )
        ) %>% 
        mutate(
          across(c("MG132_1","MG132_2","MGTG_1","MGTG_2","MOCK_1","MOCK_2"),~na_if(.x,0))
        ) %>% 
        mutate(
          across(
            c("MG132_1","MG132_2","MGTG_1","MGTG_2","MOCK_1","MOCK_2"),
            \(x) (x-mean(x,na.rm=TRUE))/sd(x,na.rm = TRUE)
          )
        ) %>%
        filter(
          Sequence %in% stripseq_arg_trinity_chymotrypsin
        ) %>%
        select(
          Sequence,nmod,Modifications,MG132_1,MG132_2,MGTG_1,MGTG_2,MOCK_1,MOCK_2
        )
    ) -> peptide_lfq_trinity
  
  
  ##GENE position
  peptide_lfq_trinity %>% 
    filter(nmod %in% c("D3Acetyl","D3AcetylArg","D3AcetylArgDeamid")) %>%
    mutate(
      nmod = case_when(nmod == "D3AcetylArgDeamid"~"Arg",
                       nmod == "D3Acetyl"~"Free",
                       nmod == "D3AcetylArg"~"Arg",
                       .default = nmod)
    ) %>% 
    pivot_wider(
      id_cols = Sequence,
      names_from = nmod,
      values_from = MG132_1:MOCK_2,
      values_fn = mean
    ) %>% 
    left_join(
      summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered_trinity %>% 
        select(Sequence,gene_position,prot_pos_reposition,Localizations,cleavage_type,cleavage_name) %>% distinct(Sequence,.keep_all = TRUE),
      join_by(Sequence)
    ) -> peptide_lfq_trinity
  
  summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered %>%
    write_tsv("./export/summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered_240715.tsv")
  
}
###peptide_lfq_trinity

###
####Chapter 4. Utilizing Arginylome

{
  ##0, volcano
  {
      
    peptide_lfq_trinity %>% 
      nest(
        MG132_free = c(MG132_1_Free,MG132_2_Free),
        MOCK_free = c(MOCK_1_Free,MOCK_2_Free),
        MGTG_free = c(MGTG_1_Free,MGTG_2_Free),
        MG132_arg = c(MG132_1_Arg,MG132_2_Arg),
        MOCK_arg = c(MOCK_1_Arg,MOCK_2_Arg),
        MGTG_arg = c(MGTG_1_Arg,MGTG_2_Arg)
      ) %>% 
      mutate(
        across(
          MG132_free:MGTG_arg,
          \(x) map(
            x,
            \(x) pivot_longer(x,everything()) %>% pull(value)
          )
        )
      ) %>%
      mutate(
        MGTG_MOCK_arg = map2(
          MGTG_arg,
          MOCK_arg,
          \(x,y) {
            if(any(!is.na(x) & !is.na(y)))
              t.test(x,y,var.equal = TRUE,paired = FALSE) %>% tidy()
            else
              NA
          }
        ),
        MGTG_MOCK_arg_mean = map_dbl(
          MGTG_MOCK_arg,
          \(x) {
            if(any(!is.na(x)))
              pull(x,estimate)
            else
              NA
          }
        ),
        MGTG_MOCK_arg_pvalue = map_dbl(
          MGTG_MOCK_arg,
          \(x) {
            if(any(!is.na(x)))
              pull(x,p.value)
            else
              NA
          }
        )
      ) %>%
      mutate(
        MG132_MOCK_arg = map2(
          MG132_arg,
          MOCK_arg,
          \(x,y) {
            if(any(!is.na(x) & !is.na(y)))
              t.test(x,y,var.equal = TRUE,paired = FALSE) %>% tidy()
            else
              NA
          }
        ),
        MG132_MOCK_arg_mean = map_dbl(
          MG132_MOCK_arg,
          \(x) {
            if(any(!is.na(x)))
              pull(x,estimate)
            else
              NA
          }
        ),
        MG132_MOCK_arg_pvalue = map_dbl(
          MG132_MOCK_arg,
          \(x) {
            if(any(!is.na(x)))
              pull(x,p.value)
            else
              NA
          }
        )
      )%>%
      mutate(
        MG132_MOCK_free = map2(
          MG132_free,
          MOCK_free,
          \(x,y) {
            if(any(!is.na(x) & !is.na(y)))
              t.test(x,y,var.equal = TRUE,paired = FALSE) %>% tidy()
            else
              NA
          }
        ),
        MG132_MOCK_free_mean = map_dbl(
          MG132_MOCK_free,
          \(x) {
            if(any(!is.na(x)))
              pull(x,estimate)
            else
              NA
          }
        ),
        MG132_MOCK_free_pvalue = map_dbl(
          MG132_MOCK_free,
          \(x) {
            if(any(!is.na(x)))
              pull(x,p.value)
            else
              NA
          }
        )
      )%>% 
      mutate(
        MGTG_MOCK_free = map2(
          MGTG_free,
          MOCK_free,
          \(x,y) {
            if(any(!is.na(x) & !is.na(y)))
              t.test(x,y,var.equal = TRUE,paired = FALSE) %>% tidy()
            else
              NA
          }
        ),
        MGTG_MOCK_free_mean = map_dbl(
          MGTG_MOCK_free,
          \(x) {
            if(any(!is.na(x)))
              pull(x,estimate)
            else
              NA
          }
        ),
        MGTG_MOCK_free_pvalue = map_dbl(
          MGTG_MOCK_free,
          \(x) {
            if(any(!is.na(x)))
              pull(x,p.value)
            else
              NA
          }
        )
      )%>%
      mutate(
        MGTG_MG132_arg = map2(
          MGTG_arg,
          MG132_arg,
          \(x,y) {
            if(any(!is.na(x) & !is.na(y)))
              t.test(x,y,var.equal = TRUE,paired = FALSE) %>% tidy()
            else
              NA
          }
        ),
        MGTG_MG132_arg_mean = map_dbl(
          MGTG_MG132_arg,
          \(x) {
            if(any(!is.na(x)))
              pull(x,estimate)
            else
              NA
          }
        ),
        MGTG_MG132_arg_pvalue = map_dbl(
          MGTG_MG132_arg,
          \(x) {
            if(any(!is.na(x)))
              pull(x,p.value)
            else
              NA
          }
        )
      )%>%
      mutate(
        MGTG_MG132_free = map2(
          MGTG_free,
          MG132_free,
          \(x,y) {
            if(any(!is.na(x) & !is.na(y)))
              t.test(x,y,var.equal = TRUE,paired = FALSE) %>% tidy()
            else
              NA
          }
        ),
        MGTG_MG132_free_mean = map_dbl(
          MGTG_MG132_free,
          \(x) {
            if(any(!is.na(x)))
              pull(x,estimate)
            else
              NA
          }
        ),
        MGTG_MG132_free_pvalue = map_dbl(
          MGTG_MG132_free,
          \(x) {
            if(any(!is.na(x)))
              pull(x,p.value)
            else
              NA
          }
        )
      ) -> peptide_lfq_trinity_for_volcano
      
    peptide_lfq_trinity_for_volcano %>% 
      mutate(
        color = case_when(
          str_detect(cleavage_type,"C14")~"Caspase",
          str_detect(cleavage_type,"SP")~"SP",
          str_detect(cleavage_type,"mTP")~"mTP",
          .default = NA_character_
        )
      ) %>% 
      mutate(
        sig = if_else(
          abs(MGTG_MOCK_arg_mean)>=1 & MGTG_MOCK_arg_pvalue<=0.05,
          gene_position,
          NA
        )
      ) %>% 
      ggplot(
        aes(
          x=MGTG_MOCK_arg_mean,
          y=-log(MGTG_MOCK_arg_pvalue,10)
        )
      )+
      geom_point(
        aes(
          color = color
        )
      )+
      geom_text_repel(
        aes(
          label = sig
        ),
        size = 2,
        force=0.3,
        min.segment.length = 0.2
      )+
      labs(
        x= expression(Log[2](MGTG/MOCK)),
        y= expression(-Log[10](p-value))
      )+
      scale_color_manual(
        values = RColorBrewer::brewer.pal(n=9,"Set1")[c(1:3,9)]
      )+
      theme_classic()+
      theme(
        axis.text = element_text(color = "black",size= 10),
        axis.title = element_text(color = "black",size= 10),
        plot.title = element_text(color = "black",size= 10),
        legend.position = "none",
        axis.ticks = element_line(color = "black"),
      )
    ggsave(
      "./figures/4_00_volcano_mgtg_mock.png",height = 8, width =8,dpi = 600,units = "cm"
    )
    RColorBrewer::brewer.pal(n=9,"Set1")
    RColorBrewer::display.brewer.all()
    
    peptide_lfq_trinity_for_volcano %>% 
      mutate(
        color = case_when(
          str_detect(cleavage_type,"C14")~"Caspase",
          str_detect(cleavage_type,"SP")~"SP",
          str_detect(cleavage_type,"mTP")~"mTP",
          .default = NA_character_
        )
      ) %>% 
      mutate(
        sig = if_else(
          abs(MG132_MOCK_arg_mean)>=1 & MG132_MOCK_arg_pvalue<=0.05,
          gene_position,
          NA
        )
      ) %>% 
      ggplot(
        aes(
          x=MG132_MOCK_arg_mean,
          y=-log(MG132_MOCK_arg_pvalue,10)
        )
      )+
      geom_point(
        aes(
          color = color
        )
      )+
      geom_text_repel(
        aes(
          label = sig
        ),
        size = 2,
        force=0.3,
        min.segment.length = 0.2
      )+
      labs(
        x= expression(Log[2](MG132/MOCK)),
        y= expression(-Log[10](p-value))
      )+
      scale_color_manual(
        values = RColorBrewer::brewer.pal(n=9,"Set1")[c(1:3,9)]
      )+
      theme_classic()+
      theme(
        axis.text = element_text(color = "black",size= 10),
        axis.title = element_text(color = "black",size= 10),
        plot.title = element_text(color = "black",size= 10),
        legend.position = "none",
        axis.ticks = element_line(color = "black"),
      )
    ggsave(
      "./figures/4_00_volcano_mg132_mock.png",height = 8, width =8,dpi = 600,units = "cm"
    )
    
    peptide_lfq_trinity_for_volcano %>% 
      mutate(
        color = case_when(
          str_detect(cleavage_type,"C14")~"Caspase",
          str_detect(cleavage_type,"SP")~"SP",
          str_detect(cleavage_type,"mTP")~"mTP",
          .default = NA_character_
        )
      ) %>% 
      mutate(
        sig = if_else(
          abs(MGTG_MG132_arg_mean)>=1 & MGTG_MG132_arg_pvalue<=0.1,
          gene_position,
          NA
        ),
        gene = str_split_i(gene_position,"\\_",1),
        upr = if_else(gene %in% keywords_upr,gene_position,NA)
      ) %>% 
      ggplot(
        aes(
          x=MGTG_MG132_arg_mean,
          y=-log(MGTG_MG132_arg_pvalue,10)
        )
      )+
      geom_point(
        aes(
          color = color
        )
      )+
      geom_text_repel(
        aes(
          label = sig
        ),
        size = 2,
        force=0.3,
        min.segment.length = 0.2
      )+
      labs(
        x= expression(Log[2](MGTG/MG132)),
        y= expression(-Log[10](p-value))
      )+
      scale_color_manual(
        values = RColorBrewer::brewer.pal(n=9,"Set1")[c(1:3,9)]
      )+
      theme_classic()+
      theme(
        axis.text = element_text(color = "black",size= 10),
        axis.title = element_text(color = "black",size= 10),
        plot.title = element_text(color = "black",size= 10),
        legend.position = "none",
        axis.ticks = element_line(color = "black"),
      )
    ggsave(
      "./figures/4_00_volcano_mgtg_mg132.png",height = 8, width =8,dpi = 600,units = "cm"
    )
    
    peptide_lfq_trinity %>% 
      pivot_longer(MG132_1_Free:MOCK_2_Arg) %>% 
      drop_na(value) %>% 
      mutate(
        experiment = str_split_i(name,"_",1),
        arg = str_split_i(name,"_",3)
      ) %>% 
      ggplot(
        aes(x= factor(experiment, levels = c("MOCK","MG132","MGTG")), y = value, fill = factor(experiment, levels = c("MOCK","MG132","MGTG")))
      )+
      geom_violin(linewidth = 0.2,scale="width")+
      geom_boxplot(width = 0.1,linewidth = 0.2,outlier.size = 0)+
      labs(x="Experiments", y="Z-score",fill="name")+
      stat_compare_means(
        comparisons = list(c("MOCK","MG132"),c("MOCK","MGTG"),c("MG132","MGTG")),
        label = "p.format",
        method="t.test",
        label.y = c(4,4.8,4.3),
        size = 2.5
      )+
      scale_fill_manual(
        values =c(wesanderson::wes_palette(5,name="Zissou1",type="discrete")[c(1,3,5)]),
        labels = c("MOCK","MG132","MGTG")
      )+
      lims(
        y=c(-3.1,5.5)
      )+
      theme_classic2()+
      theme(
        legend.position = "none",
        # legend.title = element_blank(),
        # legend.text = element_text(size = 8,color="black"),
        # legend.key.width = unit(0.2,"cm"),
        # legend.key.height = unit(0.2,"cm"),
        # legend.margin = margin(t = 0, unit="cm"),
        # legend.box.margin = margin(t = 0, unit="cm"),
        # legend.box.spacing = unit(2, "pt"),
        #axis.ticks.x = element_blank(),
        #axis.text.x = element_blank(),
        axis.text = element_text(size = 8,color="black"),
        axis.title = element_text(size = 8,color="black"),
        strip.text= element_text(size = 8,color="black"),
        axis.title.x = element_blank(),
        axis.ticks = element_line(color = "black")
      )
    ggsave(
      "./figures/4_00_all_violin.png",height = 5, width =5,dpi = 600,units = "cm"
    )
  }
  library(ggrepel)
  ggrepel::
  
  tibble(x=122,y=132) %>% pivot_longer(everything()) %>% pull(value)
  
  
  ##1. Heatmap
  {
    library(ComplexHeatmap)
    
    arginylome_validation_matrix3_sequence_distinct_trinity
    
    
    
    peptide_lfq_trinity %>% 
      select(-gene_position) %>% 
      left_join(
        table_gene_position, join_by(Sequence)
      ) %>% 
      arrange(Sequence,gene_position) %>% 
      distinct(gene_position,.keep_all = TRUE) -> peptide_lfq_trinity_distinct 
    
    peptide_lfq_trinity_distinct %>% view()
    
    peptide_lfq_trinity_distinct %>% select(
      contains("Arg")
    ) %>% 
      drop_na() %>%
      as.matrix() -> peptide_lfq_trinity_distinct_mat
    
    
    
    ###peptide_lfq_trinity_distinct <<< gene_position 
    peptide_lfq_trinity_distinct %>% 
      pivot_longer(cols = matches("Free|Arg")) %>% 
      drop_na(value) %>% 
      mutate(
        experiment = str_split_i(name,"\\_",1),
        type = str_split_i(name,"\\_",-1)
      ) %>% 
      summarise(
        .by=c(gene_position,experiment,type),
        mean = mean(value)
      ) %>% 
      pivot_wider(
        id_cols = c(gene_position,type),
        names_from = experiment,
        values_from = mean
      ) %>% 
      mutate(
        MG132_MOCK = MG132-MOCK,
        MGTG_MOCK = MGTG-MOCK
      ) -> peptide_lfq_trinity_distinct_all
    
    peptide_lfq_trinity_distinct_all %>% 
      filter(type == "Arg") -> peptide_lfq_trinity_distinct_arg
    
    peptide_lfq_trinity_distinct_all %>% 
      filter(type == "Free") -> peptide_lfq_trinity_distinct_free
    
    
    ###HEATMAP MATRIX CONSTRUCTION!
    peptide_lfq_trinity_distinct_arg %>% 
      select(MG132_MOCK:MGTG_MOCK) %>% 
      as.matrix() -> peptide_lfq_trinity_distinct_mat
    
    rownames(peptide_lfq_trinity_distinct_mat) <- peptide_lfq_trinity_distinct_arg %>% .$gene_position
    
    
    arginylome_validation_matrix3_sequence_distinct_trinity %>% 
      mutate(loc = str_split_i(Localizations,"\\|",1)) -> arginylome_validation_matrix3_sequence_distinct_trinity
    
    
    arginylome_validation_matrix3_sequence_distinct_trinity %>% 
      select(
        prot_pos_master,loc
      ) %>% distinct(loc) %>% pull(loc)
    
    
    
    ###color code for LOC
    tibble(
      loc = c("Extracellular","Cell membrane","Lysosome/Vacuole","Golgi apparatus","Endoplasmic reticulum","Cytoplasm","Nucleus","Mitochondrion"),
      color = c(paste0(wes_palette("Zissou1",7,type = "continuous")),"#4DAF4A")
    )
    
    ###site >>> loc
    tibble(
      site = rownames(peptide_lfq_trinity_distinct_mat)
    ) %>% left_join(
      arginylome_validation_matrix3_sequence_distinct_trinity %>% 
        select(gene_position,loc),
      join_by(site == gene_position)
    ) %>% pull(loc)
    
    ###HEATMAP
    set.seed(2421)
    
    ##GENE ANNOTATION
    row_ha = HeatmapAnnotation(
      Loc = tibble(
        site = rownames(peptide_lfq_trinity_distinct_mat)
      ) %>% left_join(
        arginylome_validation_matrix3_sequence_distinct_trinity %>% 
          select(gene_position,loc),
        join_by(site == gene_position)
      ) %>% pull(loc),
      col = list(Loc = c(
        "Extracellular"="#3B9AB2",
        "Cell membrane"="#63ADBE",
        "Lysosome/Vacuole"="#9EBE91",
        "Golgi apparatus"="#EBCC2A",
        "Endoplasmic reticulum"="#E4B80E",
        "Cytoplasm"="#D5D5D3",
        "Nucleus"="#F21A00",
        "Mitochondrion"="#4DAF4A"
      )),
      annotation_legend_param = list(
        Loc = list(
          title ="Localization",
          at = c("Extracellular","Cell membrane","Lysosome/Vacuole","Golgi apparatus","Endoplasmic reticulum","Cytoplasm","Nucleus","Mitochondrion"),
          labels = c("Extracellular","Cell membrane","Lysosome/Vacuole","Golgi apparatus","Endoplasmic reticulum","Cytoplasm","Nucleus","Mitochondrion")
        )
      )
    )
    png(
      paste0("./figures/heat2.png"), width = 5, height = 7, 
      units = "cm", res = 600
    )
    set.seed(2421)
    Heatmap(
      peptide_lfq_trinity_distinct_mat[,c("MG132_MOCK","MGTG_MOCK")],
      cluster_columns = FALSE,
      clustering_distance_rows = "pearson",
      clustering_method_rows = "average",
      row_km = 5,
      #bottom_annotation = row_ha,
      border = TRUE,
      column_names_gp = gpar(fontsize = 5),
      row_names_gp = gpar(fontsize = 8),
      column_names_rot = 0,
      show_row_names = FALSE,
      column_names_centered = TRUE,
      heatmap_legend_param = list(title = "Z score")
    )-> heatmap_data
    heatmap_data
    dev.off()
    
    peptide_lfq_trinity_distinct_mat %>% length()
  }
  
  ##heatmap data export
  {
    tibble(stack(row_order(heatmap_data))) %>% 
      left_join(
        tibble(
          num = 1:length(rownames(heatmap_data@matrix)),
          name = rownames(heatmap_data@matrix)
        ),
        join_by(values == num)
      ) -> heatmap_data_output_hclust
    
    peptide_lfq_trinity_distinct_arg
    peptide_lfq_trinity_distinct_free
    
    heatmap_data_output_hclust %>% 
      inner_join(
        peptide_lfq_trinity_distinct_free %>% 
          select(gene_position:MOCK),
        join_by(
          name == gene_position
        )
      )
    heatmap_data_output_hclust
    
    heatmap_data_output_hclust %>%
      mutate(
        gene = str_extract(name,"^.+(?=\\|)")
      ) -> heatmap_data_output_hclust
    
  }
  #@heatmap_data_output_hclust
  
  ##Heatmap clusters violins_arg
  {
    heatmap_data_output_hclust %>% 
      left_join(
        peptide_lfq_trinity_distinct_arg %>% 
          select(gene_position,MG132_MOCK,MGTG_MOCK),
        join_by(
          name == gene_position
        )
      ) %>%
      drop_na(gene) %>% 
      mutate(
        ind = fct(as.character(ind),levels = c("1","2","3","4","5","6","7"))
      ) %>% 
      rename(
        "gene_position" = "name"
      ) %>% 
      pivot_longer(MG132_MOCK:MGTG_MOCK) %>% 
      ggplot(
        aes(x= factor(name, levels = c("MG132_MOCK","MGTG_MOCK")), y = value, fill = factor(name, levels = c("MG132_MOCK","MGTG_MOCK")))
      )+
      geom_violin(linewidth = 0.2,scale="width")+
      geom_boxplot(width = 0.1,linewidth = 0.2,outlier.size = 0)+
      labs(x="Experiments", y="Z-score",fill="name")+
      stat_compare_means(
        comparisons = list(c("MG132_MOCK","MGTG_MOCK")),
        label = "p.signif",
        method="t.test",
        label.y = c(2.7),
        size = 3
      )+
      facet_wrap(~ind, strip.position = "top",nrow=1)+
      scale_fill_manual(
        values =c(wesanderson::wes_palette(5,name="AsteroidCity1",type="discrete")[c(1,3)]),
        labels = c("MG132/MOCK","MGTG/MOCK")
      )+
      lims(y = c(-1.5,3.2))+
      theme_classic2()+
      theme(
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 8,color="black"),
        legend.key.width = unit(0.2,"cm"),
        legend.key.height = unit(0.2,"cm"),
        legend.margin = margin(t = 0, unit="cm"),
        legend.box.margin = margin(t = 0, unit="cm"),
        legend.box.spacing = unit(2, "pt"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 8,color="black"),
        axis.title = element_text(size = 8,color="black"),
        strip.text= element_text(size = 8,color="black"),
        axis.title.x = element_blank()
      )
    ggsave("./figures/heatmap_clusters.png", width=10, height = 5, units = "cm", dpi = 600)
    
    RColorBrewer::display.brewer.all()
    wesanderson::wes_palette()
    image(volcano, col = wesanderson::wes_palette(21,name="AsteroidCity1",type="continuous"))
    wesanderson::wes_palette(5,name="AsteroidCity1",type="discrete")
  }
  #@heatmap_clusters.png
  
  ##Heatmap clusters violins_free
  {
    heatmap_data_output_hclust %>% 
      left_join(
        peptide_lfq_trinity_distinct_free %>% 
          select(gene_position,MG132_MOCK,MGTG_MOCK),
        join_by(
          name == gene_position
        )
      ) %>%
      drop_na(gene) %>% 
      mutate(
        ind = fct(as.character(ind),levels = c("1","2","3","4","5","6","7"))
      ) %>% 
      rename(
        "gene_position" = "name"
      ) %>% 
      pivot_longer(MG132_MOCK:MGTG_MOCK) %>% 
      ggplot(
        aes(x= factor(name, levels = c("MG132_MOCK","MGTG_MOCK")), y = value, fill = factor(name, levels = c("MG132_MOCK","MGTG_MOCK")))
      )+
      geom_violin(linewidth = 0.2,scale="width")+
      geom_boxplot(width = 0.1,linewidth = 0.2,outlier.size = 0)+
      labs(x="Experiments", y="Z-score",fill="name")+
      stat_compare_means(
        comparisons = list(c("MG132_MOCK","MGTG_MOCK")),
        label = "p.signif",
        method="t.test",
        label.y = c(2.7),
        size = 3
      )+
      facet_wrap(~ind, strip.position = "top",nrow=1)+
      scale_fill_manual(
        values =c(wesanderson::wes_palette(5,name="AsteroidCity1",type="discrete")[c(1,3)]),
        labels = c("MG132/MOCK","MGTG/MOCK")
      )+
      lims(y = c(-1.5,3.2))+
      theme_classic2()+
      theme(
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 8,color="black"),
        legend.key.width = unit(0.2,"cm"),
        legend.key.height = unit(0.2,"cm"),
        legend.margin = margin(t = 0, unit="cm"),
        legend.box.margin = margin(t = 0, unit="cm"),
        legend.box.spacing = unit(2, "pt"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 8,color="black"),
        axis.title = element_text(size = 8,color="black"),
        strip.text= element_text(size = 8,color="black"),
        axis.title.x = element_blank()
      )
    ggsave("./figures/heatmap_clusters_free.png", width=10, height = 5, units = "cm", dpi = 600)
  }
  #@heatmap_clusters_free.png
  
  ##heatmap for MGTG/MG132
  {
    keyref <- keys(org.Hs.eg.db, keytype = "SYMBOL")
    AnnotationDbi::select(
      org.Hs.eg.db, 
      keys=keyref,   ###gene name vector
      columns=c("SYMBOL","ENTREZID","GENENAME"),
      keytype="SYMBOL"
    ) -> gene_refs
    
    
    peptide_lfq_trinity_distinct_arg %>% 
      mutate(
        MGTG_MG132 = MGTG-MG132
      ) %>%
      left_join(
        arginylome_validation_matrix3_sequence_distinct_trinity %>% 
          select(gene_position,loc,gene),
        join_by(gene_position == gene_position)
      ) %>% 
      select(gene,MGTG_MG132) %>% 
      arrange(desc(MGTG_MG132)) %>% 
      distinct(gene,.keep_all = TRUE) %>% 
      left_join(
        gene_refs,
        join_by(gene == SYMBOL)
      ) %>% 
      drop_na()-> temp_data
    
    temp_data %>% pull(MGTG_MG132)->temp1
    temp_data %>% pull(ENTREZID)->temp2
      
    names(temp1) <- temp2
    
    temp_data %>% pull(MGTG_MG132)->temp3
    temp_data %>% pull(gene)->temp4
    names(temp3) <- temp4
    
    gseGO(
      geneList = temp3,
      OrgDb = org.Hs.eg.db,
      ont = "BP",
      keyType = "SYMBOL",
      eps = 1e-10,
      minGSSize    = 5,
      maxGSSize    = 500,
      pvalueCutoff = 1,
    ) -> gse_result_gobp
    
    gseGO(
      geneList = temp3,
      OrgDb = org.Hs.eg.db,
      ont = "CC",
      keyType = "SYMBOL",
      eps = 1e-10,
      minGSSize    = 5,
      maxGSSize    = 500,
      pvalueCutoff = 1,
    ) -> gse_result_gocc
    
    gseGO(
      geneList = temp3,
      OrgDb = org.Hs.eg.db,
      ont = "MF",
      keyType = "SYMBOL",
      eps = 1e-10,
      minGSSize    = 5,
      maxGSSize    = 500,
      pvalueCutoff = 1,
    ) -> gse_result_gomf
    
    tibble(gse_result_gomf@result) %>% 
      mutate(
        GO = "MF"
      ) %>% 
      bind_rows(
        tibble(gse_result_gocc@result) %>% 
          mutate(
            GO = "CC"
          ),
        tibble(gse_result_gobp@result) %>% 
          mutate(
            GO = "BP"
          )
      ) %>%
      filter(
        setSize >= 10
      ) -> gse_result_go
    
    calculateSimMatrix(
      gse_result_gobp$ID,
      orgdb="org.Hs.eg.db",
      ont="BP",
      method="Rel"
    ) -> simMatrix 
    
    scores <- setNames(-log10(gse_result_gobp$pvalue), gse_result_gobp$ID)
    reduceSimMatrix(
      simMatrix,
      scores,
      threshold=0.7,
      orgdb="org.Hs.eg.db",
    ) -> reducedTerms_bp
    
    calculateSimMatrix(
      gse_result_gocc$ID,
      orgdb="org.Hs.eg.db",
      ont="CC",
      method="Rel"
    ) -> simMatrix 
    
    scores <- setNames(-log10(gse_result_gocc$pvalue), gse_result_gocc$ID)
    reduceSimMatrix(
      simMatrix,
      scores,
      threshold=0.7,
      orgdb="org.Hs.eg.db",
    ) -> reducedTerms_cc
    
    calculateSimMatrix(
      gse_result_gomf$ID,
      orgdb="org.Hs.eg.db",
      ont="MF",
      method="Rel"
    ) -> simMatrix 
    
    scores <- setNames(-log10(gse_result_gomf$pvalue), gse_result_gomf$ID)
    reduceSimMatrix(
      simMatrix,
      scores,
      threshold=0.7,
      orgdb="org.Hs.eg.db",
    ) -> reducedTerms_mf
    
    
    bind_rows(
      reducedTerms_bp,
      reducedTerms_cc,
      reducedTerms_mf
    ) -> reducedTerms
    
    reducedTerms %>% 
      filter(
        termDispensability <=0
      ) %>% pull(go) -> reducedGO
    
    gse_result_go %>% 
      filter(
        ID %in% reducedGO,
        pvalue <= 0.05
      ) %>% 
      mutate(
        desc = paste0(GO,"::",Description)
      ) %>% 
      ggplot(
        aes(
          x=NES,
          y=reorder(desc,NES),
          fill = NES
        )
      )+
      geom_bar(
        stat = "identity"
      )+
      scale_fill_gradient2(
        high = muted("red"),
        mid = "white",
        low = muted("blue"),
      )+
      labs(
        x="NES (MGTG/MG132)",
        y="GO Description"
      )+
      theme_bw()+
      theme(
        axis.text.x = element_text(color = "black",size= 8),
        axis.text.y =element_text(color = "black",size= 7),
        axis.title = element_text(color = "black",size= 8),
        plot.title = element_text(color = "black",size= 8),
        legend.position = "none",
        axis.ticks = element_line(color = "black"),
        legend.title = element_text(color = "black",size= 8),
        legend.text = element_text(color = "black",size= 8)
      )
    ggsave(
      "./figures/4_12_gsea_tgovermg.png",height = 5, width =11,dpi = 600,units = "cm"
    )
    
    gseWP(
      geneList = temp1,
      organism = "Homo sapiens",
      pvalueCutoff = 1
    )
    
    gse
    
    gsePathway(
      geneList = temp1,
      pvalueCutoff = 1
    ) -> gse_result_reactome
    
    gse_result_reactome@result %>% 
      tibble() %>% 
      arrange(
        pvalue
      ) %>% 
      slice_head(n=10)%>% 
      ggplot(
        aes(
          x=NES,
          y=reorder(Description,NES),
          fill = NES
        )
      )+
      geom_bar(
        stat = "identity"
      )+
      scale_fill_gradient2(
        high = muted("red"),
        mid = "white",
        low = muted("blue"),
      )+
      labs(
        x="NES (MGTG/MG132)",
        y="Reactome"
      )+
      theme_bw()+
      theme(
        axis.text.x = element_text(color = "black",size= 8),
        axis.text.y =element_text(color = "black",size= 7),
        axis.title = element_text(color = "black",size= 8),
        plot.title = element_text(color = "black",size= 8),
        legend.position = "none",
        axis.ticks = element_line(color = "black"),
        legend.title = element_text(color = "black",size= 8),
        legend.text = element_text(color = "black",size= 8)
      )
    ggsave(
      "./figures/4_12_gsea_reactome_tgovermg.png",height = 5, width =7,dpi = 600,units = "cm"
    )
    
    gseKEGG(
      geneList = temp1,
      organism = "hsa",
      pvalueCutoff = 1
    )->gse_result
    
    gse_result@result %>% tibble() %>% view()
  }
  
  ##2.Heatmap Hclust GO REDUX
  {
    library(clusterProfiler)
    
    go_clusters <-tibble()
    for (i in 1:5){ ###5 is maximal cluster
      heatmap_data_output_hclust %>% 
        filter(ind == i) %>% 
        pull(gene) -> go_input
      
      enrichGO(
        go_input,
        'org.Hs.eg.db',
        keyType = "SYMBOL",
        ont = "ALL",
        pAdjustMethod = "fdr",
        readable = TRUE,
        pool = TRUE,
        qvalueCutoff = 0.05
      ) %>% .@result %>%  arrange(qvalue) %>% 
        #slice_head(n=10) %>% 
        mutate(
          clusters = paste0("Cluster ",i)
        )->go_clusters_part
      bind_rows(go_clusters,go_clusters_part) -> go_clusters
    }
    go_clusters %>% colnames()
    
    go_clusters %>% #filter(qvalue<0.05) %>% #loc 6, 27 
      filter(!is.na(qvalue)) %>% 
      filter(ONTOLOGY == "BP") %>% 
      nest(
        .by = clusters
      ) %>% 
      mutate(
        data = map(
          data,
          \(x) slice_head(x,n=10)
        )
      ) %>% 
      unnest(data) %>% 
      arrange(
        qvalue
      ) %>% 
      #slice_head(n=30) %>% 
      mutate(neg_log_qvalue = -log(qvalue)) %>% 
      mutate(truncated_desc = paste0(ID,": ",str_sub(Description,1,40))) %>% 
      ggplot(
        aes(
          x= clusters,
          y= reorder(truncated_desc,desc(qvalue)),
          size = neg_log_qvalue
        )
      )+geom_point(
        aes(
          fill = ONTOLOGY
        ),
        shape =21)+
      scale_fill_manual(values = c(wesanderson::wes_palette("Moonrise3")[1],wesanderson::wes_palette("Moonrise3")[3],wesanderson::wes_palette("Moonrise3")[5]))+
      labs(
        x="Localization",
        y="Gene Ontology"
      )+
      theme_classic2()+
      theme(
        axis.text.y = element_text(color = "black",size =6,hjust =0),
        axis.text.x = element_text(color = "black",size =7,angle = 90,vjust = 0.5,hjust=1),
        legend.position = "none"
      )
    ggsave(paste0("./figures/gobp_cluster.png"), width=20, height = 16, units = "cm", dpi = 600)
    
    go_clusters %>% select(ONTOLOGY,Description,ID, clusters,qvalue) %>% 
      pivot_wider(id_cols=c(ONTOLOGY,Description,ID),names_from = clusters,values_from = qvalue) %>% 
      write_tsv(
        "./export/ontology.tsv"
      )
  }
  
  ##2.Heatmap Hclust KEGG
  {
    
    
    heatmap_data_output_hclust$gene
    
    AnnotationDbi::keytypes(org.Hs.eg.db)
    
    heatmap_data_output_hclust %>% 
      left_join(
        as_tibble(AnnotationDbi::select(org.Hs.eg.db,keys=heatmap_data_output_hclust$gene, columns=c("SYMBOL","ENTREZID"),keytype="SYMBOL")),
        join_by(
          gene == SYMBOL
        ),
        multiple = "first"
      ) -> heatmap_data_output_hclust
    
    
    heatmap_data_output_hclust %>% 
      write_tsv(
        "./export/heatmap_data_output_hclust.tsv"
      )
    
    kegg_clusters <-tibble()
    for (i in 1:5){ ###5 is maximal cluster
      heatmap_data_output_hclust %>% 
        filter(ind == i) %>% 
        pull(ENTREZID) -> kegg_input
      
      enrichKEGG(
        kegg_input,
        organism = "hsa",
        keyType = "kegg",
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH"
      ) -> test_kegg
      
      test_kegg@result %>%  arrange(p.adjust) %>% 
        #slice_head(n=10) %>% 
        mutate(
          clusters = paste0("Cluster ",i)
        )->kegg_clusters_part
      
      bind_rows(kegg_clusters,kegg_clusters_part) -> kegg_clusters
    }
    
    kegg_clusters %>% write_tsv(
      "./export/kegg_cluster.tsv"
    )
  }
  
  ###ALL GENES into GO! 240213
  {
    
    arginylome_validation_matrix3_sequence_distinct_trinity %>% 
      pull(gene)->genes_arg_trinity
    
    arginylome_validation_matrix3_sequence_distinct_trinity %>% 
      left_join(
        as_tibble(AnnotationDbi::select(org.Hs.eg.db,keys=heatmap_data_output_hclust$gene, columns=c("SYMBOL","ENTREZID"),keytype="SYMBOL")),
        join_by(
          gene == SYMBOL
        ),
        multiple = "first"
      ) %>% pull(ENTREZID) -> genes_arg_trinity_entrezid
    
    enrichKEGG(
      genes_arg_trinity_entrezid,
      organism = "hsa",
      keyType = "kegg",
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH"
    ) %>% 
      setReadable(OrgDb = org.Hs.eg.db, keyType="ENTREZID")->enriched_kegg_all_arg
    
    enriched_kegg_all_arg@result %>% 
      tibble() -> enriched_kegg_all_arg
    
    enriched_kegg_all_arg %>% 
      write_tsv(
        "./export/enriched_kegg_all_arg.tsv"
      )
    
    
    BiocManager::valid()
    BiocManager::install("ReactomePA")
    library(ReactomePA)
    library(DOSE)
    get_wp_organisms()
    enrichWP(genes_arg_trinity_entrezid, organism = "Homo sapiens") %>% 
      setReadable(OrgDb = org.Hs.eg.db, keyType="ENTREZID")
    enrichDGN(genes_arg_trinity_entrezid) %>% 
      setReadable(OrgDb = org.Hs.eg.db, keyType="ENTREZID")
    enrichPathway(gene=genes_arg_trinity_entrezid, pvalueCutoff = 0.05, readable=TRUE)@result %>% 
      tibble() -> enriched_reactome_all_arg
    
    enriched_reactome_all_arg %>% 
      write_tsv(
        "./export/enriched_reactome_all_arg.tsv"
      )
    
    BiocManager::install("MeSHDbi")
    
    library(AnnotationHub)
    library(MeSHDbi)
    ah <- AnnotationHub(localHub=TRUE)
    hsa <- query(ah, c("MeSHDb", "Homo sapiens"))
    file_hsa <- hsa[[1]]
    setReadable(
      enrichWP(genes_arg_trinity_entrezid, organism = "Homo sapiens")
      , OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  }
  
  ####protease GO #240214
  {
    arginylome_validation_matrix3_sequence_distinct_trinity %>% 
      filter(
        site_type == "SP"
      ) %>% 
      distinct(gene) %>% 
      pull(gene)->genes_arg_trinity_signal
    
    enrichGO(
      genes_arg_trinity_signal,
      'org.Hs.eg.db',
      keyType = "SYMBOL",
      ont = "ALL",
      pAdjustMethod = "fdr",
      readable = TRUE,
      pool = TRUE,
      qvalueCutoff = 0.05
    ) %>% .@result %>%  arrange(qvalue) ->go_signal
    
    go_signal %>% #filter(qvalue<0.05) %>% #loc 6, 27 
      filter(!is.na(qvalue)) %>% 
      arrange(
        qvalue
      ) %>% 
      slice_head(n=10) %>% 
      mutate(neg_log_qvalue = -log(qvalue)) %>% 
      mutate(truncated_desc = paste0(ID,": ",str_sub(Description,1,50))) %>% 
      ggplot(
        aes(
          x=ONTOLOGY,
          y= reorder(truncated_desc,desc(qvalue)),
          size = neg_log_qvalue
        )
      )+geom_point(
        aes(
          fill = ONTOLOGY
        ),
        shape =21)+
      scale_fill_manual(values = c(wesanderson::wes_palette("Moonrise3")[1],wesanderson::wes_palette("Moonrise3")[3],wesanderson::wes_palette("Moonrise3")[5]))+
      labs(
        x="Localization",
        y="Gene Ontology"
      )+
      theme_classic2()+
      theme(
        axis.text.y = element_text(color = "black",size =6,hjust =0),
        axis.text.x = element_text(color = "black",size =7,angle = 90,vjust = 0.5,hjust=1),
        legend.position = "none"
      )
    ggsave(paste0("./figures/go_signal.png"), width=10, height = 10, units = "cm", dpi = 600)
    
    
    arginylome_validation_matrix3_sequence_distinct_trinity %>% 
      filter(
        str_detect(cleavage_type,"^C14")
      ) %>% 
      distinct(gene) %>% 
      pull(gene) -> genes_arg_trinity_caspase
    
    enrichGO(
      genes_arg_trinity_caspase,
      'org.Hs.eg.db',
      keyType = "SYMBOL",
      ont = "ALL",
      pAdjustMethod = "fdr",
      readable = TRUE,
      pool = TRUE,
      qvalueCutoff = 0.05
    ) %>% .@result %>%  arrange(qvalue) ->go_caspase
    
    go_caspase %>% #filter(qvalue<0.05) %>% #loc 6, 27 
      filter(!is.na(qvalue)) %>% 
      arrange(
        qvalue
      ) %>% 
      slice_head(n=10) %>% 
      mutate(neg_log_qvalue = -log(qvalue)) %>% 
      mutate(truncated_desc = paste0(ID,": ",str_sub(Description,1,50))) %>% 
      ggplot(
        aes(
          x=ONTOLOGY,
          y= reorder(truncated_desc,desc(qvalue)),
          size = neg_log_qvalue
        )
      )+geom_point(
        aes(
          fill = ONTOLOGY
        ),
        shape =21)+
      scale_fill_manual(values = c(wesanderson::wes_palette("Moonrise3")[1],wesanderson::wes_palette("Moonrise3")[3],wesanderson::wes_palette("Moonrise3")[5]))+
      labs(
        x="Localization",
        y="Gene Ontology"
      )+
      theme_classic2()+
      theme(
        axis.text.y = element_text(color = "black",size =6,hjust =0),
        axis.text.x = element_text(color = "black",size =7,angle = 90,vjust = 0.5,hjust=1),
        legend.position = "none"
      )
    ggsave(paste0("./figures/go_caspase.png"), width=10, height = 10, units = "cm", dpi = 600)
  }
  
  
  ###Localization chord diagram
  {
    library(circlize)
    
    peptide_lfq_trinity_distinct_arg %>% 
      left_join(
        arginylome_validation_matrix3_sequence_distinct_trinity %>% 
          select(gene_position,loc)
      ) %>%
      mutate(
        loc = if_else(loc == "Endoplasmic reticulum","ER",loc)
      ) %>% 
      mutate(
        loc = fct(loc,level = c("Extracellular","Cell membrane","Lysosome/Vacuole","Golgi apparatus","ER","Cytoplasm","Nucleus","Mitochondrion"))
      ) %>% 
      left_join(
        heatmap_data_output_hclust %>% select(name, ind) %>% 
          rename(
            gene_position = name,
            cluster = ind
          )
      ) %>%
      summarize(
        .by = c(loc,cluster),
        count = n()
      ) -> data_for_chord_diagram
    
    data_for_chord_diagram %>% distinct(loc)
    c(
      "Extracellular"="#3B9AB2",
      "Cell membrane"="#63ADBE",
      "Lysosome/Vacuole"="#9EBE91",
      "Golgi apparatus"="#EBCC2A",
      "Endoplasmic reticulum"="#E4B80E",
      "Cytoplasm"="#D5D5D3",
      "Nucleus"="#F21A00",
      "Mitochondrion"="#4DAF4A"
    )
    
    data_for_chord_diagram %>% 
      mutate(
        loc = str_replace(loc,"Cell membrane","CM"),
        loc = str_replace(loc,"Extracellular","EC"),
        loc = str_replace(loc,"Golgi apparatus","GA"),
        loc = str_replace(loc,"Lysosome/Vacuole","LV"),
        loc = str_replace(loc,"Mitochondrion","MT"),
        loc = str_replace(loc,"Nucleus","NU"),
        loc = str_replace(loc,"Cytoplasm","CP")
      ) %>% 
      drop_na(loc)-> data_for_chord_diagram
    
    png(
      paste0("./figures/cluster_loc.png"), width = 15, height = 15, 
      units = "cm", res = 600
    )
    data_for_chord_diagram %>% 
      mutate(
        loc = fct(loc, levels = c("NU","CP","ER","LV","CM","EC","MT","GA")),
        cluster = fct(as.character(cluster),levels = c("1","2","3","4","5"))
      ) %>% 
      arrange(
        loc,cluster) %>% 
      chordDiagram(
        grid.col = c(
          "EC"="#3B9AB2",
          "CM"="#63ADBE",
          "LV"="#9EBE91",
          "GA"="#EBCC2A",
          "ER"="#E4B80E",
          "CP"="#D5D5D3",
          "NU"="#F21A00",
          "MT"="#4DAF4A"
        ),
        annotationTrack = c("grid"),
        link.sort=TRUE,
        preAllocateTracks = list(track.height = max(strwidth(data_for_chord_diagram %>% distinct(loc) %>% pull(loc))))
      )
    circos.track(track.index = 1, panel.fun = function(x, y) {
      circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
                  facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
    }, bg.border = NA)
    dev.off()
    heatmap_data_output_hclust
    
    data_for_chord_diagram %>% 
      pivot_wider(
        id_cols = loc,
        names_from = cluster,
        values_from = count
      ) %>% 
      select(
        loc,`1`,`2`,`3`,`4`,`5`
      ) %>% 
      mutate(
        across(
          `1`:`5`,
          \(x) replace_na(x,0)
        )
      ) %>% drop_na() %>% 
      arrange(loc) %>% 
      column_to_rownames(var="loc") %>% 
      chisq.test() -> data_loc_chisq
    
    library(corrplot)
    corrplot(
      data_loc_chisq$residuals,
      is.cor = FALSE
    )
    data_loc_contrib <- 100*data_loc_chisq$residuals^2/data_loc_chisq$statistic
    png(
      paste0("./figures/cluster_loc_sq.png"), width = 20, height = 10, 
      units = "cm", res = 600
    )
    corrplot(
      t(data_loc_contrib), is.cor = FALSE,
      tl.col = "black",
      tl.srt = 30
    )
    dev.off()
    data_loc_contrib %>% 
      tibble()
    
  }
  #@cluster_loc.png
  #@cluster_loc_sq.png
  
  ##3. localization violin
  #@peptide_lfq_trinity_distinct_arg
  #@peptide_lfq_trinity_distinct_free
  #@arginylome_validation_matrix3_sequence_distinct_trinity
  {
    peptide_lfq_trinity_distinct_arg %>% 
      left_join(
        arginylome_validation_matrix3_sequence_distinct_trinity %>% 
          select(gene_position,loc)
      ) %>%
      mutate(
        loc = if_else(loc == "Endoplasmic reticulum","ER",loc)
      ) %>% 
      mutate(
        loc = fct(loc,level = c("Extracellular","Cell membrane","Lysosome/Vacuole","Golgi apparatus","ER","Cytoplasm","Nucleus","Mitochondrion"))
      ) %>%
      drop_na(loc) %>% 
      pivot_longer(MG132_MOCK:MGTG_MOCK) %>% 
      ggplot()+
      aes(x= factor(name, levels = c("MG132_MOCK","MGTG_MOCK")), y = value, fill = factor(name, levels = c("MG132_MOCK","MGTG_MOCK")))+
      geom_violin(linewidth = 0.2,scale="width")+
      geom_boxplot(width = 0.1,linewidth = 0.2,outlier.size = 0)+
      stat_compare_means(
        comparisons = list(c("MG132_MOCK","MGTG_MOCK")),
        label = "p.signif",
        method="t.test",
        #label.y = c(2,2.7,3.4),
        size = 3
      )+
      labs(x="Experiments", y="Z-score", fill = "Experiment")+
      facet_wrap(~loc, strip.position = "top",nrow=2)+
      scale_fill_manual(values =c(wesanderson::wes_palette(5,name="AsteroidCity1",type="discrete")[c(1,3)]))+
      #lims(y = c(-1.5,2.3))+
      theme_classic2()+
      theme(
        legend.position = "top",
        legend.title = element_text(size = 8,color="black"),
        legend.text = element_text(size = 8,color="black"),
        legend.key.width = unit(0.2,"cm"),
        legend.key.height = unit(0.2,"cm"),
        legend.margin = margin(t = 0, unit="cm"),
        legend.box.margin = margin(t = 0, unit="cm"),
        legend.box.spacing = unit(2, "pt"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 8,color="black"),
        axis.title = element_text(size = 8,color="black"),
        strip.text= element_text(size = 8,color="black"),
        axis.title.x = element_blank()
      )
    ggsave("./figures/4_loc_violins_arg.png", width=15, height = 15, units = "cm", dpi = 600)
    
    peptide_lfq_trinity_distinct_free %>% 
      left_join(
        arginylome_validation_matrix3_sequence_distinct_trinity %>% 
          select(gene_position,loc)
      ) %>%
      mutate(
        loc = if_else(loc == "Endoplasmic reticulum","ER",loc)
      ) %>% 
      mutate(
        loc = fct(loc,level = c("Extracellular","Cell membrane","Lysosome/Vacuole","Golgi apparatus","ER","Cytoplasm","Nucleus","Mitochondrion"))
      ) %>%
      drop_na(loc) %>% 
      pivot_longer(MG132_MOCK:MGTG_MOCK) %>% 
      ggplot()+
      aes(x= factor(name, levels = c("MG132_MOCK","MGTG_MOCK")), y = value, fill = factor(name, levels = c("MG132_MOCK","MGTG_MOCK")))+
      geom_violin(linewidth = 0.2,scale="width")+
      geom_boxplot(width = 0.1,linewidth = 0.2,outlier.size = 0)+
      stat_compare_means(
        comparisons = list(c("MG132_MOCK","MGTG_MOCK")),
        label = "p.signif",
        method="t.test",
        #label.y = c(2,2.7,3.4),
        size = 3
      )+
      labs(x="Experiments", y="Z-score", fill = "Experiment")+
      facet_wrap(~loc, strip.position = "top",nrow=2)+
      scale_fill_manual(values =c(wesanderson::wes_palette(5,name="AsteroidCity1",type="discrete")[c(1,3)]))+
      #lims(y = c(-1.5,2.3))+
      theme_classic2()+
      theme(
        legend.position = "top",
        legend.title = element_text(size = 8,color="black"),
        legend.text = element_text(size = 8,color="black"),
        legend.key.width = unit(0.2,"cm"),
        legend.key.height = unit(0.2,"cm"),
        legend.margin = margin(t = 0, unit="cm"),
        legend.box.margin = margin(t = 0, unit="cm"),
        legend.box.spacing = unit(2, "pt"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 8,color="black"),
        axis.title = element_text(size = 8,color="black"),
        strip.text= element_text(size = 8,color="black"),
        axis.title.x = element_blank()
      )
    ggsave("./figures/4_loc_violins_free.png", width=10, height = 15, units = "cm", dpi = 600)
    
    peptide_lfq_trinity_distinct_arg %>% 
      left_join(
        arginylome_validation_matrix3_sequence_distinct_trinity %>% 
          select(gene_position,loc)
      ) %>%
      mutate(
        loc = if_else(loc == "Endoplasmic reticulum","ER",loc),
        loc = if_else(loc == "Extracellular","Exc",loc),
        loc = if_else(loc == "Cell membrane","Mem",loc),
        loc = if_else(loc == "Lysosome/Vacuole","LV",loc),
        loc = if_else(loc == "Golgi apparatus","GG",loc),
        loc = if_else(loc == "Cytoplasm","Cyt",loc),
        loc = if_else(loc == "Nucleus","Nuc",loc),
        loc = if_else(loc == "Mitochondrion","Mit",loc),
      ) %>% 
      mutate(
        loc = fct(loc,level = c("Exc","Mem","LV","GG","ER","Cyt","Nuc","Mit"))
      ) %>%
      filter(
        !(loc %in% c("LV","GG"))
      ) %>% 
      drop_na(loc) %>% 
      pivot_longer(MG132:MOCK) %>% 
      ggplot()+
      aes(x= factor(name, levels = c("MOCK","MG132","MGTG")), y = value, fill = factor(name, levels = c("MOCK","MG132","MGTG")))+
      #geom_violin(linewidth = 0.2,scale="width")+
      geom_boxplot(
        #width = 0.1,linewidth = 0.2,outlier.size = 0
      )+
      stat_compare_means(
        comparisons = list(c("MOCK","MG132"),c("MG132","MGTG"),c("MOCK","MGTG")),
        label = "p.signif",
        method="t.test",
        #label.y = c(2,2.7,3.4),
        size = 3
      )+
      labs(x="Experiments", y="Z-score", fill = "Experiment")+
      facet_wrap(~loc, strip.position = "top",nrow=1)+
      scale_fill_manual(values =c(wesanderson::wes_palette(5,name="AsteroidCity1",type="discrete")[c(1,3,5)]))+
      lims(y = c(-2,4.5))+
      theme_classic2()+
      theme(
        legend.position = "bottom",
        legend.title = element_text(size = 8,color="black"),
        legend.text = element_text(size = 8,color="black"),
        legend.key.width = unit(0.2,"cm"),
        legend.key.height = unit(0.2,"cm"),
        legend.margin = margin(t = 0, unit="cm"),
        legend.box.margin = margin(t = 0, unit="cm"),
        legend.box.spacing = unit(2, "pt"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 8,color="black"),
        axis.title = element_text(size = 8,color="black"),
        strip.text= element_text(size = 8,color="black"),
        axis.title.x = element_blank()
      )
    ggsave("./figures/4_loc_violins_arg_all.png", width=8.5, height = 6, units = "cm", dpi = 600)
  }
  
  ##4. localization go
  {
    arginylome_validation_matrix3_sequence_distinct_trinity %>% 
      select(gene,loc) -> genes_loc
    
    locs =c("Extracellular","Cell membrane","Lysosome/Vacuole","Golgi apparatus","Endoplasmic reticulum","Cytoplasm","Nucleus","Mitochondrion")
    
    
      for (i in 1:length(locs)){
        genes_loc %>% 
          filter(loc == locs[i]) %>% 
          pull(gene) -> go_input_loc
        
        enrichGO(
          go_input_loc,
          'org.Hs.eg.db',
          keyType = "SYMBOL",
          ont = "ALL",
          pAdjustMethod = "fdr",
          readable = TRUE,
          pool = TRUE
        ) %>% .@result %>%  slice_head(n=5) %>% arrange(qvalue) %>%
          ggplot(
            aes(
              x= reorder(Description,qvalue), y= -log(qvalue), fill = ONTOLOGY
            )
          )+
          geom_bar(stat="identity")+
          scale_x_discrete(labels= function(x) str_sub(x, start = 1, end = 50))+
          scale_fill_manual(values = c(wesanderson::wes_palette("Moonrise3")[1],wesanderson::wes_palette("Moonrise3")[3],wesanderson::wes_palette("Moonrise3")[5]))+
          xlab("GO Description")+
          ylab(expression(-log[10](Q-value)))+
          theme_classic2()+
          theme(
            axis.text.x = element_text(size = 5,angle=90,hjust=1,vjust=0.5),
            legend.position = "none"
          )
        ggsave(paste0("./figures/cluster_loc_",i,".png"), width=6, height = 10, units = "cm", dpi = 600)
      }
    
  }
  
  ###4. GO for clusters REDUX
  {
    arginylome_validation_matrix3_sequence_distinct_trinity %>% 
      select(gene,loc) -> genes_loc
    
    locs =c("Extracellular","Cell membrane","Lysosome/Vacuole","Golgi apparatus","Endoplasmic reticulum","Cytoplasm","Nucleus","Mitochondrion")
    
    go_locs <-tibble()
    for (i in 1:length(locs)){
      genes_loc %>% 
        filter(loc == locs[i]) %>% 
        pull(gene) -> go_input_loc
      
      enrichGO(
        go_input_loc,
        'org.Hs.eg.db',
        keyType = "SYMBOL",
        ont = "ALL",
        pAdjustMethod = "fdr",
        readable = TRUE,
        pool = TRUE,
        qvalueCutoff = 0.05
      ) %>% .@result %>%  arrange(qvalue) %>% 
        #slice_head(n=10) %>% 
        mutate(
          loc = locs[i]
        )->go_locs_part
      bind_rows(go_locs,go_locs_part) -> go_locs
    }
    go_locs %>% colnames()
    
    go_locs %>% #filter(qvalue<0.05) %>% #loc 6, 27 
      filter(!is.na(qvalue)) %>% 
      arrange(
        qvalue
      ) %>% 
      slice_head(n=30) %>% 
      mutate(neg_log_qvalue = -log(qvalue)) %>% 
      mutate(truncated_desc = paste0(ID,"\n",str_sub(Description,1,30))) %>% 
      ggplot(
        aes(
          x= fct(loc,levels = c("Extracellular","Cell membrane","Lysosome/Vacuole","Golgi apparatus","Endoplasmic reticulum","Cytoplasm","Nucleus","Mitochondrion")), 
          y= reorder(truncated_desc,desc(qvalue)),
          size = neg_log_qvalue
        )
      )+geom_point(
        aes(
          fill = ONTOLOGY
        ),
        shape =21)+
      scale_fill_manual(values = c(wesanderson::wes_palette("Moonrise3")[1],wesanderson::wes_palette("Moonrise3")[3],wesanderson::wes_palette("Moonrise3")[5]))+
      labs(
        x="Localization",
        y="Gene Ontology"
      )+
      theme_classic2()+
      theme(
        axis.text.y = element_text(color = "black",size =6,hjust =0),
        axis.text.x = element_text(color = "black",size =7,angle = 90,vjust = 0.5,hjust=1),
        legend.position = "none"
      )
    ggsave(paste0("./figures/go_loc.png"), width=10, height = 20, units = "cm", dpi = 600)
  }
  
  ###4. GO for locs REDUX for main
  {
    arginylome_validation_matrix3_sequence_distinct_trinity %>% 
      select(gene,loc) -> go_input_loc
    
    locs =c("Extracellular","Cell membrane","Lysosome/Vacuole","Golgi apparatus","Endoplasmic reticulum","Cytoplasm","Nucleus","Mitochondrion")
    
    go_input_loc %>% distinct(loc)
    
    go_locs <-tibble()
    for (i in 1:length(locs)){
      genes_loc %>% 
        filter(loc == locs[i]) %>% 
        pull(gene) -> go_input_loc
      
      enrichGO(
        go_input_loc,
        'org.Hs.eg.db',
        keyType = "SYMBOL",
        ont = "ALL",
        pAdjustMethod = "fdr",
        readable = TRUE,
        pool = TRUE,
        qvalueCutoff = 0.05
      ) %>% .@result %>%  arrange(qvalue) %>% 
        slice_head(n=3) %>%  ####this is important
        mutate(
          loc = locs[i]
        )->go_locs_part
      bind_rows(go_locs,go_locs_part) -> go_locs
    }
    go_locs %>% #filter(qvalue<0.05) %>% #loc 6, 27 
      filter(!is.na(qvalue)) %>% 
      arrange(
        qvalue
      ) %>% 
      slice_head(n=30) %>% 
      mutate(neg_log_qvalue = -log(qvalue)) %>% 
      mutate(truncated_desc = paste0(ID,"\n",str_sub(Description,1,40))) %>% 
      ggplot(
        aes(
          x= fct(loc,levels = c("Extracellular","Cell membrane","Lysosome/Vacuole","Golgi apparatus","Endoplasmic reticulum","Cytoplasm","Nucleus","Mitochondrion")), 
          y= reorder(truncated_desc,desc(qvalue)),
          size = neg_log_qvalue
        )
      )+geom_point(
        aes(
          fill = ONTOLOGY
        ),
        shape =21)+
      scale_fill_manual(values = c(wesanderson::wes_palette("Moonrise3")[1],wesanderson::wes_palette("Moonrise3")[3],wesanderson::wes_palette("Moonrise3")[5]))+
      labs(
        x="Localization",
        y="Gene Ontology"
      )+
      theme_classic2()+
      theme(
        axis.text.y = element_text(color = "black",size =6,hjust =0),
        axis.text.x = element_text(color = "black",size =7,angle = 90,vjust = 0.5,hjust=1),
        legend.position = "none"
      )
    ggsave(paste0("./figures/go_loc_main.png"), width=8, height = 13, units = "cm", dpi = 600)
  }
  ##4. go for cluster! main!
  {
    heatmap_data_output_hclust %>% 
      filter(ind == i) %>% 
      pull(gene) -> go_input
    
    go_clusters <-tibble()
    for (i in 1:6){
      heatmap_data_output_hclust %>% 
        filter(ind == i) %>% 
        pull(gene) -> go_input
      
      enrichGO(
        go_input,
        'org.Hs.eg.db',
        keyType = "SYMBOL",
        ont = "ALL",
        pAdjustMethod = "fdr",
        readable = TRUE,
        pool = TRUE,
        qvalueCutoff = 0.05
      ) %>% .@result %>%  arrange(qvalue) %>% 
        #slice_head(n=10) %>% 
        mutate(
          clusters = paste0("Cluster ",i)
        )->go_clusters_part
      bind_rows(go_clusters,go_clusters_part) -> go_clusters
    }
    go_clusters %>%
      distinct(clusters)
    
    go_clusters %>% #filter(qvalue<0.05) %>% #loc 6, 27 
      filter(!is.na(qvalue)) %>% 
      arrange(
        qvalue
      ) %>% 
      slice_head(n=10) %>% 
      mutate(neg_log_qvalue = -log(qvalue)) %>% 
      mutate(truncated_desc = paste0(ID,"\n",str_sub(Description,1,40))) %>% 
      ggplot(
        aes(
          x= fct(clusters,levels = c("Cluster 1","Cluster 2","Cluster 3","Cluster 4","Cluster 5","Cluster 6")), 
          y= reorder(truncated_desc,desc(qvalue)),
          size = neg_log_qvalue
        )
      )+geom_point(
        aes(
          fill = ONTOLOGY
        ),
        shape =21)+
      scale_fill_manual(values = c(wesanderson::wes_palette("Moonrise3")[1],wesanderson::wes_palette("Moonrise3")[3],wesanderson::wes_palette("Moonrise3")[5]))+
      labs(
        x="Clusters",
        y="Gene Ontology"
      )+
      theme_classic2()+
      theme(
        axis.text.y = element_text(color = "black",size =6,hjust =0),
        axis.text.x = element_text(color = "black",size =7,angle = 90,vjust = 0.5,hjust=1),
        legend.position = "none"
      )
    ggsave(paste0("./figures/go_cluster_main.png"), width=8, height = 13, units = "cm", dpi = 600)
  }
  
  
  ##5. Network
  {
    arginylome_validation_matrix3_sequence_distinct_trinity %>% 
      pull(gene) %>% 
      writeClipboard()
    ##TO STRING
    
    arginylome_validation_matrix3_sequence_distinct_trinity %>% 
      select(gene,loc) %>% 
      write_tsv("./network/gene_loc.tsv")
    
    arginylome_validation_matrix3_sequence_distinct_trinity %>% 
      filter(loc == "Mitochondrion")
    
    arginylome_validation_matrix3_sequence_distinct_trinity %>% 
      write_tsv("./export/arginylome_validation_matrix3_sequence_distinct_trinity.tsv")
    
    peptide_lfq_trinity_distinct_arg %>% 
      write_tsv(
        "./export/peptide_lfq_trinity_distinct_arg_240122.tsv"
      )
  }
  
  ##6.HPA function
  #@peptide_lfq_trinity_distinct_arg
  #@peptide_lfq_trinity_distinct_free
  #@arginylome_validation_matrix3_sequence_distinct_trinity
  {
    peptide_lfq_trinity_distinct_arg %>% 
      left_join(
        arginylome_validation_matrix3_sequence_distinct_trinity %>% 
          select(gene_position,`Molecular function`),
        join_by(gene_position)
      ) %>%
      separate_longer_delim(
        `Molecular function`,", "
      ) %>% 
      pivot_longer(MG132:MOCK) %>% 
      ggplot()+
      aes(x= factor(name, levels = c("MOCK","MG132","MGTG")), y = value, fill = factor(name, levels = c("MOCK","MG132","MGTG")))+
      geom_violin(linewidth = 0.2,scale="width")+
      geom_boxplot(width = 0.1,linewidth = 0.2,outlier.size = 0)+
      stat_compare_means(
        comparisons = list(c("MOCK","MG132"),c("MG132","MGTG"),c("MOCK","MGTG")),
        label = "p.signif",
        method="t.test",
        #label.y = c(2,2.7,3.4),
        size = 3
      )+
      labs(x="Experiments", y="Z-score", fill = "Experiment")+
      facet_wrap(~`Molecular function`, strip.position = "top",nrow=5)+
      scale_fill_manual(values =c(wesanderson::wes_palette(5,name="AsteroidCity1",type="discrete")[c(1,2,3)]))+
      lims(y = c(-1.5,2.3))+
      theme_classic2()+
      theme(
        legend.position = "top",
        legend.title = element_text(size = 8,color="black"),
        legend.text = element_text(size = 8,color="black"),
        legend.key.width = unit(0.2,"cm"),
        legend.key.height = unit(0.2,"cm"),
        legend.margin = margin(t = 0, unit="cm"),
        legend.box.margin = margin(t = 0, unit="cm"),
        legend.box.spacing = unit(2, "pt"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 8,color="black"),
        axis.title = element_text(size = 8,color="black"),
        strip.text= element_text(size = 8,color="black"),
        axis.title.x = element_blank()
      )
    ggsave("./figures/4_6_hpa_MF.png", width=50, height = 50, units = "cm", dpi = 600)
    
    peptide_lfq_trinity_distinct_arg %>% 
      left_join(
        arginylome_validation_matrix3_sequence_distinct_trinity %>% 
          select(gene_position,`Biological process`),
        join_by(gene_position)
      ) %>%
      separate_longer_delim(
        `Biological process`,", "
      ) %>% 
      pivot_longer(MG132:MOCK) %>% 
      ggplot()+
      aes(x= factor(name, levels = c("MOCK","MG132","MGTG")), y = value, fill = factor(name, levels = c("MOCK","MG132","MGTG")))+
      geom_violin(linewidth = 0.2,scale="width")+
      geom_boxplot(width = 0.1,linewidth = 0.2,outlier.size = 0)+
      stat_compare_means(
        comparisons = list(c("MOCK","MG132"),c("MG132","MGTG"),c("MOCK","MGTG")),
        label = "p.signif",
        method="t.test",
        #label.y = c(2,2.7,3.4),
        size = 3
      )+
      labs(x="Experiments", y="Z-score", fill = "Experiment")+
      facet_wrap(~`Biological process`, strip.position = "top",nrow=5)+
      scale_fill_manual(values =c(wesanderson::wes_palette(5,name="AsteroidCity1",type="discrete")[c(1,2,3)]))+
      lims(y = c(-1.5,2.3))+
      theme_classic2()+
      theme(
        legend.position = "top",
        legend.title = element_text(size = 8,color="black"),
        legend.text = element_text(size = 8,color="black"),
        legend.key.width = unit(0.2,"cm"),
        legend.key.height = unit(0.2,"cm"),
        legend.margin = margin(t = 0, unit="cm"),
        legend.box.margin = margin(t = 0, unit="cm"),
        legend.box.spacing = unit(2, "pt"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 8,color="black"),
        axis.title = element_text(size = 8,color="black"),
        strip.text= element_text(size = 8,color="black"),
        axis.title.x = element_blank()
      )
    ggsave("./figures/4_6_hpa_BP.png", width=50, height = 50, units = "cm", dpi = 600)
    
    peptide_lfq_trinity_distinct_arg %>% 
      left_join(
        arginylome_validation_matrix3_sequence_distinct_trinity %>% 
          select(gene_position,`Biological process`),
        join_by(gene_position)
      ) %>%
      separate_longer_delim(
        `Biological process`,", "
      ) %>%
      summarize(
        .by = `Biological process`,
        count = n(),
        gene_pos = paste(gene_position,collapse = ";")
      ) %>% 
      write_tsv(
        "./export/4_6_hpa_BP.tsv"
      )
  }
  
  ###4.08.MOCK only PSM TILE
  {
    arginylome_validation_matrix3_sequence_distinct_trinity %>% 
      filter(
        MOCK >0
      ) %>% 
      select(
        gene_position,MOCK,MG132,MGTG
      ) %>% 
      arrange(
        desc(gene_position)
      ) %>% 
      pivot_longer(MOCK:MGTG) %>% 
      mutate(
        value = na_if(value, 0)
      ) %>% 
      ggplot(
        aes(y=fct(name,levels=c("MGTG","MG132","MOCK")), x= fct(gene_position,levels = unique(gene_position)), fill = value)
      )+
      geom_tile()+
      geom_text(
        aes(label = value)
      )+
      scale_fill_distiller(
        palette = "YlOrRd",
        direction = 1,
        na.value = "grey80"
      )+
      labs(
        y="Experiment",
        x="Arginylated Sites only identified in MOCK",
        fill = "#PSMs"
      )+
      theme_bw()+
      theme(
        axis.text.x = element_text(color = "black",size= 7,angle = 45,hjust = 1),
        axis.text.y =element_text(color = "black",size= 8),
        axis.title = element_text(color = "black",size= 8),
        plot.title = element_text(color = "black",size= 8),
        legend.position = "right",
        axis.ticks = element_line(color = "black"),
        legend.title = element_text(color = "black",size= 8),
        legend.text = element_text(color = "black",size= 8)
      )+
      guides(
        fill=guide_colourbar(barwidth=0.5, barheight = 4,label.position="right")
      )
    ggsave(
      "./figures/4_08_tile_mock.png",height = 5, width =12,dpi = 600,units = "cm"
    )
    
    arginylome_validation_matrix3_sequence_distinct_trinity %>% 
      filter(
        MOCK >0
      ) %>% 
      select(
        gene_position,argpair_MOCK,argpair_MG132,argpair_MGTG
      ) %>% 
      arrange(
        desc(gene_position)
      ) %>% 
      pivot_longer(cols=c(argpair_MOCK,argpair_MG132,argpair_MGTG)) %>% 
      mutate(
        value = na_if(value, 0)
      ) %>% 
      ggplot(
        aes(y=fct(name,levels=c("argpair_MGTG","argpair_MG132","argpair_MOCK")), x= fct(gene_position,levels = unique(gene_position)), fill = value)
      )+
      geom_tile()+
      geom_text(
        aes(label = value)
      )+
      scale_fill_distiller(
        palette = "YlOrRd",
        direction = 1,
        na.value = "grey80"
      )+
      scale_y_discrete(
        label = c("MGTG","MG132","MOCK")
      )+
      labs(
        y="Experiment",
        x="Arginylated Sites only identified in MOCK",
        fill = "#PSMs"
      )+
      theme_bw()+
      theme(
        axis.text.x = element_text(color = "black",size= 7,angle = 45,hjust = 1),
        axis.text.y =element_text(color = "black",size= 8),
        axis.title = element_text(color = "black",size= 8),
        plot.title = element_text(color = "black",size= 8),
        legend.position = "right",
        axis.ticks = element_line(color = "black"),
        legend.title = element_text(color = "black",size= 8),
        legend.text = element_text(color = "black",size= 8)
      )+
      guides(
        fill=guide_colourbar(barwidth=0.5, barheight = 4,label.position="right")
      )
    ggsave(
      "./figures/4_08_tile_mock_free.png",height = 5, width =12,dpi = 600,units = "cm"
    )
    
    arginylome_validation_matrix3_sequence_distinct_trinity %>% 
      filter(
        MOCK >0
      ) %>%
      distinct(gene_position) %>% 
      mutate(
        gene_position = str_replace_all(gene_position,"\\|","_")
      ) %>% 
      pull(gene_position) ->gene_pos_mock_only
    
    ### mock only violin
    peptide_lfq_trinity %>% 
      pivot_longer(MG132_1_Free:MOCK_2_Arg) %>% 
      drop_na(value) %>% 
      filter(
        gene_position %in% gene_pos_mock_only
      ) %>% 
      mutate(
        experiment = str_split_i(name,"_",1),
        arg = str_split_i(name,"_",3)
      ) %>% 
      ggplot(
        aes(x= factor(experiment, levels = c("MOCK","MG132","MGTG")), y = value, fill = factor(experiment, levels = c("MOCK","MG132","MGTG")))
      )+
      geom_violin(linewidth = 0.2,scale="width")+
      geom_boxplot(width = 0.1,linewidth = 0.2,outlier.size = 0)+
      labs(x="Experiments", y="Z-score",fill="name")+
      stat_compare_means(
        comparisons = list(c("MOCK","MG132"),c("MOCK","MGTG"),c("MG132","MGTG")),
        label = "p.format",
        method="t.test",
        label.y = c(4,4.8,4.3),
        size = 3
      )+
      scale_fill_manual(
        values =c(wesanderson::wes_palette(5,name="Zissou1",type="discrete")[c(1,3,5)]),
        labels = c("MOCK","MG132","MGTG")
      )+
      lims(
        y=c(-3.1,5.5)
      )+
      theme_classic2()+
      theme(
        legend.position = "none",
        # legend.title = element_blank(),
        # legend.text = element_text(size = 8,color="black"),
        # legend.key.width = unit(0.2,"cm"),
        # legend.key.height = unit(0.2,"cm"),
        # legend.margin = margin(t = 0, unit="cm"),
        # legend.box.margin = margin(t = 0, unit="cm"),
        # legend.box.spacing = unit(2, "pt"),
        #axis.ticks.x = element_blank(),
        #axis.text.x = element_blank(),
        axis.text = element_text(size = 8,color="black"),
        axis.title = element_text(size = 8,color="black"),
        strip.text= element_text(size = 8,color="black"),
        axis.title.x = element_blank(),
        axis.ticks = element_line(color = "black")
      )
    ggsave(
      "./figures/4_08_mock_only_violin.png",height = 5, width =5,dpi = 600,units = "cm"
    )
  }
  RColorBrewer::display.brewer.all()
  
  ##specific site boxplot
  {
    #CALR
    peptide_lfq_trinity %>% 
      filter(
        gene_position %in% c("CALR_18")#,"P4HB_18")
      ) %>% 
      select(
        gene_position,MG132_1_Free:MOCK_2_Arg
      ) %>% 
      pivot_longer(
        MG132_1_Free:MOCK_2_Arg
      ) %>% 
      mutate(
        experiment = fct(str_split_i(name,"_",1),levels = c("MOCK","MG132","MGTG")),
        arg = str_split_i(name,"_",3)
      ) %>% 
      ggplot(
        aes(x=experiment,y=value)
      )+
      geom_boxplot(
        aes(fill = experiment)
      )+
      scale_fill_manual(
        values = c(wesanderson::wes_palette(5,name="Zissou1",type="discrete")[c(1,3,5)])
      )+
      stat_compare_means(
        comparisons = list(c("MOCK","MG132"),c("MG132","MGTG"),c("MOCK","MGTG")),
        label = "p.format",
        method="t.test",
        #label.y = c(4,4.8,4.3),
        size = 2.5
      )+
      facet_wrap(
        vars(arg)
      )+
      lims(y=c(-1.5,4.2))+
      labs(
        x="Experiment",
        y="Z-score"
      )+
      theme_classic2()+
      theme(
        axis.text = element_text(color = "black",size= 8),
        axis.text.x = element_text(color = "black",size= 8,angle = 45,hjust=1),
        axis.title = element_text(color = "black",size= 8),
        plot.title = element_text(color = "black",size= 8),
        legend.position = "none",
        axis.ticks = element_line(color = "black"),
        #axis.title = element_blank(),
        axis.title.x = element_blank(),
        strip.text.x = element_blank()
      )
    ggsave(paste0("./figures/4_08_box_calr.png"),height = 5, width =6,dpi = 600,units = "cm")
    
    #p4hb
    peptide_lfq_trinity %>% 
      filter(
        gene_position %in% c("P4HB_18")
      ) %>% 
      select(
        gene_position,MG132_1_Free:MOCK_2_Arg
      ) %>% 
      pivot_longer(
        MG132_1_Free:MOCK_2_Arg
      ) %>% 
      mutate(
        experiment = fct(str_split_i(name,"_",1),levels = c("MOCK","MG132","MGTG")),
        arg = str_split_i(name,"_",3)
      ) %>% 
      ggplot(
        aes(x=experiment,y=value)
      )+
      geom_boxplot(
        aes(fill = experiment)
      )+
      scale_fill_manual(
        values = c(wesanderson::wes_palette(5,name="Zissou1",type="discrete")[c(1,3,5)])
      )+
      stat_compare_means(
        comparisons = list(c("MOCK","MG132"),c("MOCK","MGTG"),c("MG132","MGTG")),
        label = "p.format",
        method="t.test",
        #label.y = c(4,4.8,4.3),
        size = 2.5
      )+
      facet_wrap(
        vars(arg)
      )+
      labs(
        x="Experiment",
        y="Z-score"
      )+
      lims(y=c(0,5))+
      theme_classic2()+
      theme(
        axis.text = element_text(color = "black",size= 8),
        axis.text.x = element_text(color = "black",size= 8,angle = 45,hjust=1),
        axis.title = element_text(color = "black",size= 8),
        plot.title = element_text(color = "black",size= 8),
        legend.position = "none",
        axis.ticks = element_line(color = "black"),
        #axis.title = element_blank(),
        axis.title.x = element_blank(),
        strip.text.x = element_blank()
      )
    ggsave(paste0("./figures/4_08_box_p4hb.png"),height = 5, width =6,dpi = 600,units = "cm")
  }
  
  ##9. GO for each experiment
  {
    arginylome_validation_matrix3_sequence_distinct_trinity %>% 
      filter(MOCK >0) %>% pull(gene) -> go_input
    
    enrichGO(
      go_input,
      'org.Hs.eg.db',
      keyType = "SYMBOL",
      ont = "ALL",
      pAdjustMethod = "fdr",
      readable = TRUE,
      pool = TRUE,
      qvalueCutoff = 0.05
    ) %>% .@result %>%  arrange(qvalue) -> go_enriched_mock
    
    arginylome_validation_matrix3_sequence_distinct_trinity %>% 
      filter(MG132 >0) %>% pull(gene) -> go_input
    
    enrichGO(
      go_input,
      'org.Hs.eg.db',
      keyType = "SYMBOL",
      ont = "ALL",
      pAdjustMethod = "fdr",
      readable = TRUE,
      pool = TRUE,
      qvalueCutoff = 0.05
    ) %>% .@result %>%  arrange(qvalue) -> go_enriched_mg
    
    arginylome_validation_matrix3_sequence_distinct_trinity %>% 
      filter(MGTG >0) %>% pull(gene) -> go_input
    
    enrichGO(
      go_input,
      'org.Hs.eg.db',
      keyType = "SYMBOL",
      ont = "ALL",
      pAdjustMethod = "fdr",
      readable = TRUE,
      pool = TRUE,
      qvalueCutoff = 0.05
    ) %>% .@result %>%  arrange(qvalue) -> go_enriched_mgtg
    go_enriched_mock %>% 
      mutate(
        experiment = "MOCK"
      ) %>% 
      bind_rows(
        go_enriched_mg %>% 
          mutate(
            experiment = "MG132"
          )
      )%>% 
      bind_rows(
        go_enriched_mgtg %>% 
          mutate(
            experiment = "MGTG"
          )
      ) -> go_enriched_all
    
    go_enriched_all %>% 
      mutate(
        fold_enrich = (as.integer(str_split_i(GeneRatio,"\\/",1))/as.integer(str_split_i(GeneRatio,"\\/",2)))/(as.integer(str_split_i(BgRatio,"\\/",1))/as.integer(str_split_i(BgRatio,"\\/",2)))
      )-> go_enriched_all
    
    go_enriched_all %>% 
      select(ONTOLOGY,ID,Description,experiment,qvalue,fold_enrich) %>% 
      pivot_wider(
        id_cols = c(ONTOLOGY,ID,Description),
        names_from = experiment,
        values_from = qvalue
      ) %>%
      filter(
        MOCK < 0.001 | MG132 < 0.001 | MGTG <0.001
      ) %>% 
      pivot_longer(
        MOCK:MGTG
      ) %>% 
      arrange(value) -> go_enriched_all_fig
    
    go_enriched_all_fig %>%
      mutate(
        Description = str_replace(Description,"(?<=[:blank:]in)[:blank:]","\n"),
        Description = str_replace(Description,"(?<=[:blank:]to)[:blank:]","\n")
      ) ->  go_enriched_all_fig2
    go_enriched_all_fig2 %>% 
      ggplot(
        aes(x=fct(name,levels=c("MOCK","MG132","MGTG")),y=fct(Description,levels = go_enriched_all_fig2 %>% pull(Description) %>% unique() %>% rev()),fill = -log(value,10))
      )+
      geom_tile()+
      scale_fill_distiller(
        palette = "YlOrRd",
        direction = 1,
        na.value = "grey80"
      )+
      labs(
        x="Experiment",
        y="GO description",
        fill = "Q value"
      )+
      theme_bw()+
      theme(
        axis.text.x = element_text(color = "black",size= 7,angle = 45, hjust=1),
        axis.text.y =element_text(color = "black",size= 6),
        axis.title = element_text(color = "black",size= 7),
        plot.title = element_text(color = "black",size= 7),
        #legend.position = "none",
        axis.ticks = element_line(color = "black"),
        legend.title = element_text(color = "black",size= 7),
        legend.text = element_text(color = "black",size= 7)
      )
    ggsave(
      "./figures/4_09_tile_GO_all_main_legend.png",height = 5, width =5,dpi = 600,units = "cm"
    )
      
    
    
  }
  
  ##10. focal adhesion
  {
    arginylome_validation_matrix3_sequence_distinct_trinity
    
    go_enriched_all %>%
      filter(
        Description == "focal adhesion"
      ) %>% 
      select(geneID) %>% 
      separate_longer_delim(
        cols=geneID,delim=regex("\\/")
      ) %>% 
      distinct(geneID) %>% 
      pull(geneID)
    
    go_enriched_all %>%
      filter(
        Description == "endoplasmic reticulum lumen"
      ) %>% 
      select(geneID) %>% 
      separate_longer_delim(
        cols=geneID,delim=regex("\\/")
      ) %>% 
      distinct(geneID) %>% 
      pull(geneID)
    
    go_enriched_all %>%
      filter(
        Description == "protein folding in endoplasmic reticulum"
      ) %>% 
      select(geneID) %>% 
      separate_longer_delim(
        cols=geneID,delim=regex("\\/")
      ) %>% 
      distinct(geneID) %>% 
      pull(geneID)
    
    go_enriched_all %>% 
      filter(
        str_detect(geneID,"SHC")
      )
    
    ggVennDiagram(
      list(
        A=go_enriched_all %>%
          filter(
            Description == "focal adhesion"
          ) %>% 
          select(geneID) %>% 
          separate_longer_delim(
            cols=geneID,delim=regex("\\/")
          ) %>% 
          distinct(geneID) %>% 
          pull(geneID),
        B=go_enriched_all %>%
          filter(
            Description == "endoplasmic reticulum lumen"
          ) %>% 
          select(geneID) %>% 
          separate_longer_delim(
            cols=geneID,delim=regex("\\/")
          ) %>% 
          distinct(geneID) %>% 
          pull(geneID)
      ),
      category.names = c("",""),
      set_size = 3,
      label_alpha = 0,
      label_size = 3,
      edge_size = 0.5
    )+
      theme_void() +
      #scale_fill_distiller(palette = "Reds",direction=1)+
      scale_fill_gradient(low="grey90",high = "red")+
      scale_color_manual(
        values = c("black","black","black")
      )+
      coord_sf(clip="off")+
      theme(
        legend.position = "none",
        plot.margin = margin(rep(0.5,4),unit="cm")
      )
    ggsave("./figures/4_10_focal_go.png",width = 5, height = 5, units = "cm",dpi=600)
    
    go_enriched_all %>%
      filter(
        Description == "focal adhesion"
      ) %>% 
      select(geneID) %>% 
      separate_longer_delim(
        cols=geneID,delim=regex("\\/")
      ) %>% 
      distinct(geneID) %>% 
      pull(geneID)%>% 
      writeClipboard()
    
    go_enriched_all %>%
      filter(
        Description == "endoplasmic reticulum lumen"
      ) %>% 
      select(geneID) %>% 
      separate_longer_delim(
        cols=geneID,delim=regex("\\/")
      ) %>% 
      distinct(geneID) %>% 
      pull(geneID) %>% 
      writeClipboard()
  }
  
  ###11. reactome ##entrez
  {
    arginylome_validation_matrix3_sequence_distinct_trinity %>% 
      filter(MOCK >0) %>% pull(gene) -> go_input
    as_tibble(AnnotationDbi::select(org.Hs.eg.db,keys=go_input, columns=c("SYMBOL","ENTREZID"),keytype="SYMBOL")) %>% 
      pull(ENTREZID) -> pa_input
    
    enrichPathway(
      pa_input,
      readable = TRUE,
      pvalueCutoff = 0.05
    ) %>% .@result %>%  arrange(qvalue) -> pa_enriched_mock
    
    arginylome_validation_matrix3_sequence_distinct_trinity %>% 
      filter(MG132 >0) %>% pull(gene) -> go_input
    as_tibble(AnnotationDbi::select(org.Hs.eg.db,keys=go_input, columns=c("SYMBOL","ENTREZID"),keytype="SYMBOL")) %>% 
      pull(ENTREZID) -> pa_input
    
    enrichPathway(
      pa_input,
      readable = TRUE,
      pvalueCutoff = 0.05
    ) %>% .@result %>%  arrange(qvalue) -> pa_enriched_mg
    
    arginylome_validation_matrix3_sequence_distinct_trinity %>% 
      filter(MGTG >0) %>% pull(gene) -> go_input
    as_tibble(AnnotationDbi::select(org.Hs.eg.db,keys=go_input, columns=c("SYMBOL","ENTREZID"),keytype="SYMBOL")) %>% 
      pull(ENTREZID) -> pa_input
    
    enrichPathway(
      pa_input,
      readable = TRUE,
      pvalueCutoff = 0.05
    ) %>% .@result %>%  arrange(qvalue) -> pa_enriched_mgtg
    
    pa_enriched_mock %>% 
      mutate(
        experiment = "MOCK"
      ) %>% 
      bind_rows(
        pa_enriched_mg %>% 
          mutate(
            experiment = "MG132"
          )
      )%>% 
      bind_rows(
        pa_enriched_mgtg %>% 
          mutate(
            experiment = "MGTG"
          )
      ) -> pa_enriched_all
    
    
    arginylome_validation_matrix3_sequence_distinct_trinity %>% 
      filter(MOCK >0) %>% pull(gene) -> go_input
    as_tibble(AnnotationDbi::select(org.Hs.eg.db,keys=go_input, columns=c("SYMBOL","ENTREZID"),keytype="SYMBOL")) %>% 
      pull(ENTREZID) -> pa_input
    enrichKEGG(
      pa_input,
      organism = 'hsa'
    )
    
    arginylome_validation_matrix3_sequence_distinct_trinity %>% 
      filter(MOCK >0) %>% pull(gene) -> go_input
    as_tibble(AnnotationDbi::select(org.Hs.eg.db,keys=go_input, columns=c("SYMBOL","ENTREZID"),keytype="SYMBOL")) %>% 
      pull(ENTREZID) -> kegg_input
    enrichKEGG(
      kegg_input,
      pvalueCutoff = 0.05
    ) %>% .@result %>%  arrange(qvalue) -> kegg_enriched_mock
    
    arginylome_validation_matrix3_sequence_distinct_trinity %>% 
      filter(MG132 >0) %>% pull(gene) -> go_input
    as_tibble(AnnotationDbi::select(org.Hs.eg.db,keys=go_input, columns=c("SYMBOL","ENTREZID"),keytype="SYMBOL")) %>% 
      pull(ENTREZID) -> kegg_input
    enrichKEGG(
      kegg_input,
      pvalueCutoff = 0.05
    ) %>% .@result %>%  arrange(qvalue) -> kegg_enriched_mg
    
    arginylome_validation_matrix3_sequence_distinct_trinity %>% 
      filter(MGTG >0) %>% pull(gene) -> go_input
    as_tibble(AnnotationDbi::select(org.Hs.eg.db,keys=go_input, columns=c("SYMBOL","ENTREZID"),keytype="SYMBOL")) %>% 
      pull(ENTREZID) -> kegg_input
    enrichKEGG(
      kegg_input,
      pvalueCutoff = 0.05
    ) %>% .@result %>%  arrange(qvalue) -> kegg_enriched_mgtg
    
    kegg_enriched_mock %>% 
      mutate(
        experiment = "MOCK"
      ) %>% 
      bind_rows(
        kegg_enriched_mg %>% 
          mutate(
            experiment = "MG132"
          )
      )%>% 
      bind_rows(
        kegg_enriched_mgtg %>% 
          mutate(
            experiment = "MGTG"
          )
      ) -> kegg_enriched_all
    
    arginylome_validation_matrix3_sequence_distinct_trinity %>% 
      filter(MOCK >0) %>% pull(gene) -> go_input
    as_tibble(AnnotationDbi::select(org.Hs.eg.db,keys=go_input, columns=c("SYMBOL","ENTREZID"),keytype="SYMBOL")) %>% 
      pull(ENTREZID) -> kegg_input
    enrichWP(
      kegg_input,
      organism = "Homo sapiens"
    ) %>% .@result %>%  arrange(qvalue) -> wp_enriched_mock
    
    arginylome_validation_matrix3_sequence_distinct_trinity %>% 
      filter(MOCK >0) %>% pull(gene) -> go_input
    as_tibble(AnnotationDbi::select(org.Hs.eg.db,keys=go_input, columns=c("SYMBOL","ENTREZID"),keytype="SYMBOL")) %>% 
      pull(ENTREZID) -> kegg_input
    enrichWP(
      kegg_input,
      organism = "Homo sapiens"
    ) %>% .@result %>%  arrange(qvalue) -> wp_enriched_mock
    
    kegg_enriched_mock %>% 
      mutate(
        experiment = "MOCK"
      ) %>% 
      bind_rows(
        kegg_enriched_mg %>% 
          mutate(
            experiment = "MG132"
          )
      )%>% 
      bind_rows(
        kegg_enriched_mgtg %>% 
          mutate(
            experiment = "MGTG"
          )
      ) -> kegg_enriched_all
    
    pa_enriched_all %>% 
      mutate(
        Description = str_replace(Description,"(?<=.{10}[:blank:]in)[:blank:]","\n"),
        Description = str_replace(Description,"(?<=.{15}[:blank:]by)[:blank:]","\n"),
        Description = str_replace(Description,"\\(ATF6\\-alpha\\)[:blank:]",""),
      ) %>% 
      select(ID,Description,experiment,pvalue) %>% 
      pivot_wider(
        id_cols = c(ID,Description),
        names_from = experiment,
        values_from = pvalue
      ) %>% 
      mutate(
        all = MOCK*MG132*MGTG
      ) %>% 
      arrange(all) %>% 
      slice_head(n=10) %>% pull(Description) -> pa_description
    
    pa_enriched_all %>% 
      mutate(
        Description = str_replace(Description,"(?<=.{10}[:blank:]in)[:blank:]","\n"),
        Description = str_replace(Description,"(?<=.{15}[:blank:]by)[:blank:]","\n"),
        Description = str_replace(Description,"\\(ATF6\\-alpha\\)[:blank:]",""),
      ) %>%
      select(ID,Description,experiment,pvalue) %>% 
      pivot_wider(
        id_cols = c(ID,Description),
        names_from = experiment,
        values_from = pvalue
      ) %>%
      filter(
        Description %in% pa_description
      ) %>%
      pivot_longer(
        MOCK:MGTG
      ) %>% 
      arrange(value) -> pa_enriched_all_fig
    
    pa_enriched_all_fig %>% 
      ggplot(
        aes(x=fct(name,levels=c("MOCK","MG132","MGTG")),y=fct(Description,levels = rev(pa_description)),fill = -log(value,10))
      )+
      geom_tile()+
      scale_fill_distiller(
        palette = "YlOrRd",
        direction = 1,
        na.value = "grey80"
      )+
      labs(
        x="Experiment",
        y="Reactome",
        fill = "Q value"
      )+
      theme_bw()+
      theme(
        axis.text.x = element_text(color = "black",size= 7,angle = 45, hjust=1),
        axis.text.y =element_text(color = "black",size= 6),
        axis.title = element_text(color = "black",size= 7),
        plot.title = element_text(color = "black",size= 7),
        legend.position = "none",
        axis.ticks = element_line(color = "black"),
        legend.title = element_text(color = "black",size= 7),
        legend.text = element_text(color = "black",size= 7)
      )
    ggsave(
      "./figures/4_11_tile_RA_all_main_legend.png",height = 5, width =6,dpi = 600,units = "cm"
    )
    
    pa_enriched_all %>% 
      mutate(
        Description = str_replace(Description,"(?<=.{10}[:blank:]in)[:blank:]","\n"),
        Description = str_replace(Description,"(?<=.{15}[:blank:]by)[:blank:]","\n"),
        Description = str_replace(Description,"\\(ATF6\\-alpha\\)[:blank:]",""),
      ) %>% 
      select(ID,Description,experiment,pvalue) %>% 
      pivot_wider(
        id_cols = c(ID,Description),
        names_from = experiment,
        values_from = pvalue
      ) %>% 
      mutate(
        all = MOCK*MG132*MGTG
      ) %>% 
      arrange(all) %>% 
      write_tsv("./export/pa_enriched_all.tsv")
  }
  RColorBrewer::display.brewer.all()
  vignette("fully-customed", package = "ggVennDiagram")
  
  
  
  ###13. GOCC
  {
    #library(UniProtKeywords)
    arginylome_validation_matrix3_sequence_distinct_trinity
    
    keyref <- keys(org.Hs.eg.db, keytype = "SYMBOL")
    AnnotationDbi::select(
      org.Hs.eg.db, 
      keys=keyref,   ###gene name vector
      columns=c("SYMBOL","ENTREZID","GENENAME"),
      keytype="SYMBOL"
    ) -> table_genenames
    
    arginylome_validation_matrix3_sequence_distinct_trinity %>% 
      pull(gene)
    UniProtKeywords::load_keyword_genesets(
      as_table = TRUE
    ) %>% tibble() %>% 
      left_join(
        table_genenames,
        join_by(gene == ENTREZID)
      ) -> table_uniprot_keywords
    
    load_keyword_genesets(
      as_table = TRUE,
      category = "Cellular component"
    ) -> table_uniprot_keywords_cc
    
    table_uniprot_keywords_cc %>% 
      tibble() %>% 
      left_join(
        table_genenames,
        join_by(gene == ENTREZID)
      ) -> table_uniprot_keywords_cc
    
    table_uniprot_keywords %>% write_tsv("table_uniprot_keywords.tsv")
    table_uniprot_keywords_cc %>% write_tsv("table_uniprot_keywords_cc.tsv")
    
    table_uniprot_keywords_cc %>% 
      distinct(gene) #8298
    arginylome_validation_matrix3_sequence_distinct_trinity %>% 
      distinct(gene) #119
    
    table_uniprot_keywords_cc %>%
      filter(
        SYMBOL %in% dbtable$gene
      ) %>% 
      summarize(
        .by = keyword,
        count =n()
      ) %>% 
      left_join(
        table_uniprot_keywords_cc %>%
          filter(
            SYMBOL %in% arginylome_validation_matrix3_sequence_distinct_trinity$gene
          ) %>% 
          summarize(
            .by = keyword,
            count_arg =n()
          )
      ) -> table_uniprot_keywords_cc_arg
    
    table_uniprot_keywords_cc_arg %>% 
      drop_na() %>% 
      bind_cols(
        genes_all = 8298L,
        genes_interest = 119L
      ) %>% 
      mutate(
        matrix = pmap(
          list(count, count_arg, genes_all, genes_interest),
          \(a,b,c,d) matrix(
            c(b,d-b,a-b,c-a-d+b),
            nrow=2,ncol=2
          )
        ),
        fisher_p= map_dbl(
          matrix,
          \(x) fisher.test(x)$p.value
        ),
        fisher_e= map_dbl(
          matrix,
          \(x) fisher.test(x)$estimate
        )
      ) %>% 
      arrange(fisher_p)
    
    arginylome_validation_matrix3_sequence_distinct_trinity %>% 
      filter(
        MOCK >0
      ) %>% 
      distinct(gene) %>% 
      pull(gene) -> genes_mock
    
    arginylome_validation_matrix3_sequence_distinct_trinity %>% 
      filter(
        MG132 >0
      ) %>% 
      distinct(gene) %>% 
      pull(gene) -> genes_mg132
    
    arginylome_validation_matrix3_sequence_distinct_trinity %>% 
      filter(
        MGTG >0
      ) %>% 
      distinct(gene) %>% 
      pull(gene) -> genes_mgtg
    
    
    clusterProfiler::enrichGO(
      genes_mock,
      org.Hs.eg.db,
      keyType = "SYMBOL",
      ont = "CC",
      readable = TRUE
    ) %>% .@result %>% 
      mutate(
        experiment = "MOCK"
      ) %>% 
      bind_rows(
        clusterProfiler::enrichGO(
          genes_mg132,
          org.Hs.eg.db,
          keyType = "SYMBOL",
          ont = "CC",
          readable = TRUE
        ) %>% .@result %>% 
          mutate(
            experiment = "MG132"
          ),
        clusterProfiler::enrichGO(
          genes_mgtg,
          org.Hs.eg.db,
          keyType = "SYMBOL",
          ont = "CC",
          readable = TRUE
        ) %>% .@result %>% 
          mutate(
            experiment = "MGTG"
          )
      ) -> table_go_cc
    
    calculateSimMatrix(
      table_go_cc$ID,
      orgdb="org.Hs.eg.db",
      ont="CC",
      method="Rel"
    ) -> simMatrix
    
    scores <- setNames(-log10(table_go_cc$pvalue), table_go_cc$ID)
    reduceSimMatrix(
      simMatrix,
      scores,
      threshold=0.7,
      orgdb="org.Hs.eg.db",
    ) -> reducedTerms
    
    table_go_cc %>% 
      left_join(
        reducedTerms,
        join_by(
          ID == go
        )
      ) -> table_go_cc
    
    table_go_cc %>% 
      filter(
        termDispensability == 0 
      ) -> table_go_cc_reduced
    
    table_go_cc_reduced %>% 
      filter(
        qvalue <0.05
      )-> table_go_cc_reduced_fdr0.05
    
    table_go_cc_reduced_fdr0.05 %>% 
      filter(Count >2) %>% 
      ggplot(
        aes(
          x=fct(experiment,levels = c("MOCK","MG132","MGTG")),
          y=term,
          fill = -log(qvalue,10)
        )
      )+
      geom_tile()+
      scale_fill_distiller(
        palette = "YlOrRd",
        direction = 1,
        na.value = "grey80"
      )+
      labs(
        x="Experiment",
        y="GOCC Terms",
        fill = expression(-Log[10](FDR))
      )+
      theme_classic2()+
      theme(
        axis.text.x = element_text(color = "black",size= 8,angle = 45,hjust = 1),
        axis.text.y =element_text(color = "black",size= 7),
        axis.title = element_text(color = "black",size= 8,face="bold"),
        plot.title = element_text(color = "black",size= 8),
        #legend.position = "bottom",
        axis.ticks = element_line(color = "black"),
        legend.title = element_text(color = "black",size=8),
        legend.text = element_text(color = "black",size= 8)
      )
    ggsave(
      "./figures/4_13_tile_gocc.png",height = 6, width =15,dpi = 600,units = "cm"
    )
    
    table_go_cc_reduced_fdr0.05 %>% 
      ggplot(
        aes(
          x=Count,
          y=-log(qvalue,10)
        )
      )+
      geom_point()+
      geom_text_repel(
        mapping= aes(label = term),
        nudge_x = 20 -table_go_cc_reduced_fdr0.05$Count,
        force = 1,
        force_pull = 0,
        max.overlaps = 10,
        box.padding   = 0.17,
        point.padding = 0,
        size = 3,
        direction = "y",
        segment.size  = 0.1,
        xlim = c(20),
        ylim = c(0,10),
        color = "black",
        segment.alpha = 0.2
      )+
      labs(
        y= expression(-Log[10](FDR)),
        x= "Gene Count"
      )+
      lims(
        y=c(0,10)
      )+
      theme_classic2()+
      theme(
        axis.text.x = element_text(color = "black",size= 10),
        axis.text.y =element_text(color = "black",size= 10),
        axis.title = element_text(color = "black",size= 10),
        plot.title = element_text(color = "black",size= 10),
        legend.position = "none",
        axis.ticks = element_line(color = "black"),
        legend.title = element_text(color = "black",size=10),
        legend.text = element_text(color = "black",size= 10)
      )
    ggsave(
      "./figures/4_13_gocc.png",height = 18, width =18,dpi = 600,units = "cm"
    )
    
    table_go_cc_reduced_fdr0.05 %>% 
      select(
        term,Count,qvalue
      ) %>% 
      ggplot(
        aes(y=term,
            x=Count,
            fill = Count)
      )+
      geom_tile()
    
  }
  
  #@table_go_cc_reduced
  #@table_go_cc
}


###
####Chapter 5. Rcatcher validation

###
####Chapter 6. Recalling arginylome using PRM


###
#others()
{
  ###GENE_POSITION!
  arginylome_validation_matrix3_sequence_distinct_trinity %>% 
    mutate(
      gene_position = paste0(
        str_split_i(prot_pos_reposition,"\\|",1),
        "|",
        str_split_i(prot_pos_reposition,"\\|",4)
      )
    ) -> arginylome_validation_matrix3_sequence_distinct_trinity
  
  summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered_trinity %>% 
    mutate(
      gene_position = paste0(
        str_split_i(prot_pos_reposition,"\\|",1),
        "|",
        str_split_i(prot_pos_reposition,"\\|",4)
      )
    )->summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered_trinity
  
  
  arginylome_validation_matrix3_sequence_distinct_trinity %>% 
    filter(trinity == "Passed")
  
  arginylome_validation_matrix3_sequence_distinct_trinity %>% 
    mutate(
      gene_position = paste0(gene,"|",str_extract(prot_pos_reposition,"(?<=\\|)[0-9]+"))
    ) -> arginylome_validation_matrix3_sequence_distinct_trinity
  
  arginylome_validation_matrix3_sequence_distinct_trinity %>% 
    left_join(
      summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered_trinity %>% 
        select(prot_pos_master,p5p1),
      join_by(prot_pos_master==prot_pos_master),
      multiple = "first"
    ) -> arginylome_validation_matrix3_sequence_distinct_trinity
  
  arginylome_validation_matrix3_sequence_distinct_trinity %>% 
    mutate(
      p5p5prime = paste0(p5p1,str_sub(Sequence,1,5))
    ) -> arginylome_validation_matrix3_sequence_distinct_trinity
  
  summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered %>% 
    mutate(
      p5p5prime = paste0(p5p1,str_sub(Sequence,1,5))
    ) -> summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered
  
  arginylome_validation_matrix3_sequence_distinct_trinity %>% 
    mutate(
      protease_site = case_when(
        cleavage_type == "C14.001" ~ "Caspase-1",
        cleavage_type == "C14.005" ~ "Caspase-6",
        cleavage_type == "C14.003" ~ "Caspase-3",
        cleavage_type == "M12.033" ~ "LAST_MAM peptidase",
        cleavage_type == "SP" ~ "Signal Peptide",
        cleavage_type == "mTP" ~ "Transit Peptide",
        cleavage_type == "M12.004" ~ "meprin beta subunit",
        cleavage_type == "A01.009" ~ "cathepsin D",
        cleavage_type == "M12.002" ~ "matrix metallopeptidase-8",
        cleavage_type == "S01.010" ~ "granzyme B",
        .default = cleavage_type
      )
    ) ->arginylome_validation_matrix3_sequence_distinct_trinity
  
  psm_563516 %>% 
    left_join(
      file_raw_index2 %>% 
        select(-enzyme),
      join_by(file_id_series==set)
    ) -> psm_563516
  
  psm_563516 %>% 
    mutate(
      rt_norm = (`RT in min`-left_time)/(right_time-left_time)
    ) -> psm_563516
  
  psm_563516 %>% 
    mutate(
      artefact_arg = if_else(str_detect(strip_seq,"^R") & str_detect(nmod,"Acetyl\\:2H\\(3"),TRUE,FALSE)
    ) -> psm_563516
  
  psm_563516 %>%
    separate_rows(mod_sites,sep = "\\;") %>% 
    mutate(
      mod_sites = as.integer(mod_sites)
    ) %>% 
    mutate(
      mod_sites = if_else(mod_sites == 0, mod_sites,mod_sites-1)
    ) %>% 
    nest(
      mod_sites = mod_sites
    ) %>% 
    mutate(
      mod_sites = map_chr(mod_sites,\(x) pull(x,mod_sites) %>% paste0(collapse = ";"))
    ) -> psm_563516
  psm_563516 %>% 
    distinct(peak_file,.keep_all = TRUE) -> psm_563516
  
}


fdr_pcc_negative_set10 %>% 
  summarize(
    .by = moddb_result,
    count = n()
  )

fdr_pcc_negative_set10 %>% 
  filter(
    Confidence == "High"
  ) %>% 
  filter(
    moddb_result == "false"
  )

##################
##export
{
  saveRDS(
    summaryfile_psm_for_testing3_pcc_rt,
    file = "./export/summaryfile_psm_for_testing3_pcc_rt.rds"
  )
  saveRDS(
    file_index_arg6_ttest,
    file = "./export/file_index_arg6_ttest.rds"
  )
  saveRDS(
    summaryfile_psm_for_testing3_pcc_rt_hfsm,
    file = "./export/summaryfile_psm_for_testing3_pcc_rt_hfsm.rds"
  )
  saveRDS(
    arginylome_validation_matrix2,
    file = "./export/arginylome_validation_matrix2.rds"
  )
  saveRDS(
    alphafold_fetched,
    file = "./export/alphafold_fetched.rds"
  )
  saveRDS(
    summaryfile_psm2,
    file = "./export/summaryfile_psm2.rds"
  )
  saveRDS(
    summaryfile_psm_for_training2,
    file = "./export/summaryfile_psm_for_training2.rds"
  )
  saveRDS(
    psm_1015276,
    file = "./export/psm_1015276.rds"
  )
  saveRDS(
    psm_563516,
    file = "./export/psm_563516.rds"
  )
  saveRDS(
    psm_251969_chymo,
    file = "./export/psm_251969_chymo.rds"
  )
  write_rds(
    psm_311547_trypsin,
    file = "./export/psm_311547_trypsin.rds"
  )
  
  
  saveRDS(
    psm_17669,
    file = "./export/psm_17669.rds"
  )
  saveRDS(
    psm_536653,
    file = "./export/psm_536653.rds"
  )
  write_rds(
    summaryfile_psm_for_testing3_pcc_rt_hfsm2,
    file = "./export/summaryfile_psm_for_testing3_pcc_rt_hfsm2.rds"
  )
  write_rds(
    summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered,
    file = "./export/summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered.rds"
  )
  write_rds(
    summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered_trinity,
    file = "./export/summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered_trinity.rds"
  )
  write_rds(
    arginylome_validation_matrix3,
    file = "./export/arginylome_validation_matrix3.rds"
  )
  write_rds(
    arginylome_validation_matrix3_sequence_distinct,
    file = "./export/arginylome_validation_matrix3_sequence_distinct.rds"
  )
  
  write_rds(
    arginylome_validation_matrix3_sequence_distinct_trinity,
    file = "./export/arginylome_validation_matrix3_sequence_distinct_trinity.rds"
  )
  
  write_tsv(
    arginylome_validation_matrix3_sequence_distinct_trinity,
    "./export/arginylome_validation_matrix3_sequence_distinct_trinity.tsv"
  )
  summaryfile_psm_for_testing3_pcc_rt_hfsm2 %>% 
    write_tsv(
      "./export/summaryfile_psm_for_testing3_pcc_rt_hfsm2.tsv"
    )
  
  saveRDS(
    summaryfile2,
    file = "./export/summaryfile2.rds"
  )
  saveRDS(
    fdr_pcc_negative_set9,
    file = "./export/fdr_pcc_negative_set9.rds"
  )##11407
  saveRDS(
    fdr_pcc_negative_set4,
    file = "./export/fdr_pcc_negative_set4.rds"
  ) ##11667
  
  saveRDS(
    prm_training_psm,
    file = "./export/prm_training_psm.rds"
  ) ##11667
  
  remove(fdr_pcc_negative_set2_dbcheck)
  remove(prm_training_psm)
  
  read_tsv(
    "./export/arginylome_validation_matrix3_sequence_distinct_240110.tsv"
  ) -> arginylome_validation_matrix3_sequence_distinct
  
  read_tsv(
    "./export/arginylome_validation_matrix3_sequence_distinct_trinity_134_240110.tsv"
  ) -> arginylome_validation_matrix3_sequence_distinct_trinity
}

remove(summaryfile2)
remove(summaryfile_psm2)
remove(summaryfile_psm_for_training2)
remove(arginylome_validation_matrix3)
remove(alphafold_fetched)
remove(psm_563516)
remove(psm_1015276)
remove(psm_17669)
remove(psm_17669_chymo)
remove(psm_536653)
remove(psm_251969_chymo)
remove(psm_311547_trypsin)
remove(psm_299486_trypsin)
remove(psm_237167_chymotrypsin)


save.image("240130_new.RData")
save.image()

col = list(Loc = c(
  "Extracellular"="#3B9AB2",
  "Cell membrane"="#63ADBE",
  "Lysosome/Vacuole"="#9EBE91",
  "Golgi apparatus"="#EBCC2A",
  "Endoplasmic reticulum"="#E4B80E",
  "Cytoplasm"="#D5D5D3",
  "Nucleus"="#F21A00",
  "Mitochondrion"="#4DAF4A"
))


####outside
{
  read_tsv(
    "C:/Users/admin/alphapeptdeep/project/folder/alphapeptdeep/inrich/output/output_metrics_describ_for_arg_trypsin_231201.tsv"
  ) %>% 
    bind_rows(
      read_tsv(
        "C:/Users/admin/alphapeptdeep/project/folder/alphapeptdeep/inrich/output/output_metrics_describ_for_arg_chymotrypsin_231201.tsv"
      )
    ) -> temp_pcc
  
  temp_pcc %>% 
    select(PCC) %>% 
    rename(PCC2=PCC)
  
  summaryfile_psm_for_testing3_pcc %>% 
    bind_cols(
      temp_pcc %>% 
        select(PCC) %>% 
        rename(PCC2=PCC)
    ) -> temp_pcc2
  
  temp_pcc2 %>% 
    mutate(
      diff_pcc = PCC2-PCC
    ) %>% 
    ggplot(
      aes(x=diff_pcc)
    )+
    geom_histogram()
  
  
  
  temp_pcc2 %>% 
    filter(
      PCC2 >= cutoff_pcc
    )
  
  dplyr::setdiff(
    temp_pcc2 %>% 
      filter(
        PCC2 >= cutoff_pcc
      ),temp_pcc2 %>% 
      filter(
        PCC >= cutoff_pcc
      )
  ) %>% view()
    
  temp_pcc2 %>% 
    filter(
      PCC >= cutoff_pcc
    )
}


####read
{
  read_rds("./export/summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered.rds") %>% 
    filter(
      dl_score == 3
    )
  read_rds("./export/arginylome_validation_matrix3_sequence_distinct_trinity.rds")

  arginylome_validation_matrix3_sequence_distinct_trinity %>% 
    select(
      gene,loc
    ) %>% 
    write_tsv("./export/gene_loc.tsv")
}

##Export
{
  summaryfile_psm %>% write_tsv("./export/summaryfile_psm.tsv")
  arginylome_validation_matrix3_sequence_distinct_trinity %>% write_tsv("./export/arginylome_validation_matrix3_sequence_distinct_trinity.tsv")
  peptide_lfq_trinity_distinct_arg %>% write_tsv("./export/peptide_lfq_trinity_distinct_arg.tsv")
}

##UniprotR signal/transit/propep
{
  library(UniprotR)
  
  arginylome_validation_matrix3_sequence_distinct_trinity %>%
    distinct(prot_pos_master) %>% 
    mutate(
      accession= str_split_i(prot_pos_master,"\\|",1)
    ) %>% 
    pull(accession)
  
  dbtable %>% pull(
    accession
  )
  
  GetPTM_Processing(
    arginylome_validation_matrix3_sequence_distinct_trinity %>%
      distinct(prot_pos_master) %>% 
      mutate(
        accession= str_split_i(prot_pos_master,"\\|",1)
      ) %>% 
      pull(accession)
  ) -> test
  test %>% 
    rownames_to_column() -> test2
  
  arginylome_validation_matrix3_sequence_distinct_trinity %>% 
    mutate(
      accession= str_split_i(prot_pos_master,"\\|",1),
      pos= str_split_i(prot_pos_master,"\\|",2)
    ) %>% 
    left_join(
      test2,
      join_by(
        accession == rowname
      )
    ) -> arginylome_validation_matrix3_sequence_distinct_trinity2
  
  arginylome_validation_matrix3_sequence_distinct_trinity2 %>%
    mutate(
      uniprot_signal = as.integer(str_extract(Signal.peptide,"(?<=\\.\\.)[0-9]+")),
      uniprot_transit = as.integer(str_extract(Transit.peptide,"(?<=\\.\\.)[0-9]+")),
      uniprot_propep = as.integer(str_extract(Propeptide,"(?<=\\.\\.)[0-9]+")),
    )->arginylome_validation_matrix3_sequence_distinct_trinity2
  
  arginylome_validation_matrix3_sequence_distinct_trinity2 %>% 
    write_tsv(
      "./export/arginylome_validation_matrix3_sequence_distinct_trinity2.tsv"
    )
}

###
{
  read
  tibble(
    origin = read_file(
      "C:/Users/admin/Desktop/Data/DB/UniprotRefHuman/202302/UP000005640_9606.dat"
    )
  ) -> test
  
    str_extract(
      "^.."
    ) %>% unique()
}

###LFQ DATA ALL
#@summaryfile
{
  summaryfile %>% 
    select(Accessions.Single,Position.start)
  
  summaryfile %>% # peptide file
    rename(
      "MG132_1"="Abundances Normalized F12 Sample MG132 1ST trypsin",
      "MG132_2"="Abundances Normalized F6 Sample MG132 2ND trypsin",
      "MGTG_1"="Abundances Normalized F1 Sample MGTG 1ST trypsin",
      "MGTG_2"="Abundances Normalized F7 Sample MGTG 2ND trypsin",
      "MOCK_1"="Abundances Normalized F2 Sample MOCK 1ST trypsin",
      "MOCK_2"="Abundances Normalized F5 Sample MOCK 2ND trypsin"
    ) %>% 
    select(
      Sequence,nmod,gene,Accessions.Single,Position.start,Modifications,MG132_1,MG132_2,MGTG_1,MGTG_2,MOCK_1,MOCK_2
    ) %>% 
    bind_rows(
      summaryfile %>% # peptide file
        filter(
          Sequence %in% stripseq_arg_trinity_chymotrypsin
        ) %>%
        rename(
          "MG132_1"="Abundances Normalized F11 Sample MG132 1ST chymotrypsin",
          "MG132_2"="Abundances Normalized F9 Sample MG132 2ND chymotrypsin",
          "MGTG_1"="Abundances Normalized F3 Sample MGTG 1ST chymotrypsin",
          "MGTG_2"="Abundances Normalized F10 Sample MGTG 2ND chymotrypsin",
          "MOCK_1"="Abundances Normalized F4 Sample MOCK 1ST chymotrypsin",
          "MOCK_2"="Abundances Normalized F8 Sample MOCK 2ND chymotrypsin"
        ) %>% 
        select(
          Sequence,nmod,gene,Modifications,Accessions.Single,Position.start,MG132_1,MG132_2,MGTG_1,MGTG_2,MOCK_1,MOCK_2
        )
    ) -> peptide_lfq_all
  
  
  
  ##GENE position
  peptide_lfq_all %>% 
    filter(nmod %in% c("Acetyl","D3Acetyl","D3AcetylArg","D3AcetylArgDeamid")) %>%
    mutate(
      nmod = case_when(nmod == "D3AcetylArgDeamid"~"Arg",
                       nmod == "D3Acetyl"~"Free",
                       nmod == "D3AcetylArg"~"Arg",
                       nmod == "D3AcetylArg"~"Acetyl",
                       .default = nmod)
    ) %>% 
    pivot_wider(
      id_cols = c(Sequence,gene,Accessions.Single,Position.start),
      names_from = nmod,
      values_from = MG132_1:MOCK_2,
      values_fn = mean
    ) -> peptide_lfq_all
  
  peptide_lfq_all %>% mutate(
    nsite = paste0(gene,"_",Position.start)
  ) -> peptide_lfq_all
  
  peptide_lfq_all %>% 
    write_tsv(
      "./export/peptide_lfq_all.tsv"
    )
}
#@peptide_lfq_all

save.image()
save.image("240730_backup.RData")
