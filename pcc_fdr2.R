{
  library(tidyverse)
  library(pROC)
  library(tidymodels)
  tidymodels_prefer()
  library(rmarkdown)
}

##test set 
###index file 230623 Gradient time 110min, min value = 48
##from index file to test data
##nmod , increase +1 on mod_sites
{
  folder_pdspectra <- "./annotation_psm/200821_HELA_falseset_final-(2)/"
  
  #backup
  {
    
    write_rds(file_index_falseset,"./export/file_index_falseset.rds")
    write_rds(file_index_falseset3,"./export/file_index_falseset3.rds")
  }
  {
    ##read xml
    xml2::read_html(
      paste0(folder_pdspectra,"index.html")
    ) -> file_index_falseset
    #from xml object read table
    rvest::html_table(
      file_index_falseset,
      header = TRUE
    ) ->file_index_falseset
    
    ##table trimming
    file_index_falseset[[1]] -> file_index_falseset #unlist
    file_index_falseset[seq(1,nrow(file_index_falseset),by=2),seq(2,ncol(file_index_falseset),by=2)] -> file_index_falseset ## delete unnessary cols and rows
    colnames(file_index_falseset)<- as.vector(file_index_falseset[1,]) %>% unlist ## colnames
    file_index_falseset[c(-1),] -> file_index_falseset
    
    #peak list handling
    readLines(paste0(folder_pdspectra,"index.html"))-> file_index_falseset_lines
    
    ##add column of peak list
    file_index_falseset %>% 
      bind_cols(tibble(peak_file = str_extract(file_index_falseset_lines[str_detect(file_index_falseset_lines,"Peak List")],"data\\/.+\\.txt"))) -> file_index_falseset
    
    file_index_falseset %>% 
      mutate(
        peak_table = map(
          peak_file,\(x) read_tsv(paste0(folder_pdspectra,x),show_col_types = FALSE), .progress = TRUE
        )
      ) %>% 
      mutate(
        peak_table_dropna = map(
          peak_table,\(x) drop_na(x)
        )
      ) -> file_index_falseset
    #### peak_table_dropna_relative intensity!
    file_index_falseset %>% 
      mutate(
        peak_table_dropna_wo_modloss = map(peak_table_dropna,\(x) dplyr::filter(x,!str_detect(matches,"NH3|Ox|H2O")))
      ) %>% 
      mutate(
        peak_table_dropna_wo_modloss = map(peak_table_dropna_wo_modloss, \(x) mutate(x,i_relative=i/max(i)))
      )%>% 
      mutate(
        peak_table_alpha_ms2= map(peak_table_dropna_wo_modloss,\(x) mutate(x, series = str_extract(matches,"b|y")))
      ) %>% 
      mutate(
        peak_table_alpha_ms2= map(peak_table_alpha_ms2,\(x) drop_na(x))
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
        nmod = str_extract(Modifications,"(?<=N\\-Term\\()[^\\)]+")
      ) %>%
      # mutate(
      #   Modifications = str_replace(Modifications,"N\\-Term\\([^\\)]+[\\)\\;]+","") ##remove Nmod
      # ) %>% 
      mutate(
        strip_seq=str_to_upper(str_extract(Sequence,"(?<=\\.)[:alpha:]+(?=\\.)"))
      ) %>% 
      mutate( ###add R in front of strip_seq
        strip_seq = paste0("R",strip_seq)
      ) %>%
      mutate(
        peak_table_alpha_ms2 = map(peak_table_alpha_ms2, \(x) mutate(x,i_relative=i/max(i)))
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
          \(x) mutate(x,ion_num2 = if_else(series == "y",seq_len - ion_num,ion_num+1)) ##bion +1, yion stayed
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
        fragment_inten_df2 = map(seq_len,\(x) tibble(ion_num2=seq(1,x-1),b_z1 = NA_integer_, b_z2 = NA_integer_, y_z1 = NA_integer_,y_z2 = NA_integer_))
      ) %>% 
      mutate(
        fragment_inten_df3 = map2(
          fragment_inten_df,
          fragment_inten_df2,
          \(x,y) left_join(
            y,x,
            join_by(ion_num2==ion_num2),
            keep=FALSE
          )
        )
      ) %>% 
      mutate(
        fragment_inten_df3 = map(
          fragment_inten_df3,
          \(x) select(x,!ends_with("x"))
        )
      ) %>% 
      mutate(
        fragment_inten_df3 = map(
          fragment_inten_df3,
          \(x) rename_with(x,\(y) gsub("\\.y$","",y))
        )
      ) %>% 
      mutate(
        fragment_inten_df3 = map(
          fragment_inten_df3,
          \(x) select(x,c(ion_num2,b_z1,b_z2,y_z1,y_z2))
        )
      )-> file_index_falseset
    
    colnames(file_index_falseset)[16] <-"scan_num"
    
    
    file_index_falseset %>% 
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
        Modifications = str_replace_all(Modifications,"[:blank:]","")
      ) %>% 
      distinct(
        strip_seq,scan_num,charge,XCorr,
        .keep_all = TRUE
      ) -> file_index_falseset
    
    colnames(file_index_falseset)[15] <- "RT"
    file_index_falseset %>% 
      mutate(
        rt = as.numeric(RT) %>% round(digits = 3)
      ) %>% 
      mutate(
        fragment_inten_df3 = map(
          fragment_inten_df3,
          \(x) unnest(x,keep_empty = TRUE),
          .progress = TRUE
        )
      ) %>% 
      mutate(
        fragment_inten_df3 = map(
          fragment_inten_df3,
          \(x) distinct(x,ion_num2,.keep_all = TRUE)
        )
      ) -> file_index_falseset
    
    ###modification check
    file_index_falseset %>% 
      separate_rows(
        Modifications,sep = ";"
      ) %>% 
      distinct(Modifications) %>% print(n=Inf)
    
    ##index creation
    file_index_falseset %>% 
      mutate(
        index = seq(1:nrow(.))
      ) %>% 
      separate_longer_delim(Modifications, delim=";") %>% 
      mutate(
        mods_aa = case_when(
          str_detect(Modifications,"^[A-Z]")~str_extract(Modifications,"^[A-Z]"),
          str_detect(Modifications,"N\\-Term")~"Any N-term",
          .default = Modifications
        ),
        nAA = str_length(strip_seq)
      ) %>% 
      mutate( ###Check whether these mods 
        mods_name = case_when(
          str_detect(Modifications,"Acetyl\\)")~"Acetyl",
          str_detect(Modifications,"Acetyl\\:2H\\(3\\)")~"Acetyl:2H(3)",
          str_detect(Modifications,"Oxidation")~"Oxidation",
          str_detect(Modifications,"Carbamidomethyl")~"Carbamidomethyl",
          str_detect(Modifications,"Ethanolamine\\)")~"Ethanolamine",
          str_detect(Modifications,"Ethanolamine_siderxn")~"EDC",
          str_detect(Modifications,"ArgNQ$")~"",
          str_detect(Modifications,"Arg$")~"",
          str_detect(Modifications,"ArgAcetylD3")~"Acetyl:2H(3)",
          str_detect(Modifications,"ArgAcetylD3NQ")~"Acetyl:2H(3)",
          .default = Modifications
        )
      ) %>% 
      mutate(
        mod_sites = case_when(
          str_detect(Modifications,"^N\\-Term")~"0",
          str_detect(Modifications,"^C\\-Term")~as.character(nAA),
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
        nestedcol = map(nestedcol,\(x) mutate(x,mod_sites=as.integer(mod_sites)) %>% mutate(mod_sites=case_when(mod_sites == 0 ~ 0, mod_sites >0 ~ mod_sites+1)))
      )%>%
      mutate(
        mod_sites = map_chr(nestedcol,\(x) x$mod_sites %>% paste(collapse = ";"))
      )%>% 
      mutate(
        nce = as.integer(27),
        instrument = "QE",
        charge = as.integer(Charge)
      ) %>% 
      mutate(
        mods=if_else(mods == "@",NA,mods)
      ) %>% 
      mutate(
        mod_sites=if_else(mod_sites == "",NA,mod_sites)
      ) %>% 
      distinct(strip_seq,scan_num,charge,peak_file,.keep_all =TRUE) %>%
      mutate(
        rt_norm = (rt-48)/110 ###min 24.895 , >> 48.427
      ) -> file_index_falseset
    
    file_index_falseset %>% 
      mutate(
        mods = str_replace(mods,"\\@N","@Any N-term")
      ) -> file_index_falseset
    
    
    # file_index_falseset %>%
    #   filter(
    #     !(strip_seq %in% summaryfile_psm_for_training_pcc$sequence) ##remove sequences in tp
    #   ) %>%
    #   filter(
    #     !(strip_seq %in% fdr_pcc_positive_set$sequence) ##remove sequences in tp
    #   )-> file_index_falseset2
    
    file_index_falseset -> file_index_falseset2
    
    fdr_pcc_positive_set$sequence
    
    file_index_falseset2 %>% 
      mutate(
        enzyme = if_else(str_detect(peak_file,"\\_\\-1196"),"trypsin","chymotrypsin")
      ) -> file_index_falseset2
    
    file_index_falseset2 %>% 
      unnest(fragment_inten_df3) -> file_index_falseset3 ##this is frag_df!
    
    file_index_falseset3 %>% 
      filter(
        enzyme == "trypsin"
      ) %>% 
      mutate(
        across(
          .cols = c(b_z1,b_z2,y_z1,y_z2),
          .fns = \(x) if_else(is.na(x),0,x)
        )
      ) %>%
      arrange(
        index,ion_num2
      ) %>% 
      mutate(
        frag_index = seq(0,nrow(.)-1)
      ) %>% 
      select(-c(ion_num2)) %>% 
      nest(frag_df= c(frag_index,b_z1,b_z2,y_z1,y_z2)) %>% 
      mutate(
        frag_start_idx = map_dbl(frag_df,\(x) min(x$frag_index))
      ) %>% 
      mutate(
        frag_stop_idx = map_dbl(frag_df,\(x) max(x$frag_index)+1)
      ) %>% 
      select(
        strip_seq,charge,index,rt,mods,mod_sites,nAA,frag_start_idx,frag_stop_idx,nce,instrument,rt_norm,frag_df
      ) %>% 
      rename(
        sequence = strip_seq
      )-> file_index_falseset3_trypsin
    
    file_index_falseset3_trypsin %>%
      select(-frag_df) %>%
      write_tsv(
        paste0(folder_alphapeptdeep,"false_set_test_trypsin_231211.tsv")
      )
    
    file_index_falseset3_trypsin %>%
      select(frag_df) %>%
      unnest(frag_df) %>%
      select(-frag_index) %>%
      write_tsv(
        paste0(folder_alphapeptdeep,"false_set_test_fragdf_trypsin_231211.tsv")
      )
    
    file_index_falseset3 %>% 
      filter(
        enzyme == "chymotrypsin"
      ) %>% 
      mutate(
        across(
          .cols = c(b_z1,b_z2,y_z1,y_z2),
          .fns = \(x) if_else(is.na(x),0,x)
        )
      ) %>%
      arrange(
        index,ion_num2
      ) %>% 
      mutate(
        frag_index = seq(0,nrow(.)-1)
      ) %>% 
      select(-c(ion_num2)) %>% 
      nest(frag_df= c(frag_index,b_z1,b_z2,y_z1,y_z2)) %>% 
      mutate(
        frag_start_idx = map_dbl(frag_df,\(x) min(x$frag_index))
      ) %>% 
      mutate(
        frag_stop_idx = map_dbl(frag_df,\(x) max(x$frag_index)+1)
      ) %>% 
      select(
        strip_seq,charge,index,rt,mods,mod_sites,nAA,frag_start_idx,frag_stop_idx,nce,instrument,rt_norm,frag_df
      ) %>% 
      rename(
        sequence = strip_seq
      )-> file_index_falseset3_chymotrypsin
    
    file_index_falseset3_chymotrypsin %>%
      select(-frag_df) %>%
      write_tsv(
        paste0(folder_alphapeptdeep,"false_set_test_chymotrypsin_231211.tsv")
      )
    
    file_index_falseset3_chymotrypsin %>%
      select(frag_df) %>%
      unnest(frag_df) %>%
      select(-frag_index) %>%
      write_tsv(
        paste0(folder_alphapeptdeep,"false_set_test_fragdf_chymotrypsin_231211.tsv")
      )
  }
  
  ### to find true positives! (RR) ##file_index_falseset_fdr
  {
    file_index_falseset2 %>% 
      mutate(
        strip_seq_dbmatch = str_sub(strip_seq,start=2,end=-1)
      ) -> file_index_falseset_fdr
    
    colnames(file_index_falseset_fdr)[7] <- "percolator_score"
    
  }
}
#file_index_falseset
#file_index_falseset3
#file_index_falseset3_trypsin
#file_index_falseset3_chymotrypsin
#false_set_test_chymotrypsin_231211.tsv
#false_set_test_fragdf_chymotrypsin_231211.tsv


#backup
{
  
  write_rds(fdr_pcc_negative_set,"./export/fdr_pcc_negative_set_231211_old.rds")
  write_rds(file_index_falseset,"./export/file_index_falseset_231211.rds")
  write_rds(file_index_falseset3,"./export/file_index_falseset3_231211.rds")
}

##neg_set
{
  read_tsv(
    paste0(folder_alphapeptdeep,"inrich/output/false_output_metrics_finetuned_model_trypsin_231211.tsv")
  ) %>% 
    select(-...1) %>% 
    mutate(
      false_set = 1,
      enzyme = "trypsin"
    ) %>% 
    bind_rows(
      read_tsv(
        paste0(folder_alphapeptdeep,"inrich/output/false_output_metrics_finetuned_model_chymotrypsin_231211.tsv")
      ) %>% 
        select(-...1) %>% 
        mutate(
          false_set = 1,
          enzyme = "chymotrypsin"
        )
    )-> fdr_pcc_negative_set
  
  read_tsv(
    paste0(folder_alphapeptdeep,"false_set_test_fragdf_trypsin_231211.tsv")
  ) %>% 
    bind_rows(
      read_tsv(
        paste0(folder_alphapeptdeep,"false_set_test_fragdf_chymotrypsin_231211.tsv")
      )
    ) -> fdr_pcc_negative_set_fragdf
  
  read_tsv(
    paste0(folder_alphapeptdeep,"inrich/output/false_output_preddf_trypsin_231211.tsv")
  ) %>% 
    bind_rows(
      read_tsv(
        paste0(folder_alphapeptdeep,"inrich/output/false_output_preddf_chymotrypsin_231211.tsv")
      )
    ) %>%
    select(-c(...1,b_modloss_z1,b_modloss_z2,y_modloss_z1,y_modloss_z2)) %>% 
    rename(
      b_z1_pred = b_z1,
      b_z2_pred = b_z2,
      y_z1_pred = y_z1,
      y_z2_pred = y_z2
    ) -> fdr_pcc_negative_set_preddf
  
  fdr_pcc_negative_set %>% 
    mutate(
      fragindex = map2(frag_start_idx,frag_stop_idx-1,\(x,y) seq(x,y,1))
    ) %>%
    unnest(fragindex) %>% 
    bind_cols(
      fdr_pcc_negative_set_fragdf,
      fdr_pcc_negative_set_preddf
    ) %>% 
    mutate(
      fragindex2=fragindex
    ) %>% 
    nest(
      fragdf = fragindex:y_z2,
      preddf = c(fragindex2,b_z1_pred:y_z2_pred)
    ) %>% 
    mutate(
      fragdf2 = map(fragdf,\(x) mutate(x,across(.cols=c(b_z1:y_z2),.fns=~na_if(.,y=0)))),
      preddf2 = map(preddf,\(x) mutate(x,across(.cols=c(b_z1_pred:y_z2_pred),.fns=~na_if(.,y=0))))
    ) %>% 
    mutate(
      bion_count = map_dbl(fragdf2,\(x) length(na.omit(append(pull(x,b_z1),pull(x,b_z2)))))
    ) %>%
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
    ) -> fdr_pcc_negative_set
  
  
  file_index_falseset_fdr %>%
    bind_cols(
      fdr_pcc_negative_set %>% select(PCC,COS,SA,SPC,false_set,fragdf,preddf,fragdf2,preddf2,bion_count,bion_pcc,yion_pcc)
    )-> fdr_pcc_negative_set2
}

#DB true positive
{
  tibble(
    strip_seq_dbmatch = fdr_pcc_negative_set2 %>% distinct(strip_seq_dbmatch) %>% pull(strip_seq_dbmatch)
  ) %>% 
    mutate(
      p2p1 = map(
        strip_seq_dbmatch,
        \(x) dbtable %>% filter(str_detect(sequence,x)),
        .progress = TRUE
      )
    ) %>% 
    mutate(
      protein_count = map_dbl(p2p1,\(x) nrow(x))
    ) %>% 
    unnest(p2p1) %>% 
    mutate(
      pos = str_locate(sequence,strip_seq_dbmatch)[,"start"]
    ) %>% 
    select(-protein) %>% 
    mutate(
      p2p1 = str_sub(sequence,start = pos-2,end =pos-1)
    ) %>% 
    select(-sequence) %>% 
    mutate(
      tp = if_else(p2p1 %in% c("KR","RR"),1,0)
    ) %>% 
    mutate(
      seq_arg=paste0("R",strip_seq_dbmatch)
    ) -> fdr_pcc_negative_set2_dbcheck
  fdr_pcc_negative_set2_dbcheck %>%
    left_join(
      fdr_pcc_negative_set2_dbcheck %>% 
        summarise(
          .by = strip_seq_dbmatch,
          tp_avg = mean(tp)
        )
    )-> fdr_pcc_negative_set2_dbcheck
}

##tp merge fdr_pcc_negative_set2 + fdr_pcc_negative_set2_dbcheck
{
  fdr_pcc_negative_set2 %>% 
    left_join(
      fdr_pcc_negative_set2_dbcheck %>% select(strip_seq_dbmatch,p2p1,tp_avg)
    ) %>% 
    nest(p2p1=p2p1) %>% 
    mutate(
      p2p1 = map_chr(p2p1,\(x) pull(x) %>% paste0(collapse=";"))
    ) -> fdr_pcc_negative_set3
  
  fdr_pcc_negative_set3 %>% pull(p2p1)
  
  fdr_pcc_negative_set3 -> fdr_pcc_negative_set4
  
}

##FDR
{
  read_tsv("F:/PDdata/newargnrich/200821_trainingset_semi_PSMs.txt") %>% 
    pull(Sequence) %>% unique() -> peptides_tp
  
  fdr_pcc_negative_set4 %>% colnames()
  # fdr_pcc_negative_set4 %>% 
  #   # filter(Confidence=="High") %>% 
  #   mutate(
  #     true_set = if_else(strip_seq %in% peptides_tandem_training_set | strip_seq %in% peptides_training_set,"true","false")
  #   )-> fdr_pcc_negative_set4
  # 
  fdr_pcc_negative_set4 %>%
    mutate(
      true_set = if_else(p2p1 %in% c("NA"),"false","true")
    )-> fdr_pcc_negative_set4
  
  # fdr_pcc_negative_set4 %>%
  #   mutate(
  #     true_set = if_else(strip_seq %in% peptides_tandem_training_set,"true","false")
  #   )-> fdr_pcc_negative_set4
  
  # fdr_pcc_negative_set4 %>%
  #   mutate(
  #     true_set = case_when(
  #       strip_seq %in% peptides_tandem_training_set ~ "true",
  #       strip_seq %in% peptides_training_set ~ "true",
  #       str_detect(p2p1,"[RK]R") ~ "true",
  #       .default = "false"
  #     )
  #   )-> fdr_pcc_negative_set4
  
  fdr_pcc_negative_set4 %>% 
    summarize(
      .by = true_set,
      count =n()
    )
  
  fdr_pcc_negative_set4 %>% 
    select(
      Confidence,Sequence,percolator_score,PCC,bion_pcc,yion_pcc,p2p1,true_set
    ) %>% write_tsv("./export/fdr_pcc_db_modified_set4.tsv")
}

##ROC
{
  library(pROC)
  
  RColorBrewer::brewer.pal(n=9,"Set1")
  
  # fdr_pcc_negative_set4 %>% filter(
  #   !is.na(true_set)
  # ) -> fdr_pcc_negative_set5
  
  fdr_pcc_negative_set4 %>% 
    mutate(
      enzyme = case_when(
        str_detect(peak_file,"\\-1196") ~ "trypsin",
        str_detect(peak_file,"\\-1199") ~ "chymotrypsin"
      )
    ) %>% 
    # filter(
    #   enzyme == "trypsin"
    # ) %>%
    mutate(
      fdr = if_else(true_set =="true",0,1)
    ) %>% 
    # filter(
    #   p2p1 != "NA"
    # ) %>%
    mutate(
      percolator_score = as.numeric(percolator_score)
    )-> fdr_pcc_negative_set5
  
  fdr_pcc_negative_set5 %>% 
    summarize(
      .by = true_set,
      count =n()
    )
  
  fdr_pcc_negative_set5 %>% 
    ggplot(
      aes(x=fct(true_set,levels = c("true","false")),y=PCC, fill = fct(true_set,levels = c("true","false")))
    )+
    geom_violin(
    )+
    geom_boxplot(width = 0.02, outlier.size = 0.01)+
    labs(
      x=NA,
      y="PCC",
      fill = NA
    )+
    scale_y_continuous(breaks = seq(-0.2,1,0.1))+
    scale_x_discrete(
      labels = c("TRUE","FALSE")
    )+
    facet_wrap(~enzyme)+
    theme_classic()+
    theme(
      legend.position = "none",
      legend.title = element_blank(),
      text= element_text(color = "black", size = 10),
      axis.text.y = element_text(color = "black", size = 10, hjust=1),
      axis.text.x = element_text(color = "black", size = 10),
      axis.title.x = element_blank()
    )
  ggsave("./figures/fdr_distribution_pcc.png",width=10,height=10,dpi = 600,units="cm")
  
  fdr_pcc_negative_set5 %>% 
    ggplot(
      aes(x=fct(true_set,levels = c("true","false")),y=percolator_score, fill = fct(true_set,levels = c("true","false")))
    )+
    geom_violin(
    )+
    geom_boxplot(width = 0.02, outlier.size = 0.01)+
    labs(
      #x="PCC",
      y="Percolator Score"
    )+
    #scale_y_continuous(breaks = seq(-0.2,1,0.1))+
    scale_x_discrete(
      labels = c("TRUE","FALSE")
    )+
    facet_wrap(~enzyme)+
    # annotate(
    #   "text",
    #   x=1:2,
    #   y = 5,
    #   label = c("n = 11263","n = 144")
    # )+
    theme_classic()+
    theme(
      legend.position = "none",
      legend.title = element_blank(),
      text= element_text(color = "black", size = 10),
      axis.text.y = element_text(color = "black", size = 10, hjust=1),
      axis.text.x = element_text(color = "black", size = 10),
      axis.title.x = element_blank()
    )
  ggsave("./figures/fdr_distribution_percolator.png",width=10,height=10,dpi = 600,units="cm")
  
  fdr_pcc_negative_set5 %>% arrange(desc(PCC)) %>% mutate(FDR=cummean(fdr)) %>% filter(FDR<=0.01) %>% slice_tail(n=1)  %>% .$PCC -> label_pcc_fdr001
  fdr_pcc_negative_set5 %>% arrange(desc(bion_pcc)) %>% mutate(FDR=cummean(fdr)) %>% filter(FDR<=0.01) %>% slice_tail(n=1)  %>% .$bion_pcc -> label_bpcc_fdr001
  fdr_pcc_negative_set5 %>% arrange(desc(yion_pcc)) %>% mutate(FDR=cummean(fdr)) %>% filter(FDR<=0.01) %>% slice_tail(n=1)  %>% .$yion_pcc -> label_ypcc_fdr001
  
  ##fig:fdr_trypsin_pcc.png
  {
    fdr_pcc_negative_set5 %>% 
      filter(enzyme == "trypsin") %>% 
      arrange(desc(PCC)) %>% 
      mutate(
        PCC_index = 1:nrow(.)
      ) %>% 
      mutate(
        FDR1 = cummean(fdr)
      ) %>% 
      arrange(desc(bion_pcc)) %>%
      mutate(
        bion_index = 1:nrow(.)
      ) %>% 
      mutate(
        FDR2 = cummean(fdr)
      ) %>% 
      arrange(desc(percolator_score)) %>% 
      mutate(
        pecolator_index = 1:nrow(.)
      ) %>% 
      mutate(
        FDR3 = cummean(fdr)
      ) %>% 
      mutate(
        fdr = if_else(true_set =="true",0,1)
      ) %>% 
      ggplot(
      )+
      geom_point(
        aes(x=PCC_index,y=FDR1),color="#E41A1C",size=0.2
      )+
      # geom_point(
      #   aes(x=bion_index,y=FDR2),color="#377EB8",size=0.1
      # )+
      # geom_point(
      #   aes(x=pecolator_index,y=FDR3),color="#4DAF4A",size=0.1
      # )+
      geom_hline(yintercept = 0.01,linetype = 2)+
      labs(
        x="index",
        y="FDR"
      )+
      scale_y_continuous(
        labels = scales::percent,
        breaks = seq(0,0.05,0.01),
        limits = c(0,0.05)
      )+
      scale_x_reverse(
        breaks = seq(0,15000,2000)
      )+
      theme_classic()+
      theme(
        legend.position = "none",
        legend.title = element_blank(),
        text= element_text(color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10, hjust=1),
        axis.text.x = element_text(color = "black", size = 10),
        #axis.title.x = element_blank()
      )
    ggsave("./figures/fdr_trypsin_pcc.png",width=8,height=8,dpi = 600,units="cm")
  }
  
  ##fig:fdr_chymotrypsin_pcc.png
  {
    fdr_pcc_negative_set5 %>% 
      filter(enzyme == "chymotrypsin") %>% 
      arrange(desc(PCC)) %>% 
      mutate(
        PCC_index = 1:nrow(.)
      ) %>% 
      mutate(
        FDR1 = cummean(fdr)
      ) %>% 
      arrange(desc(bion_pcc)) %>%
      mutate(
        bion_index = 1:nrow(.)
      ) %>% 
      mutate(
        FDR2 = cummean(fdr)
      ) %>% 
      arrange(desc(percolator_score)) %>% 
      mutate(
        pecolator_index = 1:nrow(.)
      ) %>% 
      mutate(
        FDR3 = cummean(fdr)
      ) %>% 
      mutate(
        fdr = if_else(true_set =="true",0,1)
      ) %>% 
      ggplot(
      )+
      geom_point(
        aes(x=PCC_index,y=FDR1),color="#E41A1C",size=0.1
      )+
      # geom_point(
      #   aes(x=bion_index,y=FDR2),color="#377EB8",size=0.1
      # )+
      # geom_point(
      #   aes(x=pecolator_index,y=FDR3),color="#4DAF4A",size=0.1
      # )+
      geom_hline(yintercept = 0.01,linetype = 2)+
      labs(
        x="index",
        y="FDR"
      )+
      scale_y_continuous(
        labels = scales::percent,
        breaks = seq(0,0.05,0.01),
        limits = c(0,0.05)
      )+
      scale_x_reverse()+
      theme_classic()+
      theme(
        legend.position = "none",
        legend.title = element_blank(),
        text= element_text(color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10, hjust=1),
        axis.text.x = element_text(color = "black", size = 10),
        #axis.title.x = element_blank()
      )
    ggsave("./figures/fdr_chymotrypsin_pcc.png",width=10,height=10,dpi = 600,units="cm")
  }
  
  ##fig:fdr_trypsin_percolator.png
  {
    fdr_pcc_negative_set5 %>% 
      filter(enzyme == "trypsin") %>% 
      arrange(desc(percolator_score)) %>% 
      mutate(
        pecolator_index = 1:nrow(.)
      ) %>% 
      mutate(
        FDR3 = cummean(fdr)
      ) %>% 
      mutate(
        fdr = if_else(true_set =="true",0,1)
      ) %>% 
      ggplot(
      )+
      geom_point(
        aes(x=pecolator_index,y=FDR3),color="#4DAF4A",size=0.2
      )+
      # geom_point(
      #   aes(x=bion_index,y=FDR2),color="#377EB8",size=0.1
      # )+
      # geom_point(
      #   aes(x=pecolator_index,y=FDR3),color="#4DAF4A",size=0.1
      # )+
      geom_hline(yintercept = 0.01,linetype = 2)+
      labs(
        x="index",
        y="FDR"
      )+
      scale_y_continuous(
        labels = scales::percent,
        breaks = seq(0,0.05,0.01),
        limits = c(0,0.05)
      )+
      scale_x_reverse(
        breaks = seq(0,15000,2000)
      )+
      theme_classic()+
      theme(
        legend.position = "none",
        legend.title = element_blank(),
        text= element_text(color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10, hjust=1),
        axis.text.x = element_text(color = "black", size = 10),
        #axis.title.x = element_blank()
      )
    ggsave("./figures/fdr_trypsin_percolator.png",width=8,height=8,dpi = 600,units="cm")
  }
  
  ##fig:fdr_chymotrypsin_percolator.png
  {
    fdr_pcc_negative_set5 %>% 
      filter(enzyme == "chymotrypsin") %>% 
      arrange(desc(percolator_score)) %>% 
      mutate(
        pecolator_index = 1:nrow(.)
      ) %>% 
      mutate(
        FDR3 = cummean(fdr)
      ) %>% 
      mutate(
        fdr = if_else(true_set =="true",0,1)
      ) %>% 
      ggplot(
      )+
      geom_point(
        aes(x=pecolator_index,y=FDR3),color="#4DAF4A",size=0.1
      )+
      # geom_point(
      #   aes(x=bion_index,y=FDR2),color="#377EB8",size=0.1
      # )+
      # geom_point(
      #   aes(x=pecolator_index,y=FDR3),color="#4DAF4A",size=0.1
      # )+
      geom_hline(yintercept = 0.01,linetype = 2)+
      labs(
        x="index",
        y="FDR"
      )+
      scale_y_continuous(
        labels = scales::percent,
        breaks = seq(0,0.05,0.01),
        limits = c(0,0.05)
      )+
      scale_x_reverse()+
      theme_classic()+
      theme(
        legend.position = "none",
        legend.title = element_blank(),
        text= element_text(color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10, hjust=1),
        axis.text.x = element_text(color = "black", size = 10),
        #axis.title.x = element_blank()
      )
    ggsave("./figures/fdr_chymotrypsin_percolator.png",width=10,height=10,dpi = 600,units="cm")
  }
  
  label_pcc_fdr001 ##0.868823
  ##export raw data 
  {
    fdr_pcc_negative_set5 %>% 
      filter(enzyme == "trypsin") %>% 
      arrange(desc(PCC)) %>% 
      mutate(
        PCC_index = 1:nrow(.)
      ) %>% 
      mutate(
        FDR1 = cummean(fdr)
      ) %>% 
      arrange(desc(bion_pcc)) %>%
      mutate(
        bion_index = 1:nrow(.)
      ) %>% 
      mutate(
        FDR2 = cummean(fdr)
      ) %>% 
      arrange(desc(percolator_score)) %>% 
      mutate(
        pecolator_index = 1:nrow(.)
      ) %>% 
      mutate(
        FDR3 = cummean(fdr)
      ) %>% 
      mutate(
        fdr = if_else(true_set =="true",0,1)
      ) %>% write_tsv("./export/fdr_raw_data.tsv")
  }
  #fdr_raw_data.tsv
  
  ##spectral_FDR, spectral_FNR
  {
    fdr_pcc_negative_set5 %>% 
      mutate(
        fragdf3 = map2(
          fragdf2,
          preddf2,
          \(x,y) left_join(x,y,join_by(fragindex == fragindex2))
        )
      ) -> fdr_pcc_negative_set7
    
    fdr_pcc_negative_set7 %>% 
      mutate(
        fragdf3 = map(
          fragdf3,
          \(x) replace_na(x,list(b_z1=0,b_z2=0,y_z1=0,y_z2=0,b_z1_pred=0,b_z2_pred=0,y_z1_pred=0,y_z2_pred=0))
        )
      ) -> fdr_pcc_negative_set7
    
    fdr_pcc_negative_set7 %>%
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
      ) -> fdr_pcc_negative_set7
    
  }
  
  
  
  ##ROC
  {
    # fdr_pcc_negative_set5 %>% 
    #   filter(Confidence == "High")->fdr_pcc_negative_set6
    
    
    fdr_pcc_negative_set5 %>% 
      filter(enzyme == "trypsin") %>% 
      arrange(desc(PCC)) %>% 
      mutate(
        PCC_index = 1:nrow(.)
      ) %>% 
      mutate(
        FDR1 = cummean(fdr)
      ) %>% 
      arrange(desc(bion_pcc)) %>%
      mutate(
        bion_index = 1:nrow(.)
      ) %>% 
      mutate(
        FDR2 = cummean(fdr)
      ) %>% 
      arrange(desc(percolator_score)) %>% 
      mutate(
        pecolator_index = 1:nrow(.)
      ) %>% 
      mutate(
        FDR3 = cummean(fdr)
      ) %>% 
      arrange(desc(yion_pcc)) %>%
      mutate(
        yion_index = 1:nrow(.)
      ) %>% 
      mutate(
        FDR4 = cummean(fdr)
      ) %>% 
      mutate(
        fdr = if_else(true_set =="true",0,1)
      ) -> fdr_pcc_negative_set6
    
    fdr_pcc_negative_set7 %>% 
      filter(enzyme == "trypsin") -> fdr_pcc_negative_set7_trypsin
    ggroc(
      list(
        PCC=roc(response=fdr_pcc_negative_set7_trypsin$true_set,predictor = fdr_pcc_negative_set7_trypsin$PCC), #0.6236
        COS=roc(response=fdr_pcc_negative_set7_trypsin$true_set,predictor = fdr_pcc_negative_set7_trypsin$COS), #0.622
        SPC=roc(response=fdr_pcc_negative_set7_trypsin$true_set,predictor = fdr_pcc_negative_set7_trypsin$SPC), #0.5721
        spec_FNR=roc(response=fdr_pcc_negative_set7_trypsin$true_set,predictor = fdr_pcc_negative_set7_trypsin$spec_FNR), #0.606
        spec_FDR=roc(response=fdr_pcc_negative_set7_trypsin$true_set,predictor = fdr_pcc_negative_set7_trypsin$spec_FDR), #0.4986
        percolator=roc(response=fdr_pcc_negative_set7_trypsin$true_set,predictor = fdr_pcc_negative_set7_trypsin$percolator_score) #0.4979
      )
    )+
      scale_color_manual(
        values =c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#A65628")
      )+
      labs(
        y="Sensitivity",
        x="Specificity"
      )+
      theme_classic()+
      theme(
        legend.position = "none",
        legend.title = element_blank(),
        text= element_text(color = "black", size = 7),
        axis.text.y = element_text(color = "black", size = 7, hjust=1),
        axis.text.x = element_text(color = "black", size = 7),
        #axis.title.x = element_blank()
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm")
      )
    ggsave("./figures/fdr_roc_trypsin.png",width=4.5,height=4.5,dpi = 600,units="cm")
  }
  
  #fdr_roc_chymotrypsin.png
  {
    fdr_pcc_negative_set7 %>% 
      filter(enzyme == "chymotrypsin") -> fdr_pcc_negative_set7_chymotrypsin
    ggroc(
      list(
        PCC=roc(response=fdr_pcc_negative_set7_chymotrypsin$true_set,predictor = fdr_pcc_negative_set7_chymotrypsin$PCC), #0.4549
        COS=roc(response=fdr_pcc_negative_set7_chymotrypsin$true_set,predictor = fdr_pcc_negative_set7_chymotrypsin$COS), #0.4901
        SPC=roc(response=fdr_pcc_negative_set7_chymotrypsin$true_set,predictor = fdr_pcc_negative_set7_chymotrypsin$SPC), #0.5779
        spec_FNR=roc(response=fdr_pcc_negative_set7_chymotrypsin$true_set,predictor = fdr_pcc_negative_set7_chymotrypsin$spec_FNR), #0.7966
        spec_FDR=roc(response=fdr_pcc_negative_set7_chymotrypsin$true_set,predictor = fdr_pcc_negative_set7_chymotrypsin$spec_FDR), #0.7629
        percolator=roc(response=fdr_pcc_negative_set7_chymotrypsin$true_set,predictor = fdr_pcc_negative_set7_chymotrypsin$percolator_score) #0.688
      )
    )+
      scale_color_manual(
        values =c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#A65628")
      )+
      theme_classic()+
      theme(
        legend.position = "none",
        legend.title = element_blank(),
        text= element_text(color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10, hjust=1),
        axis.text.x = element_text(color = "black", size = 10),
        #axis.title.x = element_blank()
      )
    ggsave("./figures/fdr_roc_chymotrypsin.png",width=10,height=10,dpi = 600,units="cm")
  }
}

##FDR fdr_roc_v5_list.png
{
  tibble(
    PCC=as.numeric(auc(roc(response=fdr_pcc_negative_set7_trypsin$true_set,predictor = fdr_pcc_negative_set7_trypsin$PCC))), #0.6236
    COS=as.numeric(auc(roc(response=fdr_pcc_negative_set7_trypsin$true_set,predictor = fdr_pcc_negative_set7_trypsin$COS))), #0.622
    SPC=as.numeric(auc(roc(response=fdr_pcc_negative_set7_trypsin$true_set,predictor = fdr_pcc_negative_set7_trypsin$SPC))), #0.5721
    spec_FNR=as.numeric(auc(roc(response=fdr_pcc_negative_set7_trypsin$true_set,predictor = fdr_pcc_negative_set7_trypsin$spec_FNR))), #0.606
    spec_FDR=as.numeric(auc(roc(response=fdr_pcc_negative_set7_trypsin$true_set,predictor = fdr_pcc_negative_set7_trypsin$spec_FDR))), #0.4986
    percolator=as.numeric(auc(roc(response=fdr_pcc_negative_set7_trypsin$true_set,predictor = fdr_pcc_negative_set7_trypsin$percolator_score))) #0.4979
  ) %>% 
    pivot_longer(everything()) %>% 
    ggplot(
      aes(x=fct(name,levels=c("PCC","COS","SPC","spec_FNR","spec_FDR","percolator")),y=value)
    )+
    geom_bar(stat="identity",fill="black",width = 0.4)+
    geom_text(aes(label = round(value,3)),hjust=-0.1, size =1.5,angle =90)+
    labs(y="AUC of ROC")+
    scale_y_continuous(limits = c(0,0.7))+
    scale_x_discrete(
      labels=c("PCC","COS","SPC","Spectral FNR","Spectral FDR","Percolator")
    )+
    theme_classic2()+
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      text= element_text(color = "black", size = 5),
      axis.text.y = element_text(color = "black", size = 5, hjust=1),
      axis.text.x = element_text(color = "black", size = 5, angle=45,hjust=1),
      #axis.title.x = element_blank()
    )
  ggsave("./figures/fdr_roc_v5_list.png",width=5,height=5,dpi = 600,units="cm")
}

###
{
  fdr_pcc_negative_set7_trypsin
  
  fdr_pcc_negative_set7_trypsin %>% 
    arrange(desc(PCC)) %>% 
    mutate(
      PCC_index = 1:nrow(.)
    ) %>% 
    mutate(
      FDR1 = cummean(fdr)
    ) %>% 
    arrange(desc(percolator_score)) %>% 
    mutate(
      pecolator_index = 1:nrow(.)
    ) %>% 
    mutate(
      FDR2 = cummean(fdr)
    ) %>% 
    mutate(
      fdr = if_else(true_set =="true",0,1)
    ) %>% 
    ggplot(
    )+
    geom_point(
      aes(x=PCC_index,y=FDR1),color="#FF0000",size=0.05,alpha = 0.1
    )+
    geom_point(
      aes(x=pecolator_index,y=FDR2),color="#00A08A",size=0.05,alpha = 0.1
    )+
    geom_hline(yintercept = 0.01,linetype = 2)+
    labs(
      x="Index",
      y="FDR"
    )+
    scale_y_continuous(
      labels = scales::percent,
      breaks = seq(0,0.05,0.01),
      limits = c(0,0.05)
    )+
    scale_x_reverse()+
    theme_classic2()+
    theme(
      legend.position = "none",
      #axis.title.x = element_blank(),
      text= element_text(color = "black", size = 8),
      axis.text.y = element_text(color = "black", size = 7, hjust=1),
      axis.text.x = element_text(color = "black", size = 7),
      #axis.title.x = element_blank()
      plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm")
    )
  ggsave("./figures/2_06_new_fdrcurve.png",width=4.5,height=4.5,dpi = 600,units="cm")
  wes_palette(name = "Darjeeling1",n=2)[2]
}
RColorBrewer::brewer.pal(9,"Set1")
label_pcc_fdr001
