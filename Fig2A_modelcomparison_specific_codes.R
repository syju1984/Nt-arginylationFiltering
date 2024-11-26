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


##original codes

##chapter ms2 prediction model
###output_prediction_comparison ##9569
##output_prediction_comparison_limit
{
  ##read alphapeptdeep prediction result of trypsin 9569
  read_tsv(
    "C:/Users/admin/alphapeptdeep/project/alphapeptdeep_115/alphapeptdeep/inrich/output/prediction_17669_output_metrics_pretrained_model_trypsin_231206.tsv"
  ) %>% 
    mutate(
      model = "pretrained"
    ) %>% 
    bind_rows(
      read_tsv(
        "C:/Users/admin/alphapeptdeep/project/alphapeptdeep_115/alphapeptdeep/inrich/output/prediction_17669_output_metrics_built_model_trypsin_231206.tsv"
      ) %>% 
        mutate(
          model = "built"
        )
    ) %>% 
    bind_rows(
      read_tsv(
        "C:/Users/admin/alphapeptdeep/project/alphapeptdeep_115/alphapeptdeep/inrich/output/prediction_17669_output_metrics_finetuned_model_trypsin_231206.tsv"
      ) %>% 
        mutate(
          model = "finetuned"
        )
    ) %>% 
    bind_rows(
      read_tsv(
        "C:/Users/admin/alphapeptdeep/project/alphapeptdeep_115/alphapeptdeep/inrich/output/prediction_9569_output_metrics_311547woarg_model_trypsin_231213.tsv"
      ) %>% 
        mutate(
          model = "finetunedwoarg"
        )
    )-> output_prediction_comparison_trypsin_241009
  
  output_prediction_comparison_trypsin_241009 %>%
    mutate(
      model = fct(model,levels=c("pretrained","built","finetuned","finetunedwoarg"))
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
    theme_classic2()+
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.title.y = element_text(color = "black",size = 7),
      axis.text.x = element_blank(),
      axis.text.y = element_text(color = "black",size = 7)
    )
  ggsave("./figures/2_03_finetuning1.png",width =4.5, height = 4,dpi = 600,units = "cm")
  
  output_prediction_comparison_trypsin_241009 %>% 
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

##
{
  read_rds(
    "./export/psm_17669.rds"
  ) -> psm_17669
  
  psm_17669 %>% 
    filter(
      enzyme == "trypsin"
    ) %>% 
    mutate(
      index = 1:nrow(.)
    ) -> psm_9569
  
  psm_9569 %>% 
    select(index,frag_df) %>% 
    unnest(frag_df) %>% 
    mutate(
      across(everything(),\(x) replace_na(x,0))
    ) %>% 
    bind_cols(
      read_tsv(
        "C:/Users/admin/alphapeptdeep/project/alphapeptdeep_115/alphapeptdeep/inrich/output/prediction_17669_predict_inten_df_pretrained_model_trypsin_231206.tsv"
      ) %>% select(
        b_z1,b_z2,y_z1,y_z2
      ) %>% 
        rename_with(
          .fn = \(x) paste0(x,"_pretrained")
        ),
      read_tsv(
        "C:/Users/admin/alphapeptdeep/project/alphapeptdeep_115/alphapeptdeep/inrich/output/prediction_17669_predict_inten_df_built_model_trypsin_231206.tsv"
      ) %>% select(
        b_z1,b_z2,y_z1,y_z2
      ) %>% 
        rename_with(
          .fn = \(x) paste0(x,"_built")
        ),
      read_tsv(
        "C:/Users/admin/alphapeptdeep/project/alphapeptdeep_115/alphapeptdeep/inrich/output/prediction_17669_predict_inten_df_finetuned_model_trypsin_231206.tsv"
      ) %>% select(
        b_z1,b_z2,y_z1,y_z2
      ) %>% 
        rename_with(
          .fn = \(x) paste0(x,"_finetuned")
        ),
      read_tsv(
        "C:/Users/admin/alphapeptdeep/project/alphapeptdeep_115/alphapeptdeep/inrich/output/prediction_9569_predicted_inten_df_311547woarg_model_trypsin_231213.tsv"
      ) %>% select(
        b_z1,b_z2,y_z1,y_z2
      ) %>% 
        rename_with(
          .fn = \(x) paste0(x,"_finetunedwoarg")
        )
    ) %>% 
    select(
      -c(frag_index,ion_num2)
    ) %>% 
    nest(
      data = c(b_z1,b_z2,y_z1,y_z2),
      pretrained = contains("pretrained"),
      built = contains("built"),
      finetuned = contains("finetuned"),
      finetunedwoarg = contains("finetunedwoarg")
    ) -> psm_9569_ionspecies_test_241010
  
  psm_9569_ionspecies_test_241010 %>% #head() %>% pull(data)
    mutate(
      pretrained_pcc = map2_dbl(
        data,
        pretrained,
        \(x,y) cor(
          append(
            pull(x,b_z1),
            pull(x,b_z2)
          ) %>% append(
            pull(x,y_z1)
          ) %>% 
            append(
              pull(x,y_z2)
            )  ,
          append(
            pull(y,b_z1_pretrained),
            pull(y,b_z2_pretrained)
          ) %>% append(
            pull(y,y_z1_pretrained)
          ) %>% 
            append(
              pull(y,y_z2_pretrained)
            ),
          method = "pearson"
        )
      ),
      pretrained_pcc_b = map2_dbl(
        data,
        pretrained,
        \(x,y) cor(
          append(pull(x,b_z1),pull(x,b_z2)),
          append(pull(y,b_z1_pretrained),pull(y,b_z2_pretrained)),
          method = "pearson"
        )
      ),
      pretrained_pcc_y = map2_dbl(
        data,
        pretrained,
        \(x,y) cor(
          append(pull(x,y_z1),pull(x,y_z2)),
          append(pull(y,y_z1_pretrained),pull(y,y_z2_pretrained)),
          method = "pearson"
        )
      )
    ) %>% 
    mutate(
      built_pcc = map2_dbl(
        data,
        built,
        \(x,y) cor(
          append(
            pull(x,b_z1),
            pull(x,b_z2)
          ) %>% append(
            pull(x,y_z1)
          ) %>% 
            append(
              pull(x,y_z2)
            )  ,
          append(
            pull(y,b_z1_built),
            pull(y,b_z2_built)
          ) %>% append(
            pull(y,y_z1_built)
          ) %>% 
            append(
              pull(y,y_z2_built)
            ),
          method = "pearson"
        )
      ),
      built_pcc_b = map2_dbl(
        data,
        built,
        \(x,y) cor(
          append(pull(x,b_z1),pull(x,b_z2)),
          append(pull(y,b_z1_built),pull(y,b_z2_built)),
          method = "pearson"
        )
      ),
      built_pcc_y = map2_dbl(
        data,
        built,
        \(x,y) cor(
          append(pull(x,y_z1),pull(x,y_z2)),
          append(pull(y,y_z1_built),pull(y,y_z2_built)),
          method = "pearson"
        )
      )
    )%>% 
    mutate(
      finetuned_pcc = map2_dbl(
        data,
        finetuned,
        \(x,y) cor(
          append(
            pull(x,b_z1),
            pull(x,b_z2)
          ) %>% append(
            pull(x,y_z1)
          ) %>% 
            append(
              pull(x,y_z2)
            )  ,
          append(
            pull(y,b_z1_finetuned),
            pull(y,b_z2_finetuned)
          ) %>% append(
            pull(y,y_z1_finetuned)
          ) %>% 
            append(
              pull(y,y_z2_finetuned)
            ),
          method = "pearson"
        )
      ),
      finetuned_pcc_b = map2_dbl(
        data,
        finetuned,
        \(x,y) cor(
          append(pull(x,b_z1),pull(x,b_z2)),
          append(pull(y,b_z1_finetuned),pull(y,b_z2_finetuned)),
          method = "pearson"
        )
      ),
      finetuned_pcc_y = map2_dbl(
        data,
        finetuned,
        \(x,y) cor(
          append(pull(x,y_z1),pull(x,y_z2)),
          append(pull(y,y_z1_finetuned),pull(y,y_z2_finetuned)),
          method = "pearson"
        )
      )
    )%>% 
    mutate(
      finetunedwoarg_pcc = map2_dbl(
        data,
        finetunedwoarg,
        \(x,y) cor(
          append(
            pull(x,b_z1),
            pull(x,b_z2)
          ) %>% append(
            pull(x,y_z1)
          ) %>% 
            append(
              pull(x,y_z2)
            )  ,
          append(
            pull(y,b_z1_finetunedwoarg),
            pull(y,b_z2_finetunedwoarg)
          ) %>% append(
            pull(y,y_z1_finetunedwoarg)
          ) %>% 
            append(
              pull(y,y_z2_finetunedwoarg)
            ),
          method = "pearson"
        )
      ),
      finetunedwoarg_pcc_b = map2_dbl(
        data,
        finetunedwoarg,
        \(x,y) cor(
          append(pull(x,b_z1),pull(x,b_z2)),
          append(pull(y,b_z1_finetunedwoarg),pull(y,b_z2_finetunedwoarg)),
          method = "pearson"
        )
      ),
      finetunedwoarg_pcc_y = map2_dbl(
        data,
        finetunedwoarg,
        \(x,y) cor(
          append(pull(x,y_z1),pull(x,y_z2)),
          append(pull(y,y_z1_finetunedwoarg),pull(y,y_z2_finetunedwoarg)),
          method = "pearson"
        )
      )
    ) -> psm_9569_ionspecies_test_241010_done
  
  psm_9569_ionspecies_test_241010_done
  
  psm_9569_ionspecies_test_241010_done %>% 
    select(
      -c(data,pretrained,built,finetuned,finetunedwoarg)
    ) %>% 
    pivot_longer(
      cols = -c(index)
    ) -> psm_9569_ionspecies_test_241010_done
  
  psm_9569_ionspecies_test_241010_done %>% 
    mutate(
      model = str_split_i(name,"_",1)
    ) %>% 
    filter(
      !str_detect(name,"pcc$")
    ) %>% 
    mutate(
      frag_series = str_split_i(name,"_",3)
    ) %>% 
    mutate(
      model = factor(model,levels=c("pretrained","built","finetuned","finetunedwoarg"),labels=c("Pretrained","Built","Fine-tuned","Fine-tuned w/o \n R-starting Peptides"))
    ) %>% 
    ggplot(aes(x=frag_series,y=value))+
    geom_violin(aes(fill = frag_series),linewidth = 0.1,scale = "width")+
    labs(y="PCC")+
    geom_boxplot(width = 0.25,outlier.shape = NA,linewidth = 0.1)+
    # scale_x_discrete(labels =c("Pretrained","Built","Fine-tuned","Fine-tuned w/o\nR-starting peptides"))+
    scale_fill_manual(
      values = c(wes_palette("Zissou1")[1],wes_palette("Zissou1")[5])
    )+
    scale_x_discrete(
      labels = c("b-ion","y-ion")
    )+
    scale_y_continuous(
      limits=c(0.5,1)
    )+
    facet_wrap(
      facets = vars(model),
      nrow = 1
    )+
    theme_classic2()+
    theme(
      legend.position = "none",
      strip.text = element_text(color = "black",size = 2),
      axis.title.x = element_blank(),
      axis.title.y = element_text(color = "black",size = 7),
      axis.text.x = element_blank(),
      axis.text.y = element_text(color = "black",size = 7),
      strip.background = element_blank(),
      strip.text.x = element_blank(),
      panel.spacing = unit(0.1, "cm"),
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank(),
      plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm")
    )
  ggsave("./figures/2_03_model_yions.png",width=4.5,height=3,dpi = 600,units="cm")
  
  read_tsv(
    "C:/Users/admin/alphapeptdeep/project/alphapeptdeep_115/alphapeptdeep/inrich/output/prediction_17669_output_metrics_pretrained_model_trypsin_231206.tsv"
  ) %>% 
    mutate(
      model = "pretrained"
    ) %>% 
    
}




##chymo
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
  
  
}