##RT prediction figure ::@summaryfile_psm_for_testing_pcc_rt
{
  ##figure
  for (i in 1:12){
    summaryfile_psm_for_training_result$set %>% unique -> set_names
    set_names[i]->set_name
    summaryfile_psm_for_training_result %>% 
      filter(set == set_name) %>% 
      ggplot(aes(x=rt_norm, y=rt_pred))+
      geom_point(size=0.5,alpha=0.5) +
      geom_line(aes(y=lwr), color = "red", linetype = "dashed")+
      geom_line(aes(y=upr), color = "red", linetype = "dashed")+
      geom_smooth(method=lm, se=TRUE)+
      lims(x=c(0,1),y=c(0,1))+
      labs(
        x="Observed Normalized RT",
        y="Predicted Normalized RT"
      )+
      theme_classic()
    ggsave(
      paste0("./figure2/rt_prediction_training_set",set_name,".png"),width=10,height=10,dpi = 600,units="cm"
    )
  }
  i=1
  for (i in 1:12){
    summaryfile_psm_for_testing_pcc_rt
    set_names[i]->set_name
    summaryfile_psm_for_testing_pcc_rt %>% 
      filter(set == set_name) %>% 
      ggplot(aes(x=rt_norm, y=rt_pred))+
      geom_point(
        aes(color = fct(as.character(pcc_check),levels=c("1","0"))),
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
      theme_classic()+
      theme(
        legend.position = "top"
      )
    ggsave(
      paste0("./figure2/rt_prediction_test_set",set_name,".png"),width=10,height=10,dpi = 600,units="cm"
    )
  }
  
  for (i in 1:12){
    summaryfile_psm_for_testing_pcc_rt
    set_names[i]->set_name
    summaryfile_psm_for_testing_pcc_rt %>% 
      filter(set == set_name) %>% 
      ggplot(aes(x=rt_norm, y=rt_pred))+
      geom_point(
        aes(color = fct(as.character(diagnostic_ion_check),levels=c("1","0"))),
        size=1,alpha=0.5
      ) +
      geom_line(aes(y=lwr), color = "red", linetype = "dashed")+
      geom_line(aes(y=upr), color = "red", linetype = "dashed")+
      geom_smooth(method=lm, se=TRUE)+
      scale_color_manual(labels = c("true","false"),values = c("red","black"))+
      lims(x=c(0,1),y=c(0,1))+
      labs(
        x="Observed Normalized RT",
        y="Predicted Normalized RT",
        color = "Diagnostic Ion"
      )+
      theme_classic()+
      theme(
        legend.position = "top"
      )
    ggsave(
      paste0("./figure2/rt_prediction_test_di_set",set_name,".png"),width=10,height=10,dpi = 600,units="cm"
    )
  }
  
}

##alphafold figure :: @alphafold_fetched_tidied2
#alphafold_pcc.png
#alphafold_pcc_rt.png
#alphafold_rt.png
#alphafold_wo_pcc_rt.png
{
  alphafold_fetched_tidied2 %>% 
    dplyr::filter(pcc_check == 1 & rt_check == 1) %>% 
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
      legend.position = "none"
    )
  ggsave("./figure2/alphafold_pcc_rt.png",width=10,height=10,dpi = 600,units="cm")
  
  alphafold_fetched_tidied2 %>% 
    dplyr::filter(pcc_check == 0 & rt_check == 0) %>% 
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
      legend.position = "none"
    )
  ggsave("./figure2/alphafold_wo_pcc_rt.png",width=10,height=10,dpi = 600,units="cm")
  
  # alphafold_fetched_tidied2 %>% 
  #   dplyr::filter(diagnostic_ion_check ==1) %>% 
  #   ggplot(
  #     aes(
  #       x=fct(as.character(mod_position)),
  #       y=prediction_score,
  #       fill = fct(as.character(mod_position))
  #     )
  #   )+
  #   geom_boxplot()+
  #   lims(y=c(20,100))+
  #   labs(
  #     x = "Arginylation Position", y = "AlphaFold pLDDT score"
  #   )+
  #   theme_classic()+
  #   theme(
  #     legend.position = "none"
  #   )
  # ggsave("./figure2/alphafold_diagnosticion.png",width=10,height=10,dpi = 600,units="cm")
  
  alphafold_fetched_tidied2 %>% 
    dplyr::filter(rt_check ==1) %>% 
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
      legend.position = "none"
    )
  ggsave("./figure2/alphafold_rt.png",width=10,height=10,dpi = 600,units="cm")
  
  alphafold_fetched_tidied2 %>% 
    dplyr::filter(pcc_check ==1) %>% 
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
      legend.position = "none"
    )
  ggsave("./figure2/alphafold_pcc.png",width=10,height=10,dpi = 600,units="cm")
  
  alphafold_fetched_tidied2 %>% 
    dplyr::filter(pcc_check ==1 & rt_check ==1) %>% 
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
      legend.position = "none"
    )
  ggsave("./figure2/alphafold_pcc_rt.png",width=10,height=10,dpi = 600,units="cm")
}

##prediction performance @output_prediction_comparison
#prediction_performance.png
{
  output_prediction_comparison %>% 
    pivot_longer(
      cols = c(PCC_build:PCC_pretrained)
    ) %>%
    ggplot(
      aes(x=name,y=value,fill=name)
    )+
    geom_boxplot()+
    scale_x_discrete(
      labels = c("Built","Pretrained","Transferred")
    )+
    scale_y_continuous(labels = scales::percent,limits=c(0,1))+
    scale_fill_manual(values = wes_palette("GrandBudapest2")[1:3])+
    labs(
      x="Model",
      y="PCC"
    )+
    theme_classic()+
    theme(
      legend.position = "none"
    )
  ggsave("./figure2/prediction_performance.png",width=10,height=10,dpi = 600,units = "cm")
}

##ms2 comparison figure! @summaryfile_psm_for_testing_pcc4
{
  summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered_trinity
  ###figure generation
  #for (j in 1:nrow(summaryfile_psm_for_testing_pcc4)){
  for (j in 1:nrow(summaryfile_psm_for_testing_pcc4)){  
    summaryfile_psm_for_testing_pcc4$index[j] -> target_scan
    summaryfile_psm_for_testing_pcc4 |> 
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
    spectrum_file_plot %>% select(peak_table_forgraph) %>%
      unnest(cols=peak_table_forgraph) %>% 
      ggplot( 
        aes(x=`m/z`,y=rel_int, label = iontype)
      )+
      geom_bar(stat="identity",width = 0.2)+
      geom_col(aes(fill = iontype),width = spec_plot_xmax/1000)+
      geom_text(aes(y=if_else(rel_int<0,rel_int-0.08,rel_int),label=name,color= iontype),angle=90,hjust=-0.05,vjust=0,size=2, na.rm = TRUE)+
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
        values = c("b"=RColorBrewer::brewer.pal(n=9,"Set1")[1],"y"=RColorBrewer::brewer.pal(n=9,"Set1")[2]),
        na.value = NA
      )+
      scale_color_manual(
        values = c("b"=RColorBrewer::brewer.pal(n=9,"Set1")[1],"y"=RColorBrewer::brewer.pal(n=9,"Set1")[2]),
        na.value = NA
      )+
      scale_x_continuous(limits=c(0, spec_plot_xmax),expand=c(0,0))+ #limits=c(250,800),
      scale_y_continuous(breaks = c(1,0.5,0,-0.5,-1),labels = scales::percent(c(1,0.5,0,0.5,1)), limits=c(-1,1),expand=c(0,0))+
      theme_classic()+
      labs(x="m/z",y="Normalized Intensity", fill = "Ion Type", color = "Ion Type", title = paste0(title1,"  ",title3), subtitle = paste0("PCC = ",title2,"  index: ",subtitle1))+
      theme(
        legend.position = "none",
        plot.subtitle = element_text(size=10, color="black")
      )
    ggsave(paste0("./ms2figures/",j,".png"), width=20, height = 10, units = "cm", dpi = 600)
  }
  
  ###figure generation
  
  
  for (j in 1:nrow(summaryfile_psm_for_testing_pcc5)){
    summaryfile_psm_for_testing_pcc5$index[j] -> target_scan
    summaryfile_psm_for_testing_pcc5 |> 
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
    spectrum_file_plot %>% select(peak_table_forgraph) %>%
      unnest(cols=peak_table_forgraph) %>% 
      ggplot( 
        aes(x=`m/z`,y=rel_int, label = iontype)
      )+
      geom_bar(stat="identity",width = 0.2)+
      geom_col(aes(fill = iontype),width = spec_plot_xmax/1000)+
      geom_text(aes(y=if_else(rel_int<0,rel_int-0.08,rel_int),label=name,color= iontype),angle=90,hjust=-0.05,vjust=0,size=2, na.rm = TRUE)+
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
        values = c("b"=RColorBrewer::brewer.pal(n=9,"Set1")[1],"y"=RColorBrewer::brewer.pal(n=9,"Set1")[2]),
        na.value = NA
      )+
      scale_color_manual(
        values = c("b"=RColorBrewer::brewer.pal(n=9,"Set1")[1],"y"=RColorBrewer::brewer.pal(n=9,"Set1")[2]),
        na.value = NA
      )+
      scale_x_continuous(limits=c(0, spec_plot_xmax),expand=c(0,0))+ #limits=c(250,800),
      scale_y_continuous(breaks = c(1,0.5,0,-0.5,-1),labels = scales::percent(c(1,0.5,0,0.5,1)), limits=c(-1.2,1.2),expand=c(0,0))+
      theme_classic()+
      labs(x="m/z",y="Normalized Intensity", fill = "Ion Type", color = "Ion Type", title = paste0(title1,"  ",title3), subtitle = paste0("PCC = ",title2,"  index: ",subtitle1))+
      theme(
        legend.position = "none",
        plot.subtitle = element_text(size=10, color="black")
      )
    ggsave(paste0("./ms2figures/high/",j,".png"), width=20, height = 10, units = "cm", dpi = 600)
  }
  
  
  for (j in 1:nrow(summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered_trinity)){
    summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered_trinity$index[j] -> target_scan
    summaryfile_psm_for_testing3_pcc_rt_hfsm2_basic_filtered_trinity |> 
      dplyr::filter(
        index == target_scan
      ) -> spectrum_file_plot
    
    spectrum_file_plot$sequence_for_fragspec -> title1
    spectrum_file_plot$index -> subtitle1 ##index!
    round(spectrum_file_plot$PCC,3) -> title2 ##pcc
    spectrum_file_plot$`MHplus in Da` -> spec_plot_xmax
    spectrum_file_plot$gene_position -> title3
    
    spectrum_file_plot %>% select(peak_table_for_graph) %>%
      unnest(cols=peak_table_for_graph) %>% filter(
        !is.na(matches)
      ) %>% 
      arrange(desc(frag)) %>% 
      slice_head(n=1) %>% 
      pull(frag)->graph_frag_max
    
    spectrum_file_plot %>% select(peak_table_for_graph) %>%
      unnest(cols=peak_table_for_graph) %>% filter(
        !is.na(matches)
      ) %>% 
      arrange(desc(`m/z`)) %>% 
      slice_head(n=1) %>% 
      pull(`m/z`) -> spec_plot_xmax
    #summary_filtered_PRM_table_final_DL_library_merged_spec_DLspec_crossed$pearson_correl_coeff[j] |> round(digits = 4) -> subtitle2
    #summary_filtered_PRM_table_final_DL_library_merged_spec_DLspec_crossed$cosine_similarity [j] |> round(digits = 4) -> subtitle3
    spectrum_file_plot %>% select(peak_table_for_graph) %>%
      unnest(cols=peak_table_for_graph) %>% 
      mutate(
        frag = if_else(frag>=0,frag/graph_frag_max,frag)
      ) %>% 
      ggplot( 
        aes(x=`m/z`,y=frag, label = series)
      )+
      geom_bar(stat="identity",width = 0.2)+
      geom_col(aes(fill = series),width = spec_plot_xmax/1000)+
      geom_text(aes(y=if_else(frag<0,frag-0.28,frag),label=ion_type,color= series),angle=90,hjust=-0.05,vjust=0.5,size=1.1, na.rm = TRUE)+
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
        values = c("b"=RColorBrewer::brewer.pal(n=9,"Set1")[1],"y"=RColorBrewer::brewer.pal(n=9,"Set1")[2]),
        na.value = NA
      )+
      scale_color_manual(
        values = c("b"=RColorBrewer::brewer.pal(n=9,"Set1")[1],"y"=RColorBrewer::brewer.pal(n=9,"Set1")[2]),
        na.value = NA
      )+
      scale_x_continuous(limits=c(0, spec_plot_xmax+100),expand=c(0,0))+ #limits=c(250,800),
      scale_y_continuous(breaks = c(1,0.5,0,-0.5,-1),labels = scales::percent(c(1,0.5,0,0.5,1)), limits=c(-1.4,1.4),expand=c(0,0))+
      theme_classic()+
      labs(x="m/z",y="Normalized Intensity", fill = "Ion Type", color = "Ion Type", title = paste0(title1,"  PCC=",title2))+
      theme(
        legend.position = "none",
        plot.subtitle = element_text(size=8, color="black"),
        axis.text.x = element_text(color = "black",size= 7),
        axis.text.y =element_text(color = "black",size= 7),
        axis.title = element_text(color = "black",size= 8),
        plot.title = element_text(color = "black",size= 7),
        axis.ticks = element_line(color = "black"),
        legend.title = element_text(color = "black",size= 8),
        legend.text = element_text(color = "black",size= 8),
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm")
      )
    ggsave(paste0("./ms2figures/241023/",title3,"_",j,".png"), width=9, height = 4.5, units = "cm", dpi = 600)
  }
  ggsave(paste0("./ms2figures/241023/",j,".png"), width=9, height = 4.5, units = "cm", dpi = 600)
}

###venn diagram for arginylome validation matrix 
##venn_plot1 @arginylome_validation_matrix2
{
  arginylome_validation_matrix2[,1:3] %>% 
    dplyr::rename(
      PCC = pcc_check,
      RT = rt_check
      #DI = diagnostic_ion_check
    ) %>% 
    column_to_rownames(var = "prot_pos_master") %>% 
    eulerr::euler() ->venn_plot1
  png(file = "./figure2/venn2.png", width = 15, height = 12, units = "cm", res = 600)
  plot(
    venn_plot1,
    quantities = TRUE,
    fills = c(wes_palette("Royal1")[1],wes_palette("Royal1")[4]),
    legend = list(side = "right")
  )
  dev.off()
}

##Experiment @summaryfile_psm_for_testing_pcc_rt2 ##site_count_experiment2.png
{
  summaryfile_psm_for_testing_pcc_rt2 %>% colnames
  
  summaryfile_psm_for_testing_pcc_rt2 %>% 
    filter(pcc_check == 1 & rt_check == 1 & p1r_check == 0 & p2gv_check == 0) %>% 
    distinct(
      prot_pos_master,experiment,.keep_all = TRUE
    ) %>% 
    group_by(experiment) %>% 
    summarize(
      count = n()
    ) %>% 
    mutate(
      quality = "MS2-RT filtered"
    ) -> table_experiment_filtered
  
  summaryfile_psm_for_testing_pcc_rt2 %>% 
    #filter(dl_score == 3 & p1r_check == 0 & p2gv_check == 0) %>% 
    distinct(
      prot_pos_master,experiment,.keep_all = TRUE
    ) %>% 
    group_by(experiment) %>% 
    summarize(
      count = n()
    ) %>% 
    mutate(
      quality = "All Arg. Sites"
    )-> table_experiment_all
  bind_rows(
    table_experiment_filtered,
    table_experiment_all
  ) %>% 
    ggplot(
      aes(
        x=fct(experiment,levels =c("MOCK","MG132","MGTG")),
        fill = fct(quality,levels = c("All Arg. Sites","MS2-RT filtered")),
        y=count
      )
    )+
    geom_bar(stat="identity",position = "dodge",width = 0.5)+
    scale_fill_manual(values = c(wes_palette("Royal1")[1],wes_palette("Royal1")[2]))+
    labs(
      x = "Experiment", y = "Arginylation Site Count",
      fill = "Protease"
    )+
    theme_classic()+
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      text= element_text(color = "black", size = 15),
      axis.text = element_text(color = "black", size = 15),
      axis.title.x = element_blank()
    )
  ggsave("./figure2/site_count_experiment2.png",width=10,height=10,dpi = 600,units="cm")
}

##Localization @summaryfile_psm_for_testing_pcc_rt2
{
  summaryfile_psm_for_testing_pcc_rt2 %>%
    filter(pcc_check == 1 & rt_check == 1 & p1r_check == 0 & p2gv_check == 0) %>% 
    distinct(
      prot_pos_master,.keep_all = TRUE
    ) %>% 
    separate_rows(
      Localizations, sep = "\\|"
    ) %>% 
    summarize(
      count = n(),
      .by = Localizations
    ) %>% 
    mutate(
      percent = count/sum(count)
    ) %>% 
    mutate(
      quality = "MS2-RT filtered"
    )  %>% 
    arrange(desc(count)) -> table_localization_passed
  
  summaryfile_psm_for_testing_pcc_rt2 %>% 
    distinct(
      prot_pos_master,.keep_all = TRUE
    ) %>% 
    separate_rows(
      Localizations, sep = "\\|"
    ) %>% 
    summarize(
      count = n(),
      .by = Localizations
    ) %>% 
    mutate(
      percent = count/sum(count)
    )%>% 
    mutate(
      quality = "All Arg. Sites"
    ) %>% 
    arrange(desc(count))-> table_localization_all
  
  table_localization_all$Localizations
  bind_rows(
    table_localization_passed,
    table_localization_all
  ) %>% 
    ggplot(
      aes(
        x=fct(Localizations,levels=table_localization_all$Localizations),
        fill = fct(quality,levels = c("All Arg. Sites","MS2-RT filtered")),
        y=percent
      )
    )+
    geom_bar(
      stat="identity",position = "dodge",width = 0.5
    )+
    coord_flip()+
    scale_y_continuous(labels = scales::percent)+
    scale_fill_manual(values = c(wes_palette("Royal1")[1],wes_palette("Royal1")[2]))+
    labs(
      x = "Experiment", y = "Arginylation Site Count",
      fill = "Protease"
    )+
    theme_classic()+
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      text= element_text(color = "black", size = 10),
      axis.text.y = element_text(color = "black", size = 10, hjust=1),
      axis.text.x = element_text(color = "black", size = 10),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    )
  ggsave("./figure2/site_count_local2.png",width=10,height=10,dpi = 600,units="cm")
}

####box plot pcc @summaryfile_psm_for_testing_pcc6 ##pcc_ions
{
  summaryfile_psm_for_testing_pcc6 %>% 
    select(
      prot_pos_reposition,PCC,pcc_b,pcc_y
    ) %>% 
    pivot_longer(cols=c(pcc_b,pcc_y)) %>% 
    mutate(
      group = case_when(
        PCC >= 0.9 ~ "PCC > 0.90",
        0.9 > PCC & PCC >= 0.75 ~ "PCC > 0.75",
        0.75 > PCC & PCC >= 0.5 ~ "PCC > 0.5",
        .default = "PCC < 0.5"
      )
    ) %>% 
    summarize(
      count = n(),
      .by=group
    )
  
  summaryfile_psm_for_testing_pcc6 %>% 
    select(
      prot_pos_reposition,PCC,pcc_b,pcc_y
    ) %>% 
    pivot_longer(cols=c(pcc_b,pcc_y)) %>% 
    mutate(
      group = case_when(
        PCC >= 0.9 ~ "PCC > 0.90",
        0.9 > PCC & PCC >= 0.75 ~ "PCC > 0.75",
        0.75 > PCC & PCC >= 0.5 ~ "PCC > 0.5",
        .default = "PCC < 0.5"
      )
    ) %>% 
    ggplot(
      aes(y=value,x=fct(group,levels=c("PCC < 0.5","PCC > 0.5","PCC > 0.75","PCC > 0.90") ),fill=name)
    )+
    geom_boxplot()+
    annotate(
      "text",
      x=1:4,
      y = 1.08,
      label = paste("n = ",c(294,560,602,978))
    )+
    scale_fill_manual(labels = c("b-ion","y-ion"),values = c(wes_palette("Zissou1")[5],wes_palette("Zissou1")[1]))+
    labs(y="PCC Value")+
    theme_classic()+
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      text= element_text(color = "black", size = 10),
      axis.text.y = element_text(color = "black", size = 10, hjust=1),
      axis.text.x = element_text(color = "black", size = 10),
      axis.title.x = element_blank()
    )
  ggsave("./figure2/pcc_ions.png",width=10,height=10,dpi = 600,units="cm")
}

####SP mTP protease @arginylome_validation_matrix3  #site_perc_sp2.png
{
  arginylome_validation_matrix3 %>% 
    select(
      site_type,score_sum
    ) %>% 
    mutate(
      score_sum=fct(as.character(score_sum)),
      site_type=fct(site_type,levels= c("SP","mTP","Protease",NA))
    ) %>% 
    group_by(
      site_type,score_sum
    ) %>% 
    summarise(
      count = n()
    ) %>% 
    mutate(
      perc = count/sum(count),
      count2 = sum(count)
    ) %>% 
    ggplot(
      aes(x= site_type,y=perc,fill = score_sum)
    )+
    geom_bar(
      stat="identity"
    )+
    geom_text(
      aes(label = percent(perc,accuracy = 0.1)),
      position = position_fill(vjust = 0.5),
      size = 2
    )+
    geom_text(aes(label = paste0("n = ",count2),y=1.02),size = 2)+
    scale_fill_manual(values = rev(c(wes_palette("Royal1")[1],wes_palette("Royal2")[2:3])))+
    scale_x_discrete(labels = c("SP","mTP","Protease","None"))+
    scale_y_continuous(
      labels = scales::percent
    )+
    labs(fill = "Validation Counts")+
    theme_classic()+
    theme(
      legend.position = "top",
      #legend.title = element_blank(),
      #text= element_text(color = "black", size = 8),
      axis.text.y = element_text(color = "black", size = 10, hjust=1),
      axis.text.x = element_text(color = "black", size = 10),
      legend.text = element_text(color = "black", size = 10),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    )
  ggsave("./figure2/site_perc_sp2.png",width=10,height=10,dpi = 600,units="cm")
}

###PCC90 counts @summary_ms2prediction
{
  summary_ms2prediction %>% 
    filter(...1 == ">0.90") %>% 
    mutate(type = fct(type,levels = c("built","pretrained","transfer","arg"))) %>% 
    ggplot(
      aes(x=type,y=PCC,fill= type)
    )+
    geom_bar(
      color="black",
      stat="identity")+
    geom_text(
      aes(label = percent(PCC,accuracy = 0.1)),
      vjust = -1
    )+
    labs(
      y="PCC >= 0.9"
    )+
    scale_y_continuous(labels = scales::percent,limits=c(0,1))+
    scale_fill_manual(values = wes_palette("GrandBudapest2")[1:4])+
    scale_x_discrete(
      labels = c("Built","Pretrained","Transferred","Arginylome")
    )+
    theme_classic()+
    theme(
      legend.position = "none",
      #legend.title = element_blank(),
      #text= element_text(color = "black", size = 8),
      axis.text.y = element_text(color = "black", size = 10, hjust=1),
      axis.text.x = element_text(color = "black", size = 10),
      legend.text = element_text(color = "black", size = 10),
      axis.title.x = element_blank(),
      #axis.title.y = element_blank()
    )
  ggsave("./figure2/prediction_result.png",width=10,height=10,dpi = 600,units = "cm")
}

###rt box plot
{
  summaryfile_psm_for_testing_pcc_rt2 %>% 
    select(
      PCC,rt_norm,replicate,experiment
    ) %>% 
    mutate(
      PCC90 = if_else(PCC>= 0.9, "PCC >= 0.9","PCC < 0.9"),
      experiment = fct(experiment,levels = c("MOCK","MG132","MGTG"))
    ) %>% 
    ggplot(
      aes(x=fct(PCC90,levels=c("PCC >= 0.9","PCC < 0.9")),y=rt_norm)
    )+
    geom_boxplot(
      aes(fill = fct(PCC90,levels=c("PCC >= 0.9","PCC < 0.9")))
    )+
    facet_grid(~replicate+experiment)+
    labs(y="Observed Normalized RT", fill = "PCC Value")+
    scale_fill_manual(values = c(wes_palette("Royal1")[2],wes_palette("Royal1")[1]))+
    theme_classic()+
    theme(
      legend.position = "top",
      #legend.title = element_blank(),
      #text= element_text(color = "black", size = 8),
      axis.text.y = element_text(color = "black", size = 10, hjust=1),
      axis.text.x = element_text(color = "black", size = 10,angle = 45,hjust=1),
      legend.text = element_text(color = "black", size = 10),
      axis.title.x = element_blank(),
      #axis.title.y = element_blank()
    )
  ggsave("./figure2/rt_norm_boxplot.png",width=10,height=10,dpi = 600,units="cm")
}

###diag ion model b1+ int @file_index3 >> file_index_temp ##deprecated ##di_box.png
{ 
  file_index_temp %>%
    select(
      diag_ion_rel_int
    ) %>% 
    mutate(
      type="training"
    ) %>% 
    bind_rows(
      summaryfile_psm_for_testing_pcc4 %>% 
        select(
          PCC,diagnostic_ion_int_relative
        ) %>% 
        mutate(
          diagnostic_ion_int_relative = if_else(!is.na(diagnostic_ion_int_relative),diagnostic_ion_int_relative,0)
        ) %>% 
        mutate(
          pccgroup = if_else(PCC>=0.9,"1","0")
        ) %>% 
        select(diagnostic_ion_int_relative ,pccgroup) %>% 
        rename(
          diag_ion_rel_int = diagnostic_ion_int_relative,
          type=pccgroup
        )
    ) %>% 
    summarize(
      count = n(),
      .by = type
    )
  
  file_index3 %>% 
    mutate(
      diag_ion_rel_int = map(
        peak_table_alpha_ms2,
        \(x) filter(x, alpha_col == "b_z1" & ion_num2 == 1)
      )
    ) %>% 
    mutate(
      diag_ion_rel_int = map(
        diag_ion_rel_int,
        \(x) x$i_relative
      )
    ) %>% 
    mutate(
      diag_ion_rel_int = map_dbl(
        diag_ion_rel_int,
        \(x) ifelse(length(x)>0,x,NA)
      )
    ) %>% 
    mutate(
      diag_ion_rel_int = if_else(!is.na(diag_ion_rel_int),diag_ion_rel_int,0)
    ) -> file_index_temp
  
  file_index_temp %>% 
    mutate(
      peak_table = map(
        peak_table,
        \(x) mutate(x,relative_intensity = i/max(i))
      )
    ) -> file_index_temp
  file_index_temp %>% 
    mutate(
      diag_ion_rel_int2 = map(
        peak_table,
        \(x) filter(x, `m/z` > 202.138331-0.005 & `m/z` < 202.138331+0.005) ##decision point 0.001 7427 out 0.01>> 2403 out 0.05>> 1903 out 0.005 >>2428
      )
    ) %>% 
    mutate(
      diag_ion_rel_int2 = map(
        diag_ion_rel_int2,
        \(x) arrange(x,desc(relative_intensity))
      )
    ) %>% 
    mutate(
      diag_ion_rel_int2 = map(
        diag_ion_rel_int2,
        \(x) slice_head(x,n=1)
      )
    ) %>%
    mutate(
      diag_ion_rel_int2 = map(
        diag_ion_rel_int2,
        \(x) pull(x,var = relative_intensity)
      )
    ) %>% 
    mutate(
      diag_ion_rel_int2 = map_dbl(
        diag_ion_rel_int2,
        \(x) ifelse(length(x)>0,x,NA)
      )
    ) -> file_index_temp2
    
  
  
  file_index_temp2 %>%
    select(
      diag_ion_rel_int2
    ) %>% 
    mutate(
      type="training"
    ) %>% 
    bind_rows(
      summaryfile_psm_for_testing_pcc4 %>% 
        select(
          PCC,diagnostic_ion_int_relative
        ) %>% 
        mutate(
          diagnostic_ion_int_relative = if_else(!is.na(diagnostic_ion_int_relative),diagnostic_ion_int_relative,0)
        ) %>% 
        mutate(
          pccgroup = if_else(PCC>=0.9,"1","0")
        ) %>% 
        select(diagnostic_ion_int_relative ,pccgroup) %>% 
        rename(
          diag_ion_rel_int2 = diagnostic_ion_int_relative,
          type=pccgroup
        )
    ) %>% 
    ggplot(
      aes(x=type,y=diag_ion_rel_int2,fill=type)
    )+
    geom_boxplot()+
    annotate(
      "text",
      x = 1:3,
      y = 1.08,
      label = paste("n = ",c(728,489,20907))
    )+
    scale_y_continuous(
      labels = percent,limits = c(0,1.08),
      breaks = seq(0,1,0.2)
    )+
    scale_x_discrete(
      labels = c("PCC < 0.9","PCC >= 0.9","Trained")
    )+
    scale_fill_manual(
      values = wes_palette("Royal1")[1:3]
    )+
    stat_compare_means(
      comparisons = list(c("0","1")), label.y = 0.8
    )+
    labs(y="Relative Intensity of Diagnostic Ion")+
    theme_classic()+
    theme(
      legend.position = "none",
      legend.title = element_blank(),
      text= element_text(color = "black", size = 10),
      axis.text.y = element_text(color = "black", size = 10, hjust=1),
      axis.text.x = element_text(color = "black", size = 10),
      axis.title.x = element_blank()
    )
  ggsave("./figure2/di_box.png",width=10,height=10,dpi = 600,units="cm")
  
}

#DI cutoff decision ##deprecated
{
  file_index_temp2 %>% 
    filter(!is.na(diag_ion_rel_int2)) %>% 
    arrange(diag_ion_rel_int2) %>% 
    slice_head(prop = 0.01) %>% 
    slice_tail(n=1) %>%
    pull(diag_ion_rel_int2) -> cutoff_diagion_rel_int
    
  file_index_temp2 %>% 
    filter(!is.na(diag_ion_rel_int2)) %>%
    ggplot(
      aes(x=log(diag_ion_rel_int2,10))
    )+
    stat_ecdf(geom = "point")+
    geom_vline(xintercept = log(cutoff_diagion_rel_int,10),color="red",linetype = "dashed")+
    labs(
      x=expression(log[10]~relative~intensity~of~diagnostic~ion),
      y="Density"
      )+
    theme_classic()+
    theme(
      legend.position = "none",
      legend.title = element_blank(),
      text= element_text(color = "black", size = 10),
      axis.text.y = element_text(color = "black", size = 10, hjust=1),
      axis.text.x = element_text(color = "black", size = 10),
      axis.title = element_text(color = "black", size = 10)
    )
  ggsave("./figure2/di_cutoff.png",width=10,height=10,dpi = 600,units="cm")
}

#pie chart HPA
{
  table_hpa_proteinclass_count_total_arg %>% 
    arrange(desc(count.y)) %>% 
    slice_head(n=9) %>% 
    ggplot(
      aes(x="",y=count.y,fill=reorder(name,count.y))
    )+
    geom_bar(stat="identity")+
    coord_polar(theta = "y")+
    scale_fill_brewer(palette="Set1",direction = -1)
}

###pcc rt on vg gv
{
  summaryfile_psm_for_testing_pcc_rt2 %>% 
    mutate(
      rt_cut=upr-fit
    ) %>% pull(rt_cut)
  summaryfile_psm_for_testing_pcc_rt2 %>% 
    filter(
      p2gv_check==1
    ) %>% 
    mutate(
      rt_deviation = rt_norm-rt_pred
    ) %>% 
    ggplot(
      aes(x= p2)
    )+
    geom_boxplot(
      aes(y=rt_deviation),
      # stat="identity",
      # position = "dodge2",
      # width = 0.1
    )+
    geom_hline(yintercept = 0.03, color = "red", linetype = 2)+
    labs(y="Deviation of Normalized RT ",x="")+
    theme_classic()+
    theme(
      legend.position = "none",
      legend.title = element_blank(),
      text= element_text(color = "black", size = 10),
      axis.text.y = element_text(color = "black", size = 10, hjust=1),
      axis.text.x = element_text(color = "black", size = 10,angle = 90,vjust=0.5),
      axis.title.y = element_text(color = "black", size = 10)
    )
  ggsave("./figure2/validation_gv_rt.png",width=10,height=10,dpi = 600,units="cm")
  
  summaryfile_psm_for_testing_pcc_rt2 %>% 
    filter(
      p2gv_check==1
    ) %>% 
    ggplot(
      aes(x= p2)
    )+
    geom_boxplot(
      aes(y=PCC),
      # stat="identity",
      # position = "dodge2",
      # width = 0.1
    )+
    geom_hline(yintercept = 0.767, color = "red", linetype = 2)+
    labs(y="PCC",x="")+
    theme_classic()+
    theme(
      legend.position = "none",
      legend.title = element_blank(),
      text= element_text(color = "black", size = 10),
      axis.text.y = element_text(color = "black", size = 10, hjust=1),
      axis.text.x = element_text(color = "black", size = 10,angle = 90,vjust=0.5),
      axis.title.y = element_text(color = "black", size = 10)
    )
  ggsave("./figure2/validation_gv_pcc.png",width=10,height=10,dpi = 600,units="cm")
  
  summaryfile_psm_for_testing_pcc_rt2 %>%
    filter(
      p2gv_check==1
    ) %>%summarize(.by = p2, count = n())
}

##chord diagram
{
  
}


