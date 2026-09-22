# Helper functions for manuscript figure and table generation






#### 1. General plotting theme ####
.shared_theme <- function() {
  list(
    theme_classic(),
    theme(
      text                  = element_text(size = 10),
      axis.text.x           = element_text(angle = 45, hjust = 1),
      legend.title          = element_blank(),
      legend.margin         = margin(6, 6, 6, 6),
      legend.background     = element_rect(fill = NA, color = NA),
      strip.background      = element_blank(),
      panel.background = element_rect(fill = "transparent", color = NA),
      plot.background  = element_rect(fill = "transparent", color = NA),
      legend.position       = "none"
    )
  )
}

#### 2. Alpha diversity plotting ####
alpha_plot <- function(
    table,
    x_,
    y_ = "Shannon",
    colx_ = x_,
    ylim = 7,
    y_kruskal = 6.5,
    col = c("#E31A1C", "#1F78B4", "#4D4D4D")
) {
  
  comp_ <- as.vector(unique(table[, x_]))
  Num   <- length(comp_)
  
  base_plot <- ggplot(table, aes(x = .data[[x_]], y = .data[[y_]])) +
    geom_violin(aes(fill = .data[[x_]], color = .data[[x_]]),
                show.legend = FALSE, alpha = 0.6) +
    geom_jitter(aes(color = .data[[x_]]),
                width = 0.15, alpha = 0.3, size = 2) +
    geom_boxplot(aes(color = .data[[x_]]),
                 show.legend = FALSE, width = 0.2) +
    scale_y_continuous(limits = c(0, ylim)) +
    scale_fill_manual(values = col) +
    scale_color_manual(values = col) +
    ylab(y_) +
    theme_classic() + 
    .shared_theme()
  
  if (Num > 3) {
    
    pp <- base_plot +
      stat_compare_means(method = "kruskal.test",  tip.length=0.00 , # 0.05
                         size = 3,
                         label.y = y_kruskal, label.x =((Num+1)/2-0.4),
                         aes(label = paste0("kruskal.test, p = ", ..p.format..))
      )+
      geom_pwc(
        method     = "wilcox_test",
        label      = "p = {ifelse(p.adj < 0.001, scales::scientific(p.adj, digits = 2), p.adj)}{p.adj.signif}",
        label.size = 3,
        p.adjust.method = "fdr",
        tip.length = 0.00
      )
    
  } else if (Num == 3) {
    
    pp <- base_plot +
      stat_compare_means(
        method     = "kruskal.test",
        tip.length = 0.01,
        size       = 3,
        label.y    = y_kruskal,
        label.x    = ((Num + 1) / 2 - 0.4),
        aes(label  = paste0("kruskal.test, p = ", ..p.format..))
      ) +
      geom_pwc(
        method     = "wilcox_test",
        label      = "p = {ifelse(p.adj < 0.001, scales::scientific(p.adj, digits = 2), p.adj)}{p.adj.signif}",
        label.size = 3,
        p.adjust.method = "fdr",
        tip.length = 0.00
      )
    
  } else if (Num == 2) {
    
    pp <- base_plot +
      geom_pwc(
        method     = "wilcox_test",
        label      = "p = {ifelse(p.adj < 0.001, scales::scientific(p.adj, digits = 2), p.adj)}{p.adj.signif}",
        tip.length = 0.00,
        label.size = 3,
        p.adjust.method = "fdr",
      )
  }
  
  return(base_plot)
}



#### 3. Beta diversity plotting and statistics ####
beta_plot <- function(phyloseq, type, shap = NULL, seed = 42, plot = "PCoA",
                      SampleID = "SampleID", type_col, col_inout       = c("out"), col_right_left  = c("right"), indices = c("bray", "jaccard", "unifrac", "wunifrac")) {
  # Round p-values
  value <- function(val) {
    if (val > 0.001)
      round(val, 3)
    else
      "<0.001"
  }
  
  out_name <- function(v1)
    deparse(substitute(v1))
  
  plots       <- list()
  result_list <- list()
  
  for (index in indices) {
    set.seed(seed)
    x.dist <- phyloseq::distance(phyloseq, method = index)
    dist   <- out_name(x.dist)
    meta   <- phyloseq %>% sample_data() %>% data.frame()
    
    # PERMANOVA
    set.seed(seed)
    Perm    <- adonis2(
      as.formula(glue("{dist} ~ {type}")),
      data         = data.frame(sample_data(phyloseq)),
      permutations = 9999,
      method = index,
      strata = meta$SubjectID
    ) 
    Perm.p  <- value(Perm$`Pr(>F)`[1])
    Perm.R2 <- round(Perm$R2[1], 3)
    
    result <- data.frame(
      Statistical.test = "PERMANOVA",
      DF               = Perm$Df[1],
      Sum.Sq           = Perm$SumOfSqs[1],
      Mean.Sq          = NA,
      F.Statist        = Perm$F[1],
      R.Squared        = Perm.R2,
      P.value          = Perm$`Pr(>F)`[1]
    )
    result_list[[index]] <- result
    
    # Ordinate
    set.seed(seed)
    ord    <- ordinate(phyloseq, plot, index)
    pcoa_df <- data.frame(meta, ord$vectors[, 1:2])
    PC1 <- round(ord$values["Relative_eig"][1, ] * 100, 1)
    PC2 <- round(ord$values["Relative_eig"][2, ] * 100, 1)
    
    Title <- switch(
      index,
      bray     = "Bray-curtis",
      jaccard  = "Jaccard",
      unifrac  = "Unweighted UniFrac",
      wunifrac = "Weighted UniFrac"
    )
    
    main.plot <- pcoa_df %>%
      ggplot(aes(x = Axis.1, y = Axis.2)) +
      geom_vline(xintercept = 0, colour = "grey80") +
      geom_hline(yintercept = 0, colour = "grey80") +
      geom_point(aes_string(shape = shap, color = type),
                 alpha = 0.5,
                 size = 2) +
      stat_ellipse(aes_string(color = type)) +
      scale_color_manual(values = type_col) +
      labs(x = paste0("PCoA1 (", PC1, "%)"),
           y = paste0("PCoA2 (", PC2, "%)")) +
      annotate(
        "text",
        hjust = 0,
        vjust = -0.2,
        x = -Inf,
        y = -Inf,
        size = 3,
        # size 10pt ≈ 3 in ggplot units
        label = paste0(
          # Title,
          "\n PERMANOVA: ",
          "\n R2 = ",
          Perm.R2,
          ", p-value = ", Perm.p ) ) +
      .shared_theme() +
      theme(
        legend.position = "right",
        # beta plot은 legend 표시
        plot.caption    = element_markdown(),
        aspect.ratio    = 1,
        legend.background = element_rect(fill = scales::alpha("white", 0.5)),
        legend.key = element_rect(fill = "transparent")
      )
    
    # Legend position inside plot (optional)
    if (col_inout == "in" && col_right_left == "right") {
      main.plot <- main.plot + theme(
        legend.position      = c(1, 1),
        legend.justification = c("right", "top"),
        legend.box.just      = "right"
      )
    } else if (col_inout == "in" && col_right_left == "left") {
      main.plot <- main.plot + theme(
        legend.position      = c(0, 1),
        legend.justification = c("left", "top"),
        legend.box.just      = "left"
      )
    }
    
    assign(paste0(index, "_p"), main.plot)
    plots[[index]] <- main.plot
  }
  
  plot_list     <- lapply(indices, function(index)
    get(paste0(index, "_p")))
  plots[["Total"]] <- ggarrange(plotlist = plot_list,
                                ncol = 2,
                                nrow = 2)
  
  return(list(plots = plots, results = result_list))
}


run_beta_stats <- function(physeq, index = "bray", group_var = "Site2", seed = 42) {
  set.seed(seed)
  dist_obj <- phyloseq::distance(physeq, method = index)
  meta <- as(sample_data(physeq), "data.frame")
  
  # PERMANOVA
  permanova <- adonis2(as.formula(paste("dist_obj ~", group_var)),
                       data = meta, permutations = 9999, 
                       method = index, 
                       strata = meta[["SubjectID"]]) ############
  p_perm <- format_pval(permanova$`Pr(>F)`[1])
  r2_perm <- round(permanova$R2[1], 3)
  
  list(
    dist = dist_obj,
    permanova = list(p = p_perm, R2 = r2_perm)
    
  )
}

format_pval <- function(val) {
  if (val > 0.05) round(val, 3)
  else if (val > 0.001) round(val, 3)
  else "<0.001"
}

out_name <- function(v) {
  deparse(substitute(v))
}

extract_envfit_vectors <- function(ord_vectors, physeq, sig_level = 0.05) {
  set.seed(42)
  fit <- vegan::envfit(ord_vectors, otu_table(physeq), perm = 9999)
  arrows <- as.data.frame(fit$vectors$arrows * sqrt(fit$vectors$r))
  arrows$p.val <- fit$vectors$pvals
  arrows$r2 <- fit$vectors$r
  sig <- arrows[arrows$p.val < sig_level, , drop = FALSE]
  sig$p.adj <- p.adjust(sig$p.val, method = "BH")
  tax <- as.data.frame(tax_table(physeq))
  sig$Feature <- rownames(sig)
  merged <- merge(sig, tax[, "Species", drop = FALSE], by.x = "Feature", by.y = "row.names")
  tibble::column_to_rownames(merged, var = "Feature")
}


plot_pcoa <- function(physeq, sample_id = "SampleID", group_var = "Site2", 
                      color_var = "Skin.type", shape_var = "Site2", index = "bray",
                      levels_skin = c("Dry", "Moist", "Sebaceous"),
                      seed = 42) {
  stats <- run_beta_stats(physeq, index, group_var, seed)
  ord <- ordinate(physeq, method = "PCoA", distance = index)
  eig <- ord$values$Relative_eig
  PC1 <- round(eig[1] * 100, 1)
  PC2 <- round(eig[2] * 100, 1)
  
  mat <- as.data.frame(ord$vectors[, 1:2])
  mat[[sample_id]] <- rownames(mat)
  meta <- sample_data(physeq) |> data.frame()
  df <- merge(meta, mat, by = "row.names")
  
  df$Skin.type <- factor(df$Skin.type, levels = levels_skin)
  set.seed(seed)
  envfit_df <- extract_envfit_vectors(ord_vectors = ord$vectors[, 1:2], physeq = physeq)
  
  envfit_df2 <- envfit_df[envfit_df$p.adj <= 0.001 & envfit_df$r2 >=0.1,  ]
  p <- ggplot(df, aes(x = Axis.1, y = Axis.2)) +
    geom_vline(xintercept = 0, color = "grey80") +
    geom_hline(yintercept = 0, color = "grey80") +
    geom_point(aes_string(shape = shape_var, color = color_var),
               size = 2.5,
               alpha = 0.7) +
    scale_color_manual(values = c("#efba61", "#009999", "#E56666", "#8eb1c1")) +
    scale_shape_manual(values = c(17, 16, 17, 16, 17, 16, 3, 4, 17)) +
    labs(x = paste0("PCoA1 (", PC1, "%)"),
         y = paste0("PCoA2 (", PC2, "%)")) +
    theme_test() +
    theme(
      legend.title = element_blank(),
      aspect.ratio = 1,
      plot.caption = element_markdown(),
      legend.position = "right",
      plot.margin = unit(rep(0, 4), "points")
    ) +
    geom_segment(
      data = envfit_df2,
      aes(
        x = 0,
        xend = Axis.1 / 2,
        y = 0,
        yend = Axis.2 / 2
      ),
      arrow = arrow(length = unit(0.25, "cm")),
      color = "grey20"
    ) +
    ggrepel::geom_text_repel(data = envfit_df2,
                             aes(x = Axis.1 / 2, y = Axis.2 / 2, label = Species),
                             size = 3)
  
  
  
  out = list(
    plot = p, 
    envfit_res = envfit_df,
    permanova = stats
  )
  return(out)
}



dendrogram_group_centroid <- function(physeq, 
                                      group_var = "Site2",
                                      dist_method = "bray",
                                      hclust_method = "complete",
                                      seed = 42) {
  library(phyloseq)
  library(dplyr)
  library(vegan)
  library(ggdendro)
  library(ggplot2)
  
  set.seed(seed)
  
  # OTU/ASV 테이블 추출
  otu <- data.frame(otu_table(physeq))
  if (taxa_are_rows(physeq)) {
    otu <- t(otu)
  }
  
  # 메타데이터 병합
  meta <- sample_data(physeq) |> data.frame()
  if (!(group_var %in% colnames(meta))) {
    stop(paste0("'", group_var, "' not found in sample_data"))
  }
  otu$Group <- meta[[group_var]]
  group_levels <- unique(meta$Site2)
  
  # 그룹별 평균 산출
  group_mean <- otu %>%
    group_by(Group) %>%
    summarise(across(where(is.numeric), mean), .groups = "drop") %>%
    column_to_rownames("Group")
  
  # 거리 행렬 계산
  # dist_mat <- vegan::vegdist(group_mean, method = dist_method)
  if (dist_method %in% c("wunifrac", "unifrac")) {
    dist_mat <- phyloseq::distance(physeq, method = dist_method)
    dist_mat <- as.matrix(dist_mat)
    
    group_dist <- matrix(
      NA,
      nrow = length(group_levels),
      ncol = length(group_levels),
      dimnames = list(group_levels, group_levels)
    )
    
    for (i in group_levels) {
      for (j in group_levels) {
        samp_i <- rownames(meta[meta$Site2 == i, ])
        samp_j <- rownames(meta[meta$Site2 == j, ])
        
        group_dist[i, j] <- mean(dist_mat[samp_i, samp_j])
      }
    }
    dist_mat.out <- group_dist %>% as.dist()
    
  } else {
    dist_mat.out <- vegan::vegdist(group_mean, method = dist_method)
  }
  
  
  # 계층적 군집화 및 덴드로그램 데이터 변환
  hc <- hclust(dist_mat.out, method = hclust_method)
  dd <- dendro_data(hc)
  
  # 시각화
  p <- ggplot() +
    geom_segment(data = dd$segments,
                 aes(x = x, y = y, xend = xend, yend = yend)) +
    # geom_text(data = dd$labels,
    #           aes(x = x, y = y - 0.05 * max(dd$segments$y), label = label),
    #           angle = 90, hjust = 1, size = 3) +
    labs(# title = paste("Hierarchical Clustering (", dist_method, ")", sep = ""),
      x = NULL, 
      y = "Distance") +
    theme_minimal()
  
  
  return(out = list(plot = p, 
                    dendrogram = hc))
}









pairwise_adonis_to_df <- function(pwres) {
  # 이름이 "parent_call"인 건 제외
  result_names <- names(pwres)[-1]
  
  # 각 비교 결과를 데이터프레임으로 변환
  df_list <- lapply(result_names, function(name) {
    res <- as.data.frame(pwres[[name]])
    res$Comparison <- name
    res
  })
  
  # 전체 병합
  all_df <- do.call(rbind, df_list)
  
  # Comparison 열을 맨 앞으로 정렬
  all_df <- all_df[, c("Comparison", setdiff(colnames(all_df), "Comparison"))]
  
  # "Model" 행만 추출
  final_df <- all_df %>%
    filter(grepl("^Model", rownames(.))) %>%
    mutate(FDR = p.adjust(`Pr(>F)`, method = "fdr"))
  
  return(final_df)
}

#### 4. Taxonomic composition  ####


Abund_cal <- function(ps.glom, tax_level, group, path) {
  
  # melt
  melt <- psmelt(ps.glom)
  
  # setting
  Taxonomy <- as.vector(unique(melt[, tax_level]))
  meta <- data.frame(sample_data(ps.glom))
  meta_com <- meta[, group, drop = TRUE] %>% unique
  tax_num <- length(Taxonomy)
  
  # reset
  Total_result <- NULL
  
  ## Total abundance calculation
  Total_result <- melt %>%
    dplyr::group_by(!!rlang::sym(tax_level)) %>%
    dplyr::summarize(
      Total.Mean = mean(Abundance), # 평균
      Total.N    = n(),             # 행 개수
      Total.Sd   = sd(Abundance)    # 표준편차
    )
  
  ## Each group abundance calculation
  for (com in meta_com) {
    # 각 그룹별 계산
    melt.2 <- melt %>% filter(!!rlang::sym(group) == com)
    result <- melt.2 %>%
      dplyr::group_by(!!rlang::sym(tax_level)) %>%
      dplyr::summarize(
        !!paste0(com, ".Mean") := mean(Abundance),  # 평균
        !!paste0(com, ".N")    := n(),              # 행 개수
        !!paste0(com, ".Sd")  := sd(Abundance)      # 표준편차
      ) %>%   
      dplyr::mutate(!!paste0(com, ".se")      := !!rlang::sym(paste0(com, ".Sd"))    / sqrt(!!rlang::sym(paste0(com, ".N"))),           # 표준오차
                    !!paste0(com, ".lower")   := !!rlang::sym(paste0(com, ".Mean"))  - qnorm(0.975) * !!rlang::sym(paste0(com, ".se")), # 95% 신뢰 구간 하한
                    !!paste0(com, ".upper")   := !!rlang::sym(paste0(com, ".Mean"))  + qnorm(0.975) * !!rlang::sym(paste0(com, ".se")), # 95% 신뢰 구간 상한
                    !!paste0(com, ".CI95per") := !!rlang::sym(paste0(com, ".upper")) - !!rlang::sym(paste0(com, ".lower"))              # 95% 신뢰 구간
      )
    
    # abundance 결과 합치기
    Total_result <- bind_cols(Total_result, result[, -1])
  }
  return(Total_result)
}


taxa_plot <- function(melt, taxa, tax_otu, x_axis, phylum_or = NULL){ # 2024 12 11
  

  F1.process_data = function(melt , taxa, tax_otu) {
    
    tax_tab <-   melt[, c("OTU", "Phylum", taxa)] %>% unique
    # tax_tab
    tax_tab2 <- tax_tab[tax_tab$OTU %in% tax_otu, ]
    tax_phylum <- tax_tab2$Phylum %>% unique
    tax_index <- tax_tab2[,  taxa,  drop=T]
    
    # Reconstruct taxa classified as Others
    melt.2 <- melt # Back up
    melt.2[!melt.2[, "Phylum"] %in% tax_phylum, "Phylum"] <- "Other"
    melt.2[!melt.2[, "Phylum"] %in% tax_phylum, taxa] <- "Other"
    
    
    # Genus 와 Phylum정렬
    if (taxa != "Species") {
      for (i in tax_phylum) {
        G <-tax_tab2[tax_tab2[, "Phylum"] == i, taxa]
        
        melt.2[melt.2[, "Phylum"] == i, taxa]
        
        
        
        melt.2[melt.2[, "Phylum"] == i & !melt.2[, taxa] %in% G, taxa] <- paste0(i, "_Other")
      }
      for (i in tax_phylum) {
        G <- tax_tab2[tax_tab2[, "Phylum"] == i, taxa]
        for (g in G) {
          melt.2[melt.2[, taxa] == g, taxa] <- paste0(i, "_", g)
        }
      }
    } else {
      for (i in tax_phylum) {
        G <- tax_tab2[tax_tab2[, "Phylum"] == i, taxa]
        melt.2[melt.2[, "Phylum"] == i & !melt.2[, taxa] %in% G, taxa] <- "Other"
        melt.2[melt.2[, "Phylum"] == i & !melt.2[, taxa] %in% G, "Phylum"] <- "Other"
      }
    }
    return(list(df = melt.2,
                Phylum_list = tax_phylum))
    
  }
  # F2.Order_data
  F2.Order_data = function(processed_data, Top_p, taxa){
    # phylum level 
    table <- processed_data[processed_data[, "Phylum" ] %in% Top_p, c("Abundance", "Phylum", taxa)]
    p_order <- table %>%  .[,"Phylum" ]%>% unique
    
    processed_data[,"Phylum"]  <- factor(processed_data[,"Phylum" ],  levels = c(sort(p_order), "Other"))
    
    # Genus order
    table_2 <- table %>% 
      dplyr::group_by(Phylum, !!rlang::sym(taxa)) %>%
      dplyr::summarise(sum.Abundance=sum(Abundance), .groups = 'drop') %>%
      dplyr::arrange( -sum.Abundance) %>%
      ungroup() %>% 
      as.data.frame()
    
    g_order <- table_2 %>% 
      dplyr::arrange(Phylum) %>%
      select(all_of(taxa))  %>% .[[1]]
    
    processed_data[ ,taxa] <- factor(processed_data[,taxa], levels = c(g_order, "Other"))
    return(list(df = processed_data, summary_df = table_2))
    
  }
  # get_palette_colors
  get_palette_colors = function(palette_name, num_taxa) {
    if (num_taxa == 1) {
      colors <- rev(brewer.pal(9, palette_name)[5])
    } else if (num_taxa == 2) {
      colors <- rev(brewer.pal(9, palette_name)[c(3, 7)])
    } else if (num_taxa >= 3 & num_taxa <= 9) {
      colors <- rev(brewer.pal(num_taxa, palette_name))
    } else {
      color_list <- rev(brewer.pal(9, palette_name))
      colors <- colorRampPalette(color_list)(num_taxa)
    }
    return(colors)
  } 
  # F3.generate_colors
  F3.generate_colors = function(df, taxa) {
    
    ## arrange by abundance and phylum
    table_3 <- df %>%
      dplyr::arrange(-sum.Abundance) %>%
      dplyr::arrange(Phylum) %>%
      dplyr::select(Phylum, !!rlang::sym(taxa))
    
    ## count taxa  (2024.09.24)
    categories <- table_3 %>% 
      dplyr::group_by(Phylum) %>% 
      dplyr::summarise(Taxa = n())
    colnames(categories)[2] <- taxa
    # categories <- aggregate(as.formula(paste(taxa, "Phylum", sep = "~")), 
    #                         table_3, 
    #                         function(x) length(unique(x))) 
    
    ## First phylum 
    # P_levels <- table_3$Phylum %>% unique
    
    ## color list 
    color_list.names <- categories[, "Phylum"]
    color_list <- vector("list", length(color_list.names))
    names(color_list) <- color_list.names
    
    ## 
    phylum_color_map <- list(
      # Bacteria
      Actinobacteria = "Reds",
      Actinobacteriota = "Reds",
      Actinomycetota = "Reds",
      
      Firmicutes = "Blues",
      Bacillota = "Blues",
      Firmicutes_A = "Blues",
      Firmicutes_B = "BuGn",
      Firmicutes_C = "PuBu",
      Firmicutes_D = "YlGnBu", 
      
      Bacteroidetes = "Purples",
      Bacteroidota = "Purples",
      
      Proteobacteria = "Greens",
      Pseudomonadota = "Greens",
      
      Fusobacteria = "YlOrBr",
      Fusobacteriota = "YlOrBr",
      
      # Fungi 
      Ascomycota = "RdPu",
      Basidiomycota = "YlOrBr"
      
    )
    
    # GTDB taxonomy의 suffix가 붙은 Phylum도 동일 계열 색상을 사용한다.
    phylum_levels <- as.character(categories$Phylum)
    palette_names <- rep(NA_character_, length(phylum_levels))
    
    palette_names[grepl("^(Actinobacteria|Actinobacteriota|Actinomycetota)", phylum_levels)] <- "Reds"
    palette_names[grepl("^(Firmicutes|Bacillota)", phylum_levels)] <- "Blues"
    palette_names[grepl("^(Bacteroidetes|Bacteroidota)", phylum_levels)] <- "Purples"
    palette_names[grepl("^(Proteobacteria|Pseudomonadota)", phylum_levels)] <- "Greens"
    palette_names[grepl("^(Fusobacteria|Fusobacteriota)", phylum_levels)] <- "YlOrBr"
    palette_names[grepl("^Ascomycota", phylum_levels)] <- "RdPu"
    palette_names[grepl("^Basidiomycota", phylum_levels)] <- "YlOrBr"
    
    # 기존 함수에 명시된 세부 Phylum 색상은 우선 적용한다.
    exact_match <- phylum_levels %in% names(phylum_color_map)
    palette_names[exact_match] <- unname(unlist(
      phylum_color_map[phylum_levels[exact_match]]
    ))
    
    other_colors <- c("BrBG", "Spectral", "PuOr")
    palette_names[is.na(palette_names)] <- rep(
      other_colors,
      length.out = sum(is.na(palette_names))
    )
    
    phylum_color_map2 <- stats::setNames(
      as.list(palette_names),
      phylum_levels
    )
    
    for (phylum in names(phylum_color_map2)) {
      num_taxa <- categories[categories$Phylum == phylum, taxa]
      palette_name <- phylum_color_map2[[phylum]]
      
      if (length(num_taxa) > 0) {
        color_list[[phylum]] <- get_palette_colors(palette_name, as.numeric(num_taxa))
      }
    }
    
    
    color_vector <- unlist(color_list, use.names = FALSE)
    final_color <- c(color_vector, "#D3D3D3")
    return(final_color)
  }
  
  F5.taxa_plot <- function(df, color, taxa, x_axis){
    p <- ggplot(df, aes(x = !!rlang::sym(x_axis), y = Abundance, fill = !!rlang::sym(taxa))) +
      geom_bar(stat = "identity", position="fill") +
      labs(y = "Relative abundance") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),  # x축 라벨 각도 조정
            plot.title = element_text(hjust = 0.5),  # 제목 가운데 정렬
            legend.position = "right",  # 범례를 하단에 위치
            legend.title = element_blank()) +  # 범례 제목 제거
      scale_fill_manual(values = color)   # 색상 팔레트 변경
    
    return(p)
  }
  
  
  

  F6.sampleID_order <- function(dff, phylum_or, x_axis) {
    sample_order <- dff %>% 
      dplyr::group_by(!!rlang::sym(x_axis)) %>% 
      data.frame() %>% 
      dplyr::mutate(Abundance = Abundance / sum(Abundance)) %>%
      dplyr::filter(Phylum %in% phylum_or) %>% 
      dplyr::group_by(!!rlang::sym(x_axis)) %>% 
      dplyr::summarise(Abundance = sum(Abundance)) %>% 
      dplyr::arrange(Abundance) %>%
      pull(!!rlang::sym(x_axis)) %>% as.character()
    
    dff[, x_axis] <- factor(dff[, x_axis], levels = sample_order)
    
    return(dff)
  }
  
  
  
  
  out2 <- F1.process_data(melt = melt, taxa = taxa, tax_otu = tax_otu)
  
  out3 <- F2.Order_data(processed_data = out2$df,
                        Top_p =  out2$Phylum_list, 
                        taxa = taxa)
  
  color_code <- F3.generate_colors(df = out3$summary_df, 
                                   taxa=taxa )

    plot <- F5.taxa_plot(df = out3$df,
                       color = color_code, 
                       taxa = taxa, 
                       x_axis = x_axis)
  
  if (!is.null(phylum_or)){
    out4 <- out3
    out4$df <- F6.sampleID_order(df = out3$df, 
                                 phylum_or = phylum_or, 
                                 x_axis = x_axis)
    
    plot <- F5.taxa_plot(df = out4$df, 
                         color = color_code, 
                         taxa = taxa, 
                         x_axis = x_axis)
    
  } 
  return(out = list(plot = plot, 
                    data = out3,
                    color = color_code))
  
}









sampleID_order <- function(tax_data, ph) {
  sample_order <-  tax_data$data %>% 
    data.frame() %>%
    dplyr::group_by(SampleID) %>% 
    mutate(Abundance = Abundance / sum(Abundance)) %>%
    filter(Phylum == ph) %>% 
    dplyr::group_by(SampleID) %>% 
    dplyr::summarise(Abundance = sum(Abundance)) %>% 
    dplyr::arrange(Abundance) %>%
    pull(SampleID) %>% as.character()
  
  tax_data$data$SampleID <- factor(tax_data$data$SampleID, levels = sample_order)
  
  return(tax_data)
}



#### 5. Decontamination helper ####




get_abundance_prevalence <- function(ps_obj) {
  M <- otu_table(ps_obj) %>% t() %>%  data.frame
  M[M > 0] <- 1
  
  df <- data.frame(
    Abundance = taxa_sums(ps_obj),
    Prevalence = rowSums(M)
  )
  
  return(df)
}


#### 6. ANCOM-BC2 helper ####



make_barplots_with_p <- function(
    res,
    prefix,
    adj = "Yes",              # "Yes" = q, "No" = p
    color_mode = "direction",
    custom_colors = NULL,
    direction_labels = c("Up" = "Up", "Down" = "Down")) {
  
  
  df <- res %>% select(taxon, contains(prefix))
  
  diff_cols <- grep("^diff_", names(df), value = TRUE)
  lfc_cols  <- grep("^lfc_", names(df), value = TRUE)
  qval_cols <- grep("^q_", names(df), value = TRUE)
  pval_cols <- grep("^p_", names(df), value = TRUE)
  
  
  df$taxon <- str_replace_all(df$taxon, "Unclassified", "unclassified")
  df_lfc <- df %>%
    select(taxon, all_of(lfc_cols)) %>%
    pivot_longer(-taxon, names_to = "group", values_to = "value") %>%
    mutate(group = gsub("^lfc_", "", group))
  
  df_diff <- df %>%
    select(taxon, all_of(diff_cols)) %>%
    pivot_longer(-taxon, names_to = "group", values_to = "diff") %>%
    mutate(group = gsub("^diff_", "", group))
  
  df_q <- df %>%
    select(taxon, all_of(qval_cols)) %>%
    pivot_longer(-taxon, names_to = "group", values_to = "q") %>%
    mutate(group = gsub("^q_", "", group))
  
  df_p <- df %>%
    select(taxon, all_of(pval_cols)) %>%
    pivot_longer(-taxon, names_to = "group", values_to = "p") %>%
    mutate(group = gsub("^p_", "", group))
  
  df_final <- df_lfc %>%
    left_join(df_diff, by = c("taxon", "group")) %>%
    left_join(df_p, by = c("taxon", "group")) %>%
    left_join(df_q, by = c("taxon", "group"))
  
  
  df_final <- df_final %>%
    filter(abs(value) >= 1)
  
  if (nrow(df_final) == 0) {
    message("No taxa with |LFC| ≥ 1.")
    return(NULL)
  }
  
  df_final <- df_final %>%
    mutate(
      sig = case_when(is.na(q) ~ "", q < 0.001 ~ "***", q < 0.01  ~ "**", q < 0.05  ~ "*", TRUE ~ ""),
      alpha_val = ifelse(!is.na(q) & q < 0.05, 1, 0.3)
    )
  
  df_final <- df_final %>%
    mutate(color_var = NA_character_)
  
  if (color_mode == "direction") {
    df_final <- df_final %>%
      mutate(color_var = ifelse(value > 0, "Up", "Down")) %>%
      mutate(color_var = recode(color_var, !!!direction_labels))
    
    if (is.null(custom_colors)) {
      custom_colors <- setNames(c("#E64B35", "#4DBBD5"), unname(direction_labels))
    }
  }
  
  if (adj == "Yes") {
    df_final$plot_p <- df_final$q
    label_name <- "q"
  } else {
    df_final$plot_p <- df_final$p
    label_name <- "p"
  }
  
  
  groups <- unique(df_final$group)
  
  plots <- lapply(groups, function(g) {
    sub_df <- df_final %>% filter(group == g)
    
    sub_df <- sub_df %>%
      arrange(value) %>%   # asc
      mutate(taxon = factor(taxon, levels = taxon))
    
    p1 <- ggplot(sub_df,
                 aes(
                   x = taxon,
                   y = value,
                   fill = color_var,
                   alpha = alpha_val
                 )) +
      
      geom_col(color = "black") +
      scale_fill_manual(values = custom_colors) +
      scale_alpha_identity() +
      coord_flip() +
      labs(
        x = NULL,
        y = "LFC"
      ) +
      theme_classic() +
      # .shared_theme()+
      theme(
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.title = element_blank()
      ) #  + xlim(-5, +5)
    
    p2 <- ggplot(sub_df) +
      geom_text(aes(
        x = 1,
        y = taxon,
        label = paste0(
          ifelse(plot_p < 0.001, "<0.001", sprintf("%.3f", plot_p)),
          case_when(
            plot_p < 0.001 ~ "***",
            plot_p < 0.01  ~ "**",
            plot_p < 0.05  ~ "*",
            TRUE ~ ""
          )
        )
      )) +
      # .shared_theme() +
      theme_void() + 
      labs(x = "q-value") + 
      theme(axis.title.x = element_text())
    
    p1 + p2 + plot_layout(widths = c(4, 1))
  })
  
  final_plot <- wrap_plots(plots, nrow = 1)
  
  return(list(data = df_final, plot = final_plot))
}



make_ancombc_supp_table <- function(res,
                                    ps,
                                    term,
                                    group_var,
                                    comparison_level,
                                    comparison_name,
                                    taxonomic_level,
                                    tax_rank = NULL,
                                    reference_level = NULL) {
  
  lfc_col <- paste0("lfc_", term)
  se_col <- paste0("se_", term)
  W_col <- paste0("W_", term)
  p_col <- paste0("p_", term)
  q_col <- paste0("q_", term)
  diff_col <- paste0("diff_", term)
  diff_robust_col <- paste0("diff_robust_", term)
  
  required_cols <- c("taxon", lfc_col, se_col, W_col, p_col, q_col, diff_col)
  missing_cols <- setdiff(required_cols, colnames(res))
  
  if (length(missing_cols) > 0) {
    stop(
      paste0(
        "Missing columns in ANCOM-BC2 result: ",
        paste(missing_cols, collapse = ", ")
      )
    )
  }
  
  res_sig <- res %>%
    dplyr::transmute(
      comparison = comparison_name,
      taxonomic_level = taxonomic_level,
      taxon = as.character(taxon),
      log_fold_change = .data[[lfc_col]],
      standard_error = .data[[se_col]],
      ci_lower = .data[[lfc_col]] - 1.96 * .data[[se_col]],
      ci_upper = .data[[lfc_col]] + 1.96 * .data[[se_col]],
      test_statistic_W = .data[[W_col]],
      p_value = .data[[p_col]],
      adjusted_p_value = .data[[q_col]],
      differential_abundance = .data[[diff_col]],
      robust_differential_abundance = if (diff_robust_col %in% colnames(res)) {
        .data[[diff_robust_col]]
      } else {
        NA
      }
    )
  
  if (nrow(res_sig) == 0) {
    return(res_sig)
  }
  
  prev_df <- calculate_prevalence(
    ps = ps,
    taxa = res_sig$taxon,
    group_var = group_var,
    comparison_level = comparison_level,
    reference_level = reference_level,
    tax_rank = tax_rank
  )
  
  res_sig <- res_sig %>%
    dplyr::left_join(prev_df, by = "taxon") %>%
    dplyr::mutate(
      enriched_in = ifelse(
        log_fold_change > 0,
        comparison_level,
        reference_level
      )
    ) %>%
    dplyr::mutate(
      across(
        c(
          log_fold_change,
          standard_error,
          ci_lower,
          ci_upper,
          test_statistic_W,
          p_value,
          adjusted_p_value
        ),
        ~ round(., 4)
      )
    )
  
  return(res_sig)
}




#### 7. Rarefraction helper #### 



calculate_rarefaction_curves <- function(psdata, measures, depths) {
  
  estimate_rarified_richness <- function(psdata, measures, depth) {
    if(max(sample_sums(psdata)) < depth) return()
    psdata <- prune_samples(sample_sums(psdata) >= depth, psdata)
    
    rarified_psdata <- rarefy_even_depth(psdata, depth, verbose = FALSE)
    
    alpha_diversity <- estimate_richness(rarified_psdata, measures = measures)
    
    # as.matrix forces the use of melt.array, which includes the Sample names (rownames)
    molten_alpha_diversity <- melt(as.matrix(alpha_diversity), varnames = c('Sample', 'Measure'), value.name = 'Alpha_diversity')
    
    molten_alpha_diversity
  }
  
  names(depths) <- depths # this enables automatic addition of the Depth to the output by ldply
  rarefaction_curve_data <- plyr::ldply(depths, estimate_rarified_richness, psdata = psdata, measures = measures, .id = 'Depth', .progress = ifelse(interactive(), 'text', 'none'))
  
  # convert Depth from factor to numeric
  rarefaction_curve_data$Depth <- as.numeric(levels(rarefaction_curve_data$Depth))[rarefaction_curve_data$Depth]
  
  rarefaction_curve_data
}


#### 8. Depth summary ####
summarize_phy <- function(ps) {
  
  ps_sub <- prune_taxa(taxa_sums(ps) > 0, ps)
  
  # Read count
  reads <- sample_sums(ps_sub)
  mean_reads <- mean(reads)
  sd_reads <- sd(reads)
  
  # Taxonomic richness
  tax_table_df <- as.data.frame(tax_table(ps_sub))
  n_asv     <- ntaxa(ps_sub)
  n_phylum  <- n_distinct(na.omit(tax_table_df$Phylum))
  n_genus   <- n_distinct(na.omit(tax_table_df$Genus))
  n_species <- n_distinct(na.omit(tax_table_df$Species))
  
  # 결과 반환
  result_tbl  <- tibble(
    Total_reads = round(sum(reads), 1),
    Mean_Reads = round(mean_reads, 1),
    SD_Reads = round(sd_reads, 1),
    ASV_Count = n_asv,
    Phylum_Count = n_phylum,
    Genus_Count = n_genus,
    Species_Count = n_species)
  
  return(result_tbl)
}

summarize_phyloseq_by_site <- function(ps, site_variable, site_levels_input = NULL) {
  # site_variable: sample_data(ps) 내 그룹 변수명 (예: "Site2")
  # site_levels_input: 분석할 특정 site 부위들 (예: c("Forehead", "Cheek")), 지정하지 않으면 전체 사용
  
  # 변수 유효성 검사
  if (!site_variable %in% colnames(sample_data(ps))) {
    stop("The specified site_variable is not in the sample_data of the phyloseq object.")
  }
  
  # 전체 site 값 추출
  site_values <- as.character(sample_data(ps)[[site_variable]])
  all_site_levels <- unique(site_values)
  
  # 사용할 site levels 설정
  site_levels <- if (is.null(site_levels_input)) all_site_levels else site_levels_input
  
  # 각 site별로 요약 통계 계산
  result_list <- lapply(site_levels, function(site_level) {
    # 조건 벡터 직접 생성하여 subset
    samples_to_keep <- sample_names(ps)[sample_data(ps)[[site_variable]] == site_level]
    ps_sub <- prune_samples(samples_to_keep, ps)
    ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
    
    # Read count
    reads <- sample_sums(ps_sub)
    mean_reads <- mean(reads)
    sd_reads <- sd(reads)
    
    # Taxonomic richness
    tax_table_df <- as.data.frame(tax_table(ps_sub))
    n_asv     <- ntaxa(ps_sub)
    n_phylum  <- n_distinct(na.omit(tax_table_df$Phylum))
    n_genus   <- n_distinct(na.omit(tax_table_df$Genus))
    n_species <- n_distinct(na.omit(tax_table_df$Species))
    
    # 결과 반환
    tibble(
      Site = site_level,
      Mean_Reads = round(mean_reads,1),
      SD_Reads = round(sd_reads, 1),
      ASV_Count = n_asv,
      Phylum_Count = n_phylum,
      Genus_Count = n_genus,
      Species_Count = n_species
    )
  })
  
  # 최종 결과 결합
  result_tbl <- bind_rows(result_list)
  return(result_tbl)
}

#### 9. Reproducibility helper #### 
write_session_info <- function(path = "./sessionInfo.txt") {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  capture.output(sessionInfo(), file = path)
}



