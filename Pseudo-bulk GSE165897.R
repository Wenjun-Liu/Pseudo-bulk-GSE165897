## ----setup, echo=FALSE--------------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(
  echo = TRUE,
  message = FALSE,
  warning = FALSE,
  fig.align = "center",
  fig.width = 10
)


## ----packages-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(tidyverse)
library(yaml)
library(scales)
library(pander)
library(glue)
library(edgeR)
library(AnnotationHub)
library(ensembldb)
library(cowplot)
library(ggfortify)
library(magrittr)
library(cqn)
library(ggrepel)
library(DT)
library(Seurat)
library(corrplot)
library(plotly)
library(patchwork)
library(colorspace)
library(ggpubr)
library(ggforce)
library(MAST)
library(pheatmap)
library(BiocParallel)
library(sSNAPPY)
library(fgsea)
library(msigdbr)
library(limma)
library(goseq)


## ----options------------------------------------------------------------------------------------------------------------------------------------------------------------------------
panderOptions("table.split.table", Inf)
panderOptions("big.mark", ",")
theme_set(theme_bw())


## ----formatP------------------------------------------------------------------------------------------------------------------------------------------------------------------------
formatP <- function(p, m = 0.0001){
out <- rep("", length(p))
out[p < m] <- sprintf("%.2e", p[p<m])
out[p >= m] <- sprintf("%.4f", p[p>=m])
out
}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Sys.setenv(VROOM_CONNECTION_SIZE = 5000000L)



## ----ah-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
ah <- AnnotationHub() %>%
  subset(rdataclass == "EnsDb") %>%
  subset(str_detect(description, "101")) %>%
  subset(genome == "GRCh38")
stopifnot(length(ah) == 1)


## ----ensDb--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
ensDb <- ah[[1]]
genesGR <- genes(ensDb)
transGR <- transcripts(ensDb)


## ----addTxLen-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
mcols(transGR) <- mcols(transGR) %>%
  cbind(
    transcriptLengths(ensDb)[rownames(.), c("nexon", "tx_len")]
  )


## ----addGcLen2Genes-----------------------------------------------------------------------------------------------------------------------------------------------------------------
mcols(genesGR) <- mcols(genesGR) %>%
  as.data.frame() %>%
  dplyr::select(
    gene_id, gene_name, gene_biotype, entrezid
  ) %>%
  left_join(
    mcols(transGR) %>%
      as.data.frame() %>%
      mutate(
        tx_support_level = case_when(
          is.na(tx_support_level) ~ 1L, 
          TRUE ~ tx_support_level
        )
      ) %>%
      group_by(gene_id) %>%
      summarise(
        n_tx = n(),
        longest_tx = max(tx_len),
        ave_tx_len = mean(tx_len),
        gc_content = sum(tx_len*gc_content) / sum(tx_len)
      ) %>%
      mutate(
        bin_length = cut(
          x = ave_tx_len,
          labels = seq_len(10),
          breaks = quantile(ave_tx_len, probs = seq(0, 1, length.out = 11)),
          include.lowest = TRUE
        ),
        bin_gc = cut(
          x = gc_content,
          labels = seq_len(10),
          breaks = quantile(gc_content, probs = seq(0, 1, length.out = 11)),
          include.lowest = TRUE
        ),
        bin = paste(bin_gc, bin_length, sep = "_")
      ),
    by = "gene_id"
  ) %>%
  set_rownames(.$gene_id) %>%
  as("DataFrame")


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cellInfor <- read_tsv(gzfile(here::here("GSE165897_cellInfo_HGSOC.tsv.gz")))
counts <- read_tsv(gzfile(here::here("GSE165897_UMIcounts_HGSOC.tsv.gz")))
counts <- counts %>%
    column_to_rownames("...1")


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cellInfor %>%
    dplyr::select(cell_type, cell_subtype) %>%
    unique() %>%
    chop("cell_subtype") %>%
    mutate(
        cell_subtype = vapply(.$cell_subtype, function(x){
            paste(x, collapse = "; ")
        }, character(1))
    )


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
sample_info <- read_csv(here::here("GSE165897_SampleMeta.csv"))
sample_meta <- cellInfor %>%
    dplyr::filter(cell_type == "EOC") %>%
    mutate(name = paste(patient_id, treatment_phase, sep = "_")) %>%
    dplyr::select( name, treatment_phase, patient_id, anatomical_location) %>%
    unique() %>%
    left_join(
        sample_info
    )
#write_tsv(sample_meta, here::here("sample_meta.tsv"))
sample_info %>%
  dplyr::rename(
    Patient = patient_id
  ) %>%
  pander(
    caption = "Clinical information of patients from whom samples were derived."
  )


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
minCPM <- 1.5
minSample <- 11


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
pseudo_by <- cellInfor %>%
    dplyr::filter(cell_type == "EOC") %>%
    mutate(group = paste(patient_id, treatment_phase, sep = "_")) %>%
    split(f = .$group) %>%
    lapply(pull, cell)
pseudo_bulk <- sapply(names(pseudo_by), function(x){
    counts %>%
        .[,colnames(.) %in% pseudo_by[[x]]] %>%
        rowSums()
}, simplify = FALSE) %>%
    do.call(rbind,.) %>%
    t()
gene2keep <- pseudo_bulk %>%
  edgeR::cpm(log = TRUE) %>%
    is_greater_than(minCPM) %>%
    rowSums() %>%
    is_weakly_greater_than(minSample)
# pseudo_bulk <- pseudo_bulk[gene2keep, ]


## ----plotDensities, fig.height=5, fig.width = 14, fig.cap="*Distributions of logCPM values on pseudo-bulk counts, A) before and B) after filtering of undetectable genes.*"---------
p1 <- pseudo_bulk %>%
  edgeR::cpm(log = TRUE) %>%
  as.data.frame() %>%
  rownames_to_column("gene_id") %>%
  as_tibble() %>%
  pivot_longer(
    cols = -"gene_id",
    names_to = "name",
    values_to = "logCPM"
  ) %>%
  left_join(sample_meta) %>%
  ggplot(aes(logCPM, stat(density), colour = patient_id, linetype = treatment_phase)) +
  geom_density() +
  labs(
    y = "Density",
    colour = "Patient",
    linetype = "Treatment Phase"
  )
p2 <- pseudo_bulk[gene2keep, ] %>%
  edgeR::cpm(log = TRUE) %>%
  as.data.frame() %>%
  rownames_to_column("gene_id") %>%
  as_tibble() %>%
  pivot_longer(
    cols = -"gene_id",
    names_to = "name",
    values_to = "logCPM"
  ) %>%
  left_join(sample_meta) %>%
  ggplot(aes(logCPM, stat(density), colour = patient_id, linetype = treatment_phase)) +
  geom_density() +
  labs(
    y = "Density",
    colour = "Patient",
    linetype = "Treatment Phase"
  )
plot_grid(
    plot_grid(
    p1 + theme(legend.position = "none"), 
    p2 + theme(legend.position = "none"), 
    labels = c("A", "B"),
    nrow = 1), 
    get_legend(p1), 
    nrow = 1, 
    rel_widths = c(6,1))


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
pseudo_bulk <- pseudo_bulk[gene2keep, ]
genesGR <- genesGR %>%
     as.data.frame() %>%
    .[!duplicated(.$gene_name),]
rownames(genesGR) <- genesGR[,"gene_name"]


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
gene2keep2 <- rownames(pseudo_bulk) %in% genesGR$gene_name
pseudo_bulk <- pseudo_bulk[gene2keep2, ]


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

dge <- pseudo_bulk %>%
        DGEList(
            samples = sample_meta %>%
                set_rownames(.$name) %>%
                .[match(colnames(pseudo_bulk), rownames(.)),],
            genes = genesGR[rownames(.),]
        ) %>%
  calcNormFactors()


## ----plotLibSizes, fig.height=8, fig.width = 12, fig.cap = "*Library sizes of all samples after removal of undetectable genes. The common-use minimum library size of 10 million reads is shown as a dashed line.*"----
dge$samples %>%
    ggplot(aes(treatment_phase, lib.size, fill = treatment_phase)) +
    geom_col() +
    geom_hline(yintercept = 1e7, linetype = 2) +
    facet_wrap(~patient_id, ncol = 5) +
    scale_y_continuous(
        labels = comma, expand = expansion(c(0, 0.05))
    ) +
    scale_fill_manual(values= phase_cols) +
    labs(x = "Sample Name", y = "Library Size", fill = "Treatment")  +
    theme(legend.position = c(11/12, 1/6))


## ----pca----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
pca <- dge %>%
  edgeR::cpm(log = TRUE) %>%
  t() %>%
  prcomp() 


## ----fig.cap="*PCA plot with arrows pointing from pre-treatment to corresponding post-treatment sample. No consistent pattern could be observed>*"----------------------------------
broom::tidy(pca) %>% 
  dplyr::filter(PC %in% 1:2) %>% 
  dplyr::rename(name = row) %>% 
  left_join(dge$samples) %>% 
  pivot_wider(names_from = "PC", values_from = "value", names_prefix = "PC") %>%
  ggplot(
    aes(PC1, PC2, colour = treatment_phase)
  ) +
  geom_point(size = 2) +
  geom_segment(
    aes(x = `PC1_treatment-naive`, y = `PC2_treatment-naive`, xend = `PC1_post-NACT`, yend = `PC2_post-NACT`),
    data = . %>% 
      dplyr::select(patient_id, treatment_phase,  PC1, PC2) %>%
      pivot_wider(names_from = "treatment_phase", values_from = c("PC1", "PC2")) ,
    inherit.aes = FALSE,
    arrow = arrow(length = unit(0.02, "npc"), type = "closed"),
    colour = "grey30"
  ) +
  geom_label_repel(
    aes(label = patient_id, 
        color = Stage),
    colour = "black",
    fill = rgb(1, 1, 1, 0.5),
    size = 2,
    data = . %>% 
      group_by(patient_id, Stage) %>% 
      summarise(PC1 = mean(PC1), PC2 = mean(PC2), .groups = "drop")
  ) +
  labs(
    x = glue("PC1 ({percent(summary(pca)$importance[2, 1], 0.1)})"),
    y = glue("PC2 ({percent(summary(pca)$importance[2, 2], 0.1)})"),
    colour = "Treatment Phase"
  ) +
  scale_colour_manual(values = phase_cols) 


## ----plotPCA_GC_Length, fig.height=4,fig.width=8, fig.cap = "*Contribution of each GC/Length Bin to PC1 and PC2. Fill colours indicate the t-statistic, with tranparency denoting significance as -log10(p), using Bonferroni-adjusted p-values. The number of genes in each bin is indicated by the circle size.*"----
dge$genes %>%
  dplyr::select(gene_id, ave_tx_len, gc_content) %>%
  mutate(
    GC = cut(
      x = gc_content,
      labels = seq_len(10),
      breaks = quantile(gc_content, probs = seq(0, 1, length.out = 11)),
      include.lowest = TRUE
    ),
    Length = cut(
      x = ave_tx_len,
      labels = seq_len(10),
      breaks = quantile(ave_tx_len, probs = seq(0, 1, length.out = 11)),
      include.lowest = TRUE
    ),
    bin = paste(GC, Length, sep = "_"),
    PC1 = pca$rotation[, "PC1"],
    PC2 = pca$rotation[, "PC2"]
  ) %>%
  pivot_longer(
    cols = c("PC1", "PC2"),
    names_to = "PC",
    values_to = "value"
  ) %>%
  group_by(PC, GC, Length, bin) %>%
  summarise(
    Size = n(),
    mean = mean(value),
    sd = sd(value),
    t = t.test(value)$statistic,
    p = t.test(value)$p.value,
    adjP = p.adjust(p, method = "bonf")
  ) %>%
  ggplot(
    aes(Length, GC, colour = t, alpha = -log10(adjP), size = Size)
  ) +
  geom_point() +
  facet_wrap(~PC) +
  scale_colour_gradient2() +
  scale_size_continuous(range = c(1, 10)) +
  labs(alpha = expression(paste(-log[10], p))) +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom"
  ) 


## ----cqn----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cqn <-  with(
    dge,
    cqn(
        counts= counts,
        x = genes$gc_content,
        lengths = genes$ave_tx_len
    )
)
dge$offset <- cqn$glm.offset
logCPM <- cqn$y + cqn$offset
pcaPost <- logCPM %>%
    t() %>%
    prcomp() 



## ----eval=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------------------
## write_tsv(logCPM, here::here("logCPM.tsv"))

