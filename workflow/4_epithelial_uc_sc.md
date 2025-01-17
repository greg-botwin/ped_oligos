Epithelial by Human UC Single Cell Differntial Gene Expression Analysis
================
Translational Genomics Group
10/9/2019

Differntially expressed genes from epithelial cells identified through single cell analysis in human patients with ulcerative colitits and controls.

All genes differentially expressed with a padj &lt; 0.05. No log fold expression filter was used.

Epithelial vs Non-Epithelial Cell Genes
---------------------------------------

Table S2. Marker Genes for Cell Subsets, Lineages, and Sub-clusters in Healthy Tissue, Related to Figure 1.

1264 Unqiue Genes

``` r
epi_genes <- read_excel("data/1-s2.0-S0092867419307329-mmc2.xlsx", 
                        sheet = "Epithelial")

epi_genes <- epi_genes %>%
  dplyr::select(gene) %>%
  distinct()
```

Epithelial UC Inflammed vs Healthy Genes
----------------------------------------

Table S4. Differentially Expressed Genes for Cell Subsets and Lineages during Disease, Related to Figure 3.

1241 Unique Genes

``` r
epi_non_infl_v_health <- read_excel("data/1-s2.0-S0092867419307329-mmc4.xlsx",
                           sheet = 1)

epi_non_infl_v_health <- epi_non_infl_v_health %>%
  dplyr::select(gene) %>% 
  distinct()
```

Epithelial UC Not-Inflammed vs Healthy Genes
--------------------------------------------

1200 uique genes

``` r
epi_infl_v_health <- read_excel("data/1-s2.0-S0092867419307329-mmc4.xlsx",
                                sheet = 4)

epi_infl_v_health <- epi_infl_v_health %>%
  dplyr::select(gene) %>%
  distinct()
```

Epithelial Inflammed vs Not Inflammed Genes
-------------------------------------------

969 Unqiue Genes

``` r
epi_infl_v_non_infl <- read_excel("data/1-s2.0-S0092867419307329-mmc4.xlsx",
                                  sheet = 7)

epi_infl_v_non_infl <- epi_infl_v_non_infl %>%
  dplyr::select(gene) %>%
  distinct()
```

Epithelial Single Cell Genes
----------------------------

4674 total genes 3181 non overlapping genes

``` r
epi_sc_genes <- bind_rows(epi_genes, epi_infl_v_health, epi_infl_v_non_infl, epi_non_infl_v_health)
epi_sc_genes <- epi_sc_genes %>%
  distinct()
```

Get Uniprot ID from Gene Name Using Ensemble Biomart
----------------------------------------------------

-   3181 gene names submitted
-   5051 uniprot ids returned
-   2907 unique uniprot ids returned
-   328 of the 3181 genes missing a uniprot ID

``` r
ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", 
                      GRCh = 37)

# 5051 returned
ids <- getBM(attributes = c("external_gene_name", "uniprotswissprot"),
      filters = "external_gene_name",
      values = epi_sc_genes$gene,
      mart = ensembl)

write_tsv(ids, "data/ids.tsv")
```

``` r
ids <- read_tsv("data/ids.tsv")
```

    ## Parsed with column specification:
    ## cols(
    ##   external_gene_name = col_character(),
    ##   uniprotswissprot = col_character()
    ## )

``` r
#2907 gene id combos returned
ids <- ids %>%
  dplyr::filter(uniprotswissprot != "") %>%
  distinct()

table(epi_sc_genes$gene %in% ids$external_gene_name)
```

    ## 
    ## FALSE  TRUE 
    ##   328  2853

### manually uploaded to uniprot

-   On 08-OCT-2019 Gene Names were uploaded to Uniprot Database for Search

27 out of 328 Gene name identifiers were successfully mapped to 25 UniProtKB IDs

``` r
epi_sc_genes %>%
  filter(!gene %in% ids$external_gene_name) %>% 
  write_tsv("data/single_cell_missing_uni.tsv")

uniprot_matched <- read_tsv("data/uniprot-yourlist_M20191008216DA2B77BFBD2E6699CA9B6D1C41EB24BCB6BO-filtered-rev--.tab")
```

    ## Parsed with column specification:
    ## cols(
    ##   `yourlist:M20191008216DA2B77BFBD2E6699CA9B6D1C41EB24BCB6BO` = col_character(),
    ##   Entry = col_character(),
    ##   `Entry name` = col_character(),
    ##   Status = col_character(),
    ##   `Protein names` = col_character(),
    ##   `Gene names` = col_character(),
    ##   Organism = col_character(),
    ##   Length = col_double()
    ## )

``` r
uniprot_matched  <- uniprot_matched %>%
  rename(external_gene_name = "yourlist:M20191008216DA2B77BFBD2E6699CA9B6D1C41EB24BCB6BO",
         uniprotswissprot = Entry) %>%
  dplyr::select(external_gene_name, uniprotswissprot)
```

### Remaining Manually updated to biomart

-   ON 08-OCT\_2019 Gene Names were queried manually using Biomart

unmathced genes are combination of misc\_RNA, snoRNA, lncRNA, miRNA, mt\_tRNA, processed psudeogene, rRNA\_pseduogene, snoRNA, TEX, unprocessed\_pseduogene and transcribed proc/unpr pseduogenes

2 protein coding genes remain, only 1 with a reviewed human sequence

``` r
epi_sc_genes %>%
  filter(!gene %in% ids$external_gene_name) %>% 
  filter(!gene %in% uniprot_matched$external_gene_name) %>%
  write_tsv("data/single_cell_missing_uni_after_man.tsv")

mart_matched <- read_tsv("data/mart_export 10_08_19.txt")
```

    ## Parsed with column specification:
    ## cols(
    ##   `Gene stable ID` = col_character(),
    ##   `Gene name` = col_character(),
    ##   `Gene type` = col_character()
    ## )

``` r
mart_matched <- mart_matched %>%
  filter(`Gene name` == "AC145212.1") %>%
  rename(external_gene_name = `Gene name`) %>%
  mutate(uniprotswissprot = "Q8WZ33") %>%
  dplyr::select(external_gene_name, uniprotswissprot)
```

Totals
------

five proteins needs updates to match the swissprot database - P0CB46 needs to be removed - A8MTB1 needs to be removed - Q06430 is now Q8N0V5 - Q8NFS9 is now Q8N0V5, need to remove dup - P08107 is now P0DMV8 and P0DMV9

``` r
ids <- bind_rows(ids, uniprot_matched, mart_matched)

ids <- ids %>%
  filter(uniprotswissprot != "P0CB46") %>%
  filter(uniprotswissprot != "A8MTB1") %>%
  filter(uniprotswissprot != "Q06430") %>%
  filter(uniprotswissprot != "Q8NFS9") %>%
  filter(uniprotswissprot != "P08107") %>%
  bind_rows(., new <- data.frame("uniprotswissprot" = c("P0DMV8", "P0DMV9", "Q8N0V5"))) %>%
  distinct(uniprotswissprot, .keep_all = TRUE)
```

    ## Warning in bind_rows_(x, .id): binding character and factor vector,
    ## coercing into character vector

-2873 out of 3181 genes matched to 2928 ids

``` r
table(epi_sc_genes$gene %in% ids$external_gene_name)
```

    ## 
    ## FALSE  TRUE 
    ##   308  2873

``` r
ids %>%
  dplyr::select(uniprotswissprot) %>%
  write_tsv("data/4_epithelial_uc_sc.tsv")
```

``` r
detach("package:biomaRt", unload=TRUE)
```
