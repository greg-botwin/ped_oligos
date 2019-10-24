Genes of Special Interest
================
Translational Genomics Group
10/23/2019

IBD Genes
---------

I would use all genes in locus, unless implicated gene is known

``` r
ibd_genes <- read_excel("data/COPY__deLange_TableS2_241IBDloci.xlsx", skip = 8)
```

    ## New names:
    ## * `` -> ...6

``` r
locus_genes <- ibd_genes %>%
  filter(is.na(`Implicated gene`)) %>%
  dplyr::select(`Genes in locus`) %>%
  separate_rows('Genes in locus', sep = ",") %>%
  filter(str_detect(`Genes in locus`, "[[:alpha:]]")) %>%
  rename(gene_names = `Genes in locus`)


implicated_genes <- ibd_genes %>%
  filter(!is.na(`Implicated gene`)) %>%
  dplyr::select(`Implicated gene`) %>%
  rename(gene_names = `Implicated gene`)

ibd_genes <- bind_rows(locus_genes, implicated_genes)

ibd_genes <- distinct(ibd_genes)
```

``` r
ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", 
                      GRCh = 37)

# 701 returned
ids <- getBM(attributes = c("external_gene_name", "uniprotswissprot"),
      filters = "external_gene_name",
      values = ibd_genes$gene_names,
      mart = ensembl)

ids <- ids %>%
  filter(uniprotswissprot != "")

write_tsv(ids, "data/ids_special_interest.tsv")
```

``` r
ids <- read_tsv("data/ids_special_interest.tsv")
```

    ## Parsed with column specification:
    ## cols(
    ##   external_gene_name = col_character(),
    ##   uniprotswissprot = col_character()
    ## )

``` r
ibd_genes %>%
  filter(!gene_names %in% ids$external_gene_name) %>%
  write_tsv("data/ibd_genes_need_ids.tsv")

table(ibd_genes$gene_names %in% ids$external_gene_name)
```

    ## 
    ## FALSE  TRUE 
    ##    44   641

``` r
man_mapped <- read_tsv("data/uniprot-yourlist_M201910235C475328CEF75220C360D524E9D456CE4D8F3CH-filtered-rev--.tab")
```

    ## Parsed with column specification:
    ## cols(
    ##   `yourlist:M201910235C475328CEF75220C360D524E9D456CE4D8F3CH` = col_character(),
    ##   Entry = col_character(),
    ##   `Entry name` = col_character(),
    ##   Status = col_character(),
    ##   `Protein names` = col_character(),
    ##   `Gene names` = col_character(),
    ##   Organism = col_character(),
    ##   Length = col_double()
    ## )

``` r
man_mapped <- man_mapped %>%
  rename(external_gene_name = `yourlist:M201910235C475328CEF75220C360D524E9D456CE4D8F3CH`) %>%
  rename(uniprotswissprot = Entry) %>%
  dplyr::select(external_gene_name, uniprotswissprot)

ids <- bind_rows(ids, man_mapped)

table(ibd_genes$gene_names %in% ids$external_gene_name)
```

    ## 
    ## FALSE  TRUE 
    ##    41   644

following uniprot ids need updating P08107 remapped to P0DMV8 and P0DMV9 Q9UBA6 deleted P62158 remapped to P0DP23, P0DP24 and P0DP25 deleted Q4KN68

``` r
ids <- inner_join(ibd_genes, ids, by = c("gene_names" = "external_gene_name"))
bad_ids <- c("P08107", "Q9UBA6", "P62158", "Q4KN68")
good_ids <-  c("P0DMV8", "P0DMV9", "P0DP23", "P0DP24", "P0DP25")
ids <- ids %>%
  filter(!uniprotswissprot %in% bad_ids) %>%
  dplyr::select(uniprotswissprot) %>%
  bind_rows(add_ids <- data.frame("uniprotswissprot" = good_ids))
```

    ## Warning in bind_rows_(x, .id): binding character and factor vector,
    ## coercing into character vector

``` r
ids %>%
  distinct() %>%
  write_tsv("data/6_known_ibd_genes.tsv")
```
