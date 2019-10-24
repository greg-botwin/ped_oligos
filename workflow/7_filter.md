Selecting Proteins
================
Translational Genomics Group
10/9/2019

``` r
uniprot <- read_tsv("data/uniprot-filtered-organism__Homo+sapiens+(Human)+[9606]_+AND+review--.tab")
```

    ## Parsed with column specification:
    ## cols(
    ##   Entry = col_character(),
    ##   `Entry name` = col_character(),
    ##   Status = col_character(),
    ##   `Protein names` = col_character(),
    ##   `Gene names` = col_character(),
    ##   Organism = col_character(),
    ##   Length = col_double()
    ## )

Antigenic
---------

9613 antigenic proteins

``` r
iedb <- read_tsv("data/1_anitgenic_iedb.tsv")
```

    ## Parsed with column specification:
    ## cols(
    ##   uniprotswsissprot = col_character()
    ## )

``` r
iedb <- iedb %>%
  rename(uniprotswissprot = uniprotswsissprot)
ms <- read_tsv("data/2_antigenic_ms.tsv")
```

    ## Parsed with column specification:
    ## cols(
    ##   uniprotswsissprot = col_character()
    ## )

``` r
ms_racle <- read_tsv("data/3_ms_racle.tsv")
```

    ## Parsed with column specification:
    ## cols(
    ##   uniprotswissprot = col_character()
    ## )

``` r
antigenic <- bind_rows(iedb, ms_racle, ms_racle) 

antigenic <- antigenic %>%
  distinct()
```

``` r
table(antigenic$uniprotswissprot %in% uniprot$Entry)
```

    ## 
    ## TRUE 
    ## 9613

Special Interest
----------------

### IBD

``` r
ibd <- read_tsv("data/6_known_ibd_genes.tsv")
```

    ## Parsed with column specification:
    ## cols(
    ##   uniprotswissprot = col_character()
    ## )

``` r
ibd <- ibd %>%
  mutate(rank = 0) %>%
  dplyr::select(uniprotswissprot, rank)

table(ibd$uniprotswissprot %in% uniprot$Entry)
```

    ## 
    ## TRUE 
    ##  703

GFP
---

238 aa long P42212 (GFP\_AEQVI) ProteinGreen fluorescent protein Gene GFP Organism Aequorea victoria (Jellyfish)

``` r
gfp <- data.frame("uniprotswissprot"  = c("P42212"),
                  "Length" = 238,
                  "rank" = 0)
```

Epithelial
----------

2927 proteins epithelail sc proteins 1745 antigenic

``` r
epi_sc <- read_tsv("data/4_epithelial_uc_sc.tsv")
```

    ## Parsed with column specification:
    ## cols(
    ##   uniprotswissprot = col_character()
    ## )

``` r
table(epi_sc$uniprotswissprot %in% uniprot$Entry)
```

    ## 
    ## TRUE 
    ## 2927

``` r
table(epi_sc$uniprotswissprot %in% antigenic$uniprotswissprot)
```

    ## 
    ## FALSE  TRUE 
    ##  1182  1745

``` r
epi_sc %>%
  filter(uniprotswissprot %in% antigenic$uniprotswissprot)
```

    ## # A tibble: 1,745 x 1
    ##    uniprotswissprot
    ##    <chr>           
    ##  1 P30047          
    ##  2 O43278          
    ##  3 O43482          
    ##  4 Q9BXS6          
    ##  5 P35558          
    ##  6 O43663          
    ##  7 Q15004          
    ##  8 P00915          
    ##  9 P00918          
    ## 10 Q9NQW6          
    ## # … with 1,735 more rows

``` r
gtex <- read_tsv("data/5_epithelial_gtex.tsv")
```

    ## Parsed with column specification:
    ## cols(
    ##   tissue = col_character(),
    ##   uniprotswissprot = col_character(),
    ##   avg_log2TPM_across_ensg = col_double(),
    ##   rank = col_double()
    ## )

``` r
gtex <- gtex %>%
  dplyr::select(uniprotswissprot, rank)
```

``` r
epithelial <- bind_rows(epi_sc, gtex)
epithelial <- epithelial %>%
  mutate(rank = if_else(is.na(rank), 0.5, rank))
```

Protein to Oligo
----------------

-   Each oligo can be up to 230nt long
-   prefix and suffix 21nt long
-   need remaining nucleotide length to be divisible by 3 to stay in frame
-   unique oligo sequence will be 186nt + 21nt prefix + 21nt suffix = 230
-   unique oligo sequence will cover 62 amino acids, with 31 overlap
-   oligos needs to tile one protein is protein aa length / 31 aa per tile =, rounded up
-   we can have up to 100,000 unique oligo sequences
-   to olgio the entire protein universe requires 377,134 oligos

``` r
aa_length_2_oligos <- function(length){
  ceiling(length/31)
}

uniprot <- uniprot %>%
  mutate(n_oligos = aa_length_2_oligos(Length)) %>%
  dplyr::select(Entry, `Entry name`, n_oligos)

gfp <- gfp %>%
  mutate(n_oligos = aa_length_2_oligos(Length)) %>%
  dplyr::select(uniprotswissprot, n_oligos)
```

### Oligo Selection

``` r
antigenic_epithelial <- epithelial %>%
  filter(uniprotswissprot %in% antigenic$uniprotswissprot)

antigenic_epithelial_si <- bind_rows(antigenic_epithelial, ibd)

antigenic_epithelial_si <- antigenic_epithelial_si %>%
  left_join(., uniprot, by = c("uniprotswissprot" = "Entry"))

antigenic_epithelial_si <- bind_rows(antigenic_epithelial_si, gfp)
```

    ## Warning in bind_rows_(x, .id): binding character and factor vector,
    ## coercing into character vector

``` r
antigenic_epithelial_si %>%
  filter(rank == 0 ) %>%
  summarise(total_oligo = sum(n_oligos))
```

    ## # A tibble: 1 x 1
    ##   total_oligo
    ##         <dbl>
    ## 1       13054

``` r
antigenic_epithelial_si %>%
  filter(rank < 1 ) %>%
  summarise(total_oligo = sum(n_oligos))
```

    ## # A tibble: 1 x 1
    ##   total_oligo
    ##         <dbl>
    ## 1       52544

``` r
# arrange by rank, keep only first time protein listed
oligos <- antigenic_epithelial_si %>%
  arrange(rank) %>%
  distinct(uniprotswissprot, .keep_all = TRUE) %>%
  mutate(roll_total = cumsum(n_oligos)) %>%
  filter(roll_total < 100000)
```

Annotate
--------

``` r
uniprot_annotate <- read_tsv("data/uniprot-filtered-organism__Homo+sapiens+(Human)+[9606]_+AND+review--.tab")
```

    ## Parsed with column specification:
    ## cols(
    ##   Entry = col_character(),
    ##   `Entry name` = col_character(),
    ##   Status = col_character(),
    ##   `Protein names` = col_character(),
    ##   `Gene names` = col_character(),
    ##   Organism = col_character(),
    ##   Length = col_double()
    ## )

``` r
annotated_oligos <- uniprot_annotate %>%
  dplyr::select(Entry, `Protein names`, `Gene names`) %>%
  inner_join(., oligos, by = c("Entry" = "uniprotswissprot"))

write_csv(annotated_oligos, "data/annotated_selected_proteins.csv")
```
