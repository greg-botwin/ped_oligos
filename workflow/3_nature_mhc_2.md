HLA class II epitopes by deep motif deconvolution
================
Translational Genomics Group
10/18/2019

Protein Universe
================

![UniProt Reviewed (Swiss-Prot) - Manually annotated](data/uniprot_logo.png)

On 09-OCT-2019 20430 proteins meeting the following criteria were downloaded from Uniprot:
- reviewed:yes AND organism:"Homo sapiens (Human) \[9606\]"

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

Download Data
-------------

### Supplementary Data 1

``` bash
wget https://static-content.springer.com/esm/art%3A10.1038%2Fs41587-019-0289-6/MediaObjects/41587_2019_289_MOESM4_ESM.txt
```

### Supplementary Data 2

``` bash
wget https://static-content.springer.com/esm/art%3A10.1038%2Fs41587-019-0289-6/MediaObjects/41587_2019_289_MOESM5_ESM.txt
```

``` r
ms_racle_1 <- read_tsv("data/41587_2019_289_MOESM4_ESM.txt")
```

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_double(),
    ##   Sequence = col_character(),
    ##   Proteins = col_character(),
    ##   `Leading razor protein` = col_character(),
    ##   `Gene names` = col_character(),
    ##   `Protein names` = col_character(),
    ##   `Protein group IDs` = col_character()
    ## )

    ## See spec(...) for full column specifications.

``` r
ms_racle_2 <- read_tsv("data/41587_2019_289_MOESM5_ESM.txt")
```

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_double(),
    ##   Sequence = col_character(),
    ##   Proteins = col_character(),
    ##   `Leading razor protein` = col_character(),
    ##   `Gene names` = col_character(),
    ##   `Protein names` = col_character()
    ## )

    ## See spec(...) for full column specifications.

``` r
ms_racle_1 <- ms_racle_1 %>%
  select(`Leading razor protein`) %>%
  distinct()


ms_racle_2 <- ms_racle_2 %>%
  select(`Leading razor protein`) %>%
  distinct()

ms_racle <- bind_rows(ms_racle_1, ms_racle_2)

ms_racle <- ms_racle %>%
  distinct() %>%
  mutate(uniprotswissprot = word(`Leading razor protein`, sep = "-"))
```

Filter to Protein Reviewed
--------------------------

Missing 10 proteins not in Uniprot

``` r
table(ms_racle$uniprotswissprot %in% uniprot$Entry)
```

    ## 
    ## FALSE  TRUE 
    ##    10  9394

``` r
ms_racle <- ms_racle %>%
  filter(uniprotswissprot %in% uniprot$Entry) %>%
  select(uniprotswissprot)

table(ms_racle$uniprotswissprot %in% uniprot$Entry)
```

    ## 
    ## TRUE 
    ## 9394

``` r
write_tsv(ms_racle, "data/3_ms_racle.tsv")
```
