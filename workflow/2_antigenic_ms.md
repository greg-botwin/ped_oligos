Antigenic by Mass Spec
================
Translational Genomics Group
10/9/2019

Protein Universe
----------------

![UniProt Reviewed (Swiss-Prot) - Manually annotated](data/uniprot_logo.png)

On 21-JUN-2019 20431 proteins meeting the following criteria were downloaded from Uniprot:
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

``` r
length(unique(uniprot$Entry))
```

    ## [1] 20431

Mommeng et. Al. Mass Spec Analysis
----------------------------------

![Mass Spec Analysis](data/ms_photo.jpg)

Sampling From the Proteome to the Human Leukocyte Antigen-DR (HLA-DR) Ligandome Proceeds Via High Specificity - 1,220 unique proteins id across 3 preperation types

``` r
scx <- read_csv("data/SYSMHC00015_Mommeng_161005_Heck_Netherlands_SCX_combined.csv")
```

    ## Parsed with column specification:
    ## cols(
    ##   SysteMHC_ID = col_character(),
    ##   SampleID = col_character(),
    ##   elution_event_ID = col_character(),
    ##   search_hit = col_character(),
    ##   protein_id = col_character(),
    ##   length = col_double(),
    ##   prob = col_double(),
    ##   spectral_counts = col_double(),
    ##   numNeighbors = col_double(),
    ##   MHCClass = col_character()
    ## )

``` r
scx_hcd <- read_csv("data/SYSMHC00015_Mommeng_161005_Heck_Netherlands_HCD_SCX_combined.csv")
```

    ## Parsed with column specification:
    ## cols(
    ##   SysteMHC_ID = col_character(),
    ##   SampleID = col_character(),
    ##   elution_event_ID = col_character(),
    ##   search_hit = col_character(),
    ##   protein_id = col_character(),
    ##   length = col_double(),
    ##   prob = col_double(),
    ##   spectral_counts = col_double(),
    ##   numNeighbors = col_double(),
    ##   MHCClass = col_character()
    ## )

``` r
scx_etd <- read_csv("data/SYSMHC00015_Mommeng_161005_Heck_Netherlands_ETD_SCX_combined.csv")
```

    ## Parsed with column specification:
    ## cols(
    ##   SysteMHC_ID = col_character(),
    ##   SampleID = col_character(),
    ##   elution_event_ID = col_character(),
    ##   search_hit = col_character(),
    ##   protein_id = col_character(),
    ##   length = col_double(),
    ##   prob = col_double(),
    ##   spectral_counts = col_double(),
    ##   numNeighbors = col_double(),
    ##   MHCClass = col_character()
    ## )

``` r
ms <- bind_rows(scx, scx_etd, scx_etd, .id = "id")
ms <- ms %>%
  separate(protein_id, into = c("dnk", "uniprot_id", "uniprot_name"), sep = "\\|", 
           remove = FALSE, extra = "merge") %>%
  select(uniprot_id, uniprot_name) %>%
  distinct(uniprot_id, .keep_all = TRUE)

table(ms$uniprot_id %in% uniprot$Entry)
```

    ## 
    ## FALSE  TRUE 
    ##     5  1217

``` r
ms %>%
  filter(!uniprot_id %in% uniprot$Entry)
```

    ## # A tibble: 5 x 2
    ##   uniprot_id      uniprot_name
    ##   <chr>           <chr>       
    ## 1 ZZ_FGCZCont0112 ""          
    ## 2 P62158          CALM_HUMAN  
    ## 3 P0CW22          RS17L_HUMAN 
    ## 4 ZZ_FGCZCont0255 ""          
    ## 5 ZZ_FGCZCont0153 ""

``` r
ms <- ms %>%
  mutate(uniprot_id = if_else(uniprot_id == "P62158", "P0DP23", uniprot_id)) %>%
  mutate(uniprot_id = if_else(uniprot_id == "P0CW22", "P08708", uniprot_id)) %>% 
  add_row(uniprot_id = "P0DP24", uniprot_name = "CALM2_HUMAN")

ms_list <- ms %>%
  filter(uniprot_id %in% uniprot$Entry)

ms_list %>%
  rename(uniprotswsissprot = uniprot_id) %>%
  select(uniprotswsissprot) %>%
  write_tsv("data/2_antigenic_ms.tsv")
```

Missing ids manually searched on UniProt on 21-JUN-2019

-   ZZ is NA
-   P62158 is P0DP23 (CALM1\_HUMAN) and P0DP24 (CALM2\_HUMAN)
-   P0CW22 updated to P08708
