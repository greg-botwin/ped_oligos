Antigenic IEDB Proteins
================
Translational Genomics Group
6/21/2019

Protein Universe
----------------

![UniProt Reviewed (Swiss-Prot) - Manually annotated](data/uniprot_logo.png)

On 21-JUN-2019 20431 proteins meeting the following criteria were downloaded from Uniprot:
- reviewed:yes AND organism:"Homo sapiens (Human) \[9606\]"

``` r
uniprot <- read_tsv("data/uniprot-filtered-organism__Homo+sapiens+(Human)+[9606]_+AND+review--.tab")

length(unique(uniprot$Entry))
```

    ## [1] 20431

IEDB
----

![IEDB Search Parameters](data/iedb_logo.png)

-   Epitope: Any Epitopes
-   Antigen: Homo sapiens (human) (ID:9606, Homo sapiens)
-   Host: Human
-   Assay: Positive Assays Only, T Cell Assays, B Cell Assays, MHC Ligand Assays
-   MHC Restriction: MHC CLass II
-   Disease: Any Disease

Antigen Table Export on 21-JUN-2019

``` r
# read in data
iedb <- read_csv("data/antigen_table_export_1563304473.csv", skip = 1)
colnames(iedb) <- make.names(colnames(iedb))

# extrac tuniprot id from url
iedb <- iedb %>%
  mutate(uniprot_id = str_replace(Antigen.ID, ".*\\/", "")) %>%
  mutate(uniprot_id = str_trim(uniprot_id, side = c("both")))

# summarise
iedb %>%
  summarise(n_antigens = length(unique(uniprot_id)),
                                n_epitopes = sum(X..Epitopes))
```

    ## # A tibble: 1 x 2
    ##   n_antigens n_epitopes
    ##        <int>      <dbl>
    ## 1       4589      54078

### Number of IEDB Antigens in UniProt

``` r
table(iedb$uniprot_id %in% uniprot$Entry)
```

    ## 
    ## FALSE  TRUE 
    ##   226  4363

``` r
# create list of 226 ids to query uniprot
iedb %>%
  filter(!uniprot_id %in% uniprot$Entry) %>%
  select(uniprot_id) %>%
  write_tsv("data/iedb_ids_to_match.tsv")
```

List uploaded to <https://www.uniprot.org/uploadlists/> on 16-JUL-2019 Searching: - From: UniProtKB AC/ID - To: UniProt KB

219 out of 226 UniProtKB AC/ID identifiers were successfully mapped to 222 UniProtKB IDs

-   remove 156 entries that correspond to non-revewed proteins
-   remove 35 proteins subsequently deleted from uniprot db
-   one entry Q8WYQ7,F8W9W4 will be split into two entres
-   final list is 4376

``` r
unmatched_iedb_ids <- read_tsv("data/uniprot-yourlist_M201907166746803381A1F0E0DB47453E0216320D53922A0.tab")

# remove 156 entries that correspond to non-revewed proteins
# remove 35 proteins subsequently deleted from uniprot db
unmatched_iedb_ids <- unmatched_iedb_ids %>%
  filter(Status == "reviewed") %>%
  select(`yourlist:M201907166746803381A1F0E0DB47453E0216320D53922A0`, Entry) %>%
  rename(iedb_uniprot_id = 'yourlist:M201907166746803381A1F0E0DB47453E0216320D53922A0') %>%
  distinct()

# update names
# one entry Q8WYQ7,F8W9W4 will be split into two entres
iedb <- iedb %>%
  left_join(., unmatched_iedb_ids, by = c("uniprot_id" = "iedb_uniprot_id")) %>%
  mutate(uniprot_id = if_else(is.na(Entry), uniprot_id, Entry)) %>%
  select(-Entry) %>%
  separate_rows(uniprot_id, sep = ",", convert = FALSE) %>%
  distinct(uniprot_id, .keep_all = TRUE)

table(iedb$uniprot_id %in% uniprot$Entry)
```

    ## 
    ## FALSE  TRUE 
    ##   201  4376

``` r
# final list is 4376 
iedb_list <- iedb %>%
  filter(uniprot_id %in% uniprot$Entry) %>%
  select(uniprot_id, Antigen.Name) %>%
  rename(uniprotswsissprot = uniprot_id) %>%
  select(uniprotswsissprot) %>%
  write_tsv("data/1_anitgenic_iedb.tsv")
```
