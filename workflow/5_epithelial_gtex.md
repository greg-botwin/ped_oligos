Epithelial Proteins
================
Translational Genomics Group
6/21/2019

GTEx Tissue
-----------

downloaded on 10/9/19 V8

Median TPM Per Tissue Per Gene
------------------------------

### Get Data

``` bash
wget -O workflow/data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz

gunzip worflow/data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz
```

### Summarise

``` r
gtex_tpm <- read_tsv("data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct", 
                       comment = "#", skip = 2)

gtex_tpm_int <- gtex_tpm %>%
  gather(key = "tissue", value = "TPM", -Name, - Description) %>%
  filter(tissue %in% c("Colon - Transverse", "Colon - Sigmoid",
                       "Small Intestine - Terminal Ileum")) %>%
  mutate(Log2TPM = log2(TPM + 1)) %>%
  mutate(Name = str_replace(Name,"\\..*",""))
```

### Map GTEx ENSG to UniProt

-   Human Proteom from UniProt (20,430 proteins)
-   ENSG to Uniprot from Ensembl <ftp://ftp.ensembl.org/pub/release-98/tsv/homo_sapiens/Homo_sapiens.GRCh38.98.uniprot.tsv.gz>
-   19,426 UniProt IDs found in mapping table
-   17,7701 UniProt IDs foung in gtex intestinal tisue.

``` r
uniprot <- read_tsv("data/uniprot-filtered-organism__Homo+sapiens+(Human)+[9606]_+AND+review--.tab")

ensembl_ids_38 <- read_tsv("data/Homo_sapiens.GRCh38.98.uniprot.tsv")

ensembl_ids_38 <- ensembl_ids_38 %>%
  filter(db_name == "Uniprot/SWISSPROT") %>%
  dplyr::select(gene_stable_id, xref) %>%
  rename(uniprotswissprot = xref) 

table(uniprot$Entry %in% ensembl_ids_38$uniprotswissprot)
```

    ## 
    ## FALSE  TRUE 
    ##  1004 19426

``` r
gtex_tpm_int <- gtex_tpm_int %>%
  left_join(ensembl_ids_38, by = c("Name" = "gene_stable_id"))
```

Multiple Gene can map to same protein, treat each ENSG as unique not averaging at protein level.

For Uniprot ID with multiple ENSG ids, expression is averaged

``` r
gtex_tpm_int <- gtex_tpm_int %>%
  filter(!is.na(uniprotswissprot)) %>%
  filter(Log2TPM > 0) %>%
  distinct() %>%
  group_by(tissue, uniprotswissprot) %>%
  summarise(avg_log2TPM_across_ensg = mean(Log2TPM)) %>%
  ungroup()
```

``` r
table(uniprot$Entry %in% gtex_tpm_int$uniprotswissprot)
```

    ## 
    ## FALSE  TRUE 
    ##  2729 17701

``` r
table(gtex_tpm_int$uniprotswissprot %in% uniprot$Entry)
```

    ## 
    ##  TRUE 
    ## 52163

add rank, rank each tissue uniquely
-----------------------------------

``` r
gtex_tpm_int <- gtex_tpm_int %>%
  group_by(tissue) %>%
  mutate(rank = dense_rank(desc(avg_log2TPM_across_ensg))) %>%
  ungroup()


write_tsv(gtex_tpm_int, "data/5_epithelial_gtex.tsv")
```
