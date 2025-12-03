# CREATED: 2025-04-11
# Gene enrichment




# SET UP WORKSPACE AND LOAD PACKAGES -------------------------------------
# DNAm object
library(minfi)

# Core
library(tidyverse)

# Enrichment
library(missMethyl)
# library(IlluminaHumanMethylation450kanno.ilmn12.hg19)




# LOAD DATA ----
#should have bacon corrected results
gout <- readr::read_rds("data_objects/Gout-25cov.rds")



# HIGHLY SIG ----

gout[gout$mix.qval.bac < 0.01 & abs(gout$mix.coef) > 350,] %>% 
  as.data.frame() %>% 
  arrange(desc(abs(mix.coef))) %>% 
  as.data.frame()

ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

ann[ann$Name == "cg05343316", "UCSC_RefGene_Name"]  #TESK2; hemorrage during pregnancy
ann[ann$Name == "cg14787155", "UCSC_RefGene_Name"]  #DZIP3
ann[ann$Name == "cg01994308", "UCSC_RefGene_Name"]  #PLAG1; placental growth; interacts with IGF2; oocytes
ann[ann$Name == "cg24927841", "UCSC_RefGene_Name"]  #na
ann[ann$Name == "cg08340572", "UCSC_RefGene_Name"]  #FBXO22; prenatal onset growth restriction
ann[ann$Name == "cg10342722", "UCSC_RefGene_Name"]  #SCPEP1; vascular health
ann[ann$Name == "cg14027161", "UCSC_RefGene_Name"]  #VAPA; fetal derived; placental function via vesicle trafficking
ann[ann$Name == "cg19015611", "UCSC_RefGene_Name"]  #VAX2
ann[ann$Name == "cg26327118", "UCSC_RefGene_Name"]  #KIF6
ann[ann$Name == "cg08377570", "UCSC_RefGene_Name"]  #LPAR6; early pregnancy; endometrium
ann[ann$Name == "cg21562809", "UCSC_RefGene_Name"]  #na
ann[ann$Name == "cg26647929", "UCSC_RefGene_Name"]  #na
ann[ann$Name == "cg15459822", "UCSC_RefGene_Name"]  #PARD3B6; suggests potential link with pregnancy complications
ann[ann$Name == "cg26300461", "UCSC_RefGene_Name"]  #na

## COMPARE HITS WITH TOP HITS FROM OTHER STUDIES??







# SETUP ----
# get significant CpGs and CpG universe
table(gout$mix.qval.bac < 0.05)
subverse <- names(gout[gout$mix.qval.bac < 0.01,])

# subverse.up <- names(Gout.dmp[Gout.dmp$mix.qval < 0.01 & Gout.dmp$mix.coef > 0,])
# subverse.down <- names(Gout.dmp[Gout.dmp$mix.qval < 0.01 & Gout.dmp$mix.coef < 0,])

universe <- names(gout)


check <- getMappedEntrezIDs(sig.cpg    = subverse,
                            all.cpg    = universe,
                            array.type = "450K")

names(check)   # see ?getMappedEntrezIDs
length(check$sig.eg)
length(check$universe)
check$freq
skimr::skim(check$freq)
check$equiv

# bias correction applied; check first CpG site
check$freq[1]
check$equiv[1]


# ENRICHMENT TEST ----
go.all <- gometh(sig.cpg    = subverse,
                 all.cpg    = universe,
                 collection = "GO",
                 prior.prob = TRUE,
                 plot.bias  = TRUE)
topGSA(go.all, n = 50)


go.all.sig <- gometh(sig.cpg    = subverse[1:500],
                 all.cpg    = universe,
                 collection = "GO",
                 prior.prob = TRUE,
                 plot.bias  = TRUE,
                 sig.genes = TRUE)





kegg.all <- gometh(sig.cpg    = subverse,
                   all.cpg    = universe,
                   collection = "KEGG",
                   prior.prob = TRUE,
                   plot.bias  = FALSE,
                   sig.genes  = TRUE)
topGSA(kegg.all, n = 50)


# promoter only ----
go.promoter <- gometh(sig.cpg    = subverse,
                      all.cpg    = universe,
                      collection = "GO",
                      prior.prob = TRUE,
                      plot.bias  = FALSE,
                      genomic.features = c("TSS1500","TSS200","1stExon"))
topGSA(go.promoter, n = 150)

kegg.promoter <- gometh(sig.cpg    = subverse,
                        all.cpg    = universe,
                        collection = "KEGG",
                        prior.prob = TRUE,
                        plot.bias  = FALSE,
                        genomic.features = c("TSS1500","TSS200","1stExon"))
topGSA(kegg.promoter, n = 50)


# gene body ----
go.body <- gometh(sig.cpg    = subverse,
                  all.cpg    = universe,
                  collection = "GO",
                  prior.prob = TRUE,
                  plot.bias  = FALSE,
                  genomic.features = c("Body"))
topGSA(go.body, n = 50)


kegg.body <- gometh(sig.cpg  = subverse,
                    all.cpg    = universe,
                    collection = "KEGG",
                    prior.prob = TRUE,
                    plot.bias  = FALSE,
                    genomic.features = c("Body"))
topGSA(kegg.body, n = 20)


# UTR ----
go.utr <- gometh(sig.cpg    = subverse,
                 all.cpg    = universe,
                 collection = "GO",
                 prior.prob = TRUE,
                 plot.bias  = FALSE,
                 genomic.features = c("3'UTR", "5'UTR"))
topGSA(go.utr, n = 30)

kegg.utr <- gometh(sig.cpg    = subverse,
                   all.cpg    = universe,
                   collection = "KEGG",
                   prior.prob = TRUE,
                   plot.bias  = FALSE,
                   genomic.features = c("3'UTR", "5'UTR"))
topGSA(kegg.utr, n = 50)



# OUTPUT ----

# Output your enrichment table from your *.R file
knitr::kable(topGSA(go.all, n = 20), 
             format = "html",
             digits = 3,
             align  = c("c", "l", rep("r", 4))) %>%
  kableExtra::save_kable(file = "enrichment_1.png",
                         bs_theme = "paper")

# Import the enrichment table in your *.Rmd file (within a chunk) using:
# knitr::include_graphics("table_1.png")




# ANNOTATION ----
ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# pick <- c("cg05541267", "cg02836478", "cg04014328", "cg05230392", "cg18608055", "cg05844798", "cg14472551", "cg02869235")

pick <-
  gout[gout$mix.qval.bac < 0.05,] %>% # & abs(gout$mix.stat) > 6,] %>% 
    as.data.frame() %>% 
    rownames()

ann[ann$Name %in% pick, c("Name", "chr", "pos", "UCSC_RefGene_Name")] %>% as.data.frame()

ann[ann$UCSC_RefGene_Name == "AVP", c("Name", "chr", "pos", "UCSC_RefGene_Name")]









