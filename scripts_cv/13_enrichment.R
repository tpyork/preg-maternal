# CREATED: 2025-10-20
# Gene enrichment of VIP


## SEE OTHER VERSION IN `scripts_mwas/`


# SET UP WORKSPACE AND LOAD PACKAGES -------------------------------------
# DNAm object
library(minfi)

# Core
library(tidyverse)

# Enrichment
library(missMethyl)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)



# LOAD RESULTS ----

# Load All CpGs Considered
gout <- readr::read_rds("data_objects/dmp-timepoint_day-25cov.rds")


# Load Final Model Importance
vip <- readr::read_rds("data_objects_final/vip.rds")
vip[1:10, ]




# EXPLORATION ----
ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

ann[ann$Name == "cg22801980", "UCSC_RefGene_Name"]  # HNRNPH1; potential biomarker and therapeutic target in several diseases
ann[ann$Name == "cg16573266", "UCSC_RefGene_Name"]  # na
ann[ann$Name == "cg15465571", "UCSC_RefGene_Name"]  # ATP13A2
ann[ann$Name == "cg14847502", "UCSC_RefGene_Name"]  # MAP3K8; interacts with dexamethasone; Targeting MAP3K8 could enhance the therapeutic effects of dexamethasone
ann[ann$Name == "cg16748643", "UCSC_RefGene_Name"]  # na

ann[ann$Name == "cg00734512", "UCSC_RefGene_Name"]  # na
ann[ann$Name == "cg26684289", "UCSC_RefGene_Name"]  # na
ann[ann$Name == "cg02899473", "UCSC_RefGene_Name"]  # PARVA
ann[ann$Name == "cg12564453", "UCSC_RefGene_Name"]  # CETP
ann[ann$Name == "cg26333822", "UCSC_RefGene_Name"]  # SLC5A7

ann[ann$Name == "cg00951395", "UCSC_RefGene_Name"]  # KIAA1383; microtubule stability
ann[ann$Name == "cg27131510", "UCSC_RefGene_Name"]  # FAM24A; male fertility in mice
ann[ann$Name == "cg23894086", "UCSC_RefGene_Name"]  # na
ann[ann$Name == "cg08229615", "UCSC_RefGene_Name"]  # na
ann[ann$Name == "cg24566217", "UCSC_RefGene_Name"]  # ARHGEF12; expression in myometrium

ann[ann$Name == "cg18547574", "UCSC_RefGene_Name"]  # PDGFD; preeclampsia and HELLP biomarker
ann[ann$Name == "cg17240198", "UCSC_RefGene_Name"]  # NEK5; cancer biomarker
ann[ann$Name == "cg16990557", "UCSC_RefGene_Name"]  # na
ann[ann$Name == "cg14027161", "UCSC_RefGene_Name"]  # VAPA
ann[ann$Name == "cg20914508", "UCSC_RefGene_Name"]  # GAP43


skimr::skim(vip)
hist(vip$Importance)
table(vip$Importance > .0484) / 3000



# ENRICHMENT OF TOP IMPORTANCE CpGs ----
# get significant CpGs and CpG universe
vip2 <- vip %>% 
  filter(Importance > 0.0484) %>% 
  select(Variable) %>% 
  as.data.frame()

subverse <- vip2$Variable
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
topGSA(go.promoter, n = 50)

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










