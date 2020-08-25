# (A) 16S -----------------------------------------------------------------

# Section 1: Day-10 analysis ----------------------------------------------
setwd("P:/FMHSfiles/DIAMOND/Metadata/dada2")
library(permute)
library(lattice)
library("vegan"); packageVersion("vegan")
library("tidyverse"); packageVersion("tidyverse")
library("ggplot2"); packageVersion("ggplot2")
library(devtools)
library("phyloseq"); packageVersion("phyloseq")
library("RColorBrewer"); packageVersion("RColorBrewer")
library(dplyr)
library(ggpubr)
library(rstatix)

load("diamond_metadata_20200720.RData")

ASV_table <- read.csv("diamond_asv_table2.csv", header = TRUE)
rownames(ASV_table) <- ASV_table$Sample

#Select only Day 10
full_metadata3$timepoint %in% "Day 10"
metadata_D10 = full_metadata3[full_metadata3$timepoint %in% "Day 10" , ]
asv_table_D10 = ASV_table[(metadata_D10$Sample) , ]
all(metadata_D10$Sample == rownames(asv_table_D10))

asv_table_D10$Sample <- NULL


# 1) Adonis-General factors (without birth weight) ------------------------
a <- adonis2(asv_table_D10 ~ gestational_wk + Hospital_Name + hospital_duration + Sex + delivery_mode + plural, 
             data = metadata_D10, permutations = 10000, method = "bray", by = "margin")

a

pvalue  <- c()
pvalue <- c(pvalue, a$`Pr(>F)`[1:6])
pvalue

corrected_p_values <- p.adjust(pvalue, method = "BH")
corrected_p_values



# 2) Adonis - Maternal factors --------------------------------------------
b <- adonis2(asv_table_D10 ~ maternal_age + edu + maternal_ethnicity + SES + maternal_stress_dc_PSS + maternal_Depression_dc_EPDS + 
               M_Q5c_AntibMoL4Wk + M_Q5d_ProbMoL4W, data = metadata_D10, permutations = 10000, method = "bray", by = "margin")

metadata_D10_no_NA <- metadata_D10 %>% filter(!is.na(metadata_D10$SES) & (!is.na(metadata_D10$M_Q5c_AntibMoL4Wk)) & (!is.na(metadata_D10$maternal_stress_dc_PSS)) & (!is.na(metadata_D10$maternal_Depression_dc_EPDS)) & (!is.na(metadata_D10$M_Q5d_ProbMoL4W)))

asv_table_D10[ as.character(metadata_D10_no_NA$Sample) , ]
asv_table_D10_no_NA = asv_table_D10[(metadata_D10_no_NA$Sample) , ]

all(metadata_D10_no_NA$Sample == rownames(asv_table_D10_no_NA))

b <- adonis2(asv_table_D10_no_NA ~ maternal_age + edu + maternal_ethnicity + SES + maternal_stress_dc_PSS + maternal_Depression_dc_EPDS + 
               M_Q5c_AntibMoL4Wk + M_Q5d_ProbMoL4W, data = metadata_D10_no_NA, permutations = 10000, method = "bray", by = "margin")


b

pvalue  <- c()
pvalue <- c(pvalue, b$`Pr(>F)`[1:8])
pvalue

corrected_p_values <- p.adjust(pvalue, method = "BH")
corrected_p_values


# 3) Infants nutritional factors - intravenous nutrition, milk typ --------
c <- adonis2(asv_table_D10 ~ feeding_d10_v3 + nutrition_v2 + total_smell_taste, 
             data = metadata_D10, permutations = 10000, method = "bray", by = "margin")

c

pvalue  <- c()
pvalue <- c(pvalue, c$`Pr(>F)`[1:3])
pvalue

corrected_p_values <- p.adjust(pvalue, method = "BH")
corrected_p_values

###without smell/taste

c <- adonis2(asv_table_D10 ~ feeding_d10_v3 + nutrition_v2 +  
               G_Q5_Probio + M_Q5b_AntibBby + M_Q5c_AntibMoL4Wk + M_Q5d_ProbMoL4W, 
             data = metadata_D10, permutations = 10000, method = "bray", by = "margin")

metadata_D10_no_NA <- metadata_D10 %>% filter(!is.na(metadata_D10$M_Q5d_ProbMoL4W) & !is.na(metadata_D10$M_Q5b_AntibBby))

asv_table_D10[ as.character(metadata_D10_no_NA$Sample) , ]
asv_table_D10_no_NA = asv_table_D10[(metadata_D10_no_NA$Sample) , ]

all(metadata_D10_no_NA$Sample == rownames(asv_table_D10_no_NA))
asv_table_D10_no_NA$Sample <- NULL

metadata_D10_no_NA$Sample == rownames(metadata_D10_no_NA)

c <- adonis2(asv_table_D10_no_NA ~ feeding_d10_v3 + nutrition_v2 + 
               G_Q5_Probio + M_Q5b_AntibBby + M_Q5c_AntibMoL4Wk + M_Q5d_ProbMoL4W, 
             data = metadata_D10_no_NA, permutations = 10000, method = "bray", by = "margin")

c

pvalue  <- c()
pvalue <- c(pvalue, c$`Pr(>F)`[1:6])
pvalue

corrected_p_values <- p.adjust(pvalue, method = "BH")
corrected_p_values
```

###with smell/taste

c <- adonis2(asv_table_D10 ~ feeding_d10_v3 + nutrition_v2 + total_smell_taste + 
               G_Q5_Probio + M_Q5b_AntibBby + M_Q5c_AntibMoL4Wk + M_Q5d_ProbMoL4W, 
             data = metadata_D10, permutations = 10000, method = "bray", by = "margin")

metadata_D10_no_NA <- metadata_D10 %>% filter(!is.na(metadata_D10$M_Q5d_ProbMoL4W) & !is.na(metadata_D10$M_Q5b_AntibBby))

asv_table_D10[ as.character(metadata_D10_no_NA$Sample) , ]
asv_table_D10_no_NA = asv_table_D10[(metadata_D10_no_NA$Sample) , ]

all(metadata_D10_no_NA$Sample == rownames(asv_table_D10_no_NA))
asv_table_D10_no_NA$Sample <- NULL

metadata_D10_no_NA$Sample == rownames(metadata_D10_no_NA)

c <- adonis2(asv_table_D10_no_NA ~ feeding_d10_v3 + nutrition_v2 + total_smell_taste + 
               G_Q5_Probio + M_Q5b_AntibBby + M_Q5c_AntibMoL4Wk + M_Q5d_ProbMoL4W, 
             data = metadata_D10_no_NA, permutations = 10000, method = "bray", by = "margin")
c

pvalue  <- c()
pvalue <- c(pvalue, c$`Pr(>F)`[1:7])
pvalue

corrected_p_values <- p.adjust(pvalue, method = "BH")
corrected_p_values


# 4) Maaslin2 - all Day-10 factors ----------------------------------------
library(tidyverse)
library("devtools")
library(microbiomics)
library(Maaslin2)

set.seed(10000)

rownames(metadata_D10) <- metadata_D10$Sample

###no discharge feeding
maaslin <- 
  Maaslin2(input_data = asv_table_D10, 
           input_metadata = metadata_D10,
           output = "maaslin2_D10_20200722",
           fixed_effects = c("gestational_wk", "Hospital_Name", "hospital_duration", "Sex", "delivery_mode", "birth_wt_gm_cont", "plural", 
                             "nutrition_v2", "total_smell_taste", "G_Q5_Probio", "feeding_d10_v3",
                             "maternal_age", "edu", "maternal_ethnicity", "SES", "maternal_stress_dc_PSS", "maternal_Depression_dc_EPDS", 
                             "M_Q5c_AntibMoL4Wk", "M_Q5d_ProbMoL4W", "G_Q7_Milk_Dis", "M_Q5b_AntibBby"),
           random_effects = NULL,
           plot_scatter = F,
           min_prevalence = 0.1,
           cores = 4)


# Section 2: 4-month analysis ---------------------------------------------
setwd("P:/FMHSfiles/DIAMOND/Metadata/dada2")
library(permute)
library(lattice)
library("vegan"); packageVersion("vegan")
library("tidyverse"); packageVersion("tidyverse")
library("ggplot2"); packageVersion("ggplot2")
library(devtools)
library("phyloseq"); packageVersion("phyloseq")
library("RColorBrewer"); packageVersion("RColorBrewer")
library(dplyr)
library(ggpubr)
library(rstatix)

load("diamond_metadata_20200720.RData")

ASV_table <- read.csv("diamond_asv_table2.csv", header = TRUE)
rownames(ASV_table) <- ASV_table$Sample

###Select only 4m
full_metadata3$timepoint %in% "4 month"
metadata_4m = full_metadata3[full_metadata3$timepoint %in% "4 month" , ]
asv_table_4m = ASV_table[(metadata_4m$Sample) , ]
all(metadata_4m$Sample == rownames(asv_table_4m))

asv_table_4m$Sample <- NULL


# 1) General factors  -----------------------------------------------------
### with birth weight

a <- adonis2(asv_table_4m ~ gestational_wk + Hospital_Name + hospital_duration + Sex + delivery_mode + birth_wt_gm_cont + plural, data = metadata_4m, permutations = 10000, method = "bray", by = "margin")

a

pvalues <- c()
pvalues <- c(pvalues, a$`Pr(>F)`[1:7])
pvalues

corrected_p_values <- p.adjust(pvalues, method = "BH")
corrected_p_values

###without birth weight

a <- adonis2(asv_table_4m ~ gestational_wk + Hospital_Name + hospital_duration + Sex + delivery_mode + plural, data = metadata_4m, permutations = 10000, method = "bray", by = "margin")

a

pvalues <- c()
pvalues <- c(pvalues, a$`Pr(>F)`[1:6])
pvalues

corrected_p_values <- p.adjust(pvalues, method = "BH")
corrected_p_values


# 2) Maternal factors -----------------------------------------------------

b <- adonis2(asv_table_4m ~ maternal_age + edu + maternal_ethnicity + SES + maternal_stress_4month_PSS + maternal_Depression_4month_EPDS + 
               K_Q22d_anti_mother_4M + K_Q22e_prob_mother_4M, data = metadata_4m, permutations = 10000, method = "bray", by = "margin")
b

metadata_4m_no_NA <- metadata_4m %>% filter(!is.na(metadata_4m$SES) & (!is.na(metadata_4m$K_Q22d_anti_mother_4M)) & (!is.na(metadata_4m$K_Q22e_prob_mother_4M)))

asv_table_4m[ as.character(metadata_4m_no_NA$Sample) , ]
asv_table_4m_no_NA = asv_table_4m[(metadata_4m_no_NA$Sample) , ]

all(metadata_4m_no_NA$Sample == rownames(asv_table_4m_no_NA))

b <- adonis2(asv_table_4m_no_NA ~ maternal_age + edu + maternal_ethnicity + SES + maternal_stress_4month_PSS + maternal_Depression_4month_EPDS + 
               K_Q22d_anti_mother_4M + K_Q22e_prob_mother_4M, data = metadata_4m_no_NA, permutations = 10000, method = "bray", by = "margin")

b

pvalues <- c()
pvalues <- c(pvalues, b$`Pr(>F)`[1:8])
pvalues

corrected_p_values <- p.adjust(pvalues, method = "BH")
corrected_p_values


# 3) Nutrition factors (4m) -----------------------------------------------

c <- adonis2(asv_table_4m ~ feeding + K_Q22b_anti_baby_4M + 
               K_Q22c_prob_baby_4M + K_Q22d_anti_mother_4M + K_Q22e_prob_mother_4M, data = metadata_4m, permutations = 10000, method = "bray", by = "margin")
c

metadata_4m_no_NA <- metadata_4m %>% filter(!is.na(metadata_4m$start_solid_4month_ca) & (!is.na(metadata_4m$any_formula_4month_ca)) & (!is.na(metadata_4m$K_Q22b_anti_baby_4M)) & (!is.na(metadata_4m$K_Q22c_prob_baby_4M)) & (!is.na(metadata_4m$K_Q22d_anti_mother_4M)) & (!is.na(metadata_4m$K_Q22d_anti_mother_4M)))

asv_table_4m[ as.character(metadata_4m_no_NA$Sample) , ]
asv_table_4m_no_NA = asv_table_4m[(metadata_4m_no_NA$Sample) , ]

all(metadata_4m_no_NA$Sample == rownames(asv_table_4m_no_NA))

c <- adonis2(asv_table_4m_no_NA ~ feeding + K_Q22b_anti_baby_4M + 
               K_Q22c_prob_baby_4M + K_Q22d_anti_mother_4M + K_Q22e_prob_mother_4M, data = metadata_4m_no_NA, permutations = 10000, method = "bray", by = "margin")

c

pvalues <- c()
pvalues <- c(pvalues, c$`Pr(>F)`[1:5])
pvalues

corrected_p_values <- p.adjust(pvalues, method = "BH")
corrected_p_values

###early-life factors

a <- adonis2(asv_table_4m ~ feeding_d10_v3 + nutrition_v2 + Hospital_Name + gestational_wk + delivery_mode + G_Q5_Probio, 
             data = metadata_4m, permutations = 10000, method = "bray", by = "margin")
a

metadata_4m_no_NA <- metadata_4m %>% filter(!is.na(metadata_4m$feeding_d10_v3))

asv_table_4m[ as.character(metadata_4m_no_NA$Sample) , ]
asv_table_4m_no_NA = asv_table_4m[(metadata_4m_no_NA$Sample) , ]

all(metadata_4m_no_NA$Sample == rownames(asv_table_4m_no_NA))
asv_table_4m_no_NA$Sample <- NULL

a <- adonis2(asv_table_4m_no_NA ~ feeding_d10_v3 + nutrition_v2 + Hospital_Name + gestational_wk + delivery_mode + G_Q5_Probio, 
             data = metadata_4m_no_NA, permutations = 10000, method = "bray", by = "margin")

a

pvalues <- c()
pvalues <- c(pvalues, a$`Pr(>F)`[1:6])
pvalues

corrected_p_values <- p.adjust(pvalues, method = "BH")
corrected_p_values

###early-life nutrition without smell/taste
```{r}
a <- adonis2(asv_table_4m ~ feeding_d10_v3+ nutrition_v2 + 
               G_Q5_Probio + M_Q5b_AntibBby + M_Q5c_AntibMoL4Wk + M_Q5d_ProbMoL4W, 
               data = metadata_4m, permutations = 10000, method = "bray", by = "margin")
a

metadata_4m_no_NA <- metadata_4m %>% filter(!is.na(metadata_4m$M_Q5d_ProbMoL4W) & !is.na(metadata_4m$M_Q5b_AntibBby) & !is.na(metadata_4m$feeding_d10_v3))

asv_table_4m[ as.character(metadata_4m_no_NA$Sample) , ]
asv_table_4m_no_NA = asv_table_4m[(metadata_4m_no_NA$Sample) , ]

all(metadata_4m_no_NA$Sample == rownames(asv_table_4m_no_NA))
asv_table_4m_no_NA$Sample <- NULL

metadata_4m_no_NA$Sample == rownames(metadata_4m_no_NA)

a <- adonis2(asv_table_4m_no_NA ~ feeding_d10_v3 + nutrition_v2+ 
               G_Q5_Probio + M_Q5b_AntibBby + M_Q5c_AntibMoL4Wk + M_Q5d_ProbMoL4W, 
             data = metadata_4m_no_NA, permutations = 10000, method = "bray", by = "margin")

a

pvalues <- c()
pvalues <- c(pvalues, a$`Pr(>F)`[1:6])
pvalues

corrected_p_values <- p.adjust(pvalues, method = "BH")
corrected_p_values

###early-life with smell/taste

a <- adonis2(asv_table_4m ~ feeding_d10_v3+ nutrition_v2 + total_smell_taste + 
               G_Q5_Probio + M_Q5b_AntibBby + M_Q5c_AntibMoL4Wk + M_Q5d_ProbMoL4W, 
             data = metadata_4m, permutations = 10000, method = "bray", by = "margin")
a

metadata_4m_no_NA <- metadata_4m %>% filter(!is.na(metadata_4m$M_Q5d_ProbMoL4W) & !is.na(metadata_4m$M_Q5b_AntibBby) & !is.na(metadata_4m$feeding_d10_v3))

asv_table_4m[ as.character(metadata_4m_no_NA$Sample) , ]
asv_table_4m_no_NA = asv_table_4m[(metadata_4m_no_NA$Sample) , ]

all(metadata_4m_no_NA$Sample == rownames(asv_table_4m_no_NA))
asv_table_4m_no_NA$Sample <- NULL

a <- adonis2(asv_table_4m_no_NA ~ feeding_d10_v3 + nutrition_v2 + total_smell_taste + 
               G_Q5_Probio + M_Q5b_AntibBby + M_Q5c_AntibMoL4Wk + M_Q5d_ProbMoL4W, 
             data = metadata_4m_no_NA, permutations = 10000, method = "bray", by = "margin")

a

pvalues <- c()
pvalues <- c(pvalues, a$`Pr(>F)`[1:7])
pvalues

corrected_p_values <- p.adjust(pvalues, method = "BH")
corrected_p_values


# 4) Maaslin2 - all 4M factors --------------------------------------------

library(tidyverse)
library("devtools")
library(microbiomics)
library(Maaslin2)

set.seed(10000)

rownames(metadata_4m) <- metadata_4m$Sample

maaslin <- 
  Maaslin2(input_data = asv_table_4m, 
           input_metadata = metadata_4m,
           output = "maaslin2_4M_overall_20200722",
           fixed_effects = c("gestational_wk", "delivery_mode", "Hospital_Name", "feeding",
                             "hospital_duration", "Sex", "birth_wt_gm_cont", "plural",
                             "K_Q22b_anti_baby_4M", "K_Q22c_prob_baby_4M",  
                             "K_Q22d_anti_mother_4M", "K_Q22e_prob_mother_4M", "maternal_age", 
                             "edu", "maternal_ethnicity", "SES", "maternal_stress_4month_PSS", 
                             "maternal_Depression_4month_EPDS"),
           random_effects = NULL,
           plot_scatter = F,
           min_prevalence = 0.1,
           cores = 4)



# Section 3: RQ2: Longitudinal changes ------------------------------------
setwd("P:/FMHSfiles/DIAMOND/Metadata/dada2")
library(permute)
library(lattice)
library("vegan"); packageVersion("vegan")
library("tidyverse"); packageVersion("tidyverse")
library("ggplot2"); packageVersion("ggplot2")
library(devtools)
library("phyloseq"); packageVersion("phyloseq")
library("RColorBrewer"); packageVersion("RColorBrewer")
library(dplyr)
library(ggpubr)
library(rstatix)

load("diamond_metadata_20200720.RData")

ASV_table <- read.csv("diamond_asv_table2.csv", header = TRUE)
rownames(ASV_table) <- ASV_table$Sample
###select d10 and 4m paired

metadata_paired_d10_4m = full_metadata3 %>% group_by(Study_ID) %>% summarise(n_samples = n()) %>% filter(n_samples == 2) %>% left_join(full_metadata3)

asv_table_d10_4m = ASV_table[(metadata_paired_d10_4m$Sample) , ]
all(metadata_paired_d10_4m$Sample == rownames(asv_table_d10_4m))

asv_table_d10_4m$Sample <- NULL

a <- adonis2(asv_table_d10_4m ~ timepoint + Hospital_Name + gestational_wk + delivery_mode + G_Q5_Probio + 
               plural + hospital_duration + feeding + feeding_d10_v3 + nutrition_v2,
             data = metadata_paired_d10_4m, permutations = 10000, method = "bray", by = "margin")
a

metadata_paired_d10_4m_no_NA <- metadata_paired_d10_4m %>% filter(!is.na(metadata_paired_d10_4m$feeding))

asv_table_d10_4m[ as.character(metadata_paired_d10_4m_no_NA$Sample) , ]
asv_table_d10_4m_no_NA = asv_table_d10_4m[(metadata_paired_d10_4m_no_NA$Sample) , ]

all(metadata_paired_d10_4m_no_NA$Sample == rownames(asv_table_d10_4m_no_NA))

a <- adonis2(asv_table_d10_4m_no_NA ~ timepoint + Hospital_Name + gestational_wk + delivery_mode + G_Q5_Probio + 
               plural + hospital_duration + feeding + feeding_d10_v3 + nutrition_v2,
             data = metadata_paired_d10_4m_no_NA, permutations = 10000, method = "bray", by = "margin")
a


pvalues <- c()
pvalues <- c(pvalues, a$`Pr(>F)`[1:10])
pvalues

corrected_p_values <- p.adjust(pvalues, method = "BH")
corrected_p_values

###Maaslin
library(tidyverse)
library("devtools")
library(microbiomics)
library(Maaslin2)

set.seed(10000)

rownames(metadata_paired_d10_4m) <- metadata_paired_d10_4m$Sample
maaslin <- 
  Maaslin2(input_data = asv_table_d10_4m, 
           input_metadata = metadata_paired_d10_4m,
           output = "maaslin2_D104m_RQ2_20200722(2)",
           fixed_effects = c("timepoint", "Hospital_Name", "gestational_wk", "delivery_mode", "G_Q5_Probio", "plural", "hospital_duration", "nutrition_v2", "feeding", "feeding_d10_v3"),
           random_effects = NULL,
           plot_scatter = F,
           min_prevalence = 0.1,
           cores = 4)

#correction on 14/8/2020
rownames(metadata_paired_d10_4m) <- metadata_paired_d10_4m$Sample
maaslin <- 
  Maaslin2(input_data = asv_table_d10_4m, 
           input_metadata = metadata_paired_d10_4m,
           output = "maaslin2_D104m_RQ2_20200814(3)",
           fixed_effects = c("timepoint", "Hospital_Name", "gestational_wk", "delivery_mode", "G_Q5_Probio", "plural", "hospital_duration", "nutrition_v2", "feeding", "feeding_d10_v3"),
           random_effects = "Study_ID",
           plot_scatter = F,
           min_prevalence = 0.1,
           cores = 4)


# RQ3: Sex-specific effect ------------------------------------------------
setwd("P:/FMHSfiles/DIAMOND/Metadata/dada2")
library(permute)
library(lattice)
library("vegan"); packageVersion("vegan")
library("tidyverse"); packageVersion("tidyverse")
library("ggplot2"); packageVersion("ggplot2")
library(devtools)
library("phyloseq"); packageVersion("phyloseq")
library("RColorBrewer"); packageVersion("RColorBrewer")
library(dplyr)
library(ggpubr)
library(rstatix)

load("diamond_metadata_20200720.RData")

ASV_table <- read.csv("diamond_asv_table2.csv", header = TRUE)
rownames(ASV_table) <- ASV_table$Sample

#D10
###Select only Day 10
full_metadata3$timepoint %in% "Day 10"
metadata_D10 = full_metadata3[full_metadata3$timepoint %in% "Day 10" , ]
asv_table_D10 = ASV_table[(metadata_D10$Sample) , ]
all(metadata_D10$Sample == rownames(asv_table_D10))

asv_table_D10$Sample <- NULL

b <- adonis2(asv_table_D10 ~ Sex, 
             data = metadata_D10, permutations = 10000, method = "bray", by = "margin")


b

pvalues <- c()
pvalues <- c(pvalues, b$`Pr(>F)`[1:1])
pvalues

corrected_p_values <- p.adjust(pvalues, method = "BH")
corrected_p_values

#4M
###Select only 4m
full_metadata3$timepoint %in% "4 month"
metadata_4m = full_metadata3[full_metadata3$timepoint %in% "4 month" , ]
asv_table_4m = ASV_table[(metadata_4m$Sample) , ]
all(metadata_4m$Sample == rownames(asv_table_4m))

asv_table_4m$Sample <- NULL

c <- adonis2(asv_table_4m ~ Sex, 
             data = metadata_4m, permutations = 10000, method = "bray", by = "margin")

c

pvalues <- c()
pvalues <- c(pvalues, c$`Pr(>F)`[1:1])
pvalues

corrected_p_values <- p.adjust(pvalues, method = "BH")
corrected_p_values


# Section 4: Shannon ------------------------------------------------------
setwd("P:/FMHSfiles/DIAMOND/Metadata/dada2")
library(permute)
library(lattice)
library("vegan"); packageVersion("vegan")
library("tidyverse"); packageVersion("tidyverse")
library("ggplot2"); packageVersion("ggplot2")
library(devtools)
library("phyloseq"); packageVersion("phyloseq")
library("RColorBrewer"); packageVersion("RColorBrewer")
library(dplyr)
library(ggpubr)
library(rstatix)

load("diamond_metadata_20200720.RData")

shannon <- read_tsv("shannon_diamond_20200119.tsv")

shannon_full <- shannon %>%
  left_join(full_metadata3)


###Select only Day 10 - wilcoxon
shannon_full$timepoint %in% "Day 10"
shannon_D10 = shannon_full[shannon_full$timepoint %in% "Day 10" , ]

a <- wilcox_test(shannon_D10, formula = Shannon ~ Sex, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
a

b <- wilcox_test(shannon_D10, formula = Shannon ~ delivery_mode, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
b

c <- wilcox_test(shannon_D10, formula = Shannon ~ gestational_wk, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
c

d <- wilcox_test(shannon_D10, formula = Shannon ~ edu, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
d

e <- wilcox_test(shannon_D10, formula = Shannon ~ nutrition_v2, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
e

f <- wilcox_test(shannon_D10, formula = Shannon ~ total_smell_taste, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
f

g <- wilcox_test(shannon_D10, formula = Shannon ~ maternal_Depression_dc_EPDS, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
g

pvalues <- c()
pvalues <- c(pvalues, a$p)
pvalues <- c(pvalues, b$p)
pvalues <- c(pvalues, c$p)
pvalues <- c(pvalues, d$p)
pvalues <- c(pvalues, e$p)
pvalues <- c(pvalues, f$p)
pvalues <- c(pvalues, g$p)
pvalues

corrected_p_values <- p.adjust(pvalues, method = "BH")
corrected_p_values

h <- wilcox_test(shannon_D10, formula = Shannon ~ M_Q5c_AntibMoL4Wk, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
h

i <- wilcox_test(shannon_D10, formula = Shannon ~ M_Q5d_ProbMoL4W, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
i

j <- wilcox_test(shannon_D10, formula = Shannon ~ M_Q5b_AntibBby, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
j

k <- wilcox_test(shannon_D10, formula = Shannon ~ G_Q5_Probio, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
k

pvalues <- c()
pvalues <- c(pvalues, h$p)
pvalues <- c(pvalues, i$p)
pvalues <- c(pvalues, j$p)
pvalues <- c(pvalues, k$p)
pvalues

corrected_p_values <- p.adjust(pvalues, method = "BH")
corrected_p_values

###Kruskal & Wilcoxon 

a <- kruskal.test(Shannon ~ Hospital_Name, data = shannon_D10)
a

pvalues <- c()
pvalues <- c(pvalues, a$p.value)
pvalues

b <- wilcox_test(shannon_D10, formula = Shannon ~ Hospital_Name, comparisons = list(c("Auckland","Middlemore"),c("Auckland","NorthShore"),c("Auckland","Waitakere"),c("Middlemore","NorthShore"),c("Middlemore","Waitakere"),c("NorthShore","Waitakere")), paired = FALSE) %>%
  adjust_pvalue(method = "BH")

b

c <- kruskal.test(Shannon ~ maternal_ethnicity, data = shannon_D10)
c

pvalues <- c(pvalues, c$p.value)
pvalues

d <- wilcox_test(shannon_D10, formula = Shannon ~ maternal_ethnicity, comparisons = list(c("Asian","Caucasian"),c("Asian","Maori"),c("Asian","Pacific People"),c("Asian","Other"),c("Caucasian","Maori"),c("Pacific People","Maori"),c("Pacific People","Other"),c("Other","Maori"),c("Caucasian","Pacific People"),c("Caucasian","Other")), paired = FALSE) %>% adjust_pvalue(method = "BH")
d

e <- kruskal.test(Shannon ~ SES, data = shannon_D10)
e

pvalues <- c(pvalues, e$p.value)
pvalues

f <- wilcox_test(shannon_D10, formula = Shannon ~ SES, comparisons = list(c("Least deprived","Most deprived"),c("Least deprived","Moderate"),c("Most deprived","Moderate")), paired = FALSE) %>%
  adjust_pvalue(method = "BH")
f

g <- kruskal.test(Shannon ~ maternal_stress_dc_PSS, data = shannon_D10)
g

pvalues <- c(pvalues, g$p.value)
pvalues

h <- wilcox_test(shannon_D10, formula = Shannon ~ maternal_stress_dc_PSS, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
h

i <- kruskal.test(Shannon ~ plural, data = shannon_D10)
i

pvalues <- c(pvalues, i$p.value)
pvalues

j <- wilcox_test(shannon_D10, formula = Shannon ~ plural, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
j

k <- kruskal.test(Shannon ~ feeding_d10_v3, data = shannon_D10)
k

pvalues <- c(pvalues, k$p.value)
pvalues

corrected_p_values <- p.adjust(pvalues, method = "BH")
corrected_p_values

i <- wilcox_test(shannon_D10, formula = Shannon ~ feeding_d10_v3, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
i


###Select only 4 month
shannon_4m = shannon_full[shannon_full$timepoint %in% "4 month" , ]


a <- wilcox_test(shannon_4m, formula = Shannon ~ Sex, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
a

pvalues <- c()
pvalues <- c(pvalues, a$p)
pvalues

b <- wilcox_test(shannon_4m, formula = Shannon ~ edu, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
b

pvalues <- c(pvalues, b$p)

corrected_p_values <- p.adjust(pvalues, method = "BH")
corrected_p_values

a <- wilcox_test(shannon_4m, formula = Shannon ~ gestational_wk, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
a

pvalues <- c()
pvalues <- c(pvalues, a$p)
pvalues

b <- wilcox_test(shannon_4m, formula = Shannon ~ delivery_mode, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
b

pvalues <- c(pvalues, b$p)
pvalues

corrected_p_values <- p.adjust(pvalues, method = "BH")
corrected_p_values

c <- wilcox_test(shannon_4m, formula = Shannon ~ nutrition_v2, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
c

pvalues <- c()
pvalues <- c(pvalues, c$p)
pvalues

d <- wilcox_test(shannon_4m, formula = Shannon ~ total_smell_taste, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
d

pvalues <- c(pvalues, d$p)
pvalues

e <- wilcox_test(shannon_4m, formula = Shannon ~ maternal_Depression_4month_EPDS, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
e

pvalues <- c(pvalues, e$p)
pvalues

corrected_p_values <- p.adjust(pvalues, method = "BH")
corrected_p_values

f <- wilcox_test(shannon_4m, formula = Shannon ~ K_Q22d_anti_mother_4M, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
f

pvalues <- c()
pvalues <- c(pvalues, f$p)
pvalues

g <- wilcox_test(shannon_4m, formula = Shannon ~ K_Q22e_prob_mother_4M, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
g

pvalues <- c(pvalues, g$p)
pvalues

h <- wilcox_test(shannon_4m, formula = Shannon ~ K_Q22b_anti_baby_4M, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
h

pvalues <- c(pvalues, h$p)
pvalues

i <- wilcox_test(shannon_4m, formula = Shannon ~ K_Q22c_prob_baby_4M, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
i

pvalues <- c(pvalues, i$p)
pvalues

j <- wilcox_test(shannon_4m, formula = Shannon ~ plural, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
j

pvalues <- c(pvalues, j$p)
pvalues

corrected_p_values <- p.adjust(pvalues, method = "BH")
corrected_p_values

###Kruskal-wallis

a <- kruskal.test(Shannon ~ Hospital_Name, data = shannon_4m)
a

pvalues <- c()
pvalues <- c(pvalues, a$p.value)
pvalues

b <- wilcox_test(shannon_4m, formula = Shannon ~ Hospital_Name, comparisons = list(c("Auckland","Middlemore"),c("Auckland","NorthShore"),c("Auckland","Waitakere"),c("Middlemore","NorthShore"),c("Middlemore","Waitakere"),c("NorthShore","Waitakere")), paired = FALSE) %>%
  adjust_pvalue(method = "BH")

b

c <- kruskal.test(Shannon ~ maternal_ethnicity, data = shannon_4m)
c

pvalues <- c(pvalues, c$p.value)
pvalues

d <- wilcox_test(shannon_4m, formula = Shannon ~ maternal_ethnicity, comparisons = list(c("Asian","Caucasian"),c("Asian","Maori"),c("Asian","Pacific People"),c("Asian","Other"),c("Caucasian","Maori"),c("Pacific People","Maori"),c("Pacific People","Other"),c("Other","Maori"),c("Caucasian","Pacific People"),c("Caucasian","Other")), paired = FALSE) %>% adjust_pvalue(method = "BH")
d


e <- kruskal.test(Shannon ~ SES, data = shannon_4m)
e

pvalues <- c(pvalues, e$p.value)
pvalues

f <- wilcox_test(shannon_4m, formula = Shannon ~ SES, comparisons = list(c("Least deprived","Most deprived"),c("Least deprived","Moderate"),c("Most deprived","Moderate")), paired = FALSE) %>%
  adjust_pvalue(method = "BH")
f

g <- kruskal.test(Shannon ~ maternal_stress_4month_PSS, data = shannon_4m)
g

pvalues <- c(pvalues, g$p.value)
pvalues

h <- wilcox_test(shannon_4m, formula = Shannon ~ maternal_stress_4month_PSS, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
h

i <- kruskal.test(Shannon ~ feeding_d10_v3, data = shannon_4m)
i

pvalues <- c(pvalues, i$p.value)
pvalues

corrected_p_values <- p.adjust(pvalues, method = "BH")
corrected_p_values

j <- wilcox_test(shannon_4m, formula = Shannon ~ feeding_d10_v3, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
j

a <- kruskal.test(Shannon ~ feeding, data = shannon_4m)
a

b <- wilcox_test(shannon_4m, formula = Shannon ~ feeding, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
b

c <- wilcox_test(shannon_full, formula = Shannon ~ timepoint, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
c

###Longitudinal - Select for paired
metadata_paired_d10_4m = shannon_full %>% group_by(Study_ID) %>% summarise(n_samples = n()) %>% filter(n_samples == 2) %>% left_join(shannon_full)
```{r}
d <- wilcox_test(metadata_paired_d10_4m, Shannon ~ timepoint , paired = FALSE) %>%
  adjust_pvalue(method = "BH")
d

pvalues <- c()
pvalues <- c(pvalues, a$p.value)
pvalues <- c(pvalues, c$p)
pvalues <- c(pvalues, d$p)
pvalues

corrected_p_values <- p.adjust(pvalues, method = "BH")
corrected_p_values


# Section 5: Growth rate --------------------------------------------------

setwd("P:/FMHSfiles/DIAMOND/Metadata/dada2")
library(permute)
library(lattice)
library("vegan"); packageVersion("vegan")
library("tidyverse"); packageVersion("tidyverse")
library("ggplot2"); packageVersion("ggplot2")
library(devtools)
library("phyloseq"); packageVersion("phyloseq")
library("RColorBrewer"); packageVersion("RColorBrewer")
library(dplyr)
library(ggpubr)
library(rstatix)

load("diamond_metadata_20200720.RData")

ASV_table <- read.csv("diamond_asv_table2.csv", header = TRUE)
rownames(ASV_table) <- ASV_table$Sample

###Selecet for 4m 
metadata_4m = full_metadata3[full_metadata3$timepoint == "4 month" , ]
asv_table_4m <- ASV_table[ASV_table$Sample %in% metadata_4m$Sample , ]
all(metadata_4m$Sample == rownames(asv_table_4m))

###WT - z-score
metadata_4m <- 
  metadata_4m %>%
  mutate(wt_growth_z = (Mth4_B_Gth_WT_z - birth_wt_z)) %>% 
  full_join(metadata_4m)

###LH - z-score
metadata_4m <- 
  metadata_4m%>%
  mutate(lh_growth_z = (Mth4_B_Gth_LH_z - birth_length_z)) %>% 
  full_join(metadata_4m)

###HC - z-score
metadata_4m <- 
  metadata_4m %>%
  mutate(hc_growth_z = (Mth4_B_Gth_HC_z - birth_hc_z)) %>% 
  full_join(metadata_4m)

rownames(metadata_4m) <- metadata_4m$Sample


###WT - rate
metadata_4m <- 
  metadata_4m %>%
  mutate(wt_growth_rate = ((Mth4_B_Gth_WT_g - birth_wt_gm_cont)/stoolcol_4m_day)) %>% 
  full_join(metadata_4m)

###LH - rate
metadata_4m <- 
  metadata_4m %>%
  mutate(lh_growth_rate = ((Mth4_B_Gth_LH_cm - birth_length_cm)/stoolcol_4m_day)) %>% 
  full_join(metadata_4m)

###HC - rate
metadata_4m <- 
  metadata_4m %>%
  mutate(hc_growth_rate = ((Mth4_B_Gth_HC_cm - birth_hc_cm)/stoolcol_4m_day)) %>% 
  full_join(metadata_4m)


###seprate female
metadata_female = metadata_4m[metadata_4m$Sex == "Female" , ]
asv_table_4m_female <- asv_table_4m[asv_table_4m$Sample %in% metadata_female$Sample , ]
all(metadata_female$Sample == rownames(asv_table_4m_female))
rownames(asv_table_4m_female) <- asv_table_4m_female$Sample 
asv_table_4m_female$Sample <- NULL

###predict growht using z-score (female)

a <- adonis2(asv_table_4m_female ~ wt_growth_z + lh_growth_z + hc_growth_z, 
             data = metadata_female, permutations = 10000, method = "bray", by = "margin")

a

pvalues <- c()
pvalues <- c(pvalues, a$`Pr(>F)`[1:3])
pvalues 

corrected_p_values <- p.adjust(pvalues, method = "BH")
corrected_p_values

###predict growht using rate (female)

b <- adonis2(asv_table_4m_female ~ wt_growth_rate + lh_growth_rate + hc_growth_rate, 
             data = metadata_female, permutations = 10000, method = "bray", by = "margin")

b

pvalues <- c()
pvalues <- c(pvalues, b$`Pr(>F)`[1:3])
pvalues

corrected_p_values <- p.adjust(pvalues, method = "BH")
corrected_p_values

###seprate male
metadata_male = metadata_4m[metadata_4m$Sex == "Male" , ]
asv_table_4m_male <- asv_table_4m[asv_table_4m$Sample %in% metadata_male$Sample , ]
rownames(asv_table_4m_male) <- asv_table_4m_male$Sample
all(metadata_male$Sample == rownames(asv_table_4m_male))

###predict growht using z-score (male)

c <- adonis2(asv_table_4m_male ~ wt_growth_z + lh_growth_z + hc_growth_z, 
             data = metadata_male, permutations = 10000, method = "bray", by = "margin")
c

metadata_male_no_NA <- metadata_male %>% filter(!is.na(metadata_male$wt_growth_z) & (!is.na(metadata_male$lh_growth_z)) & (!is.na(metadata_male$hc_growth_z)))

asv_table_4m_male[as.character(metadata_male_no_NA$Sample) , ]
asv_table_4m_male_no_NA = asv_table_4m_male[(metadata_male_no_NA$Sample) , ]

rownames(asv_table_4m_male_no_NA) <- asv_table_4m_male_no_NA$Sample
asv_table_4m_male_no_NA$Sample <- NULL
all(metadata_male_no_NA$Sample == rownames(asv_table_4m_male_no_NA))


c <- adonis2(asv_table_4m_male_no_NA ~ wt_growth_z + lh_growth_z + hc_growth_z, 
             data = metadata_male_no_NA, permutations = 10000, method = "bray", by = "margin")

c

pvalues <- c()
pvalues <- c(pvalues, c$`Pr(>F)`[1:3])
pvalues

corrected_p_values <- p.adjust(pvalues, method = "BH")
corrected_p_values

###predict growth using rate (male)

a <- adonis2(asv_table_4m_male ~ wt_growth_rate + lh_growth_rate + hc_growth_rate, 
             data = metadata_male, permutations = 10000, method = "bray", by = "margin")
a

metadata_male_no_NA <- metadata_male %>% filter(!is.na(metadata_male$wt_growth_rate) & (!is.na(metadata_male$lh_growth_rate)) & (!is.na(metadata_male$hc_growth_rate)))

asv_table_4m_male[as.character(metadata_male_no_NA$Sample) , ]
asv_table_4m_male_no_NA = asv_table_4m_male[(metadata_male_no_NA$Sample) , ]

rownames(asv_table_4m_male_no_NA) <- asv_table_4m_male_no_NA$Sample
asv_table_4m_male_no_NA$Sample <- NULL
all(metadata_male_no_NA$Sample == rownames(asv_table_4m_male_no_NA))

a <- adonis2(asv_table_4m_male_no_NA ~ wt_growth_rate + lh_growth_rate + hc_growth_rate, 
             data = metadata_male_no_NA, permutations = 10000, method = "bray", by = "margin")
a

pvalues <- c()
pvalues <- c(pvalues, a$`Pr(>F)`[1:3])
pvalues

corrected_p_values <- p.adjust(pvalues, method = "BH")
corrected_p_values

###Maaslin2 in growth prediction
library(tidyverse)
library("devtools")
library(microbiomics)
library(Maaslin2)

###both male and female-z-score, because is already corrected for gestation age and sex
rownames(metadata_4m) <- metadata_4m$Sample
set.seed(10000)

asv_table_4m$Sample <- NULL
maaslin <- 
  Maaslin2(input_data = asv_table_4m, 
           input_metadata = metadata_4m,
           output = "maaslin2_4m_growth_bothmalefemale_zscore_20200722",
           fixed_effects = c("wt_growth_z", "lh_growth_z", "hc_growth_z"),
           random_effects = NULL,
           transform = "LOG", 
           plot_scatter = F,
           min_prevalence = 0.1,
           cores = 4)

###both male and female-growth rate corrected for gestation age and sex
rownames(metadata_4m) <- metadata_4m$Sample
set.seed(10000)

asv_table_4m$Sample <- NULL
maaslin <- 
  Maaslin2(input_data = asv_table_4m, 
           input_metadata = metadata_4m,
           output = "maaslin2_4m_growth_bothmalefemale2_20200722",
           fixed_effects = c("Sex", "gestational_wk", "wt_growth_rate", "lh_growth_rate", "hc_growth_rate"),
           random_effects = NULL,
           transform = "LOG", 
           plot_scatter = F,
           min_prevalence = 0.1,
           cores = 4)


#Selecet for day 10
metadata_D10 = full_metadata3[full_metadata3$timepoint == "Day 10" , ]
asv_table_D10 <- ASV_table[ASV_table$Sample %in% metadata_D10$Sample , ]
all(metadata_D10$Sample == rownames(asv_table_D10))

###WT - z-score
metadata_D10 <- 
  metadata_D10 %>%
  mutate(wt_growth_z = (Mth4_B_Gth_WT_z - birth_wt_z)) %>% 
  full_join(metadata_D10)

###LH - z-score
metadata_D10 <- 
  metadata_D10 %>%
  mutate(lh_growth_z = (Mth4_B_Gth_LH_z - birth_length_z)) %>% 
  full_join(metadata_D10)

###HC - z-score
metadata_D10 <- 
  metadata_D10 %>%
  mutate(hc_growth_z = (Mth4_B_Gth_HC_z - birth_hc_z)) %>% 
  full_join(metadata_D10)

rownames(metadata_D10) <- metadata_D10$Sample


###seprate female
metadata_female = metadata_D10[metadata_D10$Sex == "Female" , ]
asv_table_D10_female <- asv_table_D10[asv_table_D10$Sample %in% metadata_female$Sample , ]
all(metadata_female$Sample == rownames(asv_table_D10_female))
rownames(asv_table_D10_female) <- asv_table_D10_female$Sample 


###predict growht using z-score (female)
```{r}
a <- adonis2(asv_table_D10_female ~ wt_growth_z + lh_growth_z + hc_growth_z, 
             data = metadata_female, permutations = 10000, method = "bray", by = "margin")

a

metadata_female_no_NA <- metadata_female %>% filter(!is.na(metadata_female$wt_growth_z) & (!is.na(metadata_female$lh_growth_z)) & (!is.na(metadata_female$hc_growth_z)))

asv_table_D10_female[as.character(metadata_female_no_NA$Sample) , ]
asv_table_D10_female_no_NA = asv_table_D10_female[(metadata_female_no_NA$Sample) , ]

rownames(asv_table_D10_female_no_NA) <- asv_table_D10_female_no_NA$Sample
asv_table_D10_female_no_NA$Sample <- NULL
all(metadata_female_no_NA$Sample == rownames(asv_table_D10_female_no_NA))


a <- adonis2(asv_table_D10_female_no_NA ~ wt_growth_z + lh_growth_z + hc_growth_z, 
             data = metadata_female_no_NA, permutations = 10000, method = "bray", by = "margin")

a

pvalues <- c()
pvalues <- c(pvalues, a$`Pr(>F)`[1:3])
pvalues 

corrected_p_values <- p.adjust(pvalues, method = "BH")
corrected_p_values


###seprate male
metadata_male = metadata_D10[metadata_D10$Sex == "Male" , ]
asv_table_D10_male <- asv_table_D10[asv_table_D10$Sample %in% metadata_male$Sample , ]
rownames(asv_table_D10_male) <- asv_table_D10_male$Sample
all(metadata_male$Sample == rownames(asv_table_D10_male))

###predict growht using z-score (male)

c <- adonis2(asv_table_D10_male ~ wt_growth_z + lh_growth_z + hc_growth_z, 
             data = metadata_male, permutations = 10000, method = "bray", by = "margin")
c

metadata_male_no_NA <- metadata_male %>% filter(!is.na(metadata_male$wt_growth_z) & (!is.na(metadata_male$lh_growth_z)) & (!is.na(metadata_male$hc_growth_z)))

asv_table_D10_male[as.character(metadata_male_no_NA$Sample) , ]
asv_table_D10_male_no_NA = asv_table_D10_male[(metadata_male_no_NA$Sample) , ]

rownames(asv_table_D10_male_no_NA) <- asv_table_D10_male_no_NA$Sample
asv_table_D10_male_no_NA$Sample <- NULL
all(metadata_male_no_NA$Sample == rownames(asv_table_D10_male_no_NA))


c <- adonis2(asv_table_D10_male_no_NA ~ wt_growth_z + lh_growth_z + hc_growth_z, 
             data = metadata_male_no_NA, permutations = 10000, method = "bray", by = "margin")

c

pvalues <- c()
pvalues <- c(pvalues, c$`Pr(>F)`[1:3])
pvalues

corrected_p_values <- p.adjust(pvalues, method = "BH")
corrected_p_values

###Maaslin2 in growth prediction
library(tidyverse)
library("devtools")
library(microbiomics)
library(Maaslin2)

###both male and female-z-score, because is already corrected for gestation age and sex
rownames(metadata_4m) <- metadata_4m$Sample
set.seed(10000)

asv_table_4m$Sample <- NULL
maaslin <- 
  Maaslin2(input_data = asv_table_4m, 
           input_metadata = metadata_4m,
           output = "maaslin2_4m_growth_bothmalefemale_zscore_20200722",
           fixed_effects = c("wt_growth_z", "lh_growth_z", "hc_growth_z"),
           random_effects = NULL,
           transform = "LOG", 
           plot_scatter = F,
           min_prevalence = 0.1,
           cores = 4)



# (B) Metagenomic analysis  -----------------------------------------------

library(permute)
library(lattice)
library("vegan"); packageVersion("vegan")
library("tidyverse"); packageVersion("tidyverse")
library("ggplot2"); packageVersion("ggplot2")
library(devtools)
library("phyloseq"); packageVersion("phyloseq")
library("RColorBrewer"); packageVersion("RColorBrewer")
library(dplyr)
library(ggpubr)

setwd("P:/FMHSfiles/DIAMOND/mgx_analysis")
load("diamond_metadata_mgx_20200721.RData")

otu_table <- read_tsv("merged_abundance_table_species_2.tsv")

otu_table2 <- otu_table %>% 
  mutate(Sample = sapply(otu_table$ID, function(x) strsplit(x, ".", fixed = TRUE)[[1]][1]))

otu_table2$Sample <- str_replace_all(otu_table2$Sample, "H", "H-")    

otu_table2$Sample <- str_replace_all(otu_table2$Sample, "WTH", "WT")    

rownames(otu_table2) <- otu_table2$Sample                       

#save(otu_table2, file = "diamond_metaphlan_20200204_mgx.RData")
otu_table2$ID <- NULL


# Adonis - General factors ------------------------------------------------
###without birth weight

otu_table2$Sample <- NULL
a <- adonis2(otu_table2 ~ gestational_wk + Hospital_Name + hospital_duration + Sex + delivery_mode + plural, 
             data = metadata_mgx, permutations = 10000, method = "bray", by = "margin")
a

pvalues <- c()
pvalues <- c(pvalues, a$`Pr(>F)`[1:6])
pvalues

corrected_p_values <- p.adjust(pvalues, method = "BH")
corrected_p_values

###with birth weight

b <- adonis2(otu_table2 ~ gestational_wk + Hospital_Name + hospital_duration + Sex + birth_wt_gm + delivery_mode + plural, 
             data = metadata_mgx, permutations = 10000, method = "bray", by = "margin")

b

pvalues <- c()
pvalues <- c(pvalues, b$`Pr(>F)`[1:7])
pvalues

corrected_p_values <- p.adjust(pvalues, method = "BH")
corrected_p_values


# Adonis - maternal factors -----------------------------------------------

b <- adonis2(otu_table2 ~ maternal_age + edu + maternal_ethnicity + 
               SES + maternal_stress_4month_PSS + maternal_Depression_4month_EPDS, 
             data = metadata_mgx, permutations = 10000, method = "bray", by = "margin")
b

otu_table2 <- as.data.frame(otu_table2)
metadata_mgx_no_NA <- metadata_mgx %>% filter(!is.na(metadata_mgx$SES))

otu_table2[ as.character(metadata_mgx_no_NA$Study_ID) , ]
otu_table2_no_NA = otu_table2[(metadata_mgx_no_NA$Study_ID) , ]

all(metadata_mgx_no_NA$Study_ID == rownames(otu_table2_no_NA))

b <- adonis2(otu_table2_no_NA ~ maternal_age + edu + maternal_ethnicity + 
               SES + maternal_stress_4month_PSS + maternal_Depression_4month_EPDS, 
             data = metadata_mgx_no_NA, permutations = 10000, method = "bray", by = "margin")

b

pvalues <- c()
pvalues <- c(pvalues, b$`Pr(>F)`[1:6])
pvalues

corrected_p_values <- p.adjust(pvalues, method = "BH")
corrected_p_values


# Adonis - Nutritional factors --------------------------------------------

c <- adonis2(otu_table2 ~ feeding + K_Q22b_anti_baby_4M +
               K_Q22d_anti_mother_4M + K_Q22e_prob_mother_4M +  K_Q22c_prob_baby_4M, 
             data = metadata_mgx, permutations = 10000, method = "bray", by = "margin")
c

metadata_mgx_no_NA <- metadata_mgx %>% filter(!is.na(metadata_mgx$K_Q22b_anti_baby_4M) & !is.na(metadata_mgx$feeding))

otu_table2[ as.character(metadata_mgx_no_NA$Study_ID) , ]
otu_table2_no_NA = otu_table2[(metadata_mgx_no_NA$Study_ID) , ]
#rownames(full_metadata) <- full_metadata$Study_ID
all(metadata_mgx_no_NA$Study_ID == rownames(otu_table2_no_NA))

c <- adonis2(otu_table2_no_NA ~ feeding + K_Q22b_anti_baby_4M +
               K_Q22d_anti_mother_4M + K_Q22e_prob_mother_4M +  K_Q22c_prob_baby_4M, 
             data = metadata_mgx_no_NA, permutations = 10000, method = "bray", by = "margin")

c

pvalues <- c()
pvalues <- c(pvalues, c$`Pr(>F)`[1:5])
pvalues

corrected_p_values <- p.adjust(pvalues, method = "BH")
corrected_p_values

###the effect of early-life nutrition - without smell/taste

rownames(metadata_mgx) <- metadata_mgx$Study_ID
a <- adonis2(otu_table2 ~ feeding_d10_v3 + nutrition_v2, 
             data = metadata_mgx, permutations = 10000, method = "bray", by = "margin")
a

metadata_mgx_no_NA <- metadata_mgx %>% filter(!is.na(metadata_mgx$feeding_d10_v3))

otu_table2[ as.character(metadata_mgx_no_NA$Study_ID) , ]
otu_table2_no_NA = otu_table2[(metadata_mgx_no_NA$Study_ID) , ]
all(metadata_mgx_no_NA$Study_ID == rownames(otu_table2_no_NA))

a <- adonis2(otu_table2_no_NA ~ feeding_d10_v3 + nutrition_v2, 
             data = metadata_mgx_no_NA, permutations = 10000, method = "bray", by = "margin")
a

pvalues <- c()
pvalues <- c(pvalues, a$`Pr(>F)`[1:2])
pvalues

corrected_p_values <- p.adjust(pvalues, method = "BH")
corrected_p_values

###the effect of early-life nutrition - with smell/taste

a <- adonis2(otu_table2 ~ feeding_d10_v3 + nutrition_v2 + total_smell_taste, 
             data = metadata_mgx, permutations = 10000, method = "bray", by = "margin")
a

metadata_mgx_no_NA <- metadata_mgx %>% filter(!is.na(metadata_mgx$feeding_d10_v3))

otu_table2[ as.character(metadata_mgx_no_NA$Study_ID) , ]
otu_table2_no_NA = otu_table2[(metadata_mgx_no_NA$Study_ID) , ]
all(metadata_mgx_no_NA$Study_ID == rownames(otu_table2_no_NA))

a <- adonis2(otu_table2_no_NA ~ feeding_d10_v3 + nutrition_v2 + total_smell_taste, 
             data = metadata_mgx_no_NA, permutations = 10000, method = "bray", by = "margin")

a

pvalues <- c()
pvalues <- c(pvalues, a$`Pr(>F)`[1:3])
pvalues

corrected_p_values <- p.adjust(pvalues, method = "BH")
corrected_p_values

###early-life factors + nutrition 

a <- adonis2(otu_table2 ~ feeding_d10_v3 + nutrition_v2 + G_Q5_Probio + M_Q5b_AntibBby + Hospital_Name, 
             data = metadata_mgx, permutations = 10000, method = "bray", by = "margin")
a

metadata_mgx_no_NA <- metadata_mgx %>% filter(!is.na(metadata_mgx$feeding_d10_v3) & !is.na(metadata_mgx$M_Q5b_AntibBby))

otu_table2[ as.character(metadata_mgx_no_NA$Study_ID) , ]
otu_table2_no_NA = otu_table2[(metadata_mgx_no_NA$Study_ID) , ]
all(metadata_mgx_no_NA$Study_ID == rownames(otu_table2_no_NA))

a <- adonis2(otu_table2_no_NA ~ feeding_d10_v3 + nutrition_v2 + G_Q5_Probio + M_Q5b_AntibBby + Hospital_Name, 
             data = metadata_mgx_no_NA, permutations = 10000, method = "bray", by = "margin")

#n=89

a

pvalues <- c()
pvalues <- c(pvalues, a$`Pr(>F)`[1:5])
pvalues

corrected_p_values <- p.adjust(pvalues, method = "BH")
corrected_p_values

###early-life factors + nutrition (Infant & mother)

a <- adonis2(otu_table2 ~ feeding_d10_v3 + nutrition_v2 + G_Q5_Probio + M_Q5b_AntibBby + Hospital_Name + M_Q5c_AntibMoL4Wk + M_Q5d_ProbMoL4W, 
             data = metadata_mgx, permutations = 10000, method = "bray", by = "margin")
a

metadata_mgx_no_NA <- metadata_mgx %>% filter(!is.na(metadata_mgx$feeding_d10_v3) & !is.na(metadata_mgx$M_Q5b_AntibBby) & !is.na(metadata_mgx$M_Q5c_AntibMoL4Wk) & !is.na(metadata_mgx$M_Q5d_ProbMoL4W))

otu_table2[ as.character(metadata_mgx_no_NA$Study_ID) , ]
otu_table2_no_NA = otu_table2[(metadata_mgx_no_NA$Study_ID) , ]
all(metadata_mgx_no_NA$Study_ID == rownames(otu_table2_no_NA))

a <- adonis2(otu_table2_no_NA ~ feeding_d10_v3 + nutrition_v2 + G_Q5_Probio + M_Q5b_AntibBby + Hospital_Name + M_Q5c_AntibMoL4Wk + M_Q5d_ProbMoL4W, 
             data = metadata_mgx_no_NA, permutations = 10000, method = "bray", by = "margin")

a

pvalues <- c()
pvalues <- c(pvalues, a$`Pr(>F)`[1:7])
pvalues

corrected_p_values <- p.adjust(pvalues, method = "BH")
corrected_p_values


# Maaslin2 - all D10 and 4M factors  --------------------------------------

library(tidyverse)
library("devtools")
library(microbiomics)
library(Maaslin2)

setwd("P:/FMHSfiles/DIAMOND/mgx_analysis")
load("diamond_metadata_mgx_20200721.RData")

otu_table <- read_tsv("merged_abundance_table_species_2.tsv")

otu_table2 <- otu_table %>% 
  mutate(Sample = sapply(otu_table$ID, function(x) strsplit(x, ".", fixed = TRUE)[[1]][1]))

otu_table2$Sample <- str_replace_all(otu_table2$Sample, "H", "H-")    

otu_table2$Sample <- str_replace_all(otu_table2$Sample, "WTH", "WT")    

rownames(otu_table2) <- otu_table2$Sample                       

#save(otu_table2, file = "diamond_metaphlan_20200204_mgx.RData")
otu_table2$ID <- NULL

rownames(metadata_mgx) <- metadata_mgx$Study_ID

otu_table2$Sample <- NULL
set.seed(10000)

maaslin <- 
  Maaslin2(input_data = otu_table2, 
           input_metadata = metadata_mgx,
           output = "maaslin2_mgx_allfactors_20200722",
           fixed_effects = c("gestational_wk", "Hospital_Name", "hospital_duration", "Sex", "delivery_mode", "plural", 
                             "nutrition_v2", "total_smell_taste", "feeding_d10_v3", 
                             "maternal_age", "edu", "maternal_ethnicity", "SES", 
                             "G_Q7_Milk_Dis", "K_Q22b_anti_baby_4M", "K_Q22c_prob_baby_4M",  
                             "K_Q22d_anti_mother_4M", "K_Q22e_prob_mother_4M", "maternal_stress_4month_PSS", 
                             "maternal_Depression_4month_EPDS", "feeding"),
           random_effects = NULL,
           transform = "LOG", 
           plot_scatter = F,
           min_prevalence = 0.1,
           cores = 4)


# Shannon -----------------------------------------------------------------
library(rstatix)
shannon <- read_tsv("shannon_mgx_20200125.tsv")

metadata_mgx <- metadata_mgx %>% 
  rename(Sample = Study_ID)

shannon_full <- shannon %>%
  left_join(metadata_mgx)

a <- wilcox_test(shannon_full, formula = Shannon ~ Sex, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
a

pvalues <- c()
pvalues <- c(pvalues, a$p)
pvalues

b <- wilcox_test(shannon_full, formula = Shannon ~ delivery_mode, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
b

pvalues <- c(pvalues, b$p)
pvalues

c <- wilcox_test(shannon_full, formula = Shannon ~ gestational_wk, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
c

pvalues <- c(pvalues, c$p)
pvalues

d <- wilcox_test(shannon_full, formula = Shannon ~ edu, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
d

pvalues <- c(pvalues, d$p)
pvalues

e <- wilcox_test(shannon_full, formula = Shannon ~ nutrition_v2, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
e

pvalues <- c(pvalues, e$p)
pvalues

f <- wilcox_test(shannon_full, formula = Shannon ~ total_smell_taste, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
f

pvalues <- c(pvalues, f$p)
pvalues

g <- wilcox_test(shannon_full, formula = Shannon ~ maternal_Depression_4month_EPDS, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
g

pvalues <- c(pvalues, g$p)
pvalues

h <- wilcox_test(shannon_full, formula = Shannon ~ K_Q22d_anti_mother_4M, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
h

pvalues <- c(pvalues, h$p)
pvalues

i <- wilcox_test(shannon_full, formula = Shannon ~ K_Q22e_prob_mother_4M, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
i

pvalues <- c(pvalues, i$p)
pvalues

j <- wilcox_test(shannon_full, formula = Shannon ~ K_Q22b_anti_baby_4M, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
j

pvalues <- c(pvalues, j$p)
pvalues

k <- wilcox_test(shannon_full, formula = Shannon ~ K_Q22c_prob_baby_4M, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
k

pvalues <- c(pvalues, k$p)
pvalues

l <- wilcox_test(shannon_full, formula = Shannon ~ plural, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
l

pvalues <- c(pvalues, l$p)
pvalues

corrected_p_values <- p.adjust(pvalues, method = "BH")
corrected_p_values

###Kruskal-wallis

a <- kruskal.test(Shannon ~ Hospital_Name, data = shannon_full)
a

pvalues <- c()
pvalues <- c(pvalues, a$p.value)
pvalues

b <- wilcox_test(shannon_full, formula = Shannon ~ Hospital_Name, comparisons = list(c("Auckland","Middlemore"),c("Auckland","NorthShore"),c("Auckland","Waitakere"),c("Middlemore","NorthShore"),c("Middlemore","Waitakere"),c("NorthShore","Waitakere")), paired = FALSE) %>%
  adjust_pvalue(method = "BH")

b

c <- kruskal.test(Shannon ~ maternal_ethnicity, data = shannon_full)
c

pvalues <- c(pvalues, c$p.value)
pvalues

d <- wilcox_test(shannon_full, formula = Shannon ~ maternal_ethnicity, comparisons = list(c("Asian","Caucasian"),c("Asian","Maori"),c("Asian","Pacific Islander"),c("Asian","Other"),c("Caucasian","Maori"),c("Pacific Islander","Maori"),c("Pacific Islander","Other"),c("Other","Maori"),c("Caucasian","Pacific Islander"),c("Caucasian","Other")), paired = FALSE) %>% adjust_pvalue(method = "BH")
d

e <- kruskal.test(Shannon ~ SES, data = shannon_full)
e

pvalues <- c(pvalues, e$p.value)
pvalues

f <- wilcox_test(shannon_full, formula = Shannon ~ SES, comparisons = list(c("Least deprived","Most deprived"),c("Least deprived","Moderate"),c("Most deprived","Moderate")), paired = FALSE) %>%
  adjust_pvalue(method = "BH")
f

g <- kruskal.test(Shannon ~ maternal_stress_4month_PSS, data = shannon_full)
g

pvalues <- c(pvalues, g$p.value)
pvalues

h <- wilcox_test(shannon_full, formula = Shannon ~ maternal_stress_4month_PSS, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
h

i <- kruskal.test(Shannon ~ feeding, data = shannon_full)
i

pvalues <- c(pvalues, i$p.value)
pvalues

j <- wilcox_test(shannon_full, formula = Shannon ~ feeding, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
j

k <- kruskal.test(Shannon ~ feeding_d10_v3, data = shannon_full)
k

pvalues <- c(pvalues, k$p.value)
pvalues

l <- wilcox_test(shannon_full, formula = Shannon ~ feeding_d10_v3, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
l

corrected_p_values <- p.adjust(pvalues, method = "BH")
corrected_p_values


# Pathway analysis --------------------------------------------------------

setwd("P:/FMHSfiles/DIAMOND/mgx_analysis")
load("diamond_metadata_mgx_20200721.RData")

setwd("P:/FMHSfiles/DIAMOND/mgx_analysis/pathwayab_analysis")
load("pathwayab_relab_transposed_20200128.RData")
load("pathwayab_cpm_transposed_20200128.RData")

library(tidyverse)
library("devtools")
library(microbiomics)
library(Maaslin2)
rownames(pathwayab_relab_transposed) <- pathwayab_relab_transposed$Sample
pathwayab_relab_transposed$Sample <- NULL

rownames(metadata_mgx) <- metadata_mgx$Study_ID
set.seed(10000)

#cpm
rownames(pathwayab_cpm_transposed) <- pathwayab_cpm_transposed$Sample
pathwayab_cpm_transposed$Sample <- NULL

rownames(metadata_mgx) <- metadata_mgx$Study_ID
set.seed(10000)



###removed two waitakere data as only two, not repesentative
grep("WT-0005", metadata_mgx$Study_ID) #find which rows are there #22

grep("WT-0004", metadata_mgx$Study_ID) #71

metadata_mgx <- metadata_mgx[-c(22, 71) , ]
set.seed(10000)
maaslin <- 
  Maaslin2(input_data = pathwayab_cpm_transposed, 
           input_metadata = metadata_mgx,
           output = "maaslin2_mgx_cpm2_20200722(removed WT)",
           fixed_effects = c("gestational_wk", "Hospital_Name", "hospital_duration", "Sex", "delivery_mode", "plural", 
                             "nutrition_v2", "feeding_d10_v3", "M_Q5c_AntibMoL4Wk", "M_Q5d_ProbMoL4W", 
                             "maternal_age", "edu", "maternal_ethnicity", "SES", 
                             "G_Q4_FullyBF", "Mth4_N_Q4_BMilk", "G_Q5_Probio", "M_Q5b_AntibBby", 
                             "K_Q22b_anti_baby_4M", "K_Q22c_prob_baby_4M",  
                             "K_Q22d_anti_mother_4M", "K_Q22e_prob_mother_4M", "maternal_stress_4month_PSS", 
                             "maternal_Depression_4month_EPDS", "feeding"),
           random_effects = NULL,
           plot_scatter = T,
           transform = "none", #maybe we need it if the number is too big/small
           normalization = "none", #we have done it manually, so dont have to do it again
           min_prevalence = 0.1,
           cores = 4)



# Predict growth ----------------------------------------------------------

setwd("P:/FMHSfiles/DIAMOND/mgx_analysis")
library(permute)
library(lattice)
library("vegan"); packageVersion("vegan")
library("tidyverse"); packageVersion("tidyverse")
library("ggplot2"); packageVersion("ggplot2")
library(devtools)
library("phyloseq"); packageVersion("phyloseq")
library("RColorBrewer"); packageVersion("RColorBrewer")
library(dplyr); packageVersion("dplyr")
library(ggpubr)
library(rstatix)

setwd("P:/FMHSfiles/DIAMOND/mgx_analysis")
load("diamond_metadata_mgx_20200721.RData")

otu_table <- read_tsv("merged_abundance_table_species_2.tsv")

otu_table2 <- otu_table %>% 
  mutate(Sample = sapply(otu_table$ID, function(x) strsplit(x, ".", fixed = TRUE)[[1]][1]))

otu_table2$Sample <- str_replace_all(otu_table2$Sample, "H", "H-")    

otu_table2$Sample <- str_replace_all(otu_table2$Sample, "WTH", "WT")    

rownames(otu_table2) <- otu_table2$Sample                       

#save(otu_table2, file = "diamond_metaphlan_20200204_mgx.RData")
otu_table2$ID <- NULL

#WT - z-score
metadata_4m <- 
  metadata_mgx %>%
  mutate(wt_growth_z = (Mth4_B_Gth_WT_z - birth_wt_z)) %>% 
  full_join(metadata_mgx)

#LH - z-score
metadata_4m <- 
  metadata_4m%>%
  mutate(lh_growth_z = (Mth4_B_Gth_LH_z - birth_length_z)) %>% 
  full_join(metadata_4m)

#HC - z-score
metadata_4m <- 
  metadata_4m %>%
  mutate(hc_growth_z = (Mth4_B_Gth_HC_z - birth_hc_z)) %>% 
  full_join(metadata_4m)

rownames(metadata_4m) <- metadata_4m$Study_ID


#seprate female
metadata_female = metadata_4m[metadata_4m$Sex == "Female" , ]
otu_table2_female <- otu_table2[otu_table2$Sample %in% metadata_female$Study_ID , ]
all(metadata_female$Sample == rownames(otu_table2_female))
rownames(otu_table2_female) <- otu_table2_female$Sample 
otu_table2_female$Sample <- NULL

#seprate male
metadata_male = metadata_4m[metadata_4m$Sex == "Male" , ]
otu_table2_male <- otu_table2[otu_table2$Sample %in% metadata_male$Study_ID , ]
rownames(otu_table2_male) <- otu_table2_male$Sample
all(metadata_male$Study_ID == rownames(otu_table2_male))
otu_table2_male$Sample <- NULL


library(tidyverse)
library("devtools")
library(microbiomics)
library(Maaslin2)

#Maaslin2 for both male and female-z-score, because is already corrected for gestation age and sex

otu_table2$Sample <- NULL
set.seed(10000)
maaslin <- 
  Maaslin2(input_data = otu_table2, 
           input_metadata = metadata_4m,
           output = "maaslin2_mgx_growth_bothmalefemale_zscore_20200722",
           fixed_effects = c("wt_growth_z", "lh_growth_z", "hc_growth_z"),
           random_effects = NULL,
           transform = "LOG", 
           plot_scatter = F,
           min_prevalence = 0.1,
           cores = 4)


#male - Maaslin2 
set.seed(10000)

maaslin <- 
  Maaslin2(input_data = otu_table2_male, 
           input_metadata = metadata_male,
           output = "maaslin2_mgx_growth_male_zscore_20200722",
           fixed_effects = c("wt_growth_z", "lh_growth_z", "hc_growth_z"),
           random_effects = "LOG", 
           plot_scatter = F,
           min_prevalence = 0.1,
           cores = 4)

#female - Maaslin2
set.seed(10000)
maaslin <- 
  Maaslin2(input_data = otu_table2_female, 
           input_metadata = metadata_female,
           output = "maaslin2_mgx_growth_female_zscore_20200722",
           fixed_effects = c("wt_growth_z", "lh_growth_z", "hc_growth_z"),
           random_effects = "LOG", 
           plot_scatter = F,
           min_prevalence = 0.1,
           cores = 4)


#Adonis - female

a <- adonis2(otu_table2_female ~ wt_growth_z + lh_growth_z + hc_growth_z, 
             data = metadata_female, permutations = 10000, method = "bray", by = "margin")

a

pvalues <- c()
pvalues <- c(pvalues, a$`Pr(>F)`[1:3])
pvalues

corrected_p_values <- p.adjust(pvalues, method = "BH")
corrected_p_values

#Adonis - male

c <- adonis2(otu_table2_male ~ wt_growth_z + lh_growth_z + hc_growth_z, 
             data = metadata_male, permutations = 10000, method = "bray", by = "margin")
c

metadata_male_no_NA <- metadata_male %>% filter(!is.na(metadata_male$wt_growth_z) & (!is.na(metadata_male$lh_growth_z)) & (!is.na(metadata_male$hc_growth_z)))

otu_table2_male[as.character(metadata_male_no_NA$Study_ID) , ]
otu_table2_male_no_NA = otu_table2_male[(metadata_male_no_NA$Study_ID) , ]

all(metadata_male_no_NA$Sample == rownames(otu_table2_male_no_NA))


c <- adonis2(otu_table2_male_no_NA ~ wt_growth_z + lh_growth_z + hc_growth_z, 
             data = metadata_male_no_NA, permutations = 10000, method = "bray", by = "margin")

c

pvalues <- c()
pvalues <- c(pvalues, c$`Pr(>F)`[1:3])
pvalues

corrected_p_values <- p.adjust(pvalues, method = "BH")
corrected_p_values

#Adonis - both sexes

d <- adonis2(otu_table2 ~ wt_growth_z + lh_growth_z + hc_growth_z, 
             data = metadata_4m, permutations = 10000, method = "bray", by = "margin")

d

metadata_4m_no_NA <- metadata_4m %>% filter(!is.na(metadata_4m$wt_growth_z) & (!is.na(metadata_4m$lh_growth_z)) & (!is.na(metadata_4m$hc_growth_z)))

otu_table2[as.character(metadata_4m_no_NA$Study_ID) , ]
otu_table2_no_NA = otu_table2[(metadata_4m_no_NA$Study_ID) , ]
rownames(otu_table2_no_NA) <- otu_table2_no_NA$Sample
otu_table2_no_NA$Sample <- NULL
all(metadata_4m_no_NA$Sample == rownames(otu_table2_no_NA))


d <- adonis2(otu_table2_no_NA ~ wt_growth_z + lh_growth_z + hc_growth_z, 
             data = metadata_4m_no_NA, permutations = 10000, method = "bray", by = "margin")

d

pvalues <- c()
pvalues <- c(pvalues, d$`Pr(>F)`[1:3])
pvalues

corrected_p_values <- p.adjust(pvalues, method = "BH")
corrected_p_values

