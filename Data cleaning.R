#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#-#-#-#-#-#-#-#-# Vivli Team 1 Project #-#-#-#-#-#-#-#-#-#-#
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# Topic: Data exploration and cleaning for team 1 project

#clear environment
rm(list=ls())

# 0. Import packages and datasets ####

#import packages
library(readxl)
library(tidyverse)
library(AMR) #may or may not use
library(countrycode)

#set wd, could change to the OneDrive if we want
# setwd("C:/Users/jakeh/OneDrive - Nexus365/DPhil project/Vivli data challenge")
setwd("")

#import datasets
atlas <- read.csv("./2024_05_28 atlas_antibiotics.csv", header=TRUE)
gears <- read_excel("./Venatorx surveillance data_2024_06_06.xlsx")

# Summaries #

ncol(atlas) #135
nrow(atlas) #917049
summary(atlas) #135 varialbes, majority are characters except Isolate.Id and Year it looks like
#Seems all data on MIC is in character format

ncol(gears) #23
nrow(gears) #29365
summary(gears) #23 variables, character except: Isolate, Year, CAZ_MIC, FEP_MIC, MEM_MIC, GM_MIC, TSP_MIC
#Some labelled MIC are characters, some numeric, so need to explore more what these variables are in the dataset


# 1. Dataset Preliminary Exploration ####

## 1.1 Atlas Descriptives####

  ###1.1.1 Countries and entries by country ####

  table(atlas$Country, useNA = "always") #visual check as few as only 14 in Estonia
  length(unique(atlas$Country)) #83
  
  atlas_countries_year <- atlas %>%
    dplyr::select(Country, Year) %>%
    group_by(Country, Year) %>%
    mutate(number_of_entries = n()) %>%
    distinct()
  
  #Comment: 13 have only one for a year, 44 country/year pairing have under 10 in a year
  #For example with one:  United Kingdom 2005, Guatemala 2014,  Venezuela 2004, Oman  2006, Hungary 2006
  #US minimum 4991 in a year
  
  #Visualizing the differences
  # ggplot(atlas_countries_year, aes(y=Country, x=number_of_entries, color=Year)) +
    # geom_point(aes(color=Year)) +
    # xlim(0,6500) +  #uncomment this if want to truncate the US to better visualize other countries
    # xlab("Number of Entries") #+
    #scale_colour_gradient(low = "yellow", high = "red")
  
  ### 1.1.2 Other non-resistance related variables #####
  
  #state variable = US state
  table(atlas$State, useNA = "always") #all US states + empty (no "NA", empty instead)
  length(unique(atlas$State)) #47 unique
  
  #organism of infection
  table(atlas$Species, useNA = "always") #quick scan, some have 1 some have 100,000
  length(unique(atlas$Species)) #367
  
  #species family
  table(atlas$Family, useNA = "always")
  length(unique(atlas$Family)) #16
  
  #age groups
  table(atlas$Age.Group, useNA = "always") #15647 "Unknown", none listed as NA
  length(unique(atlas$Age.Group)) #7 groupings, groupings are: 0-2, 3-12, 13-18, 19-64, 65-84, 85+, Unknown
  
    #plotting the age distribution
    atlas$Age.Group <- factor(atlas$Age.Group, levels=c("0 to 2 Years", "3 to 12 Years", "13 to 18 Years", "19 to 64 Years", "65 to 84 Years", "85 and Over", "Unknown"))
    # ggplot(data.frame(atlas), aes(x=Age.Group)) +
    #   geom_bar()
  
  #speciality, seems to be where the test was taken
  table(atlas$Speciality, useNA = "always") #location (e.g. Clinic/Office, Emergency Room, etc.) 
  length(unique(atlas$Speciality)) #12
  
  #Source (sampling) location from the patient
  table(atlas$Source, useNA = "always") #where sample was taken (e.g. Appendix, Blood, etc.)
  length(unique(atlas$Source)) #97
  
  #in/out patient status
  table(atlas$In...Out.Patient, useNA = "always") #either nothing listed (229,433) or Inpatient (440977), None Given (156114), Other (1766), Outpatient (88759)
  
  #year
  table(atlas$Year, useNA = "always") #2004 to 2022, no NA; 20-40k until 2012 then goes to more 50-70k in samples
  
  table(atlas$Study, useNA = "always") #I think we were planning to include all, but need to check methodologies?
  #     Atlas  Inform SPIDAAR    TEST    <NA> 
  #   283534  219885     385  413245       0 
  
  table(atlas$Phenotype, useNA = "always")
  #            (BL Neg) (BL Pos)     ESBL     MRSA     MSSA     <NA> 
  #   671140    23323     7767    57383    68066    89370        0 
  
  
  ### 1.1.3 MIC-related variables ####
  
  #Multiple columns, so will only do some as an example
  
  #Amikacin
  table(atlas$Amikacin, useNA = "always") #several empty (360284), otherwise lists a specific value or if <= or > a value
  table(atlas$Amikacin_I, useNA = "always") #several empty (385568), Intermediate (26283), Resistant (47647), Susceptible (457551) 
  #Comment: notice that there are less empty for the S,I,R classification than those empty for the MIC value
  
  #amoxycillin clavulanate
  table(atlas$Amoxycillin.clavulanate, useNA = "always") #same type of output as Amikacin
  table(atlas$Amoxycillin.clavulanate_I, useNA = "always") #same type of output as Amikacin
  #Again, a discrepancy. There are 120,026 empty for the MIC, but 477,413 with no classification in the S,I,R
  
  #Expect other columns to have the same type of information
  #Other columns available are as follows:
  #                                    "Amikacin"                    "Amikacin_I"                  "Amoxycillin.clavulanate"    
  # [17] "Amoxycillin.clavulanate_I"   "Ampicillin"                  "Ampicillin_I"                "Azithromycin"               
  # [21] "Azithromycin_I"              "Cefepime"                    "Cefepime_I"                  "Cefoxitin"                  
  # [25] "Cefoxitin_I"                 "Ceftazidime"                 "Ceftazidime_I"               "Ceftriaxone"                
  # [29] "Ceftriaxone_I"               "Clarithromycin"              "Clarithromycin_I"            "Clindamycin"                
  # [33] "Clindamycin_I"               "Erythromycin"                "Erythromycin_I"              "Imipenem"                   
  # [37] "Imipenem_I"                  "Levofloxacin"                "Levofloxacin_I"              "Linezolid"                  
  # [41] "Linezolid_I"                 "Meropenem"                   "Meropenem_I"                 "Metronidazole"              
  # [45] "Metronidazole_I"             "Minocycline"                 "Minocycline_I"               "Penicillin"                 
  # [49] "Penicillin_I"                "Piperacillin.tazobactam"     "Piperacillin.tazobactam_I"   "Tigecycline"                
  # [53] "Tigecycline_I"               "Vancomycin"                  "Vancomycin_I"                "Ampicillin.sulbactam"       
  # [57] "Ampicillin.sulbactam_I"      "Aztreonam"                   "Aztreonam_I"                 "Aztreonam.avibactam"        
  # [61] "Aztreonam.avibactam_I"       "Cefixime"                    "Cefixime_I"                  "Ceftaroline"                
  # [65] "Ceftaroline_I"               "Ceftaroline.avibactam"       "Ceftaroline.avibactam_I"     "Ceftazidime.avibactam"      
  # [69] "Ceftazidime.avibactam_I"     "Ciprofloxacin"               "Ciprofloxacin_I"             "Colistin"                   
  # [73] "Colistin_I"                  "Daptomycin"                  "Daptomycin_I"                "Doripenem"                  
  # [77] "Doripenem_I"                 "Ertapenem"                   "Ertapenem_I"                 "Gatifloxacin"               
  # [81] "Gatifloxacin_I"              "Gentamicin"                  "Gentamicin_I"                "Moxifloxacin"               
  # [85] "Moxifloxacin_I"              "Oxacillin"                   "Oxacillin_I"                 "Quinupristin.dalfopristin"  
  # [89] "Quinupristin.dalfopristin_I" "Sulbactam"                   "Sulbactam_I"                 "Teicoplanin"                
  # [93] "Teicoplanin_I"               "Tetracycline"                "Tetracycline_I"              "Trimethoprim.sulfa"         
  # [97] "Trimethoprim.sulfa_I"        "Ceftolozane.tazobactam"      "Ceftolozane.tazobactam_I"    "Cefoperazone.sulbactam"     
  # [101] "Cefoperazone.sulbactam_I"    "Meropenem.vaborbactam"       "Meropenem.vaborbactam_I"     "Cefpodoxime"                
  # [105] "Cefpodoxime_I"               "Ceftibuten"                  "Ceftibuten_I"                "Ceftibuten.avibactam"       
  # [109] "Ceftibuten.avibactam_I"      "Tebipenem"                   "Tebipenem_I"
  
  # carbapenem or cephalosporin
  # 3rd gen cephalosporin: Cefotaxime, Ceftazidime, Ceftriaxone, Ceftibuten, Cefpodoxime, Cefoperazone.sulbactam
  # 2nd gen: cefoxitin, Ceftolozane.tazobactam, 
  # 4th gen: Cefepime, -- no 
  # Carbapenem: Tebipenem (no), Meropenem
  
  # Aislin help:
  # MRSA: cefoxitin (just staph), maybe oxacillin 
  # 3rd gen: Ceftriaxone
  # carbapenem: meropenem, imipenem
  
  ### 1.1.4 Other resistance/susceptibility variables ####
  
  table(atlas$AMPC, useNA = "always")
  #               -    NEG    POS   <NA> 
  #   903183   2039  10995    832      0 
  
  table(atlas$SHV, useNA = "always") #again negative, blank, or positive for multiple different SHV beta-lactamase
  
  #Other options we could check include:
  #                                     "TEM"                         "CTXM1"                       "CTXM2"                      
  # [117] "CTXM825"                     "CTXM9"                       "VEB"                         "PER"                        
  # [121] "GES"                         "ACC"                         "CMY1MOX"                     "CMY11"                      
  # [125] "DHA"                         "FOX"                         "ACTMIR"                      "KPC"                        
  # [129] "OXA"                         "NDM"                         "IMP"                         "VIM"                        
  # [133] "SPM"                         "GIM"
  
  
## 1.2 gears Descriptives ####
  
  ###1.2.1 Countries and entries by country ####
  
  table(gears$Country, useNA = "always") #visual check as few as only 14 in Estonia
  length(unique(gears$Country)) #59
  
  gears_countries_year <- gears %>%
    dplyr::select(Country, Year) %>%
    group_by(Country, Year) %>%
    mutate(number_of_entries = n()) %>%
    distinct()
  #Comment: 2 have only one for a year, 21 country/year pairing have under 10 in a year
  #US over 1000 every year
  
  #Visualizing the differences
  # ggplot(gears_countries_year, aes(y=Country, x=number_of_entries, color=Year)) +
  #   geom_point(aes(color=Year)) +
  #   scale_colour_gradient(low = "yellow", high = "red")
  
  ###1.2.2 Other non-MIC related variables #####
  
  #geographic region
  table(gears$Region, useNA = "always") #Africa, Asia, Europe, etc.
  length(unique(gears$Region)) #7 unique
  
  #organism of infection
  table(gears$Organism, useNA = "always") #seems similar to Species above but much less
  length(unique(gears$Organism)) #42
  
  #family of organism
  table(gears$Family, useNA = "always") #Only Enterobacteriaceae, Non-Entero, Staphylococcus spp (140)
  length(unique(gears$Family)) #3
  
  #age as continuous variable
  table(gears$Age, useNA = "always") #0 NA, but 281 NG
  is.character(gears$Age) #TRUE
  summary(as.numeric(gears$Age)) #281 NAs introduced, range is 0 to 119, Mean 57.18
  
  #visualize continuous age variable
  hist(as.numeric(gears$Age))
  
  #create age groups to compare with ATLAS dataset
  gears <- gears %>%
    mutate(Age.Group = case_when(
      Age == "NG" ~ "Unknown",
      Age <= 2 ~ "O to 2 Years",
      Age > 2 & Age <= 12 ~ "3 to 12 Years",
      Age > 12 & Age <= 18 ~ "13 to 18 Years",
      Age > 18 & Age <= 64 ~ "19 to 64 Years",
      Age > 64 & Age <= 84 ~ "65 to 84 Years",
      Age > 84 ~ "85 and Over"
    ))
  
  table(gears$Age.Group, useNA = "always")  #281 Unknown as expected and no new NA, no 3-12 or 13-18
  gears$Age.Group <- factor(gears$Age.Group, levels=c("O to 2 Years", "3 to 12 Years", "13 to 18 Years", "19 to 64 Years", "65 to 84 Years", "85 and Over", "Unknown"))
  ggplot(data.frame(gears), aes(x=Age.Group)) +
    geom_bar()
  
  #site of the test in the patient
  table(gears$BodySite, useNA = "always") #where sample was taken, different classification and less than atlas
  length(unique(gears$BodySite)) #42
  
  #facility/location where the test was done
  table(gears$Facility, useNA = "always") #location, generally seems to overlap with atlas (no Clinic or Nursing Home)
  length(unique(gears$Facility)) #10

  ###1.2.3 MIC-related variables ####
  
  summary(gears$CAZ_MIC) #0.0299 to 128.0001; Median=1, Mean=6.45
  sum(is.na(gears$CAZ_MIC)) #0
  
  table(gears$C_MIC, useNA = "always") #character because have "-" for 28,784 entries; no NA, otherwise numeric
  sum(is.na(gears$CAZ_MIC)) #0

  table(gears$CIP_MIC, useNA = "always") #character because have "-" for 21,446 entries; no NA, otherwise numeric
  

  
# 2. Frequency tables to identify drug/bugs of interest ####
  
  ## 2.1 ATLAS/GEARS Bug/Site Tables ####
  
  #species/organism frequencies for atlas and gears
  
  atlas_bug <- atlas %>%
    dplyr::select(Species) %>%
    group_by(Species) %>%
    mutate(frequency = n()) %>%
    distinct() %>%
    arrange(desc(frequency))
  # write.csv(atlas_bug, "./Frequency Tables/atlas_bug_frequency_table.csv", row.names = FALSE)
  
  gears_bug <- gears %>%
    dplyr::select(Organism) %>%
    group_by(Organism) %>%
    mutate(frequency = n()) %>%
    distinct() %>%
    arrange(desc(frequency))
  # write.csv(gears_bug, "./Frequency Tables/gears_bug_frequency_table.csv", row.names = FALSE)
  
  #source/bodysite sample location for atlas and gears
  
  atlas_bodysite <- atlas %>%
    dplyr::select(Source) %>%
    group_by(Source) %>%
    mutate(frequency = n()) %>%
    distinct() %>%
    arrange(desc(frequency))
  # write.csv(atlas_bodysite, "./Frequency Tables/atlas_bodysite_frequency_table.csv", row.names = FALSE)
  
  gears_bodysite <- gears %>%
    dplyr::select(BodySite) %>%
    group_by(BodySite) %>%
    mutate(frequency = n()) %>%
    distinct() %>%
    arrange(desc(frequency))
  # write.csv(gears_bodysite, "./Frequency Tables/gears_bodysite_frequency_table.csv", row.names = FALSE)
  
  #frequencies when including both bug and site
  
  atlas_bug_bodysite <- atlas %>%
    dplyr::select(Species,Source) %>%
    group_by(Species,Source) %>%
    mutate(frequency = n()) %>%
    distinct() %>%
    ungroup %>%
    group_by(Species) %>%
    arrange(desc(frequency), .by_group = TRUE)
  # write.csv(atlas_bug_bodysite, "./Frequency Tables/atlas_bug_bodysite_frequency_table.csv", row.names = FALSE)

  gears_bug_bodysite <- gears %>%
    dplyr::select(Organism,BodySite) %>%
    group_by(Organism,BodySite) %>%
    mutate(frequency = n()) %>%
    distinct() %>%
    ungroup %>%
    group_by(Organism) %>%
    arrange(desc(frequency), .by_group = TRUE)
  # write.csv(gears_bug_bodysite, "./Frequency Tables/gears_bug_bodysite_frequency_table.csv", row.names = FALSE)
  
  #These were what we originally wanted, but seem too long to classify 
  #Running same code but grouping by family
  
  atlas_bugfamily_bodysite <- atlas %>%
    dplyr::select(Family,Source) %>%
    group_by(Family,Source) %>%
    mutate(frequency = n()) %>%
    distinct() %>%
    ungroup %>%
    group_by(Family) %>%
    arrange(desc(frequency), .by_group = TRUE)
  # write.csv(atlas_bugfamily_bodysite, "./Frequency Tables/atlas_bugfamily_bodysite_frequency_table.csv", row.names = FALSE)
  
  gears_bugfamily_bodysite <- gears %>%
    dplyr::select(Family,BodySite) %>%
    group_by(Family,BodySite) %>%
    mutate(frequency = n()) %>%
    distinct() %>%
    ungroup %>%
    group_by(Family) %>%
    arrange(desc(frequency), .by_group = TRUE)
  # write.csv(gears_bugfamily_bodysite, "./Frequency Tables/gears_bugfamily_bodysite_frequency_table.csv", row.names = FALSE)
  
  
  ## 2.2 ATLAS/GEARS Country/Bug Tables ####
  
  # write.csv(atlas_countries_year, "./Frequency Tables/atlas_countries_year_frequency_table.csv")
  # write.csv(gears_countries_year, "./Frequency Tables/gears_countries_year_frequency_table.csv")
  
  #Bug/Country
  
  atlas_country_bug <- atlas %>%
    dplyr::select(Species,Country) %>%
    group_by(Species,Country) %>%
    mutate(frequency = n()) %>%
    distinct() %>%
    ungroup() %>%
    group_by(Country) %>%
    arrange(desc(frequency), .by_group = TRUE)
  # write.csv(atlas_country_bug, "./Frequency Tables/atlas_country_bug_frequency_table.csv", row.names = FALSE)
  
  gears_country_bug <- gears %>%
    dplyr::select(Organism,Country) %>%
    group_by(Organism,Country) %>%
    mutate(frequency = n()) %>%
    distinct() %>%
    ungroup() %>%
    group_by(Country) %>%
    arrange(desc(frequency), .by_group = TRUE)
  # write.csv(gears_country_bug, "./Frequency Tables/gears_country_bug_frequency_table.csv", row.names = FALSE)
  
  #BugFamily/Country, again because these tables are very long otherwise
  
  atlas_country_bugfamily <- atlas %>%
    dplyr::select(Family,Country) %>%
    group_by(Family,Country) %>%
    mutate(frequency = n()) %>%
    distinct() %>%
    ungroup() %>%
    group_by(Country) %>%
    arrange(desc(frequency), .by_group = TRUE)
  # write.csv(atlas_country_bugfamily, "./Frequency Tables/atlas_country_bugfamily_frequency_table.csv", row.names = FALSE)
  
  gears_country_bugfamily <- gears %>%
    dplyr::select(Family,Country) %>%
    group_by(Family,Country) %>%
    mutate(frequency = n()) %>%
    distinct() %>%
    ungroup() %>%
    group_by(Country) %>%
    arrange(desc(frequency), .by_group = TRUE)
  # write.csv(gears_country_bugfamily, "./Frequency Tables/gears_country_bugfamily_frequency_table.csv", row.names = FALSE)
  
  #Bug/Year
  
  gears_bugyear <- gears %>%
    dplyr::select(Organism,Year) %>%
    group_by(Organism,Year) %>%
    mutate(frequency = n()) %>%
    distinct() %>%
    ungroup() %>%
    group_by(Year) %>%
    arrange(desc(frequency), .by_group = TRUE)
  # write.csv(gears_bugyear, "./Frequency Tables/gears_bugyear_frequency_table.csv", row.names = FALSE)
  
  atlas_bugyear <- atlas %>%
    dplyr::select(Species,Year) %>%
    group_by(Species,Year) %>%
    mutate(frequency = n()) %>%
    distinct() %>%
    ungroup() %>%
    group_by(Year) %>%
    arrange(desc(frequency), .by_group = TRUE)
  # write.csv(atlas_bugyear, "./Frequency Tables/atlas_bugyear_frequency_table.csv", row.names = FALSE)
  
  #Bug/Year/Country
  
  gears_bugyear <- gears %>%
    dplyr::select(Organism,Year, Country) %>%
    group_by(Organism,Year, Country) %>%
    mutate(frequency = n()) %>%
    distinct() %>%
    ungroup() %>%
    group_by(Year) %>%
    arrange(desc(frequency), .by_group = TRUE)
  # write.csv(gears_bugyear, "./Frequency Tables/gears_country_bugyear_frequency_table.csv", row.names = FALSE)
  
  atlas_bugyear <- atlas %>%
    dplyr::select(Species,Year, Country) %>%
    group_by(Species,Year,Country) %>%
    mutate(frequency = n()) %>%
    distinct() %>%
    ungroup() %>%
    group_by(Year) %>%
    arrange(desc(frequency), .by_group = TRUE)
  # write.csv(atlas_bugyear, "./Frequency Tables/atlas_country_bugyear_frequency_table.csv", row.names = FALSE)

  atlas_bugyear_wider <- atlas_bugyear %>%
    pivot_wider(names_from = c(Country, Year), values_from = frequency)
  # write.csv(atlas_bugyear_wider, "./Frequency Tables/atlas_country_bugyear_frequency_table.csv", row.names = FALSE)
  
  gears_bugyear_wider <- gears_bugyear %>%
    pivot_wider(names_from = c(Country, Year), values_from = frequency)
  # write.csv(gears_bugyear_wider, "./Frequency Tables/gears_country_bugyear_frequency_table.csv", row.names = FALSE)
  
  # 2.3 Remove all but ATLAS/GEARS ####
  
  #fairly messy after all the frequency tables we wanted, so clearing environment (do not have to run)
  rm(list=setdiff(ls(),c("atlas", "gears")))

# 3. Filtering the dataset ####

  # Comments: Decision currently to focus on bloodstream infection MICs for carbapenem and 3rd generation cephalosporins 
  # against e.coli and klebsiella pneumoniae, and MRSA. 
  
  # If make a mistake in running atlas_blood or gears_blood have to rerun from section 3.1
    
  ## 3.1 Bloodstream infections only ####
  
  atlas_blood <- subset(atlas, atlas$Source=="Blood") #also blood vessels but I went with just blood
  nrow(atlas_blood) #171009
  
  gears_blood <- subset(gears, gears$BodySite=="CVS: Blood") #I think only one we want for BSI
  nrow(gears_blood) #3417
  
  
  ## 3.2 Adding country code (ISO) ####
  
  #changed to the blood datasets for now since that seems to be what we are going for, can add to overall dataset too
  
  atlas_blood$country_ISO <- countrycode(atlas_blood$Country, origin = 'country.name', destination = 'iso3c') #change iso3c if diff code needed
  table(atlas_blood$country_ISO, useNA = "always") #no NA
  
  gears_blood$country_ISO <- countrycode(gears_blood$Country, origin = 'country.name', destination = 'iso3c')
  table(gears_blood$country_ISO, useNA = "always") #no NA
  
  #especially in gears, some countries do not even have 50 entries, let alone for one of the bug-drug combos
  
  ## 3.3 Limiting abx and pathogens ####
  
  #GEARS abx of interest include: ceftazidime, imipenem, and meropenem
  #ATLAS of interest include: cefoxitin (just staph, and maybe oxacillin), Ceftriaxone, meropenem, imipenem
  ##also cefotaxime if it is in there but do not see it for now
  #if others can check we can add whichever ones I missed
  
  colnames(atlas_blood) #check which columns to keep
  atlas_blood <- atlas_blood[c(1:13,24:25,28:29,36:37,42:43,76:77,86:87,110:134,136)]
  colnames(atlas_blood) #verify columns are there
  
  colnames(gears_blood)
  gears_blood <- gears_blood[c(1:11,17,19,24:25)]
  colnames(gears_blood)
  
  ## 3.4 Frequency tables: country/year/bug/drug with filtered dataset ####
  
  #can add more bugs if we want but we said we would do the following:
  #
  
  gears_bugyear_blood <- gears_blood %>%
    dplyr::select(Organism,Year, Country) %>%
    filter(Organism == "Escherichia coli"|Organism=="Klebsiella pneumoniae")  %>% #don't see staph
    group_by(Organism,Year, Country) %>%
    mutate(frequency = n()) %>%
    distinct() %>%
    ungroup() %>%
    group_by(Year) %>%
    arrange(desc(frequency), .by_group = TRUE)
  # write.csv(gears_bugyear, "./Frequency Tables/gears_country_bugyear_frequency_table.csv", row.names = FALSE)
  gears_bugyear_blood
  #only US has over 50
  
  atlas_bugyear_blood <- atlas_blood %>%
    dplyr::select(Species,Year, Country) %>%
    filter(Species == "Escherichia coli"|Species=="Klebsiella pneumoniae"|
             Species == "Staphylococcus aureus")  %>%
    group_by(Species,Year,Country) %>%
    mutate(frequency = n()) %>%
    distinct() %>%
    ungroup() %>%
    group_by(Year) %>%
    arrange(desc(frequency), .by_group = TRUE)
  
  atlas_bugyear_blood
  #more above 50 here but still many below
  
# 4. Cleaning MIC using AMR for R ####
  
  ## 4.1 Filter GEARS ####
  
  #Apply AMR for R package to clean MICs in gears_blood
  #Realized that I did the filter tables above but have not filtered here yet
  
  gears_blood <- gears_blood %>%
    filter(Organism == "Escherichia coli"|Organism=="Klebsiella pneumoniae"|
             Organism == "Staphylococcus aureus")
  colnames(gears_blood)
  
  ### 4.1.2 GEARS CAZ, ceftazidime, no longer running but keeping for now ####
  # #remove for code that is uploaded for final
  # gears_blood$CAZ_MIC[gears_blood$CAZ_MIC=="-"] <- NA #if had a "-" I think where not tested
  # table(gears_blood$CAZ_MIC, useNA = "always")
  # gears_blood$CAZ_MIC_original <- gears_blood$CAZ_MIC
  # gears_blood$CAZ_MIC <- as.numeric(gears_blood$CAZ_MIC)
  # #you may want to change these? just what I had for now
  # gears_blood$CAZ_MIC[gears_blood$CAZ_MIC=="16.0001"] <- ">16"
  # gears_blood$CAZ_MIC[gears_blood$CAZ_MIC=="32.0001"] <- ">32"
  # gears_blood$CAZ_MIC[gears_blood$CAZ_MIC=="128.0001"] <- ">128"
  # gears_blood$CAZ_MIC[gears_blood$CAZ_MIC=="0.0299"] <- "<0.03"
  # gears_blood$CAZ_MIC[gears_blood$CAZ_MIC=="0.2499"] <- "<0.25"
  # table(gears_blood$CAZ_MIC, useNA = "always")
  # gears_blood$CAZ_MIC_clean <- as.mic(gears_blood$CAZ_MIC, na.rm = FALSE)
  # # table(gears_blood$CAZ_MIC_clean) #doesn't seem too helpful since includes all possible levels, will use dplyr instead for others
  # gears_blood %>% count(CAZ_MIC_clean) %>% filter(n > 0) 
  # sum(is.na(gears_blood$CAZ_MIC_clean)) #0
  # gears_blood$CAZ_MIC_SIR <- as.sir(x=gears_blood$CAZ_MIC_clean, mo=as.mo(gears_blood$Organism), ab="CAZ", guideline="EUCAST", include_PKPD = FALSE)
  # table(gears_blood$CAZ_MIC_SIR, useNA = "always") #0 coming from PK/PD breakpoints that are left as NA
  # 
  # #table with SIR, breakpoints, and original MIC
  # 
  # gears_CAZ_table <- gears_blood %>% 
  #   group_by(Organism) %>%
  #   count(CAZ_MIC_original, CAZ_MIC, CAZ_MIC_clean, CAZ_MIC_SIR) %>% 
  #   filter(n > 0) %>%
  #   ungroup()
  
  ### 4.1.3 GEARS MEM, meropenem ####
  #need to fix somethign for the body site here possibly
  gears_blood$MEM_MIC[gears_blood$MEM_MIC=="-"] <- NA #some had a "-" I think where not tested
  gears_blood$MEM_MIC <- as.numeric(gears_blood$MEM_MIC)
  table(gears_blood$MEM_MIC, useNA = "always") 
  gears_blood$MEM_MIC_original <- gears_blood$MEM_MIC
  #we may want to change these? just what I had for now
  gears_blood$MEM_MIC[gears_blood$MEM_MIC=="0.0599"] <- "<0.06"
  gears_blood$MEM_MIC[gears_blood$MEM_MIC=="8.0001"] <- ">8"
  gears_blood$MEM_MIC[gears_blood$MEM_MIC=="32.0001"] <- ">32"
  gears_blood$MEM_MIC[gears_blood$MEM_MIC=="64.0001"] <- ">64"
  table(gears_blood$MEM_MIC, useNA = "always")
  gears_blood$MEM_MIC_clean <- as.mic(gears_blood$MEM_MIC, na.rm = FALSE)
  # table(gears_blood$MEM_MIC_clean) #doesn't seem too helpful since includes all possible levels, will use dplyr instead for others
  gears_blood %>% count(MEM_MIC_clean) %>% filter(n > 0) #seems to match
  sum(is.na(gears_blood$MEM_MIC_clean)) #0
  gears_blood$MEM_MIC_SIR <- as.sir(x=gears_blood$MEM_MIC_clean, mo=as.mo(gears_blood$Organism), ab="MEM", guideline="EUCAST", include_PKPD = FALSE)
  table(gears_blood$MEM_MIC_SIR, useNA = "always") #0 coming from PK/PD breakpoints that are left as NA
  #assues body site "Non-meningitis"
  
  #table with SIR, breakpoints, and original MIC
  gears_MEM_table <- gears_blood %>%
    group_by(Organism) %>%
    count(MEM_MIC_original, MEM_MIC, MEM_MIC_clean, MEM_MIC_SIR) %>%
    filter(n > 0) %>%
    ungroup()
  
  
  ### 4.1.4 GEARS IPM, Imipenem #### 
  #need to fix somethign for the body site here possibly
  table(gears_blood$IPM_MIC, useNA = "always")
  gears_blood$IPM_MIC[gears_blood$IPM_MIC=="-"] <- NA #all are NA
  gears_blood$IPM_MIC_SIR <- as.sir(x=gears_blood$IPM_MIC, mo=as.mo(gears_blood$Organism), ab="IPM", guideline="EUCAST", include_PKPD = FALSE)
  table(gears_blood$IPM_MIC_SIR, useNA = "always") #all still NA, just left here if we change it
  
  
  ## 4.2 ATLAS MICs ####
  #GEARS abx of interest include: ceftazidime, imipenem, and meropenem
  #ATLAS of interest include: cefoxitin (just staph, and maybe oxacillin), Ceftriaxone, meropenem, imipenem
  ##also cefotaxime if it is in there but do not see it for now
  #if others can check we can add whichever ones I missed
  
  atlas_blood <- atlas_blood %>%
    filter(Species == "Escherichia coli"|Species=="Klebsiella pneumoniae"|
             Species == "Staphylococcus aureus")
  colnames(atlas_blood)
  
  ### 4.2.1 ATLAS Meropenem (carbapenem 1) ####
  table(atlas_blood$Meropenem, useNA = "always") #some blank
  table(atlas_blood$Meropenem_I, useNA = "always") #different numbers blank
  
  #check if MIC
  is.mic(atlas_blood$Meropenem) #FALSE
  
  #convert so will be same format as GEARS
  #atlas_blood$Meropenem[atlas_blood$Meropenem=="0.1200"] <- "0.12"
  #atlas_blood$Meropenem[atlas_blood$Meropenem=="0.2500"] <- "0.25"
  #atlas_blood$Meropenem[atlas_blood$Meropenem=="0.5000"] <- "0.5"
  #atlas_blood$Meropenem[atlas_blood$Meropenem=="1.0000"] <- "1"
  #atlas_blood$Meropenem[atlas_blood$Meropenem=="16.0000"] <- "16"
  #atlas_blood$Meropenem[atlas_blood$Meropenem=="2.0000"] <- "2"
  #atlas_blood$Meropenem[atlas_blood$Meropenem=="4.0000"] <- "4"
  #atlas_blood$Meropenem[atlas_blood$Meropenem=="8.0000"] <- "8"
  
  #create new MIC
  atlas_blood$Meropenem_MIC <- as.mic(atlas_blood$Meropenem, na.rm = FALSE)
  atlas_blood$Meropenem_SIR_EUCAST <- as.sir(x=atlas_blood$Meropenem_MIC, mo=as.mo(atlas_blood$Species), ab="MEM", guideline="EUCAST", include_PKPD = FALSE)
  atlas_blood$Meropenem_SIR_CLSI <- as.sir(x=atlas_blood$Meropenem_MIC, mo=as.mo(atlas_blood$Species), ab="MEM", guideline="CLSI", include_PKPD = FALSE)
  #assuming non-meningitis again
  
  #verify counts
  atlas_MEM_table <- atlas_blood %>% 
    group_by(Species) %>%
    count(Meropenem, Meropenem_MIC, Meropenem_I, Meropenem_SIR_EUCAST, Meropenem_SIR_CLSI) %>% 
    filter(n > 0) %>% #seems to match
    ungroup()
  
  ### 4.2.2 ATLAS Imipenem ####
  atlas_blood$Imipenem_MIC <- as.mic(atlas_blood$Imipenem, na.rm = FALSE)
  
  atlas_blood$Imipenem_SIR_EUCAST <- as.sir(x=atlas_blood$Imipenem_MIC, mo=as.mo(atlas_blood$Species), ab="IPM", guideline="EUCAST", include_PKPD = FALSE)
  atlas_blood$Imipenem_SIR_CLSI <- as.sir(x=atlas_blood$Imipenem_MIC, mo=as.mo(atlas_blood$Species), ab="IPM", guideline="CLSI", include_PKPD = FALSE)
  #assuming non-meningitis again
  
  #verify counts
  atlas_IPM_table <- atlas_blood %>% 
    group_by(Species) %>%
    count(Imipenem, Imipenem_MIC, Imipenem_I, Imipenem_SIR_EUCAST, Imipenem_SIR_CLSI) %>% 
    filter(n > 0) %>% #seems to match
    ungroup()

  ### 4.2.3 ATLAS Cefoxitin ####
  atlas_blood$Cefoxitin_MIC <- as.mic(atlas_blood$Cefoxitin, na.rm = FALSE)
  
  atlas_blood$Cefoxitin_SIR_EUCAST <- as.sir(x=atlas_blood$Cefoxitin_MIC, mo=as.mo(atlas_blood$Species), ab="CFX", guideline="EUCAST", include_PKPD = FALSE)
  atlas_blood$Cefoxitin_SIR_CLSI <- as.sir(x=atlas_blood$Cefoxitin_MIC, mo=as.mo(atlas_blood$Species), ab="CFX", guideline="CLSI", include_PKPD = FALSE)
  #assuming non-meningitis again
  
  # clinical_breakpoints %>% filter(mo == "B_STPHY_AURS", ab == "FOX") %>% print(n=50)
  
  #verify counts
  atlas_CFX_table <- atlas_blood %>%
    group_by(Species) %>%
    count(Cefoxitin, Cefoxitin_MIC, Cefoxitin_I, Cefoxitin_SIR_EUCAST, Cefoxitin_SIR_CLSI) %>%
    filter(n > 0) %>% #seems to match
    ungroup()
  # 
  ### 4.2.4 ATLAS Oxacillin ####
  atlas_blood$Oxacillin_MIC <- as.mic(atlas_blood$Oxacillin, na.rm = FALSE)
  
  atlas_blood$Oxacillin_SIR_EUCAST <- as.sir(x=atlas_blood$Oxacillin_MIC, mo=as.mo(atlas_blood$Species), ab="OXA", guideline="EUCAST", include_PKPD = FALSE)
  atlas_blood$Oxacillin_SIR_CLSI <- as.sir(x=atlas_blood$Oxacillin_MIC, mo=as.mo(atlas_blood$Species), ab="OXA", guideline="CLSI", include_PKPD = FALSE)
  #assuming non-meningitis again
  
  #verify counts
  atlas_OXA_table <- atlas_blood %>% 
    group_by(Species) %>%
    count(Oxacillin, Oxacillin_MIC, Oxacillin_I, Oxacillin_SIR_EUCAST, Oxacillin_SIR_CLSI) %>% 
    filter(n > 0) %>% #seems to match
    ungroup()
  
  ### 4.2.5 ATLAS Ceftriaxone ####
  atlas_blood$Ceftriaxone_MIC <- as.mic(atlas_blood$Ceftriaxone, na.rm = FALSE)
  
  atlas_blood$Ceftriaxone_SIR_EUCAST <- as.sir(x=atlas_blood$Ceftriaxone_MIC, mo=as.mo(atlas_blood$Species), ab="CRO", guideline="EUCAST", include_PKPD = FALSE)
  atlas_blood$Ceftriaxone_SIR_CLSI <- as.sir(x=atlas_blood$Ceftriaxone_MIC, mo=as.mo(atlas_blood$Species), ab="CRO", guideline="CLSI", include_PKPD = FALSE)
  #assuming non-meningitis again
  
  #verify counts
  atlas_CRO_table <- atlas_blood %>% 
    group_by(Species) %>%
    count(Ceftriaxone, Ceftriaxone_MIC, Ceftriaxone_I, Ceftriaxone_SIR_EUCAST, Ceftriaxone_SIR_CLSI) %>% 
    filter(n > 0) %>% #seems to match
    ungroup()
  
  ### 4.2.6 ATLAS Doripenem ####
  atlas_blood$Doripenem_MIC <- as.mic(atlas_blood$Doripenem, na.rm = FALSE)
  
  atlas_blood$Doripenem_SIR_EUCAST <- as.sir(x=atlas_blood$Doripenem_MIC, mo=as.mo(atlas_blood$Species), ab="DOR", guideline="EUCAST", include_PKPD = FALSE)
  atlas_blood$Doripenem_SIR_CLSI <- as.sir(x=atlas_blood$Doripenem_MIC, mo=as.mo(atlas_blood$Species), ab="DOR", guideline="CLSI", include_PKPD = FALSE)
  #assuming non-meningitis again
  
  #verify counts
  atlas_DOR_table <- atlas_blood %>% 
    group_by(Species) %>%
    count(Doripenem, Doripenem_MIC, Doripenem_I, Doripenem_SIR_EUCAST, Doripenem_SIR_CLSI) %>% 
    filter(n > 0) %>% #seems to match
    ungroup()
  
  ### 4.2.7 ATLAS Tebipenem ####
  atlas_blood$Tebipenem_MIC <- as.mic(atlas_blood$Tebipenem, na.rm = FALSE)
  
  atlas_blood$Tebipenem_SIR_EUCAST <- as.sir(x=atlas_blood$Tebipenem_MIC, mo=as.mo(atlas_blood$Species), ab="TBP", guideline="EUCAST", include_PKPD = FALSE)
  atlas_blood$Tebipenem_SIR_CLSI <- as.sir(x=atlas_blood$Tebipenem_MIC, mo=as.mo(atlas_blood$Species), ab="TBP", guideline="CLSI", include_PKPD = FALSE)
  #assuming non-meningitis again
  
  #verify counts
  atlas_TBP_table <- atlas_blood %>% 
    group_by(Species) %>%
    count(Tebipenem, Tebipenem_MIC, Tebipenem_I, Tebipenem_SIR_EUCAST, Tebipenem_SIR_CLSI) %>% 
    filter(n > 0) %>% #seems to match
    ungroup()
  
  ## 4.3 Export tables #### 
  #table list since forgot to as going through
  #table_list <- mget(ls(pattern = "_table"))
  
  #for (i in 1:length(table_list)) {
   # write.csv(table_list[i], paste0(names(table_list[i]),"_SIR.csv"))
 # }
  
 # cb <- clinical_breakpoints
 # write.csv(cb, "./AMR_for_R_clinical_breakpoints_export.csv")
  
  #Comment: ATLAS seems to use CLSI SIR values. So, important we use our own with the EUCAST system.
  #None of the values here have 
  
  ## 5 Appending ATLAS and GEARS and then merging with covariates
  
  # 5.1 Appending ATLAS and GEARS
  
  # 5.1.1 Checking structure and variable names of the two data frames
  glimpse(atlas_blood)
  glimpse(gears_blood)
  
  # 5.1.2 adding "dataset" column to both dataframes so that we can identify 
  # ATLAS vs GEARS in the appended dataset
  atlas_blood <- atlas_blood %>%
    mutate(dataset = rep("atlas", nrow(atlas_blood)))
  gears_blood <- gears_blood %>%
    mutate(dataset = rep("gears", nrow(gears_blood)))
  
  
  # 5.1.3 Removing redundant variables
  atlas_blood_final <- atlas_blood %>%
    dplyr::select(Isolate.Id, Species, Family, Country, Gender, Age.Group,
                  Speciality, In...Out.Patient, Year, country_ISO, Meropenem_MIC,
                  Meropenem_SIR_EUCAST, Imipenem_MIC, Imipenem_SIR_EUCAST, 
                  Cefoxitin_MIC, Cefoxitin_SIR_EUCAST, Oxacillin_MIC, Oxacillin_SIR_EUCAST,
                  Ceftriaxone_MIC, Ceftriaxone_SIR_EUCAST, Doripenem_MIC, Doripenem_SIR_EUCAST,
                  Tebipenem_MIC, Tebipenem_SIR_EUCAST, dataset)
  
  gears_blood_final <- gears_blood %>%
    dplyr::select(Isolate, Year, Organism, Family, Region, Country, Gender,
           Facility, Age.Group, country_ISO, MEM_MIC_clean, MEM_MIC_SIR,
           IPM_MIC_SIR, dataset) # keeping meropenem and imipenem SIR just in case
  
  # 5.1.4 Renaming common variables in GEARS variables to be the same as ATLAS
  gears_blood_final <- gears_blood_final %>%
    rename(Isolate.Id = Isolate,
           Species = Organism,
           Speciality = Facility,
           Meropenem_MIC = MEM_MIC_clean,
           Meropenem_SIR_EUCAST = MEM_MIC_SIR,
           Imipenem_SIR_EUCAST = IPM_MIC_SIR)
  
  # 5.1.5 Identifying and resolving differences in variable TYPES for common
  # variables
  str(atlas_blood_final)
  str(gears_blood_final)
  
  unique(atlas_blood_final$Gender)
  unique(gears_blood_final$Gender)
  
  # Year is integer in ATLAS and numeric in GEARS (probably not an issue but
  # will make it an integer in both dataframes anyway)
  
  # Gender is F/M in GEARS and Male/Female in ATLAS. Will change F/M in GEARS
  # to Male/Female
  
  # Changing Year from numeric to integer in GEARS
  gears_blood_final$Year <- as.integer(gears_blood_final$Year)
  
  # Converting "" to NA for Gender in ATLAS
  atlas_blood_final$Gender <- na_if(atlas_blood_final$Gender, "")
  
  # Converting "N" to NA for Gender in GEARS
  gears_blood_final$Gender <- na_if(gears_blood_final$Gender, "N")
  
  # Changing Gender from F/M to Female/Male in GEARS 
  gears_blood_final$Gender <- recode(gears_blood_final$Gender,
                                     "F" = "Female",
                                     "M" = "Male")
  
  # 5.1.6 Identifying and resolving differences in values for common variables
  # (ignoring the MIC/SIR variables Madison has already done a lot of checking
  # here)
  lapply(gears_blood_final, unique)
  lapply(atlas_blood_final, unique)
  
  # Year
  unique(gears_blood_final$Year)
  unique(atlas_blood_final$Year)
  # Looks fine and both are integer now
  
  # Species
  unique(gears_blood_final$Species)
  unique(atlas_blood_final$Species)
  # Looks fine (note no S.aureus isolates in GEARS)
  
  #Family
  unique(gears_blood_final$Family)
  unique(atlas_blood_final$Family)
  # Looks fine (note that we have Enterobacterales in ATLAS which is really an 
  # order of bacteria not a family 
  
  # ISO-3 (ignoring country name because all merging will be done on ISO-3. We
  # could convert ISO-3 back to country name to get a consistent set of country
  # names)
  unique(gears_blood_final$country_ISO)
  unique(atlas_blood_final$country_ISO)
  # Looks fine (trusting in the countrycode package!)
  
  # Age group
  unique(gears_blood_final$Age.Group)
  unique(atlas_blood_final$Age.Group)
  # Looks fine
  
  # Speciality (probably won't use)
  unique(gears_blood_final$Speciality)
  unique(atlas_blood_final$Speciality)
  # Looks fine
  
  # 5.1.7 Calling dplyr's bind_rows() to append ATLAS and GEARS isolates
  atlas_gears <- bind_rows(atlas_blood_final, gears_blood_final)

  
  # 5.2 Loading main covariates file and merging with GRAM antibiotic use data
  # loading main country-level covariates
  
  
  # 5.2.2 Loading GRAM data on antibiotic use by ATC3
  load("gram.did.RDATA")
  gram <- x
  rm(x)
  
  # adding lags for total antibiotic consumption, consumption of "other" 
  # beta-lactam antibiotics and consumption of Watch antibiotics
  gram <- gram %>%
    group_by(iso3_code) %>%
    arrange(year) %>%
    mutate(did_total_lag = lag(did_total, 1),
           j01d_beta_lag = lag(j01d_beta, 1),
           did_watch_lag = lag(did_watch, 1))
  
  # Loading income data
  load("incomedata.RDATA")
  income <- x
  rm(x)

  
  # 5.2.3 Joining covariates and income data
  data <- full_join(gram, income, by = c("iso3_code", "year"))
  
  # 5.3 Joining ATLAS-GEARS with covariates
  
  # 5.3.1 Renaming country code and year in atlas_gears to allow join with 
  # covariates
  atlas_gears <- atlas_gears %>%
    rename(iso3_code = country_ISO,
           year = Year)
  
  # 5.3.2 Joining atlas_gears with covariates
  df <- full_join(atlas_gears, data, by = c("iso3_code", "year"))

  # 5.3.3 Filtering for K.pneumoniae isolates from 2006-2018 since we only have 
  # clean income data for these years
  df <- df %>%
    filter(Species == "Klebsiella pneumoniae", year >= 2006, year <= 2018)
  
  # 5.3.4 Removing redundant columns and reordering identifying columns to have 
  # a more logical structure
  df <- df %>%
    select(-country_name, Isolate.Id, Country, iso3_code, year, Age.Group, Gender, 
           Speciality, In...Out.Patient, Meropenem_MIC: income)

  # getting rid of Species and Family because we are just looking at K.pneumoniae
  # there was a duplicate country name column which I have removed
  
  
  
  #------------------------------------------------------------------------------#
  #------------------------------------------------------------------------------#
  
  # 1. Data preparation ####
  
  # 1.1 Select data for analysis  
  
  # 1.2 Recode/clean variables
  
  # Set factor variables
  df$Gender <- as.factor(df$Gender)
  df$Age.Group <- as.factor(df$Age.Group)
  
  
  #Recode Meropenem_MIC  
  #Strategy: remove >8 from the sample; regroup the values below/equals to 1 to "<=1"; regroup "32" "64" ">64" to ">16"
  unique(df$Meropenem_MIC)
  df$Meropenem_MIC <- as.character(df$Meropenem_MIC)
  
  df <- df %>%
    filter(Meropenem_MIC != ">8") #obs=6358
  
  df <- df %>%
    mutate (Meropenem_MIC = case_when(
      Meropenem_MIC == "0.0039"  | Meropenem_MIC == "0.008"  | Meropenem_MIC == "0.015"  | Meropenem_MIC == "0.03" | Meropenem_MIC == "0.06" | Meropenem_MIC == "<=0.06"| Meropenem_MIC == "0.12" | Meropenem_MIC == "0.25" | Meropenem_MIC == "0.5" | Meropenem_MIC == "1"~ "<=1",
      Meropenem_MIC == "32" ~ ">16",
      Meropenem_MIC == "64" ~ ">16",
      Meropenem_MIC == ">64" ~ ">16",
      TRUE ~ Meropenem_MIC
    ))
  
  sort <- sort(unique(df$Meropenem_MIC))
  sort 
  
  table(df$Meropenem_MIC)
  
  df$Meropenem_MIC <- factor(df$Meropenem_MIC,
                             levels = c("<=1", "2", "4", "8", "16", ">16"),
                             ordered=TRUE)
  table(df$Meropenem_MIC)
  
  # 1.3 Creating two separate data frames, one with country-year pairs >= 30 
  # the other >= 50 (will probably use >= 30)
  
  # calculating country-year counts based on iso3_code
  df <- df %>%
    group_by(iso3_code, year) %>%
    mutate(count = n()) %>%
    ungroup()
  
  
  # Dropping empty low-income factor and reversing the order of the remaining
  # levels to be more intuitive (income = 1 for lower-middle, income = 2 for upper
  # middle, income = 3 for high-income)
  df$income <- factor(c(df$income), levels = c("Low income", 
                                               "Lower middle income", 
                                               "Upper middle income", 
                                               "High income"))
  
  
  # filtering for country-year pairs >30 and >50 (we are probably going to use
  # the 30 threshold)
  df_30 <- df %>%
    filter(count >= 30)
  table(df_30$iso3_code)
  
  df_50 <- df %>%
    filter(count >= 50)
  table(df_50$iso3_code)
  
  # 1.3 Create datasets for high and middle income countries
  
  df_high <- df %>%
    filter(income=="High income")
  
  df_mid <- df %>%
    filter(income=="Lower middle income" | income== "Upper middle income")
  
  df_upmid <- df %>%
    filter(income== "Upper middle income")
  
  df_lowmid <- df %>%
    filter(income=="Lower middle income")
  
  nrow(df_high) # n=14058
  nrow(df_mid) # n=5481
  nrow(df_upmid) # n=4146
  nrow(df_lowmid) # n=1335
  
  
  # saving as an R workspace file
  save(atlas_gears, df, df_30, df_50, file = "cleaned_data.RData")
  