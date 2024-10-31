This GitHub respository contains the data files and R script necessary to replicate the distributional copula analysis of the relationships between country-level income and antibiotic consumption and Meropenem resistance in Klebsiella pneumoniae bloodstream isolates, for the Vivli 2024 student data challenge submission

There are five data files
1. "cleaned_data.RData" contains the cleaned R environment used for the modelling.
2. "2024_05_28 atlas_antibiotics.7z" is a compressed version of the raw ATLAS surveillance data.
3. "Venatorx surveillance data_2024_06_06.xlsx" is the raw GEARS surveillance data
4. "incomedata.RDATA" contains publicly available country-level income (GNI per capita) data
5. "gram.did.RDATA" contains publicly available antibiotic consumption data

GNI per capita using atlas method was downloaded on 8 July 2024 for the years 1998-2023 from here: https://databank.worldbank.org/source/world-development-indicators and I applied the World Bank historical GNI per capita cutoffs for income group from here: https://datahelpdesk.worldbank.org/knowledgebase/articles/906519-world-bank-country-and-lending-groups (historical excel file) and the GRAM antibiotic usage estimates were downloaded from here: https://www.tropicalmedicine.ox.ac.uk/gram/research/visualisation-app-antibiotic-usage-and-consumption with this citation  please cite this publication as a reference: Browne AJ, Chipeta MG, Haines-Woodhouse G, et al. Global antibiotic consumption and usage in humans, 2000 to 2018: a spatial modelling study. Lancet Planetary Health 2021

There are two R script files:
1. "Data cleaning.R" contains code to harmonise and merge ATLAS, GEARS, the income data, and the antibiotic data
2. "Modelling.R" contains code to fit the distributional copulas and produce the partial effects plots

The modelling code can be run using "cleaned_data.RData". Alternatively, you can produce "cleaned_data.RData" yourself using  "Data cleaning.R" on the raw data files above. With either approach, all that should need to be done is set the working directory, and install any missing R packages.
