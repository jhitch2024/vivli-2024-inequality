This GitHub respository contains the data files and R script necessary to replicate the distributional copula analysis of the relationships between country-level income and antibiotic consumption and Meropenem resistance in Klebsiella pneumoniae bloodstream isolates, for the Vivli 2024 student data challenge submission.

There are four data files
1. "2024_05_28 atlas_antibiotics.7z" is a compressed version of the raw ATLAS surveillance data.
2. "Venatorx surveillance data_2024_06_06.xlsx" is the raw GEARS surveillance data
3. "incomedata.RDATA" is the country-level GNI per capita data (World Bank)
4. "gram.did.RDATA" is the antibiotic consumption data (GRAM)

There are two R script files:
1. "Data cleaning.R" which contains code to harmonise and merge ATLAS, GEARS, the income data, and the antibiotic data
2. "Modelling.R" which contains code to fit the distributional copulas and produce the partial effects plots

All that should need to be done is set the working directory, and install any missing R packages.
