This GitHub respository contains the data files and R script necessary to replicate the distributional copula analysis of the relationships between country-level income and antibiotic consumption and Meropenem resistance in Klebsiella pneumoniae bloodstream isolates, for the Vivli 2024 student data challenge submission

There are five data files
1. "cleaned_data.RData" contains the cleaned R environment used for the modelling.
2. "2024_05_28 atlas_antibiotics.7z" is a compressed version of the raw ATLAS surveillance data.
3. "Venatorx surveillance data_2024_06_06.xlsx" is the raw GEARS surveillance data
4. "incomedata.RDATA" is the country-level GNI per capita data (World Bank)
5. "gram.did.RDATA" is the antibiotic consumption data (GRAM)

There are two R script files:
1. "Data cleaning.R" which contains code to harmonise and merge ATLAS, GEARS, the income data, and the antibiotic data
2. "Modelling.R" which contains code to fit the distributional copulas and produce the partial effects plots

The modelling code can be run using "cleaned_data.RData". Alternatively, you can produce "cleaned_data.RData" yourself using  "Data cleaning.R" on the raw data files above. With either approach, all that should need to be done is set the working directory, and install any missing R packages.
