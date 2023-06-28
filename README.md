# UniProtExtractR
`UniProtExtractR::uniprotextract()`
### R package and Shiny App for extracting desired information from syntactically structured character strings across 9 UniProtKB categories: DNA binding, Pathway, Transmembrane, Signal peptide, Protein families, Domain [FT], Motif, Involvement in disease, and Subcellular location [CC]


## Overview
UniProtKB protein entries can be challenging to handle for data analysis, visualization, and categorization due to character string syntax. This app and R package extracts desired information from syntactically complex UniProtKB entries for the 9 UniProtKB specific categories above. Further, a user can supply a mapping file for Subcellular location [CC], as many annotated locations are quite specific. For example, if a user is studying a particular organelle, it may be useful to broadly classify all other organelles (“Cytoplasm”, “Nucleus”, etc.) while maintaining resolution within one single organelle (“Inner mitochondrial membrane”, “Mitochondrial matrix”, etc.) The resolution across organelles is user defined. An example mapping file has been provided to demonstrate utility. 

The R package `UniProtExtractR` has one function, `uniprotextract()`, which takes 3 arguments detailed below in the Upload files section. The Shiny app shares the first 2 inputs.  

## Upload files
### UniProtKB query file
The only requisite for the app and R package is a UniProtKB query file. This should be a tab separated (.TSV or .TXT). UniProtKB queries can be exported directly from UniProtKB as .TSV files, which is recommended. For only the categories to be extracted by the R package or app, those column names must exactly match those exported from UniProtKB, but do not have to be in any particular order. If no column matches the aforementioned 9 categories, then nothing in the uploaded file will be changed. Any column matching the aforementioned categories will be modified (“extracted”) and added as a new column adjacent the old column in the output. An example query file with properly named columns has been provided. There is no query file size limit for the R package, but there is a 5 GB size limit for the Shiny app (feel free to run the Shiny app locally and increase the size limit in the `shiny.maxRequestSize` argument of `server` `options`).  It is **highly recommended not to open a UniProtKB exported query file in Excel**. There are some funny things that can happen, including entries getting cut off and disrupting the structure of the entire data file (as of release 1.0.0 of UniProtExtractR, there is a character limit per Excel cell of 32,767 characters, which is shorter, for example, than the amino acid sequence of the Human protein Titin). 

### Mapping file
The mapping file is optional and should be tab separated if included (.TSV or .TXT). The mapping file is only for the category Subcellular location [CC]. The app and R package will assume no mapping file is present unless provided. An example mapping file has been provided. It should contain only 3 columns: “New”, “Exact”, and “StartsWith”. For every Subcellular location [CC] per entry in the UniProtKB query file, any terms that match exactly to “Exact” or that start with the phrases in “StartsWith” will be swapped for the value in “New” in a row-wise fashion. Each row represents paired “Exact” and “StartsWith” matches for a single replacement value “New”. The “Exact” and “StartsWith” values may have many listed phrases. For example, if we want to replace all exact matches of “Apical cell membrane” and “Lateral cell membrane” with “Cell membrane”, then one row of the mapping file would appear as:

| New | Exact | StartsWith |
| ----------- | ----------- | ----------- |
| Cell membrane| Apical cell membrane,Lateral cell membrane |   |

Importantly, the values within “Exact” and “StartsWith” columns must be separated by a single comma without a space. Notice there is only a comma and no space between “Apical cell membrane” and “Lateral cell membrane”. Ensure there are no hanging terminal commas at the end of your list of terms for “Exact” and “StartsWith”. Blank values for “Exact” and “StartsWith” columns are fine; if a blank value exists in “New”, then all the “Exact” and “StartsWith” matches will be replaced with a blank. “Exact” matches are always prioritized over “StartsWith” matches in the case multiple matching mappings exist. 

### Option `write.local` (R package only)
Default is `FALSE`, and the uniprotextract function will return the modified dataframe. If `write.local = TRUE`, then no dataframe is returned but instead written to the working directory as a .CSV file. 

## Outputs
Both the app and the R package return a dataframe containing the originally uploaded UniProtKB file with additional modified or “extracted” columns. Both the original and modified columns will be present in the resulting file for user comparison. For the categories DNA binding, Transmembrane, and Signal peptide, a binary category is created indicating whether a protein entry has that feature or not. Additionally, for DNA binding and Transmembrane categories specifically, UniProtExtractR also adds a numeric column counting the number of domains present per entry. For Pathway, the R package removes trailing details for each entry and changes the terms “biosynthesis”, “degradation”, and “metabolism” to “modification” for increased coverage at the expense of specificity. For Protein families and Domain [FT], trailing details for each entry are removed. For Involvement in disease, the name of the disease is extracted for each entry and additional details are removed. For Motif, the motif is extracted. For Subcellular location [CC], a column containing the extracted location is added, and if a mapping file is present, an additional column is adding containing the mapped term. The original Subcellular location [CC] entry, extracted, and mapped columns will allow for manual inspection. If a protein entry contains multiple values for Pathway, Protein families, Domain [FT], Motif, Involvement in disease, and Subcellular location [CC], then only the first listed element is extracted. Any isoform specific terms are removed for any category (a relatively small amount for all tested datasets during app and package development). (A user may edit some of these options by redefining the `uniprotextract` function locally; there are a few documented regions in the function that allow for expansion to multiple terms per entry if desired.)

If using the Shiny app, then the outputs are a frequency table and a modified dataframe that can be downloaded as a .CSV file. The frequency table will show how many times unique values appear for each column of the modified UniProtKB query file. This may be useful for seeing overall trends of the input UniProtKB queries (How many proteins have 7 transmembrane domains? How many proteins mitochondrial? etc.) and especially useful for refining the organelle mapping file for Subcellular location [CC] (How many proteins have I assigned to each “New” category?). 

If using `uniprotextract()` from the R package, it will return a modified or “extracted” version of the original UniProtKB query dataframe. Please see `?uniprotextract` for reproducible short examples.

If the following categories are present (and not empty), then new columns (bulleted) will appear:
1. DNA binding
- DNA.binding_binary
- DNA.binding_count
2. Pathway
- Pathway.edit
3. Transmembrane
- Transmembrane_binary
- Transmembrane_count
4. Signal peptide
- Signal.peptide_binary
5. Protein families
- Protein.families.edit
6. Domain [FT]
- Domain..FT.edit
7. Motif
- Motif.edit
8. Involvement in disease
- Involvement.in.disease.edit
- Involvement.in.disease_count
9. Subcellular location [CC]
- Subcellular.location..CC.edit
- Subcellular.location..CC.map.edit (if mapping file present)

Note that R converts spaces in column names to periods in the output.


## Troubleshooting

### Weird long strings in weird places
While testing the app among multiple users, the most common issue was opening an exported file from UniProtKB in Excel. This led to all sorts of trouble. If the app runs and there are large character string chunks in weird places that look like full paragraphs, this is the likely issue. When loading data in R locally and using `read.table` for .TSV UniProtKB query files, use `fill = TRUE` and `quote = “”` so that entries are appropriately separated. 

### Mapping not working
If a mapping file is used and all terms map to one particular value of the “New” column, then there may be a hanging comma somewhere. Conversely, if Subcellular location [CC] terms remain unchanged, then there may be a space between your “Exact” and “StartsWith” values within each tab-separated element in the mapping file. 

### Empty cells/values in original file are being “extracted” as incorrect non-empty cells/values
Check the “empty” cells in the uploaded file. If using the Shiny app, consider that some cells may look empty, but really contain a single or multiple spaces. If there are cells that are “NA” or “NULL”, they will be appropriately assigned as such if the text in those cells matches one of the following: “”, “NA”, “NULL”, “null” or “NaN”. If using the R package, similarly check empty or NA values; `NA`, `NULL`, `null` or `NaN` and empty cells are all converted to `NA`. Of note, UniProtKB queries contain many types of quotation marks; to be considered if using `read.table`.

### Edge cases
There were quite a few edge cases of UniProtKB syntax that are all addressed in this release, but if more are encountered, please let me know!

## Install
```
devtools::install_github("alex-bio/UniProtExtractR")
library(UniProtExtractR)
UniProtExtractR::uniprotextract(my.uniprot.df, map.up=NULL, write.local=FALSE)
```

Last edited Alex Panov 23-06-28.
