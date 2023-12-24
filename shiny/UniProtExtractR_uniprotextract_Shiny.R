## UniProtExtractR. AP 23-10-02.
#### Install necessary packages ####
list.of.packages <- c("stringr", "stringi", "tibble", "shiny", "DT", "BiocManager", "UniProt.ws", "bslib")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library(stringr)
library(stringi)
library(tibble)
library(shiny)
library(DT)
library(BiocManager)
library(UniProt.ws)
library(bslib)

# markdown.instructions <- source("https://github.com/alex-bio/UniProtExtractR/blob/main/README.md")
####Start of Shiny App. ####
ui <- shinyUI(navbarPage("UniProtExtractR",
  theme = bs_theme(bootswatch = "cerulean"),
  # titlePanel("UniProtExtractR: extract information from UniProtKB"),
  ## this is for making the withProgress() larger and centered
  tags$head(tags$style(HTML(".shiny-notification {
              height: 100px;
              width: 500px;
              position:fixed;
              top: calc(50% - 50px);;
              left: calc(50% - 400px);;
            }
           "
  ))),
  tabPanel("ExtractR",
  sidebarLayout(
    sidebarPanel(
      h3("Step 1: Query UniProtKB"),
      ## choices here are c(human, mouse, yeast, plant)
      helpText("Query options: (1) choose/input Taxon ID below for whole proteome, (2) list custom UniProtKB entries, or (3) upload a .TSV query downloaded from UniProtKB website. Common Taxon IDs: H. sapiens (9606), M. musculus (10090), S. cerevisiae (559292), A. thaliana (3702)"),
      hr(),
      h6("Option 1: whole proteome"),
      ## the placeholder and onInitialize arguments in options make the selectizeInput have no default value upon running the app
      selectizeInput("taxon_input", "Organism taxon ID", choices = c(9606, 10090,559292 ,3702), options=list(create=TRUE, placeholder = 'Type or select taxon ID', onInitialize = I('function() { this.setValue(""); }'))),
      selectInput("review_input","UniProtKB Entry type", c("reviewed", "not reviewed", "both reviewed and not reviewed")),
      h6("Option 2: custom UniProtKB entries"),
      textAreaInput("individual_input", "Input custom UniProtKB entries (e.g. Q16864) separated by line" ),
      helpText("If query is successful, then a \"Run ExtractR\" button will appear below Step 3."),
      actionButton("query_API", label="Query UniProtKB"),
      hr(),
      h6("Option 3: upload a .TSV query downloaded from UniProtKB website"),
      helpText("It is recommended to export from UniProtKB as a .TSV file and upload directly into this app without opening the file in another application, such as Excel (see manual for details)."),
      fileInput("file1", "Choose TSV File", accept=c(".tsv", ".txt")),
      h3("Step 2 (optional): Include organelle mapping file"),
      helpText("Optional: upload your mapping file for Subcellular location. Important: values within the Exact and StartsWith columns MUST be separated by a single comma, not a comma and space. Exact values will be prioritized over StartsWith matches if there are multiple mappings."),
      fileInput("file2", "Optional: Choose TSV File", accept=c(".tsv", ".txt")),
      hr(),
      h3("Step 3: Run ExtractR"),
      helpText("A \"Run ExtractR\" button will appear below this text upon successful UniProtKB query or upload from Step 1."),
      helpText("After running, a table and download button will appear. The table is a frequency table that shows distinct values for each column, including empty values. This may be useful especially for mapping subcellular locations."),
      hr(),
      conditionalPanel("output.files_ready",
                       actionButton("run_extract",label="Run ExtractR")),
      hr()),
    mainPanel(
      uiOutput("helpme_text"),
      uiOutput("downloader"),
      uiOutput("helpme_text2"),
      uiOutput("category_select"),
      uiOutput("category_table")
    )
  )
),
tabPanel("User Guide",
         uiOutput("instructions", style="width: 75%; margin: auto;"))
)
)

server <- function(input,output,session){
  ## 5 GB max upload file
  options(shiny.usecairo=T, shiny.maxRequestSize=5000*1024^2)
  API.df <- eventReactive(input$query_API, {
    withProgress(message='Querying UniProtKB',
                 detail='beep beep boop beep...',value=0,{
    taxon.id <- input$taxon_input
    ## wholeProteome variable is the query for queryUniProt() argument
    wholeProteome <- NULL
    wholeProteome <- if(input$taxon_input!="") {
      ## create wholeProteome if taxon_input exists
      if (input$review_input=="reviewed") {
        c(paste0("organism_id:", taxon.id), paste0("reviewed:true"))
      } else {
        if (input$review_input=="not reviewed") {
          c(paste0("organism_id:", taxon.id), paste0("reviewed:false"))
          ## and then the rest, both reviewed and not reviewed, so leave reviewed blank
        } else {
          c(paste0("organism_id:", taxon.id))
        }
      }
    } 
    else {
      if(input$individual_input!="") {
        paste0(
          ## paste accession in front to ensure it's accession/entry 
          "accession:",
          unlist(strsplit(input$individual_input, split="\n"))
        )
      }
    }
    incProgress(1/3)
    ## now craft the actual UniProtKB API query
    ## if taxon_input exists, use that for the query for wholeProteome. This is first in the ifelse statement, so will be default.
    query.df <-  if (!is.null(wholeProteome)) {
      queryUniProt(
      ## query either accession #s in quotes, like "A0JNW5", etc., or whole coverage, like "organism_id:9606", "reviewed:true". This is basically a search bar
      query = wholeProteome,
      ## which fields to return?
      fields = c("accession", "reviewed", "id", "gene_names", "organism_name", "length", "ft_dna_bind", "cc_pathway", "ft_transmem", "ft_signal", "protein_families", "ft_domain", "ft_motif", "cc_disease", "cc_subcellular_location"),
      pageSize=500L,
      ## collapse tells you the operator for query components. default is " OR "
      collapse = if(input$taxon_input!=""){
      " AND "
      } else {
        if(input$individual_input!="") {
          " OR "
        }
      }
    )
    }
    query.df
    incProgress(1)
    ## This asks if Taxon ID is valid. If there are no rows in query.df, then it will return NULL and the "Extract" button will not appear. 
    if (nrow(query.df) > 0) {
    query.df
    } else {
      NULL
    }
    ##close brackets for withProgress
      })
    ##close brackets for eventReactive
  })
  output$files_ready <- reactive({
    ## shows up when either API.df() is not null (when querying API worked), or when user inputs a file.
    return(!is.null(input$file1) | !is.null(API.df()))
  })
  outputOptions(output,"files_ready", suspendWhenHidden=FALSE)
  extracted.df <- NULL
  extracted.df <- eventReactive(input$run_extract, {
    withProgress(message='Extracting',
                 detail='May take a few minutes...',value=0,{
                   if (is.null(input$file1) & is.null(API.df()))
                     return(NULL)
                   # my.uniprot.df <- read.csv(input$file1$datapath, header=TRUE, na.strings=c("", "NA", "NaN", "NULL", "null"))
                   ## the dataframe to be used will the API.df() queried dataframe, unless it's null. If it's null, then it will be input$file1.
                   my.uniprot.df <- if(!is.null(API.df())) {
                      API.df() 
                     } else { read.table(input$file1$datapath, sep = '\t', header = TRUE, na.strings=c("", "NA", "NaN", "NULL", "null"), fill = TRUE, quote = "")
                     }
                   if (!is.null(input$file2)) {
                     # map.up<- read.csv(input$file2$datapath, header=TRUE, na.strings=c("", "NA", "NaN", "NULL", "null")) 
                     map.up <- read.table(input$file2$datapath, sep = '\t', header=TRUE, na.strings=c("", "NA", "NaN", "NULL", "null"), fill = TRUE)
                   } else {
                     map.up <- NULL
                   }
                   incProgress(1/3)
                   ## start of extracting steps
                   #### UniProt details; Read data ####
                   ## breaking down into multiple sections.
                   up <- my.uniprot.df
                   ## user should read dataframe in with na.strings = "", but just in case
                   up[up == ""] <- NA
                   ## also just in case user has modified in other program
                   up[up == "NA"] <- NA
                   up[up == "NULL"] <- NA
                   up[up == "null"] <- NA
                   up[up == "NaN"] <- NA
                   ## create a vector to show which columns have been modified. This covers 9 UniProt categories (more than 9 columns are added at end if all 9 categories are present). Change this number if you add more modified categories. This is for running UniProtExtractR locally.
                   summary.changes <- c(rep(NA,9))
                   #### Disease ####
                   if (c("Involvement.in.disease") %in% colnames(up) & !((sum(is.na(up$Involvement.in.disease))) == nrow(up))) {
                     # up$disease_count <- str_count(up$Involvement.in.disease, "DISEASE: ")
                     ## example case with 16 DISEASES and 9 THOUSAND characters!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                     ## let's count the number of diseases associated with each protein. Some of these entries start with "DISEASE: Note=" which does not list a disease. There is no obvious order I have found to which is listed first "DISEASE: Note=" vs. "DISEASE: [[disease name]]" so we will have to work that out below.               
                     ## if str_count of "DISEASE: Note=" matches "DISEASE: " then you know it's NA
                     # up$disease_count1 <- str_count(up$Involvement.in.disease, "DISEASE: ")
                     up <- add_column(up, disease_count1 = str_count(up$Involvement.in.disease, "DISEASE: "), .after="Involvement.in.disease")
                     # up$disease_count2 <- str_count(up$Involvement.in.disease, "DISEASE: Note=")
                     up <- add_column(up, disease_count2 = str_count(up$Involvement.in.disease, "DISEASE: Note="), .after="disease_count1")
                     # up$Involvement.in.disease <- ifelse(up$disease_count1 == up$disease_count2, NA, up$Involvement.in.disease )
                     up <- add_column(up, Involvement.in.disease.edit = ifelse(up$disease_count1 == up$disease_count2, NA, up$Involvement.in.disease ), .after = "Involvement.in.disease")
                     # sum(!is.na(up$Involvement.in.disease))
                     ## have removed 387 values where all disease listings are Notes. 
                     ## sometimes, all values might be NA at this point. If had all "DISEASE: Note=". Only run the rest of the code below if that is not the case.
                     if(sum(is.na(up$Involvement.in.disease.edit))!=nrow(up)) {
                     # sum(startsWith(up$Involvement.in.disease.edit, "DISEASE: Note="), na.rm=T)
                     ## how many start with "DISEASE: Note="?
                     ## while there are some entries that start with "DISEASE: Note=", remove "DISEASE: Note= ... DISEASE: " and replace with "DISEASE: " until all Disease entries do not start with "DISEASE: Note= ". Basically bringing the disease name to the front of the string and eliminating the Note= repeatedly until so. 
                     while (sum(startsWith(up$Involvement.in.disease.edit, "DISEASE: Note="), na.rm=T) > 0)
                       up$Involvement.in.disease.edit <- gsub("^DISEASE: Note=(.*?) DISEASE: ", "DISEASE: ", up$Involvement.in.disease.edit)
                     ## now this should be 0
                     # sum(startsWith(up$Involvement.in.disease.edit, "DISEASE: Note="), na.rm=T)
                     ## this should be same not NA values as before
                     # sum(!is.na(up$Involvement.in.disease.edit))
                     
                     ## then run this final command to remove any additional disease notes. (Notes are not disease names). This will not be at the beginning of the string due to repeated function above, which means first entry of Involvement.in.disease will be a disease name (technically, it will be the first disease name listed in the original unedited string).
                     ## before run command, should have n values of "DISEASE: Note = " but 0 startsWith (like above)
                     # sum(str_count(up$Involvement.in.disease.edit, "DISEASE: Note="), na.rm=T)
                     up$Involvement.in.disease.edit <- gsub("DISEASE: Note=(.*?)", "", up$Involvement.in.disease.edit)
                     
                     ## now this below should be zero, should have no matches
                     # sum(str_count(up$Involvement.in.disease.edit, "DISEASE: Note="), na.rm=T)
                     
                     ## woo!
                     ## now recount actual instances of disease names
                     up <- add_column(up, Involvement.in.disease_count = str_count(up$Involvement.in.disease.edit, "DISEASE: "), .after="disease_count2")
                     # up$disease_count <- str_count(up$Involvement.in.disease, "DISEASE: ")
                     # hist(up$disease_count)
                     ## just for fun visualization
                     
                     ## now can extract the first disease for strings that have more than 1 disease entry for further cleaning below.
                     up$Involvement.in.disease.edit <- ifelse(up$Involvement.in.disease_count > 1, str_extract(up$Involvement.in.disease.edit, "DISEASE:(.*?) DISEASE: "), up$Involvement.in.disease.edit)
                     ## is there " [MIM:" in each string?
                     # sum(str_count(up$Involvement.in.disease.edit, " \\[MIM"), na.rm=T)
                     # sum(!is.na(up$Involvement.in.disease.edit))
                     ## THere is one entry with 2 " [MIM" in it. Hmm
                     # up$temp <- str_count(up$Involvement.in.disease, " \\[MIM")
                     ## It's just another annotation. We're good to go. Extract just the first disease from all entries. 
                     up$Involvement.in.disease.edit <- str_extract(up$Involvement.in.disease.edit, "DISEASE:(.*?)\\[MIM")
                     ## check that each string ends in [MIM
                     # sum(grepl("\\[MIM$", up$Involvement.in.disease.edit))
                     ## take off the first "DISEASE: " 9 characters 
                     up$Involvement.in.disease.edit <- str_sub(up$Involvement.in.disease.edit, 10, end = nchar(up$Involvement.in.disease.edit))
                     ## and take off last 5 characters " [MIM"
                     up$Involvement.in.disease.edit <- substr(up$Involvement.in.disease.edit, 1, nchar(up$Involvement.in.disease.edit)-5)
                     ## Does every disease have parentheses at the end?
                     # sum(str_count(up$Involvement.in.disease.edit, "\\((.*?)\\)$"), na.rm=T)
                     # sum(!is.na(up$Involvement.in.disease.edit))
                     ## yes, these two are the same
                     ## get rid of anything after a comma and then floating integers at end.
                     ## This is comma command. Space is crucial! Otherwise things like "46,XX" get removed.
                     up$Involvement.in.disease.edit <- gsub("\\, (.*?)$", "", up$Involvement.in.disease.edit)
                     ## get rid of parentheses so can get at floating integers
                     up$Involvement.in.disease.edit <- gsub(" \\((.*?)\\)$", "", up$Involvement.in.disease.edit)
                     ## floating integers and alphas at ends of strings, many cases
                     ## case 1: digits (x, xx, and xxx)
                     up$Involvement.in.disease.edit <- gsub(" [[:digit:]]$| [[:digit:]][[:digit:]]$| [[:digit:]][[:digit:]][[:digit:]]$", "", up$Involvement.in.disease.edit)
                     ## case 2: usual suspects
                     up$Involvement.in.disease.edit <- gsub(" [[:digit:]][[:alpha:]]$| [[:digit:]][[:digit:]][[:alpha:]]$| [[:digit:]][[:alpha:]][[:digit:]]$", "", up$Involvement.in.disease.edit)
                     ## case 3: believe it or not!
                     up$Involvement.in.disease.edit <- gsub(" [[:digit:]][[:alpha:]][[:digit:]][[:alpha:]]$| [[:digit:]][[:alpha:]][[:alpha:]]$", "", up$Involvement.in.disease.edit)
                     ## finally, remove the disease_count1 and 2 columns, have served their purpose
                     ## also change any final disease_count == 0 values to NA since there was only a Note= originally in the disease category.
                     up$Involvement.in.disease_count[up$Involvement.in.disease_count== 0] <- NA
                     }
                     up$disease_count1 <- NULL
                     up$disease_count2 <- NULL
                     
                     ## BAM!
                     # View(table(up$Involvement.in.disease))
                     summary.changes[1] <- "Involvement.in.disease"
                   }
                   
                   #### Binary variables: DNA.binding, Transmembrane, Signal.peptide ####
                   if (c("DNA.binding") %in% colnames(up) & !((sum(is.na(up$DNA.binding))) == nrow(up))) {
                     # View(table(up$DNA.binding))
                     # length(table(up$DNA.binding))
                     ## 518 terms
                     ## 612 total protein entries with DNA binding terms...change to just DNA binding
                     # sum(!is.na(up$DNA.binding))
                     # up$DNA.binding_binary <- ifelse(!is.na(up$DNA.binding), "DNA binding", NA)
                     up <- add_column(up, DNA.binding_binary = ifelse(!is.na(up$DNA.binding), "DNA binding", NA), .after="DNA.binding")
                     ## count number of DNA binding regions
                     # up$DNA.binding_count <- str_count(up$DNA.binding, "DNA_BIND")
                     up <- add_column(up, DNA.binding_count = str_count(up$DNA.binding, "DNA_BIND"), .after = "DNA.binding_binary")
                     summary.changes[2] <- "DNA.binding"
                   }
                   if (c("Transmembrane") %in% colnames(up) & !((sum(is.na(up$Transmembrane))) == nrow(up))) {
                     # View(table(up$Transmembrane))
                     # sum(table(up$Transmembrane))
                     ## 5208 proteins that are transmembrane
                     # up$Transmembrane_binary <- ifelse(!is.na(up$Transmembrane), "Transmembrane", NA)
                     up <- add_column(up, Transmembrane_binary = ifelse(!is.na(up$Transmembrane), "Transmembrane", NA), .after="Transmembrane")
                     ## count number of transmembrane domains
                     # up$Transmembrane_count <- str_count(up$Transmembrane, "TRANSMEM ")
                     up <- add_column(up, Transmembrane_count = str_count(up$Transmembrane, "TRANSMEM "), .after="Transmembrane_binary")
                     summary.changes[3] <- "Transmembrane"
                   }
                   if (c("Signal.peptide") %in% colnames(up) & !((sum(is.na(up$Signal.peptide))) == nrow(up))) {
                     # View(table(up$Signal.peptide))
                     # length(table(up$Signal.peptide))
                     # sum(table(up$Signal.peptide))
                     # up$Signal.peptide_binary <- ifelse(!is.na(up$Signal.peptide), "Signal peptide", NA)
                     up <- add_column(up, Signal.peptide_binary = ifelse(!is.na(up$Signal.peptide), "Signal peptide", NA), .after="Signal.peptide")
                     # sum(!is.na(up$Signal.peptide))
                     summary.changes[4] <- "Signal.peptide"
                   }
                   #### Protein.families ####
                   if (c("Protein.families") %in% colnames(up) & !((sum(is.na(up$Protein.families))) == nrow(up))) {
                     ## This step just removes anything after the comma in the entries.
                     # up$Protein.families.short <- str_split_fixed(up$Protein.families, ", ", n=2)[,1]
                     up <- add_column(up, Protein.families.edit = str_split_fixed(up$Protein.families, ", |; ", n=2)[,1], .after="Protein.families")
                     # length(table(up$Protein.families.short))
                     # sum(!is.na(up$Protein.families.short))
                     summary.changes[5] <- "Protein.families"
                   }
                   ## Following commented out code is to make a graph of x = potential arbitrary cutoff (# proteins/term) and y = ratio of term coverage; are these terms super specific?
                   # prot.fam.n <- as.data.frame(table(up$Protein.families.short))
                   # ## make a loop to make a plot to see where is good arbitrary cutoff
                   # max.per <- max(prot.fam.n$Freq)
                   # min.per <- min(prot.fam.n$Freq)
                   # ## table of a table?? lol. It tells you the numbers of times a protein family term is counted, then a vector of those times for the loop purpose
                   # percent.vec <- names(table(prot.fam.n$Freq))
                   # percent.vec <- as.numeric(percent.vec)
                   # percents <- vector(length = length(percent.vec))
                   # 
                   # for (i in 1:length(percent.vec)) {
                   #   x <- (subset(prot.fam.n, prot.fam.n$Freq >= percent.vec[i]))
                   #   ## when i = 1, the sum of these frequencies should be the same as sum(!is.na(up$Protein.families))
                   #   ## now fraction of total proteins that are annotated, what percent is your arbitrary cutoff covering?
                   #   percents[i] <- sum(x$Freq)/sum(!is.na(up$Protein.families))
                   # }
                   # prot.fam.n.df <- data.frame(potential_cutoff = percent.vec, ratio_represent_prot_fam = percents)
                   # 
                   # ggplot(data = prot.fam.n.df, aes(x = potential_cutoff, y = ratio_represent_prot_fam)) +
                   #   geom_point() +
                   #   xlim(0, 50)
                   # ## Thought this would have more coverage, but nope. This is showing 50% of coverage around 8 terms which is
                   # nrow(subset(prot.fam.n, prot.fam.n$Freq >= 8))
                   # ## 294 terms :(
                   #### Pathway ####
                   if (c("Pathway") %in% colnames(up) & !((sum(is.na(up$Pathway))) == nrow(up))) {
                     ## This is very similar to Subceullular location
                     ## do all entries start with "PATHWAY: "?
                     # sum(!is.na(up$Pathway))
                     # sum(is.na(up$Pathway))
                     # nrow(up[grep("PATHWAY: ", up$Pathway), ])
                     ## all do contain "PATHWAY: "
                     ## are all instances counted above at the beginning of the strings? Because we know entries have it starting in multiple positions within single string.
                     # sum(startsWith(up$Pathway, prefix="PATHWAY: "), na.rm=T)
                     ## also 1344, nice
                     ## yes all start with "PATHWAY: "
                     ## take off the first 9 (start at character 10) characters which corresponds to "PATHWAY: "
                     # up$Pathway.edit <- str_sub(up$Pathway, 10, end = nchar(up$Pathway))
                     up <- add_column(up, Pathway.edit = str_sub(up$Pathway, 10, end = nchar(up$Pathway)), .after="Pathway")
                     ## no regex for this function below
                     # sum(startsWith(up$Pathway.edit, prefix="[Isoform"), na.rm=T)
                     ## 4 start with "[Isoform x]"
                     ## eliminate "Isoform " and all text after. "[]" are included in punct. This will eliminate proteins that only have subcellular location information for specific isoforms. 
                     ## HOW MANY ELIMINATED? count NAs first
                     # sum(is.na(up$Pathway.edit))
                     #19078
                     up$Pathway.edit <- gsub("\\[Isoform .*", "", up$Pathway.edit)
                     # sum(is.na(up$Pathway.edit))
                     ## still 19078
                     ## now replace any empty values in Subcell Location column with NA; should be 3422 + 412
                     up$Pathway.edit[up$Pathway.edit == ""] <- NA
                     # sum(is.na(up$Pathway.edit))
                     ## 19082, nice! Successfully dropped isoform specific entries.
                     #  sum(!is.na(up$Pathway.edit))
                     ## totaling 1340 entries
                     ## add a "}" at the beginning for pattern matching (see loop later)
                     stri_sub(up$Pathway.edit, 1, 0) <- "}"
                     ## this only puts it in front of non-NA strings, just fyi
                     # sum(grepl("^\\}", up$Pathway.edit))
                     ## should match
                     # sum(!is.na(up$Pathway.edit))
                     ## it does
                     up$Pathway.edit <- paste0(up$Pathway.edit, "{")
                     ## so I guess just make them NA again
                     up$Pathway.edit[up$Pathway.edit == "NA{"] <- NA
                     # sum(is.na(up$Pathway.edit))
                     ## sanity check it's still 19082, same as before
                     
                     ## transient.df can take on as many Pathways annotated per protein as you want. Commented out code below shows an example for first 3 entries, code that runs is for just the first entry. Kept the loop in case of >1 subcellular locations are desired. 
                     ## create a dataframe to fill each row with Protein X's first 3 unlisted Pathways
                     # transient.df <- data.frame(a = rep(0, nrow(df)), b = rep(0, nrow(df)), c = rep(0, nrow(df)))
                     ## create a dataframe to fill each row with Protein X's first unlisted Pathway
                     transient.df <- data.frame(Pathway.edit = rep(0, nrow(up)))
                     ## loop through the rows
                     for (i in 1:nrow(up)) {
                       ## first actually extract the text you want
                       ## this step...just know a lot of thinking happened here
                       ## regex brain power
                       test <- str_match_all(up[i,c("Pathway.edit")], "\\}(.*?)\\{")[[1]][,2]
                       ## get rid of everything after the semicolon for the first term of test
                       test2 <- gsub(";.*", "", test[1])
                       ## remove random punctuation at end of strings
                       test3 <- gsub("\\. $|\\.$", "", test2)
                       ## replace these terms with "regulation" for increasing coverage per term at term specificity expense; can comment this out and change return value to test3 if want to keep specificity.
                       test4 <- gsub("biosynthesis|degradation|metabolism", "regulation", test3)
                       ## take only the first elements for each. 
                       transient.df[i,] <- test4
                     }
                     ## okay will it blend?
                     # sum(!is.na(up$Pathway.edit)) - sum(!is.na(transient.df$Pathway.edit))
                     #  sum(!is.na(up$Pathway.edit))
                     ## BAM! It blends!!!!!!! 0 values different :)
                     ## 1340 annotations, wow
                     # View(table(transient.df$Pathway.edit))
                     up$Pathway.edit <- transient.df$Pathway.edit
                     summary.changes[6] <- "Pathway"
                   }
                   #### Motif ####
                   if (c("Motif") %in% colnames(up) & !((sum(is.na(up$Motif))) == nrow(up))) {
                     # View(table(up$Motif))
                     ## how many Proteins have >1 motif? str_count
                     # hist(str_count(up$Motif, "note="), breaks=20, xlim=c(0,10), ylim=c(0,200))
                     ## ~150 proteins have 2 motifs of 2369 total entries, 6%. ~50 proteins have 3, 25 have 4...negligible after 5 really, though looks like max is 20? wow
                     transient.df <- data.frame(a = rep(0, nrow(up)))
                     ## don't need a loop if just taking first 
                     ## Apparently, UniProt .TSV export files have quotes around Motifs and Domains. If you query the API directly, there are no quotes. This if statement should take care of that.
                     transient.df$Motif.edit <- str_match(up$Motif, "note=\"(.*?)\"")[,2]
                     if(sum(is.na(transient.df$Motif.edit))== nrow(transient.df)) {
                       transient.df$Motif.edit <- str_match(up$Motif, "note=(.*?)$")[,2]
                     }
                     ## eliminate anything after a semicolon
                     transient.df$Motif.edit <- gsub(";.*", "", transient.df$Motif.edit)
                     ## str_trim eliminates spaces at start and ends of strings
                     transient.df$Motif.edit <- str_trim(transient.df$Motif.edit, side = c("both"))
                     ## There are so many instances of "(blah) 1" or "1 (blah)" or "(blah)," at the ends of these strings
                     ## case 1. 
                     # grep("\\,$", transient.df$Motif.edit)
                     transient.df$Motif.edit <- gsub("\\,$", "", transient.df$Motif.edit)
                     ## case 2. 
                     # grep(" \\((.*?)\\) [[:digit:]]$", transient.df$Motif.edit)
                     transient.df$Motif.edit <- gsub(" \\((.*?)\\) [[:digit:]]$", "", transient.df$Motif.edit)
                     ## case 3. 
                     # grep(" \\((.*?)\\)$", transient.df$Motif.edit)
                     transient.df$Motif.edit <- gsub(" \\((.*?)\\)$", "", transient.df$Motif.edit)
                     ## case 4. Floating integers. Srsly.
                     # grep(" [[:digit:]]$", transient.df$Motif.edit)
                     transient.df$Motif.edit <- gsub(" [[:digit:]]$", "", transient.df$Motif.edit)
                     ## sanity check interlude: spaces at end?
                     # grep(" $", transient.df$Motif.edit)
                     ## excellent news, it's zero.
                     ## case 5. Get rid of "motif" at end. Many instances.
                     # grep(" motif$", transient.df$Motif.edit)
                     transient.df$Motif.edit <- gsub(" motif$", "", transient.df$Motif.edit)
                     ## quality control to make sure you changed what you thought
                     # sum(!is.na(up$Motif))
                     # sum(!is.na(transient.df$Motif.edit))
                     # length(table(transient.df$Motif.edit))
                     # sum(table(transient.df$Motif.edit))
                     # View(table(transient.df$Motif.edit))
                     # up$Motif.edit <- transient.df$Motif.edit
                     up <- add_column(up, Motif.edit = transient.df$Motif.edit, .after="Motif")
                     summary.changes[7] <- "Motif"
                   }
                   #### Domain..FT.####
                   if (c("Domain..FT.") %in% colnames(up) & !((sum(is.na(up$Domain..FT.))) == nrow(up))){
                     # sum(!is.na(up$Domain..FT.))
                     # View(table(up$Domain..FT.))
                     ## Ah, our old friend. Looks just like Motif.
                     # hist(str_count(up$Domain..FT., "note="), breaks=2000, xlim=c(0, 20))
                     ## out of 8643 annotated proteins, ~2000 have 2 annotated domains. 2000/8643 ~ 23%. Hmm. We'll just take the first domain.  
                     # summary(str_count(up$Domain..FT., "note="))
                     ## holy smokes, a protein has 285 annotated domains??? lol
                     transient.df <- data.frame(a = rep(0, nrow(up)))
                     ## don't need a loop if just taking first  
                     ## Apparently, UniProt .TSV export files have quotes around Motifs and Domains. If you query the API directly, there are no quotes. This if statement should take care of that.
                     transient.df$Domain..FT.edit <- str_match(up$Domain..FT., "note=\"(.*?)\"")[,2]
                     if(sum(is.na(transient.df$Domain..FT.edit))== nrow(transient.df)) {
                       transient.df$Domain..FT.edit <- str_match(up$Motif, "note=(.*?);")[,2]
                     }
                     ## floating integers strike again!!!!
                     # View(table(transient.df$Domain..FT.edit))
                     ## case 1. Floating integers. Srsly.
                     # grep(" [[:digit:]]$", transient.df$Domain..FT.edit)
                     transient.df$Domain..FT.edit <- gsub(" [[:digit:]]$", "", transient.df$Domain..FT.edit)
                     ## case 2. Everything that comes after semicolon.
                     # transient.df[grep("; ", transient.df$Domain..FT.edit), c("Domain..FT.edit")]
                     transient.df$Domain..FT.edit <- gsub(";.*", "", transient.df$Domain..FT.edit)
                     ## nice. 692 domains.
                     # up$Domain..FT.edit <- transient.df$Domain..FT.edit
                     up <- add_column(up, Domain..FT.edit = transient.df$Domain..FT.edit, .after="Domain..FT.")
                     summary.changes[8] <- "Domain..FT."
                   }
                   #### Through the Fire and Flames on Expert: Subcellular Location ####
                   if (c("Subcellular.location..CC.") %in% colnames(up) & !((sum(is.na(up$Subcellular.location..CC.))) == nrow(up))) {
                     ## do all entries start with "SUBCELLULAR LOCATION: "? yes
                     # sum(!is.na(up$Subcellular.location..CC.))
                     # sum(is.na(up$Subcellular.location..CC.))
                     # nrow(up[grep("SUBCELLULAR LOCATION: ", up$Subcellular.location..CC.), ])
                     ## are all "SUBCELLULAR LOCATION: " instances counted above at the beginning of the strings? Because we know entries have it starting in multiple positions within single string.
                     # sum(startsWith(up$Subcellular.location..CC., prefix="SUBCELLULAR LOCATION: "), na.rm=T)
                     ## also 16982, nice
                     ## yes all start with "SUBCELLULAR LOCATION: "
                     ## take off the first 22 characters which corresponds to "SUBCELLULAR LOCATION: "
                     # up$Subcellular.location..CC. <- str_sub(up$Subcellular.location..CC., 23, end = nchar(up$Subcellular.location..CC.))
                     up <- add_column(up, Subcellular.location..CC.edit = str_sub(up$Subcellular.location..CC., 23, end = nchar(up$Subcellular.location..CC.)), .after="Subcellular.location..CC.")
                     ## no regex for this function below
                     # sum(startsWith(up$Subcellular.location..CC.edit, prefix="[Isoform"), na.rm=T)
                     ##412 start with "[Isoform x]"
                     ## eliminate "Isoform " and all text after. "[]" are included in punct. This will eliminate proteins that only have subcellular location information for specific isoforms. 
                     ## HOW MANY ELIMINATED?
                     ## count NAs first
                     # sum(is.na(up$Subcellular.location..CC.edit))
                     #3422
                     up$Subcellular.location..CC.edit <- gsub("\\[Isoform .*", "", up$Subcellular.location..CC.edit)
                     # sum(is.na(up$Subcellular.location..CC.edit))
                     ## still 3422
                     ## now replace any empty values in Subcell Location column with NA; should be 3422 + 412
                     up$Subcellular.location..CC.edit[up$Subcellular.location..CC.edit == ""] <- NA
                     # sum(is.na(up$Subcellular.location..CC.edit))
                     ## 3834, nice! Successfully dropped isoform specific entries.
                     # sum(!is.na(up$Subcellular.location..CC.edit))
                     ## totaling 16570 entries
                     ## add a "}" at the beginning for pattern matching (see lapply later)
                     stri_sub(up$Subcellular.location..CC.edit, 1, 0) <- "}"
                     ## this only puts it in front of non-NA strings, just fyi
                     # sum(grepl("^\\}", up$Subcellular.location..CC.edit))
                     ## should match
                     # sum(!is.na(up$Subcellular.location..CC.edit))
                     ## it does
                     #### Subcellular Location: Isoforms interlude; this does not run ####
                     ## Bingo, 412 proteins that start with "Isoform". Subcellular location is Isoform specific.
                     ## before you just eliminate anything after the first "Note=", what is max isoform number?
                     ## omg wow there are tons! make histogram
                     ## holy smokes, there is one protein with 15 isoforms
                     # iso.vec <- vector(length = 15)
                     # for (i in 1:15) {
                     #   # name.var <- paste("Isoform", i, sep=" ")
                     #   ## default of paste() is sep=" "!!!!!!!! 
                     #   iso.vec[i] <- nrow(up[grep(paste("Isoform", i, sep= " "), up$Subcellular.location..CC.edit),])
                     # }
                     ## another way to vectorize grep from our favorite website
                     ## https://stackoverflow.com/questions/20935910/use-the-grep-function-in-a-loop-in-r
                     # vGrep <- Vectorize(grep, vectorize.args = "pattern", SIMPLIFY = FALSE)
                     # vGrep(iso.vec, df$Subcellular.location..CC.edit)
                     ## data for ggplot has to be data.frame, so convert the vector
                     # iso.df <- data.frame(Isoform.count = 1:15, freq = iso.vec)
                     # ggplot(data = iso.df) +
                     #   geom_col(aes(x = Isoform.count, y = freq)) +
                     #   theme_bw()
                     ##bingo
                     ## wow this is beautiful S curve
                     ## okay but why would nrow(Isoform 1) =/= nrow(Isoform 2)?
                     ## One example is B0YJ81; HACD1_HUMAN. Has just Isoform 1 listed. 
                     ## okay how many are there that are like this? 
                     # iso.df[1, 2] - iso.df[2, 2]
                     ## 18 total proteins...drop for now.
                     # nrow(df[grep("Isoform 15", df$Subcellular.location..CC.edit),])
                     #### Subcellular Location: Back to our main programming ####
                     ## isoforms are gone, but some proteins start with "}[Some characters] Secreted" etc. how many?
                     # sum(startsWith(up$Subcellular.location..CC.edit, prefix="}["), na.rm=T)
                     ## 99!!!!! another edge case, excellent
                     ## how many matches for this regex? "[]" (not StartsWith)
                     # sum(grepl("\\[(.*?)\\]", up$Subcellular.location..CC.edit))
                     ## 192
                     # sum(grepl("^\\{\\[(.*?)\\]", up$Subcellular.location..CC.edit))
                     ## count how many strings start with "}[" or "}[words strings blah]"
                     # sum(grepl("^\\}\\[", up$Subcellular.location..CC.edit))
                     ## 99
                     # sum(grepl("^\\}\\[(.*?)\\]", up$Subcellular.location..CC.edit))
                     ## 99 
                     ## same as above, but not at beginning of string?
                     # sum(grepl("\\}\\[(.*?)\\]", up$Subcellular.location..CC.edit))
                     ## 99, we're good
                     ## "}[Ribosomal protein X]: Nucleus" is what we are dealing with
                     up$Subcellular.location..CC.edit <- str_replace(up[,c("Subcellular.location..CC.edit")], "\\}\\[(.*?)\\]: ", "}")
                     ## let's make sure only startsWith are altered
                     ## Nice! Triple check: 
                     # sum(grepl("\\}\\[(.*?)\\]: ", up$Subcellular.location..CC.edit))
                     ## 0 
                     # sum(grepl("\\}\\[", up$Subcellular.location..CC.edit))
                     ## 0
                     ## also NAs should be the same
                     # sum(!is.na(up$Subcellular.location..CC.edit))
                     ## nice
                     ## this test case works
                     # str_replace("}[try me]: Secreted {", "\\}\\[(.*?)\\]: ", "}")
                     ## eliminate "Note=" and all text after. Notes exist at end of string, unless protein has isoforms listed. In that case it's isoform1-location-Note; isoform2-location-Note; for n isoforms. If eliminate Notes, then assume all locations are assigned to Isoform1 in UniProt.
                     ## do any entries start with "Note="?
                     # sum(grepl("^}Note=", up$Subcellular.location..CC.edit))
                     # up$Subcellular.location..CC.edit[grep("^}Note=", up$Subcellular.location..CC.edit)]
                     ## 4 entries. Okay this is where those random 4 empty entries were coming from. 
                     up$Subcellular.location..CC.edit <- gsub("^}Note=", NA, up$Subcellular.location..CC.edit)
                     ## check these 4 values are gone now
                     # sum(grepl("^}Note=", up$Subcellular.location..CC.edit))
                     ## check number of NA
                     # sum(is.na(up$Subcellular.location..CC.edit))
                     ## added 4, nice :)
                     ## Now sub the rest of Notes with blanks
                     up$Subcellular.location..CC.edit <- gsub("Note=.*", "", up$Subcellular.location..CC.edit)
                     ## are there entries that only contain a Note?
                     ## Any of these strings would actually be "}" and not "" because of na.strings = "" in the beginning. These would have been captured already.
                     # up$Subcellular.location..CC.edit[up$Subcellular.location..CC.edit == ""] <- NA
                     # sum(is.na(up$Subcellular.location..CC.edit))
                     ## no, still 3834 NA values here
                     ## Add "{" at the end of each string to create the "}{" sandwich
                     up$Subcellular.location..CC.edit <- paste0(up$Subcellular.location..CC.edit, "{")
                     ## make the now "NA{" values NA again
                     up$Subcellular.location..CC.edit[up$Subcellular.location..CC.edit == "NA{"] <- NA
                     # sum(is.na(up$Subcellular.location..CC.edit))
                     ## sanity check it's still 3834, same as before
                     ## create a dataframe to fill each row with Protein X's n unlisted subcellular localization
                     ## this dataframe does only the first value (1 column), which does not required a loop (can change to str_match instead of str_match_all), but keeping in case user wants it.
                     # transient.up <- data.frame(a = rep(0, nrow(up)))
                     ## this dataframe commented out below does the first 3 values. There are proteins with 15+ locations. 
                     # transient.up <- data.frame(a = rep(0, nrow(up)), b = rep(0, nrow(up)), c = rep(0, nrow(up)))
                     ## loop through the rows
                     # for (i in 1:nrow(up)){
                     #   ## first actually extract the text you want
                     #   ## this step...just know a lot of thinking happened here
                     #   ## regex brain power
                     #   test <- str_match_all(up[i,c("Subcellular.location..CC.edit")], "\\}(.*?)\\{")[[1]][,2]
                     #   ## get rid of punctuation. This includes [:punct:]: punctuation characters, ! " # $ % & ’ ( ) * + , - . / : ; < = > ? @ [  ] ^ _ ` { | } ~.
                     #   test2 <- gsub('[[:punct:]]','',test)
                     #   ## squish gets rid of leading or trailing spaces
                     #   test3 <- str_squish(test2)
                     #   ## take only the first element for each. Some have n localizations. You can change this to n localizations if you want.
                     #   transient.up[i,] <- test3[1]
                     #   # string.vec[[i]] <- str_squish(test3)
                     # }
                     
                     sub.fx.1 <- function(x) {
                       ## first actually extract the text you want
                       ## this step...just know a lot of thinking happened here
                       ## regex brain power
                       if(!is.na(x)) {
                         test <- str_match_all(x, "\\}(.*?)\\{")[[1]][,2]
                         ## get rid of punctuation. This includes [:punct:]: punctuation characters, ! " # $ % & ’ ( ) * + , - . / : ; < = > ? @ [  ] ^ _ ` { | } ~.
                         test2 <- gsub('[[:punct:]]','',test)
                         ## There is one edge case I found where Subcellular location is listed twice before the "}{" sandwich. This should take care of that. Rare case.
                         test2.5 <- gsub("SUBCELLULAR LOCATION.*", "", test2)
                         ## squish gets rid of leading or trailing spaces
                         test3 <- str_squish(test2.5)
                         # test3 <- str_squish(test2)
                         ## take only the first element for each. Some have n localizations. You can change this to n localizations if you want.
                         
                         return(test3[1])
                       } else {
                         return(NA)
                       }
                     }
                     ## apply this function to the scratch work column
                     transient.up <- data.frame(a = unlist(lapply(up[,c("Subcellular.location..CC.edit")], sub.fx.1)))
                     up <- add_column(up, Subcellular.location_short = transient.up$a, .after="Subcellular.location..CC.edit")
                     # sum(is.na(transient.up$a))
                     ## okay will it blend?
                     # sum(!is.na(up$Subcellular.location..CC.edit)) - sum(!is.na(transient.up$a))
                     # sum(!is.na(up$Subcellular.location..CC.edit))
                     ## BAM! It blends!!!!!!! 0 values different :)
                     ## 16570 annotations, wow
                     #### Subcellular Location: Rounding up specificity to organelle with a custom Map ####
                     ## Reducing down the specificity.
                     ## first pass; what are we dealing with?
                     # length(table(transient.up$a))
                     ## 865 annotations...... lol
                     # View(as.data.frame(table(transient.up$a)))
                     ## sorting this table by Var1 is useful for seeing which UniProt term is more prominent for each organelle 
                     ## if a map.up dataframe is provided, then carry out the mapping function. If not, then extracted Subcellular.location terms will remain the first extracted term from above.
                     if (!is.null(map.up)) {
                       map.up <- map.up
                       ## make these regex statements for Exact matches
                       map.up$Exact <- gsub(",","$|^",map.up$Exact)
                       stri_sub(map.up$Exact, 1, 0) <- "^"
                       stri_sub(map.up$Exact, (nchar(map.up$Exact)+1), (nchar(map.up$Exact)+2)) <- "$"
                       ## make these StartsWith regex statements
                       map.up$StartsWith <- gsub(",", ".*|^",map.up$StartsWith)
                       stri_sub(map.up$StartsWith, 1, 0) <- "^"
                       stri_sub(map.up$StartsWith, (nchar(map.up$StartsWith)+1), (nchar(map.up$StartsWith)+2)) <- ".*"
                       ## this will help with later regex statements
                       map.up[is.na(map.up)] <- c("ImpossibleREGEX0000")
                       ## now replace all the Exact matches
                       up$Subcellular.location..CC.map.edit1 <- stri_replace_all_regex(up[,c("Subcellular.location_short")], map.up[,c("Exact")], map.up[,c("New")], vectorize_all = FALSE)
                       ## now replace all the StartsWith matches
                       up$Subcellular.location..CC.map.edit2 <- stri_replace_all_regex(up[,c("Subcellular.location..CC.map.edit1")], map.up[,c("StartsWith")], map.up[,c("New")], vectorize_all = FALSE)
                       ## Want to make sure Exact matches are prioritized before being potentially wiped out by StartsWith matches. Assuming the value would change between original (_short) and Exact changes (edit1), then keep the changed version, else keep the StartsWith (edit2)
                       # up$Subcellular.location..CC.map.edit3 <- ifelse(up$Subcellular.location_short==up$Subcellular.location..CC.map.edit1, up$Subcellular.location..CC.map.edit2, up$Subcellular.location..CC.map.edit1)
                       up <- add_column(up, Subcellular.location..CC.map.edit3 = ifelse(up$Subcellular.location_short==up$Subcellular.location..CC.map.edit1, up$Subcellular.location..CC.map.edit2, up$Subcellular.location..CC.map.edit1), .after="Subcellular.location_short")
                       ## after Subcellular.location..CC.map.edit3 is created, drop the other two columns
                       up[,c("Subcellular.location..CC.map.edit1", "Subcellular.location..CC.map.edit2")] <- NULL
                       ## and rename the existing column
                       colnames(up)[colnames(up)=="Subcellular.location..CC.map.edit3"] <- "Subcellular.location..CC.map.edit"
                       ## end of if (!is.null(map.up)) section
                     }
                     ## remove the edited column that is now scratchwork
                     up$Subcellular.location..CC.edit <- NULL
                     ## change column back to maintain naming consistency
                     colnames(up)[colnames(up)=="Subcellular.location_short"] <- "Subcellular.location..CC.edit" 
                     ## BAM, amazing!
                     summary.changes[9] <- "Subcellular.location..CC."
                   }
                   #### Write to .csv file and summarize changes (local version) ####
                   # ## This section is for running local version, will not be run in Shiny version.
                   # filenameme <- paste0("UniProt_edited_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv", sep = "")
                   # ## writing to working directory
                   # write.csv(up, filenameme)
                   # ## cat prints the strings!
                   # ## This prints the final message a user sees:
                   # cat(paste0(
                   #   "Edits completed. Additional edited columns have been added for the following categories:", "\n",
                   #   paste0(summary.changes[!is.na(summary.changes)], "\n", collapse=""), "\n",
                   #   "Thanks for flying with us.", "\n"
                   # )
                   # )
                   # 
                   #### Return a dataframe, up, for the Shiny ####
                   ## This section is for running Shiny version, will not be run in local version.
                   up
                   ## end of extracting steps
                 }
    )
    ##end of eventReactive({}) for Shiny
  })
  output$helpme_text <-  renderUI({
    if(!is.null(extracted.df())) {
        HTML(paste("Your UniProtKB query was successfully extracted. Your original query and added \"extracted\" data columns are available to download using the \"Download\" button below.",  "", sep= "<br/>"))
    }
      })
  
  output$downloader <- renderUI({
    if(!is.null(extracted.df())) 
      downloadButton('OutputFile', "Download extracted query file")
  })
  output$OutputFile <- downloadHandler(
    filename=function() {
      paste0("UniProt_modified","_",Sys.time(),".csv")
    },
    content=function(file){
      write.csv(extracted.df(),file)
    }
  )
  output$helpme_text2 <-  renderUI({
    if(!is.null(extracted.df())) {
      HTML(paste("This frequency table below is to browse a summary of each category of the UniProtKB query.", "", sep= "<br/>"))
    }
  })
  output$category_select <- renderUI({
    if(!is.null(extracted.df()))
      selectizeInput("category","Frequency table category:",
                     choices=colnames(extracted.df()),
                     options=list(
                       placeholder="Type column name",
                       ## may the force be with you if you have 20000 columns
                       maxOptions = 20000),
                     selected=NULL,
                     multiple=FALSE)
  })
  
  output$category_table <- renderUI({
    if(!is.null(extracted.df()))
      output$actual_table <- DT::renderDataTable(as.data.frame(table(extracted.df()[,input$category], useNA="always")) %>% setNames(c("Unique value", "Frequency count")), server = FALSE, selection = "single")
    DT::dataTableOutput("actual_table")
  })
  output$instructions <- renderUI({
    HTML( ##user guide summary
      "
      <h1> UniProtExtractR User Guide Summary </h1> 
      <h5> Shiny App for extracting desired information from syntactically structured character strings across 9 UniProtKB categories: DNA binding, Pathway, Transmembrane, Signal peptide, Protein families, Domain [FT], Motif, Involvement in disease, and Subcellular location [CC]. A full manual and example files may be found at the <a href='https://github.com/alex-bio/UniProtExtractR'>UniProtExtractR GitHub page.</a></h3>
      <h3> Overview </h3>
      <p>UniProtKB protein entries can be challenging to handle for data analysis, visualization, and categorization due to character string syntax. This app extracts desired information from syntactically complex UniProtKB entries for the 9 UniProtKB specific categories above. Further, a user can supply an optional mapping file for Subcellular location [CC], as many annotated locations are quite specific.</p>
      <h3> Step 1: Query UniProtKB </h3>
      <p>The only requisite for the app is a UniProtKB query (size limit is 5GB). UniProtKB can be directly queried, which is recommended. Else, a user can query UniProtKB on the <a href='https://uniprot.org'>UniProtKB website</a> and download a .TSV file (after searching, click 'Download' and select .TSV as file format). This should be a tab separated (.TSV or .TXT). For only the categories to be extracted by the R package or app, those column names must exactly match those exported from UniProtKB, but do not have to be in any particular order.<br>
      When querying UniProtKB directly, a loading bar will appear until the query is finished loading into the app. Upon successful completion, a 'Run ExtractR' button will appear at the bottom of the left-hand panel.</p>
      <h3> Step 2 (optional): Include organelle mapping file </h3>
      <p> The mapping file is optional and should be tab separated if included (.TSV or .TXT). The mapping file is only for the category Subcellular location [CC]. It should contain only 3 columns: 'New', 'Exact', and 'StartsWith' (see example at the <a href='https://github.com/alex-bio/UniProtExtractR'>UniProtExtractR GitHub page</a>). For every Subcellular location [CC] per entry in the UniProtKB query file, any terms that match exactly to 'Exact' or that start with the phrases in 'StartsWith' will be swapped for the value in 'New' in a row-wise fashion. Each row represents paired 'Exact' and 'StartsWith' matches for a single replacement value 'New'. The 'Exact' and 'StartsWith' values may have many listed phrases, and 'Exact' takes priority. Importantly, the values within 'Exact' and 'StartsWith' columns must be separated by a single comma without a space.</p>
      <h3> Step 3: Run the app </h3>
      <p>If UniProtKB was successfully queried or a file was uploaded in Step 1, then a 'RunExtractR' button will appear at the bottom of the left-hand panel. It may take some time for the app to process the UniProtKB query, depending on size.</p>
      <h3> Outputs </h3>
      <p>If the app successfully ran, a message will appear stating success at the top of the right-hand side. Additionally, the app will return two things: (1) a .CSV file containing the originally uploaded UniProtKB file with additional modified or 'extracted' columns, and (2) an in-app browsable frequency table of the extracted UniProtKB query. The .CSV file containing the uploaded query and 'extracted' columns is available to download by using the download button at the top of the main panel. The frequncy table that appears on the screen summarizes the counts of all unique values within each column of the UniProtKB query and 'extracted' columns. This may be useful for summarizing the UniProtKB query: as one example, the number of entries in the query with 1, 2, 3, ..., etc. transmembrane domains. The frequency table also shows a blank value for each category to demonstrate category coverage, i.e. the number of empty values within that category. If the following categories are present (and not empty), then new columns (bulleted) will appear immediately adjacent:
      <ol>
        <li>DNA Binding</li>
          <ul>
              <li>DNA.binding_binary</li>
              <li>DNA.binding_count</li>
          </ul>
        <li>Pathway</li>
          <ul>
              <li>Pathway.edit</li>
          </ul>
        <li>Transmembrane</li>
        <ul>
              <li>Transmembrane_binary</li>
              <li>Transmembrane_count</li>
          </ul>
        <li>Signal peptide</li>
        <ul>
              <li>Signal.peptide_binary</li>
          </ul>
        <li>Protein families</li>
        <ul>
              <li>Protein.families.edit</li>
          </ul>
        <li>Domain [FT]</li>
        <ul>
              <li>Domain..FT.edit</li>
          </ul>
        <li>Motif</li>
        <ul>
              <li>Motif.edit</li>
          </ul>
        <li>Involvement in disease</li>
        <ul>
              <li>Involvement.in.disease.edit</li>
              <li>Involvement.in.disease_count</li>
          </ul>
        <li>Subcellular location [CC]</li>
        <ul>
              <li>Subcellular.location..CC.edit</li>
              <li>Subcellular.location..CC.map.edit (if mapping file present)</li>
          </ul>
      
      
      </ol>
      
      <p> For troubleshooting, further instructions, or example files, please refer to the full manual at the <a href='https://github.com/alex-bio/UniProtExtractR'>UniProtExtractR GitHub page</a>.
      <br>
      <br>
      <br>
      <br>
      <p style='color:gray;'> Page last updated 2023-10-18.</p>
      "
    )
  })
}

shinyApp(ui,server)

####End of Shiny App.
