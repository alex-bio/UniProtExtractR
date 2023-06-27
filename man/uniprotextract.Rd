\name{uniprotextract}
\alias{uniprotextract}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
uniprotextract
}
\description{
Extracts desired information from syntactically structured character strings across 9 UniProtKB categories: DNA binding, Pathway, Transmembrane, Signal peptide, Protein families, Domain [FT], Motif, Involvement in disease, and Subcellular location [CC]. A mapping file for Subcellular location [CC] can be provided to further reduce specificity.
}
\usage{
uniprotextract(my.uniprot.df, map.up = NULL, write.local = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{my.uniprot.df}{
A dataframe exported from UniProt
}
  \item{map.up}{
A dataframe with 3 columns: "New", "Exact", and "StartsWith"
}
  \item{write.local}{
TRUE or FALSE; default is FALSE and function will return the modified dataframe; if TRUE, then it writes to the local directory
}
}
\details{
The dataframe \code{my.uniprot.df} should be exported as directly from UniProtKB as possible. It can contain as many entries as desired; whole proteomes will take longer than a few dozen or hundred entries.

For \code{Subcellular.location..CC.}, entries can be further mapped at user-determined resolution with \code{map.up}; for example, if the user is interested in mitochondria specific protein locations, the mapping file can assign broad, organelle-level terms to every organelle except mitochondria (e.g. "Cytoplasm", "Nucleus"). Mitochondria-specific terms can be mapped to a finer resolution (e.g. "Inner membrane", "Matrix").

The option \code{write.local} is \code{FALSE} by default. If \code{TRUE}, then the function will return a message indicating completion and write a .CSV file of the modified dataframe to the working directory.

Warning: opening exported files from UniProtKB in Excel can yield unexpected results.

Please see the manual for further detail.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## upload your dataframe from UniProt as directly as possible,
## unlike the handmade dataframe in the lines below

example.df <- data.frame(Entry = "ENTRY1",
  Gene.Names = "ProteinX",
  Transmembrane = 'TRANSMEM 3..23; /note="Helical; Name=1"; /evidence="ev"; TRANSMEM 302..322; /note="Helical; Name=2"; /evidence="ev"',
  Subcellular.location..CC. = "SUBCELLULAR LOCATION: Membrane {ev}; Multi-pass membrane protein {ev}")

uniprotextract(example.df)
## will extract "Membrane" for Subcellular.location..CC. category, and will count number of Transmembrane domains

## if provide mapping file, will assign "Membrane" to whichever category you want
mapping.df <- data.frame(New = "Cell membrane",
  Exact = "Membrane",
  StartsWith = NA)

uniprotextract(example.df, mapping.df)
## Now "Membrane" has been assigned to user-defined "Cell membrane"


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
