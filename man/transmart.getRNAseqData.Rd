\name{transmart.getRNAseqData}
\alias{transmart.getRNAseqData}
%- Also NEED an "\alias" for EACH other topic documented here.
\title{
This function gets RNAseq data from the tranSMART database.
}
\description{
This function will retrieve RNAseq data based on filters. Options are available to pivot and aggregate the data.
}
\usage{
transmart.getRNAseqData(study.list = NA, 
					gene.list = NA, 
					pathway = NA, 
					signature = NA, 
					patient.list = NA, 
					sample.list = NA,
					sample.type.list = NA, 
					tissue.type.list = NA, 
					timepoint.list = NA, 
					platform.list = NULL, 
					transcript.list = NULL, 
					platform.removeOnOverlap = NULL, 
					show.genes = FALSE, 
					print.statement = FALSE, 
					data.pivot = TRUE, 
					data.pivot.aggregate = NULL,
					data.pivot.patient_id = FALSE,
					data.pivot.sample = FALSE)
}
%- maybe also "usage" for other objects documented here.
\arguments{
  \item{study.list}{
	A list of studies to limit the text search to. This lookup is case insensitive.
}
  \item{gene.list}{
	A list of genes to filter the RNAseq by.
}
  \item{pathway}{
	A pathway to filter the RNAseq by (exact text match, case insensitive).
}
  \item{signature}{
	The name of a gene signature to filter the RNAseq by.
}
  \item{patient.list}{
	A list of patient IDs to filter the data by. This is the patient identifier generated within tranSMART.
}
 \item{sample.list}{
	A list of sample codes to filter the data by. This is the sample code used to load data into the database.
}
  \item{sample.type.list}{
	A list of sample types to filter the data by. The available sample types can be retrieved using the transmart.listHDDAttributes function.
}
  \item{tissue.type.list}{
	A list of tissue types to filter the data by. The available tissue types can be retrieved using the transmart.listHDDAttributes function.
}
  \item{timepoint.list}{
	A list of timepoints to filter the data by. The available timepoints can be retrieved using the transmart.listHDDAttributes function.
}
  \item{platform.list}{
	This is a list of platforms to filter on.
}
  \item{transcript.list}{
	A list of transcript ids to filter the results by.
}
  \item{platform.removeOnOverlap}{
	This is a list of platforms that you want to be *overlapped* if there are transcripts in multiple platforms. The script will check for transcripts in multiple platforms and remove the records that have a platform in the provided list.
}
  \item{show.genes}{
	If this is set to true a gene column will be added to data output.
}
  \item{print.statement}{
	If this is set to true the function will only print the SQL statement to retrieve the RNAseq data instead of running it.
}
  \item{data.pivot}{
	Flag indicating whether the RNAseq data should be pivoted or not.
}
  \item{data.pivot.aggregate}{
	Function to use when aggregating data during a pivot. Use this only if you know what you are aggregating.
}
  \item{data.pivot.patient_id}{
	This will use the tranSMART internal identifier which is guarenteed to be unique per study (barring any data issues).
}
  \item{data.pivot.sample}{
	This will use the sample code instead of the subject identifier during the pivot.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of "comp1"}
%%  \item{comp2 }{Description of "comp2"}
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
	#---
	#Getting RNAseq data by Study and a list of genes. Get the data in
	# a non-pivoted format, and with a column for gene names
	RNAseqData <- transmart.getRNAseqData(study.list = c("TEST_RNA_SEQ"),
								c("BRCA1", "ERBB2"), data.pivot=FALSE, show.gene=TRUE)
	#---
	
	#---
	#Getting RNAseq data by Study, a list of patients and a pathway
	RNAseqData <- transmart.getRNAseqData(study.list = c("TEST_RNA_SEQ"),
								patient.list=c("10001", "10002")), pathway="Apoptosis")
	#---

	#---
	#Getting RNAseq data by Study, a list of patiennts and list of transcripts
	RNAseqData <- transmart.getRNAseqData(study.list = c("TEST_RNA_SEQ"),
								patient.list=c("10001", "10002"))transcript.list=c("ENSG00000006025", "ENSG00000000003"))
	#---

	#---
	#Getting RNAseq data by Study and gene signature. Get the data in
	# a non-pivoted format
	RNAseqData <- transmart.getRNAseqData(study.list = c("TEST_RNA_SEQ"),
								signature=c("Gene signature"), data.pivot=FALSE)
	#---

}
% Add one or more standard keywords, see file "KEYWORDS" in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
