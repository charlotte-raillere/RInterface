\name{transmart.getMiRNAData}
\alias{transmart.getMiRNAData}
%- Also NEED an "\alias" for EACH other topic documented here.
\title{
This function gets miRNA data from the tranSMART database.
}
\description{
This function will retrieve miRNA data based on filters. Options are available to pivot and aggregate the data.
}
\usage{
transmart.getMiRNAData(mirna.type = "",
					study.list = NA, 
					mirna.list = NA, 
					patient.list = NA, 
					sample.list = NA,
					sample.type.list = NA, 
					tissue.type.list = NA, 
					timepoint.list = NA, 
					platform.list = NULL, 
					probe.list = NULL, 
					platform.removeOnOverlap = NULL, 
					show.mirna = FALSE, 
					print.statement = FALSE, 
					data.pivot = TRUE, 
					data.pivot.aggregate = NULL,
					data.pivot.patient_id = FALSE,
					data.pivot.sample = FALSE)
}
%- maybe also "usage" for other objects documented here.
\arguments{
  \item{mirna.type}{
	The type of the miRNA. This can be one of "qPCR miRNA" or "miRNAseq" (case insensitive)
}
  \item{study.list}{
	A list of studies to limit the text search to. This lookup is case insensitive.
}
  \item{mirna.list}{
	A list of miRNA to filter the miRNA by.
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
  \item{probe.list}{
	A list of probe identifiers to filter the results by.
}
  \item{platform.removeOnOverlap}{
	This is a list of platforms that you want to be *overlapped* if there are probes in multiple platforms. The script will check for probes in multiple platforms and remove the records that have a platform in the provided list.
}
  \item{show.mirna}{
	If this is set to true a mirRNA column will be added to data output.
}
  \item{print.statement}{
	If this is set to true the function will only print the SQL statement to retrieve the miRNA data instead of running it.
}
  \item{data.pivot}{
	Flag indicating whether the miRNA data should be pivoted or not.
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
	#Getting miRNAseq data by Study.
	miRNAData <- transmart.getMiRNAData(study.list=c("TEST_MIRNA_SEQ"),
										mirna.type = "miRNAseq")
	#---

	#---
	#Getting qPCR miRNA data by Study.
	miRNAData <- transmart.getMiRNAData(study.list=c('TEST_MIRNA_QPCR'),
										mirna.type = "qPCR miRNA")
	#---
	
	#---
	#Getting qPCR miRNA data by Study and probe list. Get non pivoted data.
	miRNAData <- transmart.getMiRNAData(study.list=c('TEST_MIRNA_QPCR'),
										data.pivot=FALSE, show.mirna=TRUE, probe.list=c('3','6','9'),
										mirna.type = "qPCR miRNA")
	#---
	
	#---
	#Getting qPCR miRNA data by Study and miRNA list
	miRNAData <- transmart.getMiRNAData(study.list=c('TEST_MIRNA_QPCR'),
										mirna.list=c('mmu-let-7d','mmu-let-7i','mmu-mir-10b'),
									    mirna.type = "qPCR miRNA")
	#---	

	#---
	#Getting qPCR miRNA data by Study, sample code and miRNA list. Get non pivoted data.
	miRNAData <- transmart.getMiRNAData(study.list=c('TEST_MIRNA_QPCR'),
										sample.list=c('GSM001', 'GSM002'), data.pivot=FALSE,
										show.mirna=TRUE, mirna.list=c('mmu-let-7d','mmu-let-7i','mmu-mir-10b'),
										mirna.type = "qPCR miRNA")
	#---	
}
% Add one or more standard keywords, see file "KEYWORDS" in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
