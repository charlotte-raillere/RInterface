\name{transmart.getMetabolomicsData}
\alias{transmart.getMetabolomicsData}
%- Also NEED an "\alias" for EACH other topic documented here.
\title{
This function gets Metabolomics data from the tranSMART database.
}
\description{
This function will retrieve Metabolomics data based on filters. Options are available to pivot and aggregate the data.
}
\usage{
transmart.getMetabolomicsData(study.list = NA, 
					hmdb.list = NA,  
					patient.list = NA, 
					sample.list = NA,
					sample.type.list = NA, 
					tissue.type.list = NA, 
					timepoint.list = NA, 
					platform.list = NULL, 
					biochemical.list = NULL, 
					platform.removeOnOverlap = NULL, 
					show.hmdb = FALSE, 
					print.statement = FALSE, 
					data.pivot = TRUE, 
					data.pivot.aggregate = NULL,
					data.pivot.patient_id = FALSE,
					data.pivot.sample = FALSE,
					superpathway.list = NULL,
					subpathway.list = NULL)
}
%- maybe also "usage" for other objects documented here.
\arguments{
  \item{study.list}{
	A list of studies to limit the text search to. This lookup is case insensitive.
}
  \item{hmdb.list}{
	A list of HMDB identifiers to filter the Metabolomics by.
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
  \item{biochemical.list}{
	A list of antigen names to filter the results by.
}
  \item{platform.removeOnOverlap}{
	This is a list of platforms that you want to be *overlapped* if there are probes in multiple platforms. The script will check for probes in multiple platforms and remove the records that have a platform in the provided list.
}
  \item{show.hmdb}{
	If this is set to true a protein column will be added to data output.
}
  \item{print.statement}{
	If this is set to true the function will only print the SQL statement to retrieve the Metabolomics data instead of running it.
}
  \item{data.pivot}{
	Flag indicating whether the Metabolomics data should be pivoted or not.
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
  \item{superpathway.list}{
	This is a list of superpathway names to filter on.
}
  \item{subpathway.list}{
	This is a list of subpathway names to filter on.
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
	#Getting Metabolomics data by Study and a list of biochemical.
	MetabolomicsData <- transmart.getMetabolomicsData(study.list = c("TEST_METABOLOMICS"),
										biochemical.list=c("glycerate", "citrate"))
	#---

	#---
	#Getting Metabolomics data by Study and a list of HMDB id.
	MetabolomicsData <- transmart.getMetabolomicsData(study.list = c("TEST_METABOLOMICS"),
										hmdb.list=c("HMDB00684", "HMDB00714"))
	#---
	
	#---
	#Getting Metabolomics data by study and a superpathway
	MetabolomicsData <- transmart.getMetabolomicsData(study.list = c("TEST_METABOLOMICS"),
										superpathway=c("Lipid"))
	#---

	#---
	#Getting Metabolomics data by study and a subpathway
	MetabolomicsData <- transmart.getMetabolomicsData(study.list = c("TEST_METABOLOMICS"),
										subpathway=c("Lysolipid"))
	#---
	

	
}
% Add one or more standard keywords, see file "KEYWORDS" in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
