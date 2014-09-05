
# ********************************************************************************
#   Copyright 2012   Recombinant Data Corp.
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
#
# ********************************************************************************

#-----------------------
#Retrieve GEX data based on some passed in parameters.
#-----------------------
transmart.getMiRNAData <- function(
  mirna.type      ='',
  study.list      =NA,
  mirna.list				=NA,
  patient.list			=NA,
  sample.list      =NA,
  sample.type.list 		=NA,
  tissue.type.list 		=NA,
  timepoint.list 		= NA,
  platform.list			= NULL,
  probe.list			= NULL,
  platform.removeOnOverlap = NULL,
  show.mirna			= FALSE,
  print.statement 		= FALSE,
  data.pivot			= TRUE,
  data.pivot.aggregate	= NULL,
  data.pivot.patient_id = FALSE,
  data.pivot.sample = FALSE
)
{
  #Create the connection to the oracle DB.
  tranSMART.DB.connection <- tranSMART.DB.establishConnection()
  
  if(any(is.na(study.list)) && any(is.na(patient.list)) && any(is.na(sample.list))) stop("You must provide either a study list or a patient list to run a query.")
  
  if(toupper(mirna.type) != toupper('miRNAseq') && toupper(mirna.type) != toupper('qPCR miRNA')) stop("You must provide the miRNA type: mirna.type must be miRNAseq or qPCR miRNA")
  
  #This is the base SQL statement to get microarray data, Column section only.
  baseSQLSelectStatement <- gsub("\n", "", "
                                 a.PATIENT_ID, 
                                 ssm.SAMPLE_CD,
                                 a.RAW_INTENSITY, 
                                 a.ZSCORE, 
                                 a.LOG_INTENSITY, 
                                 a.assay_id,
                                 ssm.subject_id,
                                 ssm.sample_type,
                                 ssm.timepoint,
                                 ssm.tissue_type,
                                 ssm.trial_name,
                                 ssm.GPL_ID, 
                                 b.id_ref")
  
  #This is the base of the table clauses.
  baseSQLTableStatement <- gsub("\n", "", "
                                FROM deapp.de_subject_mirna_data a 
                                INNER JOIN deapp.de_subject_sample_mapping ssm ON ssm.assay_id = A.assay_id")
  
  baseSQLWhereStatement <- ""
  
  whereStatementBegan <- FALSE
  
  #We either filter by a patient list or a study list.
  if(any(!is.na(study.list)))
  {
    baseSQLWhereStatement <- gsub("\n", "", " 
                                  WHERE SSM.trial_name IN (?1) ")
    
    studyListString <- paste("UPPER('",study.list,"')",sep="",collapse=",")
    
    baseSQLWhereStatement <- gsub("\\?1",studyListString,baseSQLWhereStatement)
    
    whereStatementBegan <- TRUE
  }
  if(any(!is.na(patient.list))){
    if(whereStatementBegan){
      baseSQLWhereStatement <- paste(baseSQLWhereStatement,"  AND SSM.patient_id IN (?) ",sep="")
    }else{
      baseSQLWhereStatement <- " WHERE SSM.patient_id IN (?) "  
    }
    
    patientListString <- paste("'",patient.list,"'",sep="",collapse=",")
    
    baseSQLWhereStatement <- gsub("\\?",patientListString,baseSQLWhereStatement)
    
    whereStatementBegan <- TRUE
  }
  if(any(!is.na(sample.list))){
    if(whereStatementBegan){
      baseSQLWhereStatement <- paste(baseSQLWhereStatement,"  AND SSM.sample_cd IN (?) ",sep="")
    }else{
      baseSQLWhereStatement <- " WHERE SSM.sample_cd IN (?) "  
    }
    
    patientListString <- paste("'",sample.list,"'",sep="",collapse=",")
    
    baseSQLWhereStatement <- gsub("\\?",patientListString,baseSQLWhereStatement)
  }
  
  if(toupper(mirna.type) == toupper('miRNAseq')){
    baseSQLWhereStatement <- paste(baseSQLWhereStatement," AND ssm.platform ='MIRNA_SEQ' ",sep="")
  }else{
    baseSQLWhereStatement <- paste(baseSQLWhereStatement," AND ssm.platform ='MIRNA_QPCR' ",sep="")
  }
  
  
  if(show.mirna == TRUE)
  {
    #Add the mirna columns to the select.
    baseSQLSelectStatement <- paste(baseSQLSelectStatement,", b.mirna_id ",sep= " ")
  }
  
  #If a pathway was included, append more info to the query here.
  if(any(!is.na(mirna.list)) || any(!is.null(probe.list)))
  {
    #Add a SELECT DISTINCT to the select statement.
    baseSQLSelectStatement <- paste("SELECT DISTINCT ",baseSQLSelectStatement,sep=" ")
    
    #Add the table joins.
    newTableJoins <- gsub("\n", "", "
                          INNER JOIN deapp.de_qpcr_mirna_annotation b ON a.probeset_id = b.probeset_id 
                          INNER JOIN biomart.bio_marker bm ON bm.PRIMARY_EXTERNAL_ID = to_char(b.mirna_id)
                          INNER JOIN biomart.bio_marker_correl_mv sbm ON sbm.asso_bio_marker_id = bm.bio_marker_id
                          INNER JOIN searchapp.search_keyword sk ON sk.bio_data_id = sbm.bio_marker_id
                          INNER JOIN searchapp.SEARCH_KEYWORD_TERM skt ON sk.SEARCH_KEYWORD_ID = skt.SEARCH_KEYWORD_ID ")
    
    baseSQLTableStatement <- paste(baseSQLTableStatement,newTableJoins,sep="")
    
    #Create the string of mirna
    if(any(!is.na(mirna.list)))
    {
      searchWords <- paste("UPPER('",mirna.list,"')",sep="",collapse=",")
      baseSQLWhereStatement <- paste(baseSQLWhereStatement," AND skt.KEYWORD_TERM IN (",searchWords,")",sep="")
    }else if(any(!is.null(probe.list)))
    {
      searchWords <- paste("'",probe.list,"'",sep="",collapse=",")
      baseSQLWhereStatement <- paste(baseSQLWhereStatement," AND b.id_ref IN (",searchWords,")",sep="")
    }
    
  }else
  {
    #Add a SELECT DISTINCT to the select statement.
    baseSQLSelectStatement <- paste("SELECT DISTINCT ",baseSQLSelectStatement,sep=" ")	
    newTableJoins <- " INNER JOIN deapp.de_qpcr_mirna_annotation b ON a.probeset_id = b.probeset_id "
    baseSQLTableStatement <- paste(baseSQLTableStatement,newTableJoins,sep="")
  }
  
  baseSQLWhereStatement <- paste(baseSQLWhereStatement," AND ssm.gpl_id = b.gpl_id",sep="")
  
  #We can add the 3 SSM filters at the end of the query here if they were passed in.
  if(any(!is.na(sample.type.list)))
  {
    sampleList <- paste("'",sample.type.list,"'",sep="",collapse=",")
    
    baseSQLWhereStatement <- paste(baseSQLWhereStatement," AND ssm.sample_type IN (",sampleList,")",sep="")
  }
  
  if(any(!is.na(tissue.type.list)))
  {
    tissueList <- paste("'",tissue.type.list,"'",sep="",collapse=",")
    
    baseSQLWhereStatement <- paste(baseSQLWhereStatement," AND ssm.tissue_type IN (",tissueList,")",sep="")
  }	
  
  if(any(!is.na(timepoint.list)))
  {
    timepointList <- paste("'",timepoint.list,"'",sep="",collapse=",")
    
    baseSQLWhereStatement <- paste(baseSQLWhereStatement," AND ssm.timepoint IN (",timepointList,")",sep="")
  }	
  
  if(any(!is.null(platform.list)))
  {
    platformlist <- paste("'",platform.list,"'",sep="",collapse=",")
    
    baseSQLWhereStatement <- paste(baseSQLWhereStatement," AND ssm.GPL_ID IN (",platformlist,")",sep="")	
  }
  
  #Put together the final query.
  finalQuery <- paste(baseSQLSelectStatement,baseSQLTableStatement,baseSQLWhereStatement,sep=" ")
  
  #Print the SQL statement if flag was set.
  if(print.statement == TRUE)
  {
    print("Statement to be run : ")
    print(finalQuery)
    return()
  }
  
  print("Sending miRNA Query.")
  
  #Send the query to the server.
  rs <- dbSendQuery(tranSMART.DB.connection, finalQuery)
  
  print("Retrieving miRNA Records.")
  
  #Retrieve the entire data set.
  dataToReturn <- fetch(rs)
  
  #Clear the results object.
  dbClearResult(rs)
  
  #Disconnect from the database.
  dbDisconnect(tranSMART.DB.connection)	
  
  #Warn the user if two platforms were returned.
  if(length(unique(dataToReturn$GPL_ID)) > 1)
  {
    print("Warning! Multiple platforms being returned.")
  }
  
  #If the user wants the data pivoted, do so here.
  if(data.pivot == TRUE)
  {
    library(reshape2)
    
    if(data.pivot.patient_id == TRUE)
    {
      pivotPatient = "PATIENT_ID"
    }else if(data.pivot.sample == TRUE)
    {
      pivotPatient = "SAMPLE_CD"
    }else
    {
      pivotPatient = "SUBJECT_ID"
    }
    
    print("Trimming columns.")
    #These are the column we pull from the GEX dataframe.
    relevantColumns <- c(pivotPatient,'LOG_INTENSITY','ID_REF','GPL_ID')
    
    #This is the formula we use to pivot the GEX data.
    pivotFormula <- paste("GPL_ID + ID_REF ~ ",pivotPatient)
    
    #If the user wants a miRNA column add the column to the list of columns to pull and the UNIPROT_ID to the pivot.
    if(show.mirna == TRUE)
    {
      relevantColumns <- c(relevantColumns,'UNIPROT_ID')
      pivotFormula <- paste("MIRNA_ID + GPL_ID + ID_REF ~ ",pivotPatient)
    }
    
    #We need to Pivot the data so we have probe id, chip, probe id, then each patient's GEX data.
    #Cut down the results to the relevant columns.
    dataToReturn <- dataToReturn[relevantColumns] 
    
    print("Pivoting data.")
    
    if(!is.null(data.pivot.aggregate))
    {
      print("WARNING! About to use an aggregation when pivoting data. Do not use this option unless you are aware of the reason for duplicate information.")
      dataToReturn <- dcast(dataToReturn, as.formula(pivotFormula), value.var = 'LOG_INTENSITY',fun.aggregate = data.pivot.aggregate) 			
    }
    else
    {
      #Pivot the trimmed data.
      dataToReturn <- dcast(dataToReturn, as.formula(pivotFormula), value.var = 'LOG_INTENSITY') 			
    }
  }
  
  #If we want to remove platforms on probe overlap, do it here.
  if(any(!is.null(platform.removeOnOverlap)))
  {
    #For any duplicate probes on different platforms we need to always use one value over another.
    #Find the probes that are in both platforms.
    probeData <- data.frame(table(dataToReturn$ID_REF))
    
    colnames(probeData) <- c('ID_REF','Freq')
    
    #Delete the records from this table with only 1 record.
    probeData <- probeData[which(probeData$Freq > 1),]
    
    #In our final data we remove any records that have a probe in the probeData table, and are the GPL we want to remove.
    dataToReturn <- dataToReturn[which(!(dataToReturn$ID_REF %in% probeData$ID_REF & dataToReturn$GPL_ID %in% platform.removeOnOverlap)),]
    
  }
  
  dataToReturn
}

