
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
#Retrieve Gene signature data
#-----------------------
transmart.getGeneSignature <- function(
  name      =''
)
{
  #Create the connection to the oracle DB.
  tranSMART.DB.connection <- tranSMART.DB.establishConnection()
  
  if(name == '') stop("You must provide the gene signature name")
  
  #This is the base SQL statement to get gene signature
  query <- paste("select case when d.bio_marker_name is null then to_char
(c.gene_symbol) else to_char(d.bio_marker_name) end as GENE_SYMBOL, c.probe_id as PROBE_ID,
b.fold_chg_metric as FOLD_CHANGE_METRIC
from searchapp.search_gene_signature a
inner join searchapp.search_gene_signature_item b on a.search_gene_signature_id = b.search_gene_signature_id
left outer join deapp.de_mrna_annotation c on c.probeset_id = b.probeset_id
left outer join biomart.bio_marker d on d.bio_marker_id = b.bio_marker_id
where deleted_flag=0 and name = '",name, "'", sep='')
  
  #Send the query to the server.
  rs <- dbSendQuery(tranSMART.DB.connection, query)
  
  print("Retrieving Gene signature.")
  
  #Retrieve the entire data set.
  dataToReturn <- fetch(rs)
  
  #Clear the results object.
  dbClearResult(rs)
  
  #Disconnect from the database.
  dbDisconnect(tranSMART.DB.connection)	
  
  dataToReturn
}
