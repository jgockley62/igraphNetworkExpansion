% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{rm_dups}
\alias{rm_dups}
\title{Remove Duplicate Gene Entries
The purpose. of this function is to pull info from a synapse table into a
dataframe. Removal can be of all entries with a lower score in a given column
or entries which lack the desired ensg}
\usage{
rm_dups(df, feature_col, feature, type, type_spec, ensg_keep = NULL)
}
\arguments{
\item{df}{the data frame of gene names}

\item{feature_col}{the column of the feature to treat as rows ie.'GName'}

\item{feature}{the values of the feature to pull ie. NPC1}

\item{type}{highest value or ensg ie 'value' or 'ensg'}

\item{type_spec}{column name of type ie 'ENSG' or 'Overall'}

\item{ensg_keep}{the ensg to retain if method is ensg default = NULL}
}
\value{
a dataframe object
}
\description{
Remove Duplicate Gene Entries
The purpose. of this function is to pull info from a synapse table into a
dataframe. Removal can be of all entries with a lower score in a given column
or entries which lack the desired ensg
}
\examples{
\dontrun{
 syn <- log_into_synapse()
 omics_scores <- dplyr::left_join(
  table_pull(
    syn_id ='syn25575156',
    feature_name = 'GeneName',
    features = names(igraph::V(net)),
    column_names = c('ENSG', 'OmicsScore', 'GeneticsScore', 'Logsdon'),
    synap_import = syn$synapse
  ), 
  table_pull(
    syn_id ='syn22758536',
    feature_name = 'GName',
    features = names(igraph::V(net)),
    column_names = c('ENSG', 'GName', 'RNA_TE', 'Pro_TE'),
    synap_import = syn$synapse
  ), 
  by = 'ENSG'
 )
 colnames(omics_scores)[ colnames(omics_scores) == 'Logsdon' ] <- 'Overall'
 omics_scores <- rm_dups(
   df = omics_scores,
   feature_col = 'GName',
   feature = 'POLR2J3',
   type = 'value',
   type_spec = 'Overall' ,
   ensg_keep = NULL
 )
 omics_scores <- rm_dups(
    df = omics_scores,
    feature_col = 'GName',
    feature = "FCGBP",
    type = 'ensg',
    type_spec = 'ENSG' ,
    ensg_keep = 'ENSG00000281123'
 )
}
}
