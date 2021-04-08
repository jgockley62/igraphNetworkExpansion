#' Sample simple network
#'
#' A simple example network without weights
#'
#' @format An igraph network object
#' \describe{
#'   \item{verticies}{genes}
#'   \item{edges}{protien interactions}
#'   ...
#' }
"sample_network"

#' Test and example network: implicit loading
#'
#' A Full Template network to use for testing and examples. All protein-protein
#' interactions in pathway commons for genes exressed in at least one tissue 
#' from AMP-AD transcriptomics data.
#'
#' @format An igraph network object
#' \describe{
#'   \item{verticies}{genes}
#'   \item{edges}{protien interactions}
#'   \item{igraph::V(slim_net$name)}{The Vertex Gene Name}
#'   \item{igraph::V(slim_net$community)}{A Sample igraph Community designation 
#'   for the vertex}
#'   \item{igraph::V(slim_net$weight)}{A Sample Weight Value}
#'   \item{igraph::V(slim_net$RNA_EffectScore)}{A Sample RNA Phenotype 
#'   Effect Score}
#'   \item{igraph::V(slim_net$Pro_EffectScore)}{A Sample Protien Phenotype
#'   Effect Score}
#'   \item{igraph::V(slim_net$RNA_Cor_EffectScore)}{A Sample RNA co-expression 
#'   value}
#'   \item{igraph::V(slim_net$Pro_Cor_EffectScore)}{A Sample Protien co-expression 
#'   value}
#'   \item{igraph::E(slim_net$interaction)}{The Edge Interaction Type}
#'   \item{igraph::E(slim_net$UniqCol)}{The From_vert:To_vert-Interaction_type
#'   unique identifier}
#'   \item{igraph::E(slim_net$Edge)}{The From_vertex:To_vertex identifier}
#'   \item{igraph::E(slim_net$pathway)}{The Pathway Commons data base(s) where 
#'   this edge is found}
#'   \item{igraph::E(slim_net$EdgeRep)}{How many interaction types are there 
#'   between the two vertices across pathway commons databases}
#'   \item{igraph::E(slim_net$Occurance)}{How many Data Bases across pathway
#'   commons record an interaction of anytype between the two vertices}
#'   \item{igraph::E(slim_net$SumOccurence)}{Total amount of times these to 
#'   vertices interact accross all pathway commons database and interaction 
#'   types}
#'   ...
#' }
"slim_net"

#' Sample Genelists to Run Network Functions On
#'
#' Simple lists of genes with associated biological function. Target and 
#' sentinal gene lists contain the same 15 lists of genes stored as character
#' vectors in lists named targets and sentinals.
#' 
#' @format An igraph network object
#' \describe{
#'   \item{targets}{target gene lists}
#'   \item{sentinals}{sentinal gene lists}
#'   ...
#' }
"genelists"