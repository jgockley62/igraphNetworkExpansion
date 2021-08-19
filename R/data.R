#' Sample  network
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

#' Basic Network Object 
#'
#' The basic Pathway Commons network. This is all SIF Pathway commons 
#' interactions loaded into a network. The CHEMI interactions are present. The
#' igraph network object consisting of 30,910 genes and 1,902,605 
#'  interactions. basic_network$graph is a list of networks. Each 
#' network is every interaction in Pathway Commons consisting of a specific type
#' of interaction.
#'
#' @format A list object
#' \describe{
#'   \item{igraph::V(basic_network$network)}{Vertex names}
#'   \item{igraph::E(basic_network$network)}{Network edges}
#'   \item{igraph::E(basic_network$network)$interaction}{Type of interaction}
#'   \item{igraph::E(basic_network$network)$Occurance}{How many Pathway Commons 
#'   databases record an interaction between to and from vertex}
#'   \item{igraph::E(basic_network$network)$UniqCol}{The Unique Code of the 
#'   entire interaction}
#'   \item{igraph::E(basic_network$network)$pathway}{The Pathway(s) in Pathway 
#'   commons this interaction is documented in as a list object}
#'   \item{igraph::E(basic_network$network)$EdgeRep}{How many interaction types 
#'   are there between to and from vertacies}
#'   \item{igraph::E(basic_network$network)$Edge}{The to:from name desigantion 
#'   of the edge}
#'   \item{igraph::E(basic_network$network)$SumOccurancel}{how many times is to
#'    vertex connected to from vertex across all pathway commons data bases and 
#'    interaction types}
#'   ...
#' }
"basic_network"

#' Basic Network Table 
#'
#' The basic Pathway Commons table that comprises the foundation of the network
#' This is all SIF Pathway commons interactions loaded into a network. The 
#' CHEMI interactions are present. The igraph network object consisting of 
#' 30,910 genes and 1,902,605 interactions. 
#'
#' @format A data table
#' \describe{
#'   \item{basic_network_table$to}{Origin Vertex}
#'   \item{basic_network_table$from}{Target Vertex}
#'   \item{basic_network_table$interaction}{Type of interaction}
#'   \item{basic_network_table$Occurance}{how many Pathway Commons databases record an
#'   interaction between to and from vertex}
#'   \item{basic_network_table$UniqCol}{The Unique Code of the entire interaction}
#'   \item{basic_network_table$pathway}{The Pathway(s) in Pathway commons this 
#'   interaction is documented in as a list object}
#'   \item{basic_network_table$EdgeRep}{How many interaction types are there between 
#'   to and from vertacies }
#'   \item{basic_network_table$Edge}{The to:from name desigantion of the edge}
#'   \item{basic_network_table$SumOccurancel}{how many times is to vertex connected to
#'    from vertex across all pathway commons data bases and interaction types}
#'   ...
#' }
"basic_network_table"

#' Basic Network Interactions
#'
#' This Resources is a list of networks. Each network is every interaction in 
#' Pathway Commons consisting of a specific type of interaction.
#'
#' @format A list object
#' \describe{
#'   \item{basic_network_interactions}{A list of igraph network objects by interaction 
#'   type, all edge and vertex attributes are the same as basic_network$network}
#'   ...
#' }
"basic_network_interactions"
