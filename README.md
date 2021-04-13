
<!-- badges: start -->
  [![R-CMD-check](https://github.com/jgockley62/igraphNetworkExpansion/workflows/R-CMD-check/badge.svg)](https://github.com/jgockley62/igraphNetworkExpansion/actions)
[![R-CMD-check](https://github.com/jgockley62/igraphNetworkExpansion/workflows/pkgdown/badge.svg)](https://github.com/jgockley62/igraphNetworkExpansion/actions)
<!-- badges: end -->

# Pathway Tracing Package

## Compute Environment and Sandbox Infrastructure
Establishing linked docker containers with an RStudio Sandbox environment and
a linux VIM hosted cytoscape application is supported in separate repository. 
This can be leveraged on and AWS EC2 instance using docker-compose.
`https://github.com/jgockley62/networktracing`

## Installation
`remotes::install_github("jgockley62/igraphNetworkExpansion")`

## Building a template tetwork

## Tracing a template network
The set of functions currently imports 2 gene lists, a targets list and a 
sentinel list. Tracing paths within the template network occurs pairwise within 
the targets list and pairwise between each target and each sentinel. Tracing is 
currently not preformed pairwise from sentinel to sentinel. 

* Tracing: 
The igraph tracing method currently employed is `igraph::get.all.shortest.paths()`. 

* Cutoff Limit: 
The path cutoff limit is currently defined from the whichever is greater between
 the median and mean of traces resulting from the target gene to the sentinel
genes. If there are no trace paths in this object the target gene to other 
target genes trace object is used. 



