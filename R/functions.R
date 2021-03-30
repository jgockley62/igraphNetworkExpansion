#' Log into Synapse
#'
#' Log into Synapse. Assumes credentials are stored.
#' User can specify UserID and Password as
#' usr and pass respectivly
#' @export
#' @param usr A single value character vector of the users Synapse ID
#' @param pass A single value character vector of the users Synapse Password
#' @return Synapse login object from
#' @examples
#' log_into_synapse()
#' log_into_synapse(
#'    usr = '<UserName>',
#'    pass = '<UserPassword>'
#'  )
#'  
log_into_synapse <- function(usr=NULL, pass=NULL) {
  #install 
  reticulate::conda_create("r-reticulate")
  reticulate::conda_install( channel = 'bioconda', packages="synapseclient")
  synapseclient <- reticulate::import("synapseclient")
  syn_temp <- synapseclient$Synapse()
  if (usr==NULL & pass==NULL){
    syn_temp$login( )
  }else{
    if(usr==NULL | pass==NULL){
      stop("Must Specify Both User ID and Password or 
           Leave Blank to use Synapse Credentials"
          )
    }else{
      eval(parse(text=paste0(
        'syn_temp$login(\'',
        usr,
        '\', \'',
        pass,
        '\')'
      )))
    }
  }
  
}

#' Loads a Gene List From Synapse
#'
#' Loads a gene list. If is_syn is TRUE file_path is interperated as a synapse
#' ID. Otherwise file_path is interperated as a file path. 
#'
#' @export
#' @param file_path the igraph network to push to synapse eg. net
#' @param network igraph network object to use to filter the vertices From
#' @param is_syn is the list a synapse ID default=FALSE
#' @return genes from the input present within the network
#' @examples 
#' all_goterms <- list(
#'    c('syn25185319', "APP_Metabolism", "APP Metabolism"),
#'    c('syn25185320', "Endolysosomal", "Endolysosomal")
#' )
#' list_load(
#'    all_goterms[[1]][1],
#'    network=<USER SPECIFIED igraph Network>,
#'    is_syn = TRUE
#' )
list_load <- function (file_path, network, is_syn = FALSE) {
  if (isTRUE(is_syn)) {
    genes <- read.table(
      file=syn_temp$get(file_path)$path, header=F, sep='\n', stringsAsFactors = F
    )[,1]
  }else{
    genes <- read.table(
      file=file_path,header=F, sep='\n', stringsAsFactors = F
    )[,1]
  }
  #Genes in network
  total <- length(genes)
  perc <- length(genes[ genes %in% names(igraph::V(network)) ]) / length(genes)
  numb <- length(genes[ genes %in% names(igraph::V(network)) ])
  
  genes <- genes[ genes %in% names( igraph::V(network) ) ]
  
  writeLines(paste0( numb, 
                     ' of ', 
                     total, 
                     ' (', signif(perc, digits = 4)*100,
                     '%) genes from the list appear in the primary network' )
  )
  return(genes)
}

#' Pull the names of a PathTrace
#' 
#' Pulls the names from an igraph get.all.shortest.paths() object
#'
#' @export
#' @param char is a list entry from get.all.shortest.paths()
#' @return names of the genes in the input trace path
#' @examples 
#' example_path <- list(
#'    c('319', "49", "23")
#' )
#' names(example_path[[1]]) <- c("GeneA", "GeneZ", "GeneAlpha")
#' name_pull(example_path[[1]])
name_pull <- function( char ){
  return( names( char ) )
}

#' Calculate the OMICS Scores across paths
#' 
#' This functions takes a vector path object from a list of path objects, 
#' usually output from an igraph path trace. This vector is comprised of
#' vertex numbered numeric vector with gene names comprising the vector
#' namespace. These names are matched to a weights vector comprised of gene 
#' weights and a namespace of gene names. The frist and last gene of the trace
#' are droped as they are the start and end of all traces in the trace object.
#' the mean score of the path is returned unless the mean score is equal to or 
#' less than zero or NA.
#' 
#' @export
#' @param indv_path an indv path
#' @param weights the omics scores as a gene-named vector
#' @return mean score of the path trace input without the starting and ending
#' gene
#' @examples 
#' example_path <- list(
#'    c('319', "49", "23", "86", "690", "238", "102")
#' )
#' names(example_path[[1]]) <- c(
#'    "GeneA","GeneZ", "GeneAlpha",
#'    "GeneB", "GeneX", "GeneOmega"
#')
#'
#' sampweights <- c(1.45, 2.45, 0.89, .003, 1.3, 2.1)
#' names(sampweights) <- c(
#'    "GeneA","GeneZ", "GeneAlpha",
#'    "GeneB", "GeneX", "GeneOmega"
#')
#' path_calc(example_path[[1]], sampweights)
path_calc <-  function (indv_path, weights) {
  scores <- weights[names(indv_path)]
  #Trim off the search genes....
  scores <- scores[1:(length(scores) -1)][-1]
  scores <- mean(scores[ !(is.na(scores)) & scores > 0 ])
  return(scores)
}

#' Filters Paths Baised on OMICS Scores
#' 
#' Fiters a path baised on an input limit value and  returns a list object.
#' If the mean path score is above the limit the path will be returned as a 
#' kept gene list and passed score will be 1. Otherwise kept value will be NA 
#' and passed value will be 0. If the mean score is zero path used vaule is
#' equal to zero instead of 1 and the keep and passed vaules resemble a failed
#' path. All genes in a path regaurdless of a path's score are returned in the 
#' genes value.
#' 
#' @export
#' @param weights named vector of genes weights
#' @param indv_path the path of genes
#' @param lim the value to out the path
#' @return A list object reporting all the genes for every path, the genes to
#' keep 
#' @return$Keep If the mean path score is greater than the limit score this is 
#' a vector of the path's gene names otherwise Keep sis NA
#' @return$Genes The names of the  genes in the path
#' @return$Used Was the path used (ie mean path score > 0): Yes = 1, No = 0
#' @return$Passed Was the path's mean score greater than lim: Yes = 1, No = 0
#' @examples 
#' example_path <- list(
#'    c('319', "49", "23", "86", "690", "238", "102")
#' )
#' names(example_path[[1]]) <- c(
#'    "GeneA","GeneZ", "GeneAlpha",
#'    "GeneB", "GeneX", "GeneOmega"
#' )
#'
#' sampweights <- c(1.45, 2.45, 0.89, .003, 1.3, 2.1)
#' names(sampweights) <- c(
#'    "GeneA","GeneZ", "GeneAlpha",
#'    "GeneB", "GeneX", "GeneOmega"
#' )
#' path_filter(example_path[[1]], sampweights, 1)
#' path_filter(example_path[[1]], sampweights, 1.5)
path_filter <-  function (indv_path, weights, lim) {
  scores <- weights[names(indv_path)]
  
  #Filter out origin and terminus
  scores <- scores[1:(length(scores) -1)][-1]
  if(length(scores) == 0){
    return(
      list( Keep = names(indv_path),
            Genes = names(indv_path),
            Used = 1,
            Passed = 1
      )
    )
  }
  scores <- scores[ !(is.na(scores)) ]
  if(length(scores) > 0){
    used <- 1
  }else{ used <- 0 }
  if (mean(scores) > lim){
    return(
      list( Keep = names(scores),
            Genes = names(scores),
            Used = used,
            Passed = 1
      )
    )
  }else{
    return(
      list( Keep = NA,
            Genes = names(scores),
            Used = used,
            Passed = 0
      )
    )
  }
}

#' Finds the limit cutoff when target and sentinal paths are given
#'
#' Using all paths from an igraph path trace this will filter all the paths in
#' the trace by the lim argument. Implementation deploys path_filter() in 
#' parallel over the number of cores specified by the cores argument.
#'
#' @export
#' @param vertices the vertices vector from find_limit
#' @param path_obj pathtrace object from short_paths()
#' @param lim the filter limit value
#' @param path_name the name of the path for message (target or sentinal)
#' @param weights the gene OMICS scores or score metric by verex to weight paths
#' @param cores the cores to run the path calulation over. default = 1
#' @return list of vertex names from paths that path the filter limit
#' @examples 
#' example_path <- list()
#' example_path$res <- list(
#'    c('319', "49", "23", "86", "690", "238"),
#'    c('422', "899", "37", "240", "970", "28")
#' )
#' names(example_path$res[[1]]) <- c(
#'    "GeneA","GeneZ", "GeneAlpha",
#'    "GeneB", "GeneX", "GeneOmega"
#' )
#' names(example_path$res[[2]]) <- c(
#'    "Gene1","Gene2", "GeneUno",
#'    "Gene12", "Gene13", "GeneOcho"
#' )
#'
#' sampweights <- c(1.45, 2.45, 0.89, .003, 1.3, 2.1, 0.02, 0, 0, 0.2, 0.6, .70)
#' names(sampweights) <- c(
#'    "GeneA","GeneZ", "GeneAlpha",
#'    "GeneB", "GeneX", "GeneOmega",
#'    "Gene1","Gene2", "GeneUno",
#'    "Gene12", "Gene13", "GeneOcho"
#' )
#' path_obj_filter(
#'    path_obj = example_path,
#'    path_name =  "test",
#'    vertices = c(0,0,0),
#'    weights = sampweights,
#'    lim = 1,
#'    cores = 2
#' )
path_obj_filter <- function( 
  path_obj, vertices, lim, path_name, weights, cores=1 
) {
  ## Target path filtering
  if ((length(vertices) == 1) & 
      is.na(vertices)[1]) {
    message('Target Path Trace Returned No Paths')
  }else{
    # Filter the target Paths
    fitlered_paths <- parallel::mclapply( 
      path_obj$res,  
      path_filter, 
      lim=lim, 
      weights=weights,
      mc.cores = cores
    ) 
    
    filter_summary <- do.call(Map, c(c, fitlered_paths))
    
    vertices <- filter_summary$Keep[ !duplicated(filter_summary$Keep)]
    vertices <- vertices[!(is.na(vertices))]
    
    all_vertices <- filter_summary$Genes[ !duplicated(filter_summary$Genes)]
    all_vertices <- all_vertices[!(is.na(all_vertices))]
    
    #Paths Kept - ISSUE Here
    message(paste0(
      sum(filter_summary$Passed),
      " ",
      path_name,
      " paths out of ",
      length(path_obj$res), 
      " ( ",
      signif(
        sum(filter_summary$Passed) / length(path_obj$res),
        4
      ) * 100,
      "% ) kept"
    ))
    
    #Genes Filtered Out
    message(paste0(
      length(all_vertices) - length(vertices),
      ' of ',
      length(all_vertices),
      ' ', 
      path_name,
      ' genes ( ', 
      signif(
        (length(all_vertices) - length(vertices)) / length(all_vertices), 
        4
      ) * 100,
      '% ) filtered out'
    ))
  }
  
  return( vertices )
}

#' Finds the limit cutoff when target and sentinal paths are given
#'
#' Limit is based on whichever value is greater between; the mean or median 
#' scores of the Target-traced-to-sentinal paths. If there are no target to
#'  sentinal paths the pairwise within target traces are used.
#'
#' @export
#' @param s_path the sentinel path object
#' @param t_path the target path object
#' @param weights OMICS/vertex scores as named vector
#' @param cores the cores to run the path calulation over. default = 1
#' @return List object with named attributes of limit (the path score cutoff 
#' limit) and both t_verts and s_verts objects used as placeholders to determine
#' if the paths were scorable. 
#' @examples 
#' 
#' example_path <- list()
#' example_path$res <- list(
#'    c('319', "49", "23", "86", "690", "238"),
#'    c('422', "899", "37", "240", "970", "28")
#' )
#' names(example_path$res[[1]]) <- c(
#'    "GeneA","GeneZ", "GeneAlpha",
#'    "GeneB", "GeneX", "GeneOmega"
#' )
#' names(example_path$res[[2]]) <- c(
#'    "Gene1","Gene2", "GeneUno",
#'    "Gene12", "Gene13", "GeneOcho"
#' )
#' example_path_b <-  example_path 
#' sampweights <- c(1.45, 2.45, 0.89, .003, 1.3, 2.1, 0.02, 0, 0, 0.2, 0.6, .70)
#' names(sampweights) <- c(
#'    "GeneA","GeneZ", "GeneAlpha",
#'    "GeneB", "GeneX", "GeneOmega",
#'    "Gene1","Gene2", "GeneUno",
#'    "Gene12", "Gene13", "GeneOcho"
#' )
#' find_limit(
#'    s_path = example_path,
#'    t_path =  example_path_b,
#'    weights = sampweights,
#'    cores = 2
#' )
find_limit <- function ( s_path, t_path, weights, cores=1) {
  #if there are no sentinel paths return empty and look at target
  if(length(s_path$res)==0) {
    s_path$res <- 1
  }
  if(length(t_path$res)==0) {
    t_path$res <- 1
  }
  if ((length(s_path$res) == 1) & (length(s_path$res[[1]]) == 1)) {
    sent_keep_vertices <- NA
    
    if ((length(t_path$res) == 1) & (length(t_path$res[[1]]) == 1)) {
      target_vertices <- NA
      limit <- NA
    }else{
      target_vertices <- c(0,0,0)
      #set limit based on targets
      scores <- do.call(c, parallel::mclapply( 
        t_path$res,  
        path_calc, 
        weights=weights,
        mc.cores = cores
      ))
      
      # Look at limit of highest between mean or median non-zero sentinel paths
      target_summary <- summary(
        scores[(scores > 0) & (is.na(scores) == F)]
      )
      if (target_summary['Mean'] > target_summary['Median']) {
        limit <- target_summary['Mean']
      }else{
        limit <- target_summary['Median']
      }
    }
  }else{
    target_vertices <- c(0,0,0)
    sent_keep_vertices <- c(0,0,0)
    sent_scores <- do.call(c, parallel::mclapply( 
      s_path$res,  
      path_calc, 
      weights=weights,
      mc.cores = cores
    ))
    
    # Look at limit of highest between mean or median non-zero sentinal paths
    sentinal_summary <- summary(
      sent_scores[(sent_scores > 0) & (is.na(sent_scores) == F)]
    )
    if (sentinal_summary['Mean'] > sentinal_summary['Median']) {
      limit <- sentinal_summary['Mean']
    }else{
      limit <- sentinal_summary['Median']
    }
  }
  return(list(
    cutoff = limit,
    t_verts = target_vertices,
    s_verts = sent_keep_vertices
  ))
}  

#' Traces the shortest paths of a gene to a vector of target and a vector of 
#' sentinal genes 
#' 
#' Traces the the shortest paths of target gene paiwise to the target gene list.
#' then traces the shortest path of the gene to a list of sentinal genes. This
#' function returns the list of genes in a path object which score over the an
#' emperically defined limit using the find_limit() function.
#' 
#' @export 
#' @param tnet igraph network (Main entire network) eg. net/net_undirected/JS_net_undirected
#' @param target the from gene target eg Genes[1]
#' @param targets List of the total list of targets in the User set eg. Genes
#' @param sentinals List of the sentinal genes to trace to eg. Sentinal
#' @param cores the number of cores to use in paralellizable functions
#' @return List object of Inter genes from target path traces and Sentinal genes
#' from sentinal gene traces
#' @examples 
#' 
short_paths <- function( tnet, target, targets, sentinals, cores = 1 ){
  
  message( paste0( 'Working on: ', target))
  # Pull paths that have median OMICS Score. ( Need to integrate a Genetics+Genomics Measure )
  omics_scores <- setNames(igraph::V(tnet)$weight, igraph::V(tnet)$name)
  
  # All Shortest paths from target to Target Genes directed
  snet <- igraph::simplify(
    tnet,
    remove.multiple = TRUE,
    remove.loops = FALSE,
  )
  paths <- igraph::get.all.shortest.paths(
    snet,
    from = target,
    to = igraph::V(snet)[ names(igraph::V(snet)) %in% targets ],
    mode = c("all")
  ) 
  
  # All Shortest paths from target to Sentinel Genes
  sent_paths <- igraph::get.all.shortest.paths(
    snet,
    from = target,
    to = igraph::V(snet)[ names(igraph::V(snet)) %in% sentinals ],
    mode = c("all")
  )
  
  ## Find the limit cut off
  cutoff_obj <- find_limit(
    s_path = sent_paths,
    t_path = paths,
    weights = omics_scores,
    cores = cores
  )
  
  ## Target path filtering
  t_vertices <- path_obj_filter(
    path_obj = paths,
    vertices = cutoff_obj$t_verts,
    lim = cutoff_obj$cutoff,
    path_name ='target',
    weights = omics_scores,
    cores = cores
  )
  
  ## Sentinal path filtering
  s_vertices <- path_obj_filter(
    path_obj = sent_paths,
    vertices = cutoff_obj$s_verts,
    lim = cutoff_obj$cutoff,
    path_name ='sentinel',
    weights = omics_scores,
    cores = cores
  )
  return(list( Inter = t_vertices, Sentinal=s_vertices))
}

#' Process a path trace list
#'
#' This function takes a path trace object from short_paths() and transforms it
#' into the list of genes to keep for the sub network generation.
#' 
#' @export 
#' @param path_obj path trace object from short_paths()
#' @return a list of genes from the path trace to use for filtering the main
#' network
#' @examples 
#' obj <- list(list( 
#'     Inter = c("GeneA","GeneZ", "GeneAlpha",
#'        "GeneB", "GeneX", "GeneOmega",
#'        "Gene1","Gene2", "GeneUno",
#'        "Gene12", "Gene13", "GeneOcho"),
#'     Sentinal = c("GeneEh","GeneZee", "GeneAlpha",
#'        "GeneB", "GeneXray", "GeneOmega",
#'        "Gene11","Gene21", "GeneUno",
#'        "Gene1", "Gene3", "GeneOcho")
#' ))
#' trace_filter(obj)
trace_filter <- function (path_obj) {
  #collapse Pairwise Pass genes and Sentinal Path Genes
  list_tar <- NULL
  sentinal_tar <- NULL
  len_lts <- NULL
  len_sts <- NULL
  for( i in 1:length(path_obj) ){
    len_lts <- c( len_lts, length( path_obj[[i]]$Inter ) )
    len_sts <- c( len_sts, length( path_obj[[i]]$Sentinal ) )
    list_tar <- c(list_tar, path_obj[[i]]$Inter)
    sentinal_tar <- c(sentinal_tar, path_obj[[i]]$Sentinal)
  }
  
  length( list_tar[!duplicated(list_tar)] )
  length( sentinal_tar[!duplicated(sentinal_tar)] )
  
  table( list_tar[!duplicated(list_tar)] %in% sentinal_tar[!duplicated(sentinal_tar)] )
  
  gene_list <- sentinal_tar[!duplicated(sentinal_tar)][ 
    sentinal_tar[ !duplicated(sentinal_tar) ] %in% 
      list_tar[ !duplicated(list_tar) ]
    ]
  return(gene_list)
}

#' Push Network to Synapse
#'
#'This function takes a network object and pushes it to synapse
#'
#' @export
#' @param network the igraph network to push to synapse eg. net
#' @param net_filename the file name of the network without file extension
#' @param net_synname the desplay name of the network in synapse
#' @param p_id the parent synapse ID of the network destination
#' @param folder the name of the storage folder in the parent synapse ID to 
#' store the net
#' @param act_name the name of the syn activity object to 
#' @param act_desc the description of the syn activity object to 
#' @param code the path of the code which generated the network for the 
#' provenance (optional)
#' @param repo the repo which generated the network for the 
#' provenance (optional)
#' @param syn_used character vector of synIDs to seed the provenance (optional)
#' @param subset An vector of vertex names to filter the network for (optional) 
#' eg. test
#' @param prov_object A pre made github code provenance object or vector of 
#' objects
#' @return a synapse entity of a .graphml subnetwork object stored in synapse
#' @examples 
#' 
store_net <- function (network, net_filename, net_synname,
                       p_id, folder, act_name, act_desc,
                       code=NULL, repo=NULL,
                       syn_used=NULL, subset=NULL,
                       prov_object = NULL) {
  #Set Activity
  activity <- syn_temp$store(synapseclient$Folder(
    name = folder, 
    parentId = p_id
  ))
  
  #Set annotations
  all.annotations = list(
    dataType = 'Network',
    summaryLevel = 'gene',
    assay	 = 'RNAseq',
    tissueTypeAbrv	= c('IFG', 'STG', 'FP', 'PHG', 'TCX', 'DLFPC'), 
    study = c( 'MSBB', 'ROSMAP', 'Mayo' ), 
    organism = 'HomoSapiens',
    consortium	= 'TreatAD',
    genomeAssemblyID = 'GRCh38'
  )
  
  #Subset the network if there is a vertex vector given 
  if (!is.null(subset)) {
    network <- igraph::induced_subgraph(
      network, v=igraph::V(network)[ names(igraph::V(network)) %in% subset ]
    )
  }
  #eg. IGRAPH ff4b668 DN-- 486 6119 -- 
  sub_net_simple <- igraph::simplify(
    network,
    remove.multiple = TRUE,
    remove.loops = FALSE,
    edge.attr.comb = list( interaction = "concat", 
                           Occurance = "concat",
                           UniqCol = "concat",
                           pathway = "concat", 
                           EdgeRep = "mean",
                           Edge = "random",
                           SumOccurence = "mean",
                           DLPFC_CE = "mean",
                           CBE_CE = "mean",
                           FP_CE = "mean",
                           IFG_CE = "mean",
                           PHG_CE = "mean",
                           STG_CE = "mean",
                           TCX_CE = "mean",
                           Avg_Cortex_CE = "mean",
                           Avg_All_CE = "mean"
    )
  )
  
  # Github link - "jgockley62/igraph_Network_Expansion" 
  if (is.null(prov_object)) {
    if (!is.null(repo) | !is.null(code)) {
      this_repo <- githubr::getRepo(
        repository = repo,
        ref="branch",
        refName='master'
      )
      this_file <- githubr::getPermlink(
        repository = this_repo,
        repositoryPath = code
      )
    }else{
      this_file <- NULL
    }
  }else{
    this_file <- prov_object
  }
  
  # write file
  igraph::write_graph(
    network,
    paste0( '~/igraph_Network_Expansion/', net_filename,'.graphml'),
    format = "graphml"
  )
  # push file
  enrich_obj <-  syn_temp$store(
    synapseclient$File(
      path=paste0( '~/igraph_Network_Expansion/', net_filename,'.graphml'),
      name = net_synname,
      parentId=activity$properties$id ),
    used = syn_used,
    executed = this_file,
    activityName = act_name,
    activityDescription = act_desc
  )
  
}

#' Re-Load the  networks and calc meterics
#' 
#' This Function Loads a network object from synapse using its synapse ID and
#' the type of network, ie graphml.
#' 
#' @export
#' @param syn_id the networks synID
#' @param form the format of the netwrok file eg. "graphml"
#' @return an igraph network object
#' @examples 
#' 
network_load <- function (syn_id, form) {
  import_net <- igraph::read_graph(
    file = syn_temp$get(syn_id)$path,
    format = form
  )
  return(import_net)
}
