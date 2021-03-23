#_# This script loads and filters Greg Cary's process enrichments
#_# Input are two syn IDs for the GeneSet process and Gene assignments for the processes

setwd('~/igraph_Network_Expansion/')

if ("fitdistrplus" %in% rownames(installed.packages()) == FALSE) {
  install.packages("fitdistrplus")
}
##Insert User Credentials here. DO NOT Save DO NOT push to git or docker hub
user <- '< USER_ID >'
pass <- '< Password >'

reticulate::use_python("/home/jgockley/.local/share/r-miniconda/envs/r-reticulate/bin/python", required = TRUE)
synapseclient <- reticulate::import("synapseclient")
syn_temp <- synapseclient$Synapse()
syn_temp$login( user, pass, '--rememberMe' )

rm(user)
rm(pass)

#Annotation of Genes to Biological Domains 
genes <- 'syn24827928' 
processes <- 'syn24827958'

genes <- readRDS(file = syn_temp$get(genes)$path)
processes <- readRDS(file = syn_temp$get(processes)$path)

sig_terms <- processes[processes$padj < 0.05, ]$pathway

#' @param proc a list object of processes eg. genes
#' @param term a go term eg. sig_terms[1] 
biodomain_enumerator <- function (proc, term) {
  return(
    length(genes[ proc$GOterm_Name == term, ]$Biodomain)
  )
}
#' @param proc a list object of processes eg. genes
#' @param term a go term eg. sig_terms[1] 
biodomain_translator <- function (proc, term) {
  biods <- proc[ proc$GOterm_Name == term, ]$Biodomain
  output <- cbind(
    c(rep(term, length(biods))),
    biods
    )
  return(output)
}

# Mapping of GoTerms to Biodomain
table(sapply(sig_terms, biodomain_enumerator, proc=genes))

# All GoTerms for domains
table(sapply(processes$pathway, biodomain_enumerator, proc=genes))

# DF of GoTerm to Biodomains
translator <- as.data.frame(
  do.call(
    rbind,
    sapply(sig_terms, biodomain_translator, proc=genes)
  )
)

# DF of GoTerm to Biodomains
translator_All <- as.data.frame(
  do.call(
    rbind,
    sapply(processes$pathway, biodomain_translator, proc=genes)
  )
)
leadingedge_genes <- list()
all_genes <- list()

biodomains <- names(table(translator$biods))
biodomains_All <- names(table(translator_All$biods))

######################################################################
## All genes in GoTerms

# Leading Edge finder
#' @param proc a list object of processes eg. genes
#' @param bd a biodomain term eg. biodomains[1]
#' @param tran translation of go to biodomaing 
all_numerator <- function (proc, tran, bd) {

  pull_goterms <- tran[tran[,2]==bd,][,1]
  proc[proc$GOterm_Name %in% pull_goterms,]$ensembl_id
  ensgs <- do.call( c, proc[proc$GOterm_Name %in% pull_goterms,]$ensembl_id)
  
  return(ensgs)
}

# Plot leading edge reccurence
par(mfrow=c(4,4))
for(go in biodomains) {
  message(go)
  hist(
    table(
      all_numerator(
        genes,
        translator,
        go)
      ),
    main = go,
    xlab='Gene Recurences')
}

# Leading edge Unique Gene number
uniq_all <- as.list(biodomains)
names(uniq_all) <- biodomains
for(go in biodomains) {
  message(go)
  enns <- all_numerator(
    genes,
    translator,
    go
  )
  message(length(enns[!duplicated(enns)]))
  uniq_all[[go]] <- enns[!duplicated(enns)]
}

# Plot leading edge reccurence
par(mfrow=c(4,4))
for(go in biodomains) {
  message(go)
  hist(
    table(
      all_numerator(
        genes,
        translator,
        go)
    ),
    main = go,
    xlab='Gene Recurences')
}

# All GosxBiodomain
all_biodomain <- as.list(biodomains_All)
names(all_biodomain) <- biodomains_All
for(go in biodomains_All) {
  message(go)
  enns <- all_numerator(
    genes,
    translator_All,
    go
  )
  message(length(enns[!duplicated(enns)]))
  all_biodomain[[go]] <- enns[!duplicated(enns)]
}
# Plot leading edge reccurence
par(mfrow=c(4,4))
for(go in biodomains_All) {
  message(go)
  hist(
    table(
      all_numerator(
        genes,
        translator_All,
        go)
    ),
    main = go,
    xlab='Gene Recurences')
}
############################################
## Greg's leading edges


# Leading Edge finder
#' @param proc a list object of processes eg. processes
#' @param bd a biodomain term eg. biodomains[1]
#' @param tran translation of go to biodomaing 
leadedge_numerator <- function (proc, tran, bd) {
  
  pull_goterms <- tran[tran[,2]==bd,][,1]
  ensgs <- do.call(c, proc[proc$pathway %in% pull_goterms,]$leadingEdge)
  
  return(ensgs)
}

# Plot leading edge reoccurence
par(mfrow=c(4,4))
for(go in biodomains) {
  message(go)
  hist(
    table(
      leadedge_numerator(
        processes,
        translator,
        go)
    ),
    main = go,
    xlab='Gene Recurences')
}

# Leading edge Unique Gene number
uniq_lead <- as.list(biodomains)
names(uniq_lead) <- biodomains
for(go in biodomains) {
  message(go)
  enns <- leadedge_numerator(
    processes,
    translator,
    go
  )
  message(length(enns[!duplicated(enns)]))
  uniq_lead[[go]] <- enns[!duplicated(enns)]
}

# Pull all genes for the complete BioDomains
domain_lead <- as.list(biodomains)
names(domain_lead) <- biodomains
for(go in biodomains) {
  message(go)
  enns <- leadedge_numerator(
    processes,
    translator_All,
    go
  )
  message(length(enns[!duplicated(enns)]))
  domain_lead[[go]] <- enns[!duplicated(enns)]
}

################################################################################
## Pull Logsdon Scores

## Function to pull gene scores from the overall score table
#'@param genevec a vector of ENSGs to pull Logsdon scores ie. uniq_lead$Myelination
#'

# SELECT Overall FROM syn24168007 WHERE ENSG = 'ENSG00000149269'
# SELECT ENSG,GeneName,OmicsScore,Overall FROM syn24168007 WHERE ENSG IN ('ENSG00000149269')
pull_lead <- as.list(biodomains)
names(pull_lead) <- biodomains

pull_all <- as.list(biodomains)
names(pull_all) <- biodomains

pull_biodomain <- as.list(biodomains_All)
names(pull_biodomain) <- biodomains_All

for (bd in biodomains) { 
  #Leading Edge Genes
  df <- read.csv(
    syn_temp$tableQuery(
      query = paste0(
        'SELECT ENSG,GeneName,OmicsScore,Overall FROM syn24168007 WHERE ENSG IN (\'',
        paste0( uniq_lead[[bd]], collapse='\', \'' ),
        '\')'),
      resultsAs = 'csv')$filepath
    )[, c('ENSG', 'GeneName', 'OmicsScore', 'Overall')]
  df$logOmics <- log2(df$OmicsScore)
  df$logOverall <- log2(df$Overall)
  pull_lead[[bd]] <- df
  
  #All Genes
  df <- read.csv(
    syn_temp$tableQuery(
      query = paste0(
        'SELECT ENSG,GeneName,OmicsScore,Overall FROM syn24168007 WHERE ENSG IN (\'',
        paste0( uniq_all[[bd]], collapse='\', \'' ),
        '\')'),
      resultsAs = 'csv')$filepath
  )[, c('ENSG', 'GeneName', 'OmicsScore', 'Overall')]
  df$logOmics <- log2(df$OmicsScore)
  df$logOverall <- log2(df$Overall)
  pull_all[[bd]] <- df
  
  #Complete Biodomains
  df <- read.csv(
    syn_temp$tableQuery(
      query = paste0(
        'SELECT ENSG,GeneName,OmicsScore,Overall FROM syn24168007 WHERE ENSG IN (\'',
        paste0( domain_lead[[bd]], collapse='\', \'' ),
        '\')'),
      resultsAs = 'csv')$filepath
  )[, c('ENSG', 'GeneName', 'OmicsScore', 'Overall')]
  df$logOmics <- log2(df$OmicsScore)
  df$logOverall <- log2(df$Overall)
  pull_biodomain[[bd]] <- df
  
}
######################################################################
library(fitdistrplus)

all_scores <- read.csv(
  syn_temp$tableQuery(
    query = paste0(
      'SELECT ENSG,GeneName,OmicsScore,Overall FROM syn24168007'),
    resultsAs = 'csv')$filepath
)[, c('ENSG', 'GeneName', 'OmicsScore', 'Overall')]


#### Functionaize...

lapply( pull_lead,dim )
lapply( pull_all,dim )
lapply( pull_biodomain,dim )

data <- pull_lead[[bd]]
look_dist <- pull_lead[[bd]]$Overall

#' @param dat list object of dataframes eg pull_all or pull_lead
#' @param val name of a data frame in list object dat eg. biodomains[1]
#' @param Max max list size eg. 50
#' @param Zstart the z value to start the pruning eg. 0.5 - 4 in increments
#' @param SDstep the value to step increase until geneset is lesst than max
list_generator <- function(dat, val, Max = 50, Zstart=1, Zstep = 0.5) {
  
  #Extract the dataframe
  data <- dat[[val]]
  # Extract the overall scores
  look_dist <- data$Overall
  
  #Return if value is less than user threshold
  if (length(look_dist) < 50){
    return(data$GeneName)
  }
  # Test the 4 fit models
  fitg <- summary( fitdistrplus::fitdist( look_dist, "gamma" ) )
  fitln <- summary( fitdistrplus::fitdist( look_dist, "lnorm" ) )
  fitW <- summary( fitdistrplus::fitdist( look_dist, "weibull" ) )
  fitn <- fitdistrplus::fitdist(look_dist,"norm")
  
  SumStat <- fitdistrplus::gofstat(
    list(fitW, fitg, fitln, fitn),
    fitnames=c("Weibull", "gamma", "lognormal", "normal")
  )
  
  # Examine if the 4 models converge
  if (!(
    names(SumStat$bic[SumStat$bic == min(SumStat$bic) ]) == 
    names(SumStat$aic[ SumStat$aic == min(SumStat$aic) ])
  )) {
    message(paste0( val, ": MODEL DO NOT CONVERGE"))
  }
  
  bic_stats <- as.data.frame(
    cbind( 
      BIC = SumStat$bic,
      delBIC = SumStat$bic - min(SumStat$bic),
      relLik. = exp(-0.5 * (SumStat$bic - min(SumStat$bic))),
      bicweight =
        exp(-0.5 * (SumStat$bic-min(SumStat$bic))) / 
        sum(exp(-0.5 * (SumStat$bic-min(SumStat$bic))))
    )
  )
  message(paste0(paste0(val, ":Overall Scores Best Predicted by ", 
                  row.names(bic_stats[bic_stats$bicweight == max(bic_stats$bicweight), ]),
                  " model measured by BIC"))
  )
  aic_stats <- as.data.frame(
    cbind( 
      AIC = SumStat$aic,
      delAIC = SumStat$aic-min(SumStat$aic),
      relLik = exp(-0.5 * (SumStat$aic - min(SumStat$aic))),
      aicweight = 
        exp(-0.5 * (SumStat$aic - min(SumStat$aic))) / 
        sum(exp(-0.5 * (SumStat$aic - min(SumStat$aic))))
         )
  )
  
  model_name <- row.names(
    bic_stats[bic_stats$bicweight == max(bic_stats$bicweight), ]
  )
  
  message(paste0(val, ": Overall Scores Best Predicted by ", 
                  row.names(aic_stats[aic_stats$aicweight == max(aic_stats$aicweight), ]),
                  " model measured by AIC")
  )
  
  # Examine the liklehood that the model is best
  if ( !(
    aic_stats[ model_name , ]$aicweight > .80 &
    bic_stats[ model_name , ]$bicweight > .80 
  )) {
    message(paste0( val, ": Consenus Model Less than 80% Likely!"))
  }
  
  # Assign Probability if Best Model is Normal
  if(model_name == 'normal') {
    data$probability <- stats::pnorm( 
      look_dist,
      mean = mean(look_dist), 
      sd = sd(look_dist),
      lower.tail = FALSE, 
      log.p = FALSE
    )
  }
  # Assign Probability if Best Model is gamma
  if(model_name == 'gamma') {
    data$probability <- stats::pgamma( 
      look_dist,
      shape = fitg$sd['shape'], 
      rate = fitg$sd['rate'],
      lower.tail = FALSE, 
      log.p = FALSE
    )
  }
  # Assign Probability if Best Model is Weibull
  if(model_name == 'Weibull') {
    data$probability <- stats::pweibull( 
      look_dist,
      shape = fitW$estimate['shape'], 
      scale = fitW$estimate['scale'],
      lower.tail = FALSE, 
      log.p = FALSE
    )
  }
  # Assign Probability if Best Model is Log Normal
  if(model_name == 'lognormal') {
    data$probability <- stats::plnorm( 
      look_dist,
      meanlog = mean(log(look_dist)), 
      sdlog = sd(log(look_dist)),
      lower.tail = FALSE, 
      log.p = FALSE
    )
  }
  
  # Order the geneset based on assigned probability
  data <- data[ order(data$probability),  ]
  
  #Step the number of sds ontop of the Zstart till the remaining geneset is < 50
  size <- dim(data)[1]
  iteration <- 0
  while (size > 50) {
    cutoff <- 1 - pnorm(iteration * Zstep + Zstart)
    output <- data[ data$probability < cutoff, ]
    GeneList <- output$GeneName
    
    if(dim(output)[1] < 50){
      message(paste0( val, " Gene Set for Tracing is ", length(GeneList), " genes"))
    }
    size <- dim(output)[1]
    iteration <- iteration+1
  }
  return(GeneList)
}

genelist_lead <- sapply(
  biodomains, 
  list_generator,
  dat=pull_lead,
  Max = 50,
  Zstart=.5,
  Zstep = 0.5
)  

genelist_all <- sapply(
  biodomains, 
  list_generator,
  dat=pull_all,
  Max = 50,
  Zstart=.5,
  Zstep = 0.5
)  

genelist_biodomain <- sapply(
  biodomains, 
  list_generator,
  dat=pull_biodomain,
  Max = 50,
  Zstart=.5,
  Zstep = 0.15
)  

for (i in 1:length(genelist_lead)) {
  lead <- as.numeric(
    table(
      genelist_lead[[i]] %in% genelist_all[[i]]
      )['TRUE']
  )
  all <- as.numeric(
    table(
      genelist_all[[i]] %in% genelist_lead[[i]]
      )['TRUE']
  )
  if (all == lead) {
    message(paste0(
      "For ",
      biodomains[i],
      ": Total Overlap = ",
      all,
      " Percent All = ",
      all / length(genelist_all[[i]]),
      " Percent Lead = ",
      lead / length(genelist_lead[[i]])
      
    ))
  }else{
    warning(paste0( "Issue matching genes for: ", biodomains[i] ))
  }
  
}


### Comp the All Biodomains
for (i in 1:length(genelist_lead)) {
  lead <- as.numeric(
    table(
      genelist_biodomain[[i]] %in% genelist_all[[i]]
    )['TRUE']
  )
  all <- as.numeric(
    table(
      genelist_all[[i]] %in% genelist_biodomain[[i]]
    )['TRUE']
  )
  if (all == lead) {
    message(paste0(
      "For ",
      biodomains[i],
      ": Total Overlap = ",
      all,
      " Percent All = ",
      all / length(genelist_all[[i]]),
      " Percent Biodomain = ",
      lead / length(genelist_biodomain[[i]])
      
    ))
  }else{
    warning(paste0( "Issue matching genes for: ", biodomains[i] ))
  }
  
}

for (i in 1:length(genelist_lead)) {
  lead <- as.numeric(
    table(
      genelist_lead[[i]] %in% genelist_biodomain[[i]]
    )['TRUE']
  )
  all <- as.numeric(
    table(
      genelist_biodomain[[i]] %in% genelist_lead[[i]]
    )['TRUE']
  )
  if (all == lead) {
    message(paste0(
      "For ",
      biodomains[i],
      ": Total Overlap = ",
      all,
      " Percent Biodomain = ",
      all / length(genelist_biodomain[[i]]),
      " Percent Lead = ",
      lead / length(genelist_lead[[i]])
      
    ))
  }else{
    warning(paste0( "Issue matching genes for: ", biodomains[i] ))
  }
  
}

#Write Lists to File: 
##test for destination paths:
rec_paths <- c('InputList/BiodomainLists',
               'InputList/BiodomainLists/LeadingEdge',
               'InputList/BiodomainLists/AllGoTerm',
               'InputList/BiodomainLists/Biodomains',
               'InputList/BiodomainLists/Biodomains_Total')
for (path in rec_paths) {
  if (dir.exists(path)){
  }else{
    dir.create(path)
  }
}

##write to file:
for (bd in biodomains) {
  data.table::fwrite(
    as.list(genelist_all[[bd]]),
    file = paste0( 'InputList/BiodomainLists/AllGoTerm/',
                   gsub( ' ', '_',  bd),
                   '.txt'
                  ),
    quote = FALSE,
    sep = '\n')
  data.table::fwrite(
    as.list(genelist_lead[[bd]]),
    file = paste0( 'InputList/BiodomainLists/LeadingEdge/',
                   gsub( ' ', '_',  bd),
                   '.txt'
    ),
    quote = FALSE,
    sep = '\n')
}
for (bd in biodomains) { 
  data.table::fwrite(
    as.list(genelist_biodomain[[bd]]),
    file = paste0( 'InputList/BiodomainLists/Biodomains/',
                   gsub( ' ', '_',  bd),
                   '.txt'
    ),
    quote = FALSE,
    sep = '\n')
  
  data.table::fwrite(
    as.list(pull_biodomain[[bd]]$GeneName),
    file = paste0( 'InputList/BiodomainLists/Biodomains_Total/',
                   gsub( ' ', '_',  bd),
                   '.txt'
    ),
    quote = FALSE,
    sep = '\n')
}

#Push to synapse: 
parent_id <- 'syn25185168'
activity_name = 'Gene Lists';
activity_description = 'List of genes for pathway tracing';

this_filename <- 'Pull_Gregs_BioDomainAnalysis.R'
this_repo <- githubr::getRepo(repository = "jgockley62/igraph_Network_Expansion", 
                             ref="branch",
                             refName='master'
                            )
this_file <- githubr::getPermlink(repository = this_repo,
                                 repositoryPath=paste0('code/',
                                                       this_filename
                                                       )
                                 )

# Set annotations
all_annotations = list(
  dataType = 'Gene Lists',
  dataSubType = 'HGNC Gene Names',
  summaryLevel = 'gene',
  assay	 = 'AD Gene Ranks',
  tissueTypeAbrv	= NULL, 
  study = 'TREAT-AD', 
  organism = 'HomoSapiens',
  consortium	= 'TREAT-AD'
)

syns_used <- c('syn24827928', 'syn24827958')

### Store files in synapse
## Push Leading Edge
code <- syn_temp$store(synapseclient$Folder(name = "All Leading Edge Lists",
                                            parentId = parent_id))

for (fil_n in list.files('InputList/BiodomainLists/LeadingEdge/')) {
  enrich_obj <- syn_temp$store(
    synapseclient$File(
      path=paste0('InputList/BiodomainLists/LeadingEdge/',fil_n),
      name = gsub(
        '_', ' ', gsub('.txt',  '', fil_n)
      ), parentId=code$properties$id ),
    used = syns_used,
    activityName = activity_name,
    executed = this_file,
    activityDescription = activity_description
  )
  syn_temp$setAnnotations(enrich_obj, annotations = all_annotations)
} 

## Push All GoTerm
code <- syn_temp$store(synapseclient$Folder(name = "All GoTerm Edge Lists",
                                            parentId = parent_id))

for (fil_n in list.files('InputList/BiodomainLists/AllGoTerm/')) {
  enrich_obj <- syn_temp$store(
    synapseclient$File(
      path=paste0('InputList/BiodomainLists/AllGoTerm/',fil_n),
      name = gsub(
        '_', ' ', gsub('.txt',  '', fil_n)
      ), parentId=code$properties$id ),
    used = syns_used,
    activityName = activity_name,
    executed = this_file,
    activityDescription = activity_description
  )
  syn_temp$setAnnotations(enrich_obj, annotations = all_annotations)
} 


#### Biodomain filtered list
## Push All GoTerm
code <- syn_temp$store(synapseclient$Folder(name = "Biodomain Filtered Genes",
                                            parentId = parent_id))

for (fil_n in list.files('InputList/BiodomainLists/Biodomains/')) {
  enrich_obj <- syn_temp$store(
    synapseclient$File(
      path=paste0('InputList/BiodomainLists/Biodomains/',fil_n),
      name = gsub(
        '_', ' ', gsub('.txt',  '', fil_n)
      ), parentId=code$properties$id ),
    used = syns_used,
    activityName = activity_name,
    executed = this_file,
    activityDescription = activity_description
  )
  syn_temp$setAnnotations(enrich_obj, annotations = all_annotations)
} 

## Push All GoTerm
code <- syn_temp$store(synapseclient$Folder(name = "Biodomain All Genes By Domain",
                                            parentId = parent_id))

# Biodomain total list
for (fil_n in list.files('InputList/BiodomainLists/Biodomains_Total/')) {
  enrich_obj <- syn_temp$store(
    synapseclient$File(
      path=paste0('InputList/BiodomainLists/Biodomains_Total/',fil_n),
      name = gsub(
        '_', ' ', gsub('.txt',  '', fil_n)
      ), parentId=code$properties$id ),
    used = syns_used,
    activityName = activity_name,
    executed = this_file,
    activityDescription = activity_description
  )
  syn_temp$setAnnotations(enrich_obj, annotations = all_annotations)
} 
