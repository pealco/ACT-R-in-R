########################################################################################################
###
###  Simple Memory Retrieval Model of Agreement Phenomena (based on ACT-R theory)
###    With HTML trace output and PDF activation plots
###        
###
###  Rick Lewis & William Badecker (rickl@umich.edu)
###  Version 3.0
###  7 Feb 2007
###
###
#########################################################################################################


## act-r library
source("actr.r")


## SET PARAMETERS function: sets the parameters as global variables

set.parameters <- function(p) {
    ## latency factor
    cat.penalty <<- p$cat.penalty  
    
    ## latency factor
    F <<- p$F
    
    ## total source activation
    G <<- p$G
    
    ## activation noise parameter for logistic distribution
    ans <<- p$ans
    
    ## fan parameter
    mas <<- p$mas
    
    ## base level decay parameter
    d <<- p$d
    
    ## match penalty
    match.penalty <<- p$match.penalty
    
    ## do VAR cues impose a mismatch penalty?
    var.mismatch.penalty <<- p$var.mismatch.penalty
    
    ## additional fan associated with VAR retrieval cues
    VAR.fan <<- p$VAR.fan
    
    ## distinctiveness parameters
    modulate.by.distinct <<- p$modulate.by.distinct
    distinctiveness <<- p$distinctiveness
}

## RUN THE MODEL given the global parameter settings

run.model <- function(quiet=TRUE) {
    trials <<- default.trials  
    
    ## Read in the item definitions: each item is a feature vector
    items <<- read.delim(file=item.file,header=FALSE,colClasses="character")
    num.columns <<- length(items)
    num.features <<- length(items[,1])-2
    
    creation.moment <<- c()
    for (cm in items[2,2:num.columns]) {
        creation.moment <<- c(creation.moment,as.integer(cm))
    }
    
    item.name <<- c()
    for (it in items[1,2:num.columns]) {
        item.name <<- c(item.name,it)
    }
    
    num.items <<- length(item.name)
    item.features <<- t(items[3:(num.features+2),2:num.columns])
    
    ## create initial history matrix to be the creation moments
    history <<- matrix(1:num.items,nrow=trials, ncol=num.items, byrow=TRUE)
    moments <<- creation.moment
        
    ## Read in the schedule of retrievals
    retrievals <- read.delim(file=retrieval.file,header=FALSE)
    cue.names <- retrievals[2:(num.features+1),1]
    retrievals <- t(retrievals[,2:length(retrievals)])
    num.retrievals <- length(retrievals[,1])
    retrieval.cue.list <- retrievals[,2:(num.features+1)]
        
    ## Now do each retrieval, incrementally updating the history matrix.
    complete.results = NULL
    
    for (rr in 1:num.retrievals) {
        
        cues = retrieval.cue.list[rr,]
        moment = as.integer(retrievals[rr, 1])
        
        result = retrieve(cue.names, cues, moment)
        
        complete.results$summary = append(complete.results$summary, list(result$summary))
        complete.results$latencies[[rr]] = result$winner.latency
        
        ## update the history with what just happened at this moment
        moments <<- c(moments, moment)
        history <<- cbind(history, result$winner)          
    }
    return(complete.results)
}