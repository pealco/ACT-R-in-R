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
###  Pedro Alcocer
###  Modifications: 8 Mar 2011
###
#########################################################################################################

library(xtable)
library(doMC)

registerDoMC()

## act-r library
source("actr.r")


## SET PARAMETERS function: sets the parameters as global variables

set.parameters = function(p) {
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
    modulate.by.distinct <<-   p$modulate.by.distinct
    distinctiveness <<-   p$distinctiveness
}



## RUN THE MODEL given the global parameter settings

run.model = function(quiet=TRUE) {
    if (!(is.na(num.experimental.items) | is.na(num.experimental.subjects))) {
        trials <<- simulated.experiments * num.experimental.items * num.experimental.subjects
    } else {
        trials <<- default.trials
    }
  
    #print("Starting to run model....")
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
    
    if (!quiet) {
    
        ## Print out header
        write(paste("<HR SIZE=8 NOSHADE><h2>Item file: ",item.file,"</h2>"),file=output.file)
        write(paste("<HR SIZE=8 NOSHADE><h2>Retrieval file: ",retrieval.file,"</h2>"),file=output.file,append=TRUE)
        
        ## Print out the set of items
        write("<HR SIZE=8 NOSHADE><h2>Memory items (feature bundles) and their creation times</h2>",file=output.file,
            append=TRUE)
        
        fb = items[2:(num.features+2),2:num.columns]
        colnames(fb) = item.name
        rownames(fb) = items[2:(num.features+2),1]
        ft = xtable(fb)
        print(file=output.file, ft, append=TRUE,type="html")
        
        write("<HR SIZE=8 NOSHADE><h2>Retrieval History</h2>",file=output.file,append=TRUE) 
    }
    
    ## create initial history matrix to be the creation moments
    history <<- matrix(1:num.items,nrow=trials, ncol=num.items, byrow=TRUE)
    moments <<- creation.moment
    
    ## Read in the schedule of retrievals
    retrievals = read.delim(file=retrieval.file,header=FALSE)
    cue.names = retrievals[2:(num.features+1),1]
    retrievals = t(retrievals[,2:length(retrievals)])
    num.retrievals = length(retrievals[,1])
    retrieval.cue.list = retrievals[,2:(num.features+1)]
    
    ## Now do each retrieval, incrementally updating the history matrix.
    complete.results = NULL
    
    for (r in 1:num.retrievals) {
        cues = retrieval.cue.list[r,]
        moment = as.integer(retrievals[r,1])
    
        if (!quiet) {
            print("",quote=FALSE)
            print ("=========================================================================================",quote=FALSE)
            print(paste("Retrieving at",moment,"ms, with cues:"),quote=FALSE)
            print ("=========================================================================================",quote=FALSE)
            
            write(file=output.file,paste("<h3>Retrieval at ",moment,"ms</h3>",sep=""),append=TRUE) 
            
            
            print.cues = as.matrix(cues)
            colnames(print.cues) = c("value")
            rownames(print.cues) = items[3:(num.features+2),1]
            print(print.cues,quote=FALSE)
            cues.tab = xtable(print.cues,caption = paste("Cues at ",moment,"ms",sep=""))
            
            write(file=output.file,"<TABLE border=0 cellspacing=30><TD>",append=TRUE)
            print(file=output.file, cues.tab, caption.placement="top",append=TRUE,type="html")
            write(file=output.file,"</TD><TD>",append=TRUE)
        }
      
        result = retrieve(cue.names, cues, moment)
        summary = result$summary
        complete.results = append(complete.results, list(summary))
          
        if (!quiet) {
            print(summary,quote=FALSE)
            
            rtable = xtable(summary,caption=paste("Result at ",moment,"ms",sep=""))
            print(file=output.file, rtable, caption.placement="top",append=TRUE,type="html")
            write(file=output.file,"</TD></TABLE><HR size=5 noshade>",append=TRUE)
            
            params = rbind(c("Latency factor", "Total source activation", "Activation noise", "Fan parameter", "Base level decay"), 
                            c(F, G, ans, mas, d))
            params.table = xtable(params)
            print(file=output.file, params.table, append=TRUE, type="html")
        }      
        ## update the history with what just happened at this moment
        moments <<- c(moments, moment)
        history <<- cbind(history, result$winner)
          
        if (!quiet) {
            for (c in 1:num.items) {
                hist(main = paste("Retrieval time distribution for item",item.name[c]),
                ##           result$final.latency[c, result$final.latency[c,]<1200])
                result$final.latency)
            }
        }
    }
    
    if (!quiet) {
        dev.off()
    }
    return(complete.results)
}


run.model.quietly = function() {
  run.model(quiet=TRUE)
}



plot.activation = function(moments, history, correct.item, distractor, experiment, condition) {
    min.time = 0
    max.time = moments[length(moments)] + 200

    ticks = seq(min.time,max.time,10)    
    base.activations = foreach(t = ticks, .combine="cbind") %dopar% {
        base.levels = compute.base.levels(t)
        exists = matrix(creation.moment <  t, ncol=trials, nrow=num.items)
        activation  = base.levels*exists + 0*!exists
        rowMeans(activation)
        #if (round(t/100)==t/100) {print(t)}
    }

    # Plot only correct item and distractor item.
    plotting.items = c(correct.item, distractor)
    clrs = c("black", "red")
    
    #print(paste("Plotting exp",experiment,"condition:",condition,"with items:",correct.item,distractor))
    #print(paste(" ... and",num.items,"total items"))
    
    matplot(ticks, t(base.activations[plotting.items,]), type="l", 
        main=paste("Mean activation of items,", trials,"trial", "Exp:", experiment, "Condition:", condition),
        ylab="Activation", xlab="Time", lwd=4, lty=1)
    
    maxb = max(base.activations, na.rm=TRUE)
    minb = min(base.activations, na.rm=TRUE)
    
    for (m in creation.moment) {
      lines(x=c(m, m), y=c(minb -0.1, minb-0.5), lend=2, lwd=2, col="darkgreen")
    }

    for (m in setdiff(moments, creation.moment)) {
      lines(x=c(m, m), y=c(minb -0.1, minb-0.5), lend=2, lwd=4, col="red")
    }

    legend("top", c("Head NP","Distractor NP"), 
        lty=1, lwd=4, bty="n", cex=1, col = clrs[1:2],)
    
    # Plot all items.
    clrs = c("black","purple","green","blue","orange", "yellow")
    clrs[correct.item] = "black"
    clrs[distractor] = "red"

    matplot(ticks, t(base.activations), type="l", 
        main=paste("Mean activation of items,", trials,"trial", "Exp:", experiment, "Condition:", condition),
        ylab="Activation", xlab="Time", lwd=4, lty=1, col=clrs)

    legend("top", item.name, 
        lty=1, lwd=4, bty="n", cex=1, col = clrs[1:num.items],)
}

#plot.activation.profiles = function(moments, history, min.time, max.time, increment=10,
#                                     creation.moments, item.names=item.name) {
#    time.span = seq(min.time, max.time, increment)
#    base.activations = matrix(nrow=num.items, ncol=length(time.span))
#
##    print("Computing complete history of activation values at times....")
#
#    #  First compute the history of activation values at each time point
#    j = 1
#    for (t in time.span) {
##    if (round(t/100)==t/100) {print(t)}
#      base.levels = compute.base.levels(t)
#
#    # make items that don't exist yet have activation of NA
#    exists = matrix(creation.moment <  t, ncol=trials, nrow=num.items)
#    activation  = base.levels*exists + 0*!exists
#        activation[activation==0] = NA
#
#    # take mean activation over all the monte carlo trials
#    base.activations[, j] = rowMeans(activation)
#        j = j + 1
#    }
#
#    maxb = max(base.activations, na.rm=TRUE)
#    minb = min(base.activations, na.rm=TRUE)
#
#    plot(base.activations[1,] ~ time.span,
#               type="l", lwd=1.5, col=clrs[1],
#               main=paste("Mean activation of items over time (", trials," runs)", sep=""),
#               sub="Green bars indicate initial encoding points, red bars indicate retrieval points",
#               ylab="Activation", xlab="Time",
#               ylim=c(minb-0.5, maxb+0.5))
#
#    lines(x=c(min(time.span), max(time.span)), y=c(0, 0), lty=3)
#
#    for (c in 1:num.items) {
#      lines(base.activations[c,] ~ time.span, type="l", lwd=1.5, col=clrs[c])
#    }
#
#    # add markers for the creation and retrieval moments at the bottom
#    for (m in creation.moment) {
#      lines(x=c(m, m), y=c(minb -0.1, minb-0.5), lend=2, lwd=2, col="darkgreen")
#    }
#
#    for (m in setdiff(moments, creation.moment)) {
#      lines(x=c(m, m), y=c(minb -0.1, minb-0.5), lend=2, lwd=4, col="red")
#    }
#
#    # add a legend
#    width = max.time - min.time
#    height = maxb
#    legend(0.8*width+min.time, height+0.5, item.name, lty=1, lwd=1.5, bty="n", col = clrs[1:num.items])
#  }
#