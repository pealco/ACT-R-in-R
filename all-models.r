########################################################################################################
###
###  Simple Memory Retrieval Model of Agreement Phenomena (based on ACT-R theory)
###    With HTML trace output and PDF activation plots
###        
###
###  Rick Lewis & William Badecker (rickl@umich.edu)
###  Version 3.0
###  2 Feb 2007
###
###
########################################################################################################

library(reshape)
library(gplots)

source("run-model.r")

history = NULL

output = function(filename) {
    output_dir = "./output/"
    return(paste(output_dir, filename, sep=""))
}

## Experiment definitions

bp_embed = list(
    name = "BPNS_embed",
    description ="Brazilian Portuguese Null Subjects",
    conditions = list(
        list(
            condition = "GOOD",
            num.experimental.subjects = 45,
            num.experimental.items = 18,
            retrievals = "NS_with_embed/bp-good-retrievals.txt",
            items      = "NS_with_embed/bp-good-items.txt",
            data = 0.011,
            measure="percent error",
            correct.item = 1,
            distractor.item = 2,
            critical.retrieval = 2),   # second retrieval is critical

        list(
            condition = "BAD",
            num.experimental.subjects = 45,
            num.experimental.items = 18,                       
            retrievals = "NS_with_embed/bp-bad-retrievals.txt",
            items      = "NS_with_embed/bp-bad-items.txt",
            data = 0.012,
            measure="percent error",
            correct.item = 1,
            distractor.item = 2,
            critical.retrieval = 2),   

        list(
            condition = "INTERFERER",
            num.experimental.subjects = 45,
            num.experimental.items = 18,
            retrievals = "NS_with_embed/bp-interferer-retrievals.txt",
            items      = "NS_with_embed/bp-interferer-items.txt",
            data = 0.013,
            measure="percent error",
            correct.item = 1,
            distractor.item = 2,
            critical.retrieval = 2)
        )
    )

bp_no_embed = list(
    name = "BPNS_no_embed",
    description ="Brazilian Portuguese Null Subjects",
    conditions = list(
        list(
            condition = "GOOD",
            num.experimental.subjects = 45,
            num.experimental.items = 18,
            retrievals = "NS_without_embed/bp-good-retrievals.txt",
            items      = "NS_without_embed/bp-good-items.txt",
            data = 0.011,
            measure="percent error",
            correct.item = 1,
            distractor.item = 2,
            critical.retrieval = 2),   # second retrieval is critical

        list(
            condition = "BAD",
            num.experimental.subjects = 45,
            num.experimental.items = 18,                       
            retrievals = "NS_without_embed/bp-bad-retrievals.txt",
            items      = "NS_without_embed/bp-bad-items.txt",
            data = 0.012,
            measure="percent error",
            correct.item = 1,
            distractor.item = 2,
            critical.retrieval = 2),   

        list(
            condition = "INTERFERER",
            num.experimental.subjects = 45,
            num.experimental.items = 18,
            retrievals = "NS_without_embed/bp-interferer-retrievals.txt",
            items      = "NS_without_embed/bp-interferer-items.txt",
            data = 0.013,
            measure="percent error",
            correct.item = 1,
            distractor.item = 2,
            critical.retrieval = 2)
        )
    )
# Entries for data cannot all be equal.


## Complete list of experiments
experiments = list(bp_embed, bp_no_embed)
num.experiments = length(experiments)


## number of monte carlo trials per experiment
default.trials = 5000 #1000

## Discrete parameter space to search.  Each parameter now contains a list
## of possible values rather than a single value.  

## Latency factor
#F = c(1)
F = c(0.85)


## Extra category penalty
cat.penalty = c(-999)


## Total source activation
#G = c(1.0)
G = c(1.0)

## Activation noise parameter for logistic distribution
#ans  = c(0.15)
ans  = c(0.45, 0.6, 0.1)

## Fan parameter
#mas = c(1.5)
mas = c(1.5, 2.0)

## Base level decay parameter
d = c(0.5)

## Match penalty
#match.penalty = c(-1, -0.7, -0.4)
match.penalty = seq(-1, -0.4, 0.1)
## match.penalty = c(0, -0.1)

## Following are non-ACT-R modifications that we turn off by default
var.mismatch.penalty = c(FALSE)
VAR.fan = c(0)    # additional fan associated with VAR retrieval cues
modulate.by.distinct = c(FALSE)
distinctiveness = c(0)                  # Note that the quantitative parameter has no
                                        # effect when modulate.by.distinct is FALSE. So this is some  wasted
                                        # effort in the simple code below.

## Generate matrix defining combinatoric space of parameter values.  Each
## column is a parameter each row is a distinct comabination of parameters
## defining a model simulation to run.

parameters = list(cat.penalty, F, G, ans, mas, d, match.penalty, VAR.fan, var.mismatch.penalty, modulate.by.distinct, distinctiveness)
num.parameters = length(parameters)


## The total number of combinations is the product of the number of values
## for each parameter
num.combinations = prod(as.vector(lapply(parameters, length), mode="integer"))


## Set up matrix of parameter combinations.  Rows are model experiments,
## columns are parameters.
num.params = length(parameters)
p.matrix = matrix(nrow=num.combinations, ncol=num.params)

cumulative.num.combs = 1
for (p in 1:num.params) {
    p.matrix[,p] = rep(parameters[p][[1]], each=cumulative.num.combs, length.out=num.combinations)
    cumulative.num.combs = cumulative.num.combs * length(parameters[p][[1]])
}

##  Now set up matrix of unique model runs.
count.conditions = function(e) {
    length(e$conditions)
}

total.conditions = sum(as.vector(lapply(experiments, count.conditions), mode="integer"))
model.runs = data.frame()

for (e in 1:num.experiments) {
    exp.name = experiments[[e]]$name
    for (c in 1:length(experiments[[e]]$conditions)) {
        cond = experiments[[e]]$conditions[[c]]
        
        # placeholder for model result
        model.runs = rbind(
            model.runs,
            data.frame(
                experiment = rep(exp.name,num.combinations),
                condition = rep(cond$condition, num.combinations),
                retrievals = rep(cond$retrievals, num.combinations),
                items = rep(cond$items, num.combinations),
                num.experimental.items = rep(cond$num.experimental.items, num.combinations),
                num.experimental.subjects = rep(cond$num.experimental.subjects, num.combinations),
                measure = rep(cond$measure, num.combinations),
                critical.retrieval = rep(cond$critical.retrieval, num.combinations),
                correct.item = rep(cond$correct.item, num.combinations),
                distractor.item = rep(cond$distractor.item, num.combinations),
                data = rep(cond$data, num.combinations),
#               data.lower = rep(cond$data - cond$data.se, num.combinations),
#               data.upper = rep(cond$data + cond$data.se, num.combinations),
                model = rep(NA, num.combinations),
                model.lower = rep(NA, num.combinations),
                model.upper = rep(NA, num.combinations)
            )
        )
    }
}



# Duplicate the parameter matrix "total.conditions" number of times.
total.runs = total.conditions * num.combinations

full.parameter.matrix = matrix(data=t(p.matrix), nrow=total.runs,
                                ncol=num.params, byrow=TRUE)
colnames(full.parameter.matrix) = c("cat.penalty", "F", "G", "ans", "mas", "d", "match.penalty", "VAR.fan",
                                     "var.mismatch.penalty", "modulate.by.distinct", "distinctiveness")


## Finally, form the complete model run matrix.
all.runs = as.data.frame(cbind(full.parameter.matrix, model.runs))

## Loop over all runs and run the models
start = Sys.time()
for (r in 1:total.runs) {
    print(paste("Executing run #",r,"of",total.runs))
    
    ## select out row corresponding to this run
    this.run = all.runs[r,]
    
    ## set output file name
    filename.prefix = paste(this.run$experiment, "-", this.run$condition, "-", r, sep="")
    output.file = output(paste(filename.prefix, "-output.html", sep=""))
    
    ## now set the model parameters according to this combination of values
    set.parameters(this.run[1:num.parameters])
    
    ## and run the model
    item.file = as.character(this.run$items)
    retrieval.file = as.character(this.run$retrievals)
    num.experimental.items = this.run$num.experimental.items
    num.experimental.subjects = this.run$num.experimental.subjects
    
    results = run.model.quietly()
    #results = run.model(quiet=FALSE)
    
    ## plot the activation profiles for the critical and distractor items
    clrs = c("black", "green","blue","orange", "brown", "red", "purple")
    
    

    pdf(file = output(paste(filename.prefix, "-activation-plots.pdf", sep="")), width=11, height=5)
    plot.activation(moments, history, this.run$correct.item, this.run$distractor.item, this.run$experiment, this.run$condition)
    dev.off()

    ## now extract the relevant measure
    
    if (this.run$measure=="percent error") {
        crit.ret = results[[this.run$critical.retrieval]]
        model.result = crit.ret$retrieval.prob[this.run$distractor.item] * 100
        model.result.lower = crit.ret$retrieval.prob.lower[this.run$distractor.item] * 100
        model.result.upper = crit.ret$retrieval.prob.upper[this.run$distractor.item] * 100        
    }
    else {
        model.result = NA
        model.result.lower = NA
        model.result.upper = NA        
        print(paste("The", this.run$measure, "measure is not yet implemented."))
    }
    all.runs[r,]$model = model.result
    all.runs[r,]$model.lower = model.result.lower
    all.runs[r,]$model.upper = model.result.upper
    
    
    print(Sys.time() - start)
    start = Sys.time()
}

## Compute MSE and R^2 for each unique combination of parameter settings

param.results = data.frame(experiment = rep(NA,num.combinations*num.experiments),
                            combo     = rep(NA,num.combinations*num.experiments),
                            r2        = rep(NA,num.combinations*num.experiments),
                            mse       = rep(NA,num.combinations*num.experiments),
                            smse      = rep(NA,num.combinations*num.experiments)
                            )

print("Computing aggregate model fits for each unique parameter setting....")

i = 1
for (e in 1:num.experiments) {
    exp = experiments[[e]]$name
  
    for (j in 1:num.combinations) {
        print(paste("    parameter setting #",j,"of",num.combinations))
        index = seq(from=j, to=total.runs, by=num.combinations)
        runs = all.runs[index,]
        runs.exp = runs[runs$experiment==exp,]
        param.results$combo[i] = j
        param.results$experiment[i] = exp
        
        d = runs.exp$data
        m = runs.exp$model
        
        param.results$r2[i] = cor(d, m)^2
        param.results$spearman[i] = cor(d, m, method="spearman")
        param.results$mse[i] = mean((d - m)^2)
        param.results$smse[i] = sqrt(param.results$mse[i])
        i = i + 1
    }
}

aggregate.parameter.matrix = matrix(data=t(p.matrix), nrow=num.combinations*num.experiments,
                                ncol=num.params, byrow=TRUE)
colnames(aggregate.parameter.matrix) = c("cat.penalty", "F", "G", "ans", "mas", "d", "match.penalty", "VAR.fan",
                                          "var.mismatch.penalty", "modulate.by.distinct", "distinctiveness")

param.results = cbind(aggregate.parameter.matrix, param.results)

r.melt = melt(param.results, measure.var=c("r2","mse","smse","spearman"), variable="variable")

r2.summary       = cast(r.melt, combo ~ ., mean, subset = (variable == "r2"))
mse.summary      = cast(r.melt, combo ~ ., mean, subset = (variable == "mse"))
spearman.summary = cast(r.melt, combo ~ ., mean, subset = (variable == "spearman"))

combo.summary = cbind(p.matrix, r2.summary[2], mse.summary[2], spearman.summary[2])

colnames(combo.summary) = c("cat.penalty", "F", "G", "ans", "mas", "d", "match.penalty", "VAR.fan",
                             "var.mismatch.penalty", "modulate.by.distinct", "distinctiveness","r2","mse","spearman")

combo.summary = as.data.frame(combo.summary)

combo.summary$combo = seq(1:length(combo.summary[,1]))

plot.experiment = function(exp, combo) {
    index = seq(from=combo, to=total.runs, by=num.combinations)
    runs = all.runs[index,]
    runs.exp = runs[runs$experiment==exp,]
    ht = as.matrix(runs.exp$model)

    ht.lower =as.matrix(runs.exp$model.lower)
    ht.upper =as.matrix(runs.exp$model.upper)

    param.names = colnames(all.runs)[1:num.parameters]
    param.values = runs[1,1:num.parameters]
  
    barplot2(
        height=ht,
        bty="n",
        ci.l = ht.lower,
        ci.u = ht.upper,
        ci.color="gray40",
        plot.ci=TRUE,
        col=rev(grey.colors(length(runs.exp$condition))),
        ##  legend.text=runs.exp$condition,
        names=c("Retrieval Interference Model"),
        beside=TRUE,
        space = c(0.1, 1),
        ylim = c(0,50),
        axes=TRUE,   
        xlab="", ylab=as.character(runs.exp$measure[1]),
        cex.lab=1.0,
        main=get.description(exp),
        sub=sprintf("noise=%f, mas=%f, mp=%f, d=%f", param.values$ans,param.values$mas,
        param.values$match.penalty,param.values$d)
        )
  
  axis(1,labels=FALSE,tcl=0)  ## zero tick length
  legend(x="topleft",
         legend=runs.exp$condition,
         fill=rev(grey.colors(length(runs.exp$condition))),
         bty="n")

  
  return(ht)
  
# text(x=2,y=seq(from=38,to=(38-(num.parameters*2)+2),by=-2),labels=param.names,cex=1.0)
# text(x=4,y=seq(from=38,to=(38-(num.parameters*2)+2),by=-2),labels=as.character(param.values),cex=1.0)
}



get.description = function(exp) {
    for (i in 1:num.experiments) {
        if (experiments[[i]]$name == exp) {
            return(experiments[[i]]$description)
        }
    }
    return("No description")
}



plot.best.overall.no.decay = function() {
    subs = combo.summary[combo.summary$d==0.001,]
    c = subs$combo[which.min(subs$mse)]
  
    pdf(file=output("best-overall-no-decay.pdf"), height=11, width=8.5)
#   par(mfrow=c(3,2))
    par(pin=c(6.5,3.5))

    model.points = c()
    data.points = c()

    for (e in 1:num.experiments) {
        exp = experiments[[e]]$name
        this.exp = param.results[(param.results$experiment==exp) & (param.results$d==0.001),]
        plotted.points = plot.experiment(exp, c)
        data.points = c(data.points, plotted.points[,1])
        model.points = c(model.points, plotted.points[,2])
        par(pin=c(3.5,3.5))    
    }

}



plot.best.overall.no.decay.no.mp = function() {
    subs = combo.summary[combo.summary$d==0.001 & combo.summary$match.penalty==0,]
    c = subs$combo[which.min(subs$mse)]
  
    pdf(file=output("best-overall-no-decay-no-mp.pdf"),height=11,width=8.5)
#   par(mfrow=c(3,2))
    par(pin=c(6.5,3.5))

    model.points = c()
    data.points = c()

    for (e in 1:num.experiments) {
        exp = experiments[[e]]$name
        this.exp = param.results[(param.results$experiment==exp) & (param.results$d==0.001),]
        plotted.points = plot.experiment(exp, c)
        data.points = c(data.points, plotted.points[,1])
        model.points = c(model.points, plotted.points[,2])
        par(pin=c(3.5,3.5))    
    }
 
    plot(model.points ~ data.points,
        xlim=c(0,30),ylim=c(0,30),
        main="Model vs. Data Across 14 Conditions",
        xlab ="Data", ylab="Model")
    lines(x=c(0,30), y=c(0,30),lty=3)
    dev.off()
}

plot.individual.no.decay = function() {
    pdf(file=output("best-individual-no-decay.pdf"),height=11,width=8.5)
#   par(mfrow=c(3,2))
#   par(pin=c(6.5,3.5))   # this is a hack to make the Slovak graph wider
    model.points = c()
    data.points = c()
  
    for (e in 1:num.experiments) {
        exp = experiments[[e]]$name
        this.exp = param.results[(param.results$experiment==exp) & (param.results$d==0.001),]
        m = which.min(this.exp$mse)
        c = this.exp$combo[m]
        plotted.points = plot.experiment(exp, c)
        data.points = c(data.points, plotted.points[,1])
        model.points = c(model.points, plotted.points[,2])
        par(pin=c(3.5,3.5))        
    }

    plot(model.points ~ data.points,
        main="Model vs. Data Across 14 Conditions", 
        xlim=c(0,30),ylim=c(0,30), 
        xlab ="Data", ylab="Model")
    lines(x=c(0,30), y=c(0,30),lty=3)
    
    dev.off()
}


  
plot.individual.no.decay.no.mp = function() {
    pdf(file=output("best-individual-no-decay-no-mp.pdf"),height=11,width=8.5)
#   par(mfrow=c(3,2))
    par(pin=c(6.5,3.5))   # this is a hack to make the Slovak graph wider
    model.points = c()
    data.points = c()
  
    for (e in 1:num.experiments) {
        exp = experiments[[e]]$name
        this.exp = param.results[(param.results$experiment==exp) & (param.results$d==0.001) & (param.results$match.penalty==0),]
        m = which.min(this.exp$mse)
        c = this.exp$combo[m]
        plotted.points = plot.experiment(exp, c)
        data.points = c(data.points, plotted.points[,1])
        model.points = c(model.points, plotted.points[,2])
        par(pin=c(3.5,3.5))        
    }

    plot(model.points ~ data.points,
        main="Model vs. Data Across 14 Conditions",
        xlim=c(0,30),ylim=c(0,30),
        xlab ="Data", ylab="Model")
    lines(x=c(0,30), y=c(0,30),lty=3)
    dev.off()
}

plot.best.overall.decay = function() {
    subs = combo.summary[combo.summary$d==0.5,]
    c = subs$combo[which.min(subs$mse)]
  
    pdf(file=output("best-overall-decay.pdf"),height=11,width=8.5)
#   par(mfrow=c(3,2))
    par(pin=c(6.5,3.5))

    model.points = c()
    data.points = c()

    for (e in 1:num.experiments) {
        exp = experiments[[e]]$name
        this.exp = param.results[(param.results$experiment==exp) & (param.results$d==0.5),]
        plotted.points = plot.experiment(exp, c)
        model.points = c(model.points, plotted.points)
        par(pin=c(3.5,3.5))    
    }
 
    #plot(model.points ~ data.points,
    #    main="Model vs. Data Across 14 Conditions",
    #    xlim=c(0,30),ylim=c(0,30),
    #    xlab ="Data", ylab="Model")
    #lines(x=c(0,30), y=c(0,30),lty=3)
    
    dev.off()
}

plot.best.overall.decay.no.mp = function() {
    subs = combo.summary[combo.summary$d==0.5 & combo.summary$match.penalty==0,]
    c = subs$combo[which.min(subs$mse)]
  
    pdf(file=output("best-overall-decay-no-mp.pdf"),height=11,width=8.5)
#   par(mfrow=c(3,2))
    par(pin=c(6.5,3.5))

    model.points = c()
    data.points = c()

    for (e in 1:num.experiments) {
        exp = experiments[[e]]$name
        this.exp = param.results[(param.results$experiment==exp) & (param.results$d==0.5),]
        plotted.points = plot.experiment(exp, c)
        data.points = c(data.points, plotted.points[,1])
        model.points = c(model.points, plotted.points[,2])
        par(pin=c(3.5,3.5))    
    }
 
    #plot(model.points ~ data.points,
    #    xlim=c(0,30),ylim=c(0,30),
    #    main="Model vs. Data Across 14 Conditions",
    #    xlab ="Data", ylab="Model")
    #lines(x=c(0,30), y=c(0,30),lty=3)
    
    dev.off()
}

plot.individual.decay = function() {
    pdf(file=output("best-individual-decay.pdf"),height=11,width=8.5)
#   par(mfrow=c(3,2))
    par(pin=c(6.5,3.5))   # this is a hack to make the Slovak graph wider
    model.points = c()
    data.points = c()
  
    for (e in 1:num.experiments) {
        exp = experiments[[e]]$name
        this.exp = param.results[(param.results$experiment==exp) & (param.results$d==0.5),]
        m = which.min(this.exp$mse)
        c = this.exp$combo[m]
        plotted.points = plot.experiment(exp, c)
        model.points = c(model.points, plotted.points)
        par(pin=c(3.5,3.5))        
    }

    #plot(model.points ~ data.points,
    #    xlim=c(0,30),ylim=c(0,30),
    #    main="Model vs. Data Across 14 Conditions",
    #    xlab ="Data", ylab="Model")
    #lines(x=c(0,30), y=c(0,30),lty=3)
    
    dev.off()
}
  
plot.individual.decay.no.mp = function() {
    pdf(file=output("best-individual-decay-no-mp.pdf"),height=11,width=8.5)
#   par(mfrow=c(3,2))
    par(pin=c(6.5,3.5))   # this is a hack to make the Slovak graph wider
    model.points = c()
    data.points = c()
  
    for (e in 1:num.experiments) {
        exp = experiments[[e]]$name
        this.exp = param.results[(param.results$experiment==exp) & (param.results$d==0.5) & (param.results$match.penalty==0),]
        m = which.min(this.exp$mse)
        c = this.exp$combo[m]
        plotted.points = plot.experiment(exp, c)
        data.points = c(data.points, plotted.points[,1])
        model.points = c(model.points, plotted.points[,2])
        par(pin=c(3.5,3.5))        
    }

    #plot(model.points ~ data.points,
    #    xlim=c(0,30),ylim=c(0,30),
    #    main="Model vs. Data Across 14 Conditions",
    #    xlab ="Data", ylab="Model")
    #lines(x=c(0,30), y=c(0,30),lty=3)
    
    dev.off()
}

plot.experiment.all.combos = function(exp,decay=TRUE) {
    ## set up plot
    if (decay) {
        runs.exp = all.runs[(all.runs$experiment==exp & all.runs$d==0.5),]
    } else {
        runs.exp = all.runs[(all.runs$experiment==exp & all.runs$d==0.001),]
    }
      
    num.runs = length(runs.exp[,1])
    num.conds = length(unique(runs.exp$condition))
    ## only works if even (e.g. both decay and not decay)
    num.combs = num.runs/num.conds

    dummy = data.frame(condition = 1:num.conds,
                      m = rep(0,num.conds))
    plot(m ~ condition,
        data = dummy,
        type="n",
        xlab="",
        ylab="",
        main = get.description(exp),
        ylim=c(0,40),
        xaxt="n"
        )

    for (combo in 1:num.combs) {
        index = seq(from=combo, to=num.runs,by=num.combs)
        r = runs.exp[index,]
        lines(x=1:num.conds, y=r$model,type="l",lty=1,col="grey50",pch=c(20))
    }

    lines(x=1:num.conds, y=r$data,type="o",lty=1,lwd=2,col="black",pch=c(20))
  
    axis(1,labels=FALSE,tcl=0)  ## zero tick length
}

plot.full.range.decay = function() {
    pdf(file=output("full-range-decay.pdf"),height=10,width=10)
#   par(mfrow=c(3,3))
    par(bty="n",xpd=NA)
    par(mgp=c(2.8,1,0))
    par(cex.sub=1.1,font.sub=3)
    par(cex.main=0.9,font.sub=3)  
    par(cex.axis=1.05,cex.lab=1.1)
    
    for (e in 1:num.experiments) {
        exp = experiments[[e]]$name
        plot.experiment.all.combos(exp,decay=TRUE)
    }
    dev.off()
}


plot.full.range.no.decay = function() {
    pdf(file=output("full-range-no-decay.pdf"),height=10,width=10)
#   par(mfrow=c(3,3))
    par(bty="n",xpd=NA)
    par(mgp=c(2.8,1,0))
    par(cex.sub=1.1,font.sub=3)
    par(cex.main=0.8,font.sub=3)  
    par(cex.axis=1.05,cex.lab=1.1)
    
    for (e in 1:num.experiments) {
        exp = experiments[[e]]$name
        plot.experiment.all.combos(exp,decay=FALSE)
    }
    dev.off()
}

## plot.best.overall.no.decay()
## plot.individual.no.decay()
## plot.best.overall.no.decay.no.mp()
## plot.individual.no.decay.no.mp()

plot.best.overall.decay()
plot.individual.decay()
#plot.best.overall.decay.no.mp()
#plot.individual.decay.no.mp()

##plot.full.range.no.decay()
plot.full.range.decay()
