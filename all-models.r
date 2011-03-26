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
    description ="Brazilian Portuguese Null Subjects w/ embed as cue",
    conditions = list(
        list(
            condition = "GOOD",
            retrievals = "NS_with_embed/bp-good-retrievals.txt",
            items      = "NS_with_embed/bp-good-items.txt",
            measure="percent error",
            correct.item = 1,
            distractor.item = 2,
            critical.retrieval = 2),   # second retrieval is critical

        list(
            condition = "BAD",                 
            retrievals = "NS_with_embed/bp-bad-retrievals.txt",
            items      = "NS_with_embed/bp-bad-items.txt",
            measure="percent error",
            correct.item = 1,
            distractor.item = 2,
            critical.retrieval = 2),   

        list(
            condition = "INTERFERER",
            retrievals = "NS_with_embed/bp-interferer-retrievals.txt",
            items      = "NS_with_embed/bp-interferer-items.txt",
            measure="percent error",
            correct.item = 1,
            distractor.item = 2,
            critical.retrieval = 2)
        )
    )

bp_no_embed = list(
    name = "BPNS_no_embed",
    description ="Brazilian Portuguese Null Subjects w/o embed as cue",
    conditions = list(
        list(
            condition = "GOOD",
            retrievals = "NS_without_embed/bp-good-retrievals.txt",
            items      = "NS_without_embed/bp-good-items.txt",
            measure="percent error",
            correct.item = 1,
            distractor.item = 2,
            critical.retrieval = 2),   # second retrieval is critical

        list(
            condition = "BAD",                       
            retrievals = "NS_without_embed/bp-bad-retrievals.txt",
            items      = "NS_without_embed/bp-bad-items.txt",
            measure="percent error",
            correct.item = 1,
            distractor.item = 2,
            critical.retrieval = 2),   

        list(
            condition = "INTERFERER",
            retrievals = "NS_without_embed/bp-interferer-retrievals.txt",
            items      = "NS_without_embed/bp-interferer-items.txt",
            measure="percent error",
            correct.item = 1,
            distractor.item = 2,
            critical.retrieval = 2)
        )
    )
# Entries for data cannot all be equal.

agr_embed = list(
    name = "agr_embed",
    description ="Brazilian Portuguese Subject-Verb Agreement w/ embed as cue",
    conditions = list(
        list(
            condition = "GP",
            retrievals = "agr_with_embed/agr-GP-retrievals.txt",
            items      = "agr_with_embed/agr-GP-items.txt",
            measure="percent error",
            correct.item = 2,
            distractor.item = 1,
            critical.retrieval = 1),

        list(
            condition = "GS",
            retrievals = "agr_with_embed/agr-GS-retrievals.txt",
            items      = "agr_with_embed/agr-GS-items.txt",
            measure="percent error",
            correct.item = 2,
            distractor.item = 1,
            critical.retrieval = 1),

        list(
            condition = "UP",
            retrievals = "agr_with_embed/agr-UP-retrievals.txt",
            items      = "agr_with_embed/agr-UP-items.txt",
            measure="percent error",
            correct.item = 2,
            distractor.item = 1,
            critical.retrieval = 1),
            
        list(
            condition = "US",
            retrievals = "agr_with_embed/agr-US-retrievals.txt",
            items      = "agr_with_embed/agr-US-items.txt",
            measure="percent error",
            correct.item = 2,
            distractor.item = 1,
            critical.retrieval = 1)
        )
    )


agr_no_embed = list(
    name = "agr_no_embed",
    description ="Brazilian Portuguese Subject-Verb Agreement w/o embed as cue",
    conditions = list(
        list(
            condition = "GP",
            retrievals = "agr_without_embed/agr-GP-retrievals.txt",
            items      = "agr_without_embed/agr-GP-items.txt",
            measure="percent error",
            correct.item = 2,
            distractor.item = 1,
            critical.retrieval = 1),

        list(
            condition = "GS",
            retrievals = "agr_without_embed/agr-GS-retrievals.txt",
            items      = "agr_without_embed/agr-GS-items.txt",
            measure="percent error",
            correct.item = 2,
            distractor.item = 1,
            critical.retrieval = 1),

        list(
            condition = "UP",
            retrievals = "agr_without_embed/agr-UP-retrievals.txt",
            items      = "agr_without_embed/agr-UP-items.txt",
            measure="percent error",
            correct.item = 2,
            distractor.item = 1,
            critical.retrieval = 1),
            
        list(
            condition = "US",
            retrievals = "agr_without_embed/agr-US-retrievals.txt",
            items      = "agr_without_embed/agr-US-items.txt",
            measure="percent error",
            correct.item = 2,
            distractor.item = 1,
            critical.retrieval = 1)
        )
    )


## Complete list of experiments
experiments = list(bp_embed, bp_no_embed, agr_embed, agr_no_embed)
num.experiments = length(experiments)


## number of monte carlo trials per experiment
default.trials = 5000

## Discrete parameter space to search.  Each parameter now contains a list
## of possible values rather than a single value.  

## Latency factor
F = c(0.85)

## Extra category penalty
cat.penalty = c(-999);

## Total source activation
G = c(0.7, 0.8, 1.0);
#G = c(1.0);


## Activation noise parameter for logistic distribution
ans  = seq(0.45, 0.6, 0.05)
#ans  = c(0.45)

## Fan parameter
mas = c(1.5, 2.0)
#mas = c(1.5)


## Base level decay parameter
d = c(0.000001, 0.5)
#d = c(0.5)

## Match penalty
match.penalty = seq(-1.0, -0.4, 0.2)
#match.penalty = c(-1.0)


## Following are non-ACT-R modifications that we turn off by default
var.mismatch.penalty = c(FALSE);
VAR.fan = c(0);    # additional fan associated with VAR retrieval cues
modulate.by.distinct = c(FALSE);
distinctiveness = c(0);                # Note that the quantitative parameter has no
                                        # effect when modulate.by.distinct is FALSE. So this is some  wasted
                                        # effort in the simple code below.



## Generate matrix defining combinatoric space of parameter values.  Each
## column is a parameter; each row is a distinct comabination of parameters
## defining a model simulation to run.

parameters = list(cat.penalty, F, G, ans, mas, d, match.penalty, VAR.fan, var.mismatch.penalty,
                   modulate.by.distinct, distinctiveness);
num.parameters = length(parameters);


## The total number of combinations is the product of the number of values
## for each parameter
num.combinations = prod(as.vector(lapply(parameters, length), mode="integer"));


## Set up matrix of parameter combinations.  Rows are model experiments,
## columns are parameters.
num.params = length(parameters);
p.matrix = matrix(nrow=num.combinations, ncol=num.params);

cumulative.num.combs = 1;
for (p in 1:num.params) {
  p.matrix[,p] = rep(parameters[p][[1]], each=cumulative.num.combs, length.out=num.combinations);
  cumulative.num.combs = cumulative.num.combs * length(parameters[p][[1]]);
}



##  Now set up matrix of unique model runs.
count.conditions = function(e) {
  length(e$conditions);
}

total.conditions = sum(as.vector(lapply(experiments, count.conditions),
                                  mode="integer"));
model.runs = data.frame();

for (e in 1:num.experiments) {
  exp.name = experiments[[e]]$name;
  for (c in 1:length(experiments[[e]]$conditions)) {
    cond = experiments[[e]]$conditions[[c]];
    
    model.runs = rbind(model.runs,
        data.frame(experiment                = rep(exp.name,num.combinations),
                   condition                 = rep(cond$condition, num.combinations),
                   retrievals                = rep(cond$retrievals, num.combinations),
                   items                     = rep(cond$items, num.combinations),
                   measure                   = rep(cond$measure, num.combinations),
                   critical.retrieval        = rep(cond$critical.retrieval, num.combinations),
                   correct.item              = rep(cond$correct.item, num.combinations),
                   distractor.item           = rep(cond$distractor.item, num.combinations),
                   model                     = rep(NA, num.combinations)
                   )
                ) # placeholder for
                  # model result
  }
};



# Duplicate the parameter matrix "total.conditions" number of times.
total.runs = total.conditions * num.combinations;

full.parameter.matrix = matrix(data=t(p.matrix), nrow=total.runs,
                                ncol=num.params, byrow=TRUE);
colnames(full.parameter.matrix) = c("cat.penalty", "F", "G", "ans", "mas", "d", "match.penalty", "VAR.fan",
                                     "var.mismatch.penalty", "modulate.by.distinct", "distinctiveness");


## Finally, form the complete model run matrix.
all.runs = as.data.frame(cbind(full.parameter.matrix, model.runs));
final.samples = matrix(data=NA, nrow=default.trials, ncol=total.runs)

## Loop over all runs and run the models
for (r in 1:total.runs) {
    print(paste("Executing run #",r,"of",total.runs));
    
    ## select out row corresponding to this run
    this.run = all.runs[r,];      
    
    ## now set the model parameters according to this combination of values
    set.parameters(this.run[1:num.parameters]);
    
    ## and run the model
    item.file = as.character(this.run$items);
    retrieval.file = as.character(this.run$retrievals);
    
    results = run.model()
        
    ## now extract the relevant measures.
    
    # Percent error.
    crit.ret = results$summary[[this.run$critical.retrieval]]
    model.result = crit.ret$retrieval.prob[this.run$distractor.item] * 100
    
    # Latency.
    final.samples[,r] = results$latencies[[this.run$critical.retrieval]]
    
    all.runs[r,]$model = model.result
}

write.csv(all.runs, "all.runs.txt")
write.csv(final.samples, "final.samples.txt")