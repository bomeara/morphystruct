NormalToMultinomial <- function(mean.continuous, sd.continuous, min.discrete, max.discrete) {
	values <- seq(from=min.discrete, to=max.discrete, by=1)
	results <- rep(0, length(values))
	for (i in sequence(length(values))) {
		results[i] <- pnorm(values[i]+.5, mean=mean.continuous, sd=sd.continuous) - pnorm(values[i]-.5, mean=mean.continuous, sd=sd.continuous)
	}
	results<-results/sum(results) #make sum to one
	names(results) <- as.character(values)
	return(results)
}

ReturnLogLikelihoodMeristicSinglePop <- function(p, observations, min.discrete, max.discrete, return.neglnL=FALSE) {
	mean.continuous <- p[1]
	sd.continuous <- exp(p[2]) #note that sd is optimized in log space b/c it must be nonnegative
	prob <- NormalToMultinomial(mean.continuous, sd.continuous, min.discrete, max.discrete)
	observation.counts <- tabulate(bin=observations, nbins=max.discrete)[min.discrete:max.discrete]
	lnL <- dmultinom(observation.counts, prob=prob, log=TRUE)
	ifelse(return.neglnL, return(-lnL), return(lnL))
}

OptimizeLogLikelihoodMeristicSinglePop  <- function(observations, min.discrete=min(observations), max.discrete=max(observations)) {
	p <- c(mean(observations), log(sd(observations))) #starting guess
	results <- optim(par=p, fn= ReturnLogLikelihoodMeristicSinglePop, observations=observations, min.discrete=min.discrete, max.discrete=max.discrete, return.neglnL=TRUE) #return.neglnL so we minimize
	return(list(mean.continuous = results$par[1], sd.continuous = exp(results$par[2]), lnL=-results$value))
}

AllValuesPresent <- function(x) {
	ifelse(length(unique(x))==max(x), return(TRUE), return(FALSE))	
}

OptimizeLogLikelihoodMeristicGivenAssignment <- function(assignment, populations.list, min.discrete=min(unlist(populations.list)), max.discrete=max(unlist(populations.list))) {
	lnL <- 0
	for (assignment.index in sequence(max(assignment))) {
		lnL <- lnL + 	OptimizeLogLikelihoodMeristicSinglePop(observations=unlist(populations.list[which(assignment==assignment.index)]), min.discrete=min.discrete, max.discrete=max.discrete)$lnL
	}
	return(lnL)
}

CompareAllAssignmentsMeristic <- function(populations.list) {
	min.discrete <- min(unlist(populations.list))
	max.discrete <- max(unlist(populations.list))
	possibilities <- lapply(sequence(length(populations.list)), sequence)	
	assignments <- expand.grid(possibilities)
	assignments <- assignments[apply(assignments, 1, AllValuesPresent),]
	all.likelihoods <- apply(assignments, 1, OptimizeLogLikelihoodMeristicGivenAssignment, populations.list=populations.list, min.discrete=min.discrete, max.discrete=max.discrete)
	all.AIC <- 2*apply(assignments, 1, max) - 2 * all.likelihoods
	delta.AIC <- all.AIC - min(all.AIC)
	relative.likelihoods <- exp(-0.5 * delta.AIC)
	Akaike.weight <- relative.likelihoods / sum(relative.likelihoods)
	result <- data.frame(lnL=all.likelihoods, AIC=all.AIC, delta.AIC=delta.AIC, Akaike.weight=Akaike.weight)
	assignments.df <- data.frame(assignments)
	colnames(assignments.df) <- paste("Pop", sequence(length(populations.list)), sep=".")
	result <- cbind(result, assignments.df)
	
	return(result)
}

PlotPopulationsMeristic <- function(populations.list) {
	plot(x=range(unlist(populations.list)), y=c(0,1+length(populations.list)), xlab="Trait Value", yaxt="n", type="n", ylab="Population", bty="n")
	axis(side=2, at=sequence(length(populations.list))) 
	for (i in sequence(length(populations.list))) {
		counts <- table(populations.list[[i]])
		for (j in sequence(length(counts))) {
			text(as.numeric(names(counts)[j]), i, counts[j])	
		}
	}
}