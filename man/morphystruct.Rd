\name{morphystruct}
\alias{morphystruct}
\alias{PlotPopulationsMeristic}
\alias{CompareAllAssignmentsMeristic}
\alias{OptimizeLogLikelihoodMeristicGivenAssignment}
\alias{AllValuesPresent}
\alias{OptimizeLogLikelihoodMeristicSinglePop}
\alias{ReturnLogLikelihoodMeristicSinglePop}
\alias{NormalToMultinomial}
\title{Structure based on morphology}
\description{Using morphology (at first, just a single meristic character, though this will be extended), to lump groups. Initially, individuals must be passed in as groups (in the same way bpp needs initial assignments) though this will be changed eventually.}
\usage{
CompareAllAssignmentsMeristic(populations.list)
}       
\arguments{
\item{populations.list}{a list, where each element is a vector of observations of a meristic character (i.e., scale number) from a population}
}
\details{
Initially, with meristic traits, it takes initial populations and tries exhaustively all groupings, calculating the likelihood and AIC for each grouping. The model used assumes that within a population you expect a normal distribution of trait values, which is discretized to meristic data. For all the observations assigned to one population, it attempts to optimize the mean and sd for the underlying continuous trait leading to the observed discrete data. The likelihood is usually best by maximally splitting, but that is not always the best model under AIC. Model weights are calculated. It's possible to do multiple traits by running CompareAllAssignmentsMeristic() for each trait and combining the AIC for the corresponding assignments.
}
\value{
\code{CompareAllAssignmentsMeristic} returns a data.frame with the log likelihood (not negative log likelihood), AIC, delta AIC, AIC weight, and assignments. For example, an assigment of 1 2 1 3 means populations 1 and 3 are assigned to one group (this is shown by 1 * 1 * in the assignments vector), population 2 (* 2 * *) is assigned to another group, and population 4 (* * * 3) is assigned to a third group.
}
\examples{
#Imagine our empirical example are scale counts from three different rivers, with 20 samples from the first river, 15 from the second, and 12 from the third:
fishies <- list(round(rnorm(20, mean=40, sd=4)), round(rnorm(15, mean=60, sd=6)), round(runif(12, min=42, max=47)))
#Look at the distribution of traits:
PlotPopulationsMeristic(fishies)

#Now look at all possible assignments:
results <- CompareAllAssignmentsMeristic(fishies)
print(results)

#See which one has the best AIC and the weight of that model
}
\author{Brian C. O'Meara}
\keywords{plot}
