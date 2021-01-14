#  File R/zzz.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2020 Statnet Commons
#######################################################################
.onAttach <- function(lib, pkg){
  #' @importFrom statnet.common statnetStartupMessage
  sm <- statnetStartupMessage("ergm", c("statnet","ergm.count","tergm"), TRUE)
  if(!is.null(sm)){
    packageStartupMessage(sm)
    packageStartupMessage(paste(c(strwrap(paste("NOTE: Versions before 3.6.1 had a bug in the implementation of the bd() constraint which distorted the sampled distribution somewhat. In addition, Sampson's Monks datasets had mislabeled vertices. See the NEWS and the documentation for more details.",sep="")),""),collapse="\n"))
    packageStartupMessage(paste(c(strwrap(paste("NOTE: Some common term arguments pertaining to vertex attribute and level selection have changed in 3.10.0. See terms help for more details. Use ",sQuote("options(ergm.term=list(version=\"3.9.4\"))")," to use old behavior.",sep="")),""),collapse="\n"))
  }
}

.onLoad <- function(libname, pkgname){
  # . is used as a placeholder by stantet.common::NVL3().
  utils::globalVariables(".")

  default_options(ergm.eval.loglik=TRUE,
                  ergm.loglik.warn_dyads=TRUE,
                  ergm.cluster.retries=5)

  eval(COLLATE_ALL_MY_CONTROLS_EXPR)

  .RegisterProposals()
}

#' @name ergm-reexports
#' @title Reexports from other packages
#' @param ... Arguments to reexported functions.
NULL

#' @describeIn ergm-reexports See [statnet.common::sctrl].
#' @param drop,init,init.method,main.method,force.main,main.hessian,checkpoint,resume,MPLE.max.dyad.types,MPLE.samplesize,init.MPLE.samplesize,MPLE.type,MPLE.nonident,MPLE.nonident.tol,MCMC.prop.weights,MCMC.prop.args,MCMC.interval,MCMC.burnin,MCMC.samplesize,MCMC.effectiveSize,MCMC.effectiveSize.damp,MCMC.effectiveSize.maxruns,MCMC.effectiveSize.base,MCMC.effectiveSize.points,MCMC.effectiveSize.order,MCMC.return.stats,MCMC.runtime.traceplot,MCMC.init.maxedges,MCMC.max.maxedges,MCMC.addto.se,MCMC.compress,MCMC.packagenames,SAN.maxit,SAN.nsteps.times,SAN,MCMLE.termination,MCMLE.maxit,MCMLE.conv.min.pval,MCMLE.NR.maxit,MCMLE.NR.reltol,obs.MCMC.samplesize,obs.MCMC.interval,obs.MCMC.burnin,obs.MCMC.burnin.min,obs.MCMC.prop.weights,obs.MCMC.prop.args,obs.MCMC.impute.min_informative,obs.MCMC.impute.default_density,MCMLE.MCMC.precision,MCMLE.MCMC.max.ESS.frac,MCMLE.metric,MCMLE.method,MCMLE.trustregion,MCMLE.dampening,MCMLE.dampening.min.ess,MCMLE.dampening.level,MCMLE.steplength.margin,MCMLE.steplength,MCMLE.steplength.parallel,MCMLE.adaptive.trustregion,MCMLE.sequential,MCMLE.density.guard.min,MCMLE.density.guard,MCMLE.effectiveSize,MCMLE.last.boost,MCMLE.steplength.esteq,MCMLE.steplength.miss.sample,MCMLE.steplength.maxit,MCMLE.steplength.min,MCMLE.effectiveSize.interval_drop,MCMLE.save_intermediates,MCMLE.nonident,MCMLE.nonident.tol,SA.phase1_n,SA.initial_gain,SA.nsubphases,SA.niterations,SA.phase3_n,SA.trustregion,RM.phase1n_base,RM.phase2n_base,RM.phase2sub,RM.init_gain,RM.phase3n,Step.MCMC.samplesize,Step.maxit,Step.gridsize,CD.nsteps,CD.multiplicity,CD.nsteps.obs,CD.multiplicity.obs,CD.maxit,CD.conv.min.pval,CD.NR.maxit,CD.NR.reltol,CD.metric,CD.method,CD.trustregion,CD.dampening,CD.dampening.min.ess,CD.dampening.level,CD.steplength.margin,CD.steplength,CD.adaptive.trustregion,CD.adaptive.epsilon,CD.steplength.esteq,CD.steplength.miss.sample,CD.steplength.maxit,CD.steplength.min,CD.steplength.parallel,loglik,term.options,seed,parallel,parallel.type,parallel.version.check,nsteps,GF.init.maxedges.mul,nsim,network.output,warn.dyads,SAN.tau,SAN.invcov,SAN.invcov.diag,SAN.nsteps.alloc,SAN.nsteps,SAN.samplesize,SAN.init.maxedges,SAN.max.maxedges,SAN.prop.weights,SAN.prop.args,SAN.packagenames,SAN.ignore.finite.offsets Arguments to the `control` functions.
#' @export
sctrl <- statnet.common::sctrl

eval(UPDATE_MY_SCTRL_EXPR)

.RegisterProposals <- function(){
  ergm_proposal_table("c", "Bernoulli", c("", "bd"),  0, "random", "randomtoggle")
  ergm_proposal_table("c", "Bernoulli", c("", "bd"),  1, "TNT", "TNT")
  ergm_proposal_table("c", "Bernoulli", c(".dyads",".dyads+bd"),  -2, "random", "RLE")
  ergm_proposal_table("c", "Bernoulli", c(".dyads",".dyads+bd"),  -1, "TNT", "RLETNT")
  ergm_proposal_table("c", "Bernoulli", "", -100, "TNT10", "TNT10")
  ergm_proposal_table("c", "Bernoulli", "degrees",  0, "random", "CondDegree")
  ergm_proposal_table("c", "Bernoulli", "degreesmix",  0, "random", "CondDegreeMix")
  ergm_proposal_table("c", "Bernoulli", c("idegrees+odegrees","b1degrees+b2degrees"),  0, "random", "CondDegree")
  ergm_proposal_table("c", "Bernoulli", "odegrees",  0, "random", "CondOutDegree")
  ergm_proposal_table("c", "Bernoulli", "idegrees",  0, "random", "CondInDegree")
  ergm_proposal_table("c", "Bernoulli", "b1degrees",  0, "random", "CondB1Degree")
  ergm_proposal_table("c", "Bernoulli", "b2degrees",  0, "random", "CondB2Degree")
  ergm_proposal_table("c", "Bernoulli", "degreedist",  0, "random", "CondDegreeDist")
  ergm_proposal_table("c", "Bernoulli", "idegreedist",  0, "random", "CondInDegreeDist")
  ergm_proposal_table("c", "Bernoulli", "odegreedist",  0, "random", "CondOutDegreeDist")
  ergm_proposal_table("c", "Bernoulli", c("bd+edges","edges"),  0, "random", "ConstantEdges")
  ergm_proposal_table("c", "Bernoulli", "edges+hamming",  0, "random", "HammingConstantEdges")
  ergm_proposal_table("c", "Bernoulli", "hamming",  0, "random", "HammingTNT")
  ergm_proposal_table("c", "Bernoulli", c("bd+observed","observed"),  0, "random", "randomtoggleNonObserved")
  ergm_proposal_table("c", "Bernoulli", c("bd+observed","observed"),  1, "TNT", "NonObservedTNT")
  ergm_proposal_table("c", "Bernoulli", c("blockdiag","bd+blockdiag"), 0, "random", "blockdiag")
  ergm_proposal_table("c", "Bernoulli", c("blockdiag","bd+blockdiag"), 1, "TNT", "blockdiagTNT")
  ergm_proposal_table("c", "Bernoulli", c("blockdiag+observed","bd+blockdiag+observed"),  0, "random", "blockdiagNonObserved")
  ergm_proposal_table("c", "Bernoulli", c("blockdiag+observed","bd+blockdiag+observed"),  1, "TNT", "blockdiagNonObservedTNT")
  ergm_proposal_table("c", "Bernoulli", "fixedas",  0, "random", "fixedas")
  ergm_proposal_table("c", "Bernoulli", "fixedas",  1, "TNT", "fixedasTNT")
  ergm_proposal_table("c", "Bernoulli", "fixallbut",  0, "random", "fixallbut")
  ergm_proposal_table("c", "Bernoulli", "fixallbut",  1, "TNT", "fixallbutTNT")

  
  ergm_proposal_table("c", "StdNormal", "",  0, "random", "StdNormal")
  ergm_proposal_table("c", "StdNormal", ".dyads",  0, "random", "DistRLE")

  ergm_proposal_table("c", "Unif", "",  0, "random", "Unif")
  ergm_proposal_table("c", "Unif", "observed",  0, "random", "UnifNonObserved")
  ergm_proposal_table("c", "Unif", ".dyads",  0, "random", "DistRLE")
  
  ergm_proposal_table("c", "DiscUnif", "",  0, "random", "DiscUnif")
  ergm_proposal_table("c", "DiscUnif", "observed",  0, "random", "DiscUnifNonObserved")  

  ergm_proposal_table("c", c("Unif","DiscUnif","StdNormal","Poisson","Binomial"), ".dyads",  -3, "random", "DistRLE")
}
