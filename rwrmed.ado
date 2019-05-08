*!TITLE: RWRMED - causal mediation analysis using regression-with-residuals	
*!AUTHOR: Geoffrey T. Wodtke, Department of Sociology, University of Toronto
*!
*! version 0.9 - 4/25/2019 - GT Wodtke
*!

program define rwrmed, eclass

	version 13	

	syntax varname(numeric), avar(varname numeric) mvar(varname numeric) ///
			[cvars(varlist numeric)] [lvars(varlist numeric)] ///
			a(real) astar(real) m(real) mreg(string) [NOINTERaction] ///
			[reps(integer 200) strata(varname numeric) cluster(varname numeric) level(cilevel) seed(passthru)]

	local mregtypes linear logistic poisson
	local nmreg : list posof "`mreg'" in mregtypes
	if !`nmreg' {
		display as error "Error: mreg must be chosen from: `mregtypes'."
		error 198		
		}
	
	***************************************
	///block 1: for mreg=linear
	if ("`mreg'"=="linear") {
		
		if ("`lvars'"!="") {
	
			bootstrap CDE=r(cde) RNDE=r(rnde) RNIE=r(rnie) RATE=r(rate), ///
				reps(`reps') strata(`strata') cluster(`cluster') level(`level') `seed' noheader: ///
				rwrmedbs `varlist', ///
				avar(`avar') mvar(`mvar') cvars(`cvars') lvars(`lvars') ///
				a(`a') astar(`astar') m(`m') mreg(`mreg') `nointeraction'
	
			di "CDE: controlled direct effect at m=`m'"
			di "RNDE: randomized intervention analogue of the natural direct effect"
			di "RNIE: randomized intervention analogue of the natural indirect effect"
			di "RATE: randomized intervention analogue of the total effect"
			}
		
		if ("`lvars'"=="") {
	
			bootstrap CDE=r(cde) NDE=r(nde) NIE=r(nie) ATE=r(ate), ///
				reps(`reps') strata(`strata') cluster(`cluster') level(`level') `seed' noheader: ///
				rwrmedbs `varlist', ///
				avar(`avar') mvar(`mvar') cvars(`cvars') lvars(`lvars') ///
				a(`a') astar(`astar') m(`m') mreg(`mreg') `nointeraction'
	
			di "CDE: controlled direct effect at m=`m'"
			di "NDE: natural direct effect"
			di "NIE: natural indirect effect"
			di "ATE: average total effect"
			}
		}
	///end block 1
	***************************************
	
	***************************************
	///block 2: for mreg=logistic
	if ("`mreg'"=="logistic") {
		
		if (("`cvars'"!="") & ("`lvars'"!="")) {
	
			bootstrap CDE=r(cde) RNDEc=r(rnde_c) RNIEc=r(rnie_c) RATEc=r(rate_c), ///
				reps(`reps') strata(`strata') cluster(`cluster') level(`level') `seed' noheader: ///
				rwrmedbs `varlist', ///
				avar(`avar') mvar(`mvar') cvars(`cvars') lvars(`lvars') ///
				a(`a') astar(`astar') m(`m') mreg(`mreg') `nointeraction'
	
			di "CDE: controlled direct effect at m=`m'"
			di "RNDEc: randomized intervention analogue of the natural direct effect at sample means of cvars"
			di "RNIEc: randomized intervention analogue of the natural indirect effect at sample means of cvars"
			di "RATEc: randomized intervention analogue of the total effect at sample means of cvars"
			}
		
		if (("`cvars'"!="") & ("`lvars'"=="")) {
	
			bootstrap CDE=r(cde) NDEc=r(nde_c) NIEc=r(nie_c) ATEc=r(ate_c), ///
				reps(`reps') strata(`strata') cluster(`cluster') level(`level') `seed' noheader: ///
				rwrmedbs `varlist', ///
				avar(`avar') mvar(`mvar') cvars(`cvars') lvars(`lvars') ///
				a(`a') astar(`astar') m(`m') mreg(`mreg') `nointeraction'
	
			di "CDE: controlled direct effect at m=`m'"
			di "NDEc: natural direct effect at sample means of cvars"
			di "NIEc: natural indirect effect at sample means of cvars"
			di "ATEc: total effect at sample means of cvars"
			}
		
		if (("`cvars'"=="") & ("`lvars'"!="")) {
	
			bootstrap CDE=r(cde) RNDE=r(rnde) RNIE=r(rnie) RATE=r(rate), ///
				reps(`reps') strata(`strata') cluster(`cluster') level(`level') `seed' noheader: ///
				rwrmedbs `varlist', ///
				avar(`avar') mvar(`mvar') cvars(`cvars') lvars(`lvars') ///
				a(`a') astar(`astar') m(`m') mreg(`mreg') `nointeraction'
	
			di "CDE: controlled direct effect at m=`m'"
			di "RNDE: randomized intervention analogue of the natural direct effect"
			di "RNIE: randomized intervention analogue of the natural indirect effect"
			di "RATE: randomized intervention analogue of the total effect"
			}
		
		if (("`cvars'"=="") & ("`lvars'"=="")) {
	
			bootstrap CDE=r(cde) NDE=r(nde) NIE=r(nie) ATE=r(ate), ///
				reps(`reps') strata(`strata') cluster(`cluster') level(`level') `seed' noheader: ///
				rwrmedbs `varlist', ///
				avar(`avar') mvar(`mvar') cvars(`cvars') lvars(`lvars') ///
				a(`a') astar(`astar') m(`m') mreg(`mreg') `nointeraction'
	
			di "CDE: controlled direct effect at m=`m'"
			di "NDE: natural direct effect"
			di "NIE: natural indirect effect"
			di "ATE: average total effect"
			}
		}
	///end block 2
	***************************************

	***************************************
	///block 3: for mreg=poisson
	if ("`mreg'"=="poisson") {
		
		if (("`cvars'"!="") & ("`lvars'"!="")) {
	
			bootstrap CDE=r(cde) RNDEc=r(rnde_c) RNIEc=r(rnie_c) RATEc=r(rate_c), ///
				reps(`reps') strata(`strata') cluster(`cluster') level(`level') `seed' noheader: ///
				rwrmedbs `varlist', ///
				avar(`avar') mvar(`mvar') cvars(`cvars') lvars(`lvars') ///
				a(`a') astar(`astar') m(`m') mreg(`mreg') `nointeraction'
	
			di "CDE: controlled direct effect at m=`m'"
			di "RNDEc: randomized intervention analogue of the natural direct effect at sample means of cvars"
			di "RNIEc: randomized intervention analogue of the natural indirect effect at sample means of cvars"
			di "RATEc: randomized intervention analogue of the total effect at sample means of cvars"
			}
		
		if (("`cvars'"!="") & ("`lvars'"=="")) {
	
			bootstrap CDE=r(cde) NDEc=r(nde_c) NIEc=r(nie_c) ATEc=r(ate_c), ///
				reps(`reps') strata(`strata') cluster(`cluster') level(`level') `seed' noheader: ///
				rwrmedbs `varlist', ///
				avar(`avar') mvar(`mvar') cvars(`cvars') lvars(`lvars') ///
				a(`a') astar(`astar') m(`m') mreg(`mreg') `nointeraction'
	
			di "CDE: controlled direct effect at m=`m'"
			di "NDEc: natural direct effect at sample means of cvars"
			di "NIEc: natural indirect effect at sample means of cvars"
			di "ATEc: total effect at sample means of cvars"
			}
		
		if (("`cvars'"=="") & ("`lvars'"!="")) {
	
			bootstrap CDE=r(cde) RNDE=r(rnde) RNIE=r(rnie) RATE=r(rate), ///
				reps(`reps') strata(`strata') cluster(`cluster') level(`level') `seed' noheader: ///
				rwrmedbs `varlist', ///
				avar(`avar') mvar(`mvar') cvars(`cvars') lvars(`lvars') ///
				a(`a') astar(`astar') m(`m') mreg(`mreg') `nointeraction'
	
			di "CDE: controlled direct effect at m=`m'"
			di "RNDE: randomized intervention analogue of the natural direct effect"
			di "RNIE: randomized intervention analogue of the natural indirect effect"
			di "RATE: randomized intervention analogue of the total effect"
			}
		
		if (("`cvars'"=="") & ("`lvars'"=="")) {
	
			bootstrap CDE=r(cde) NDE=r(nde) NIE=r(nie) ATE=r(ate), ///
				reps(`reps') strata(`strata') cluster(`cluster') level(`level') `seed' noheader: ///
				rwrmedbs `varlist', ///
				avar(`avar') mvar(`mvar') cvars(`cvars') lvars(`lvars') ///
				a(`a') astar(`astar') m(`m') mreg(`mreg') `nointeraction'
	
			di "CDE: controlled direct effect at m=`m'"
			di "NDE: natural direct effect"
			di "NIE: natural indirect effect"
			di "ATE: average total effect"
			}
		}
	///end block 3
	***************************************

end rwrmed
