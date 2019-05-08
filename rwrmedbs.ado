*!TITLE: RWRMED - causal mediation analysis using regression-with-residuals	
*!AUTHOR: Geoffrey T. Wodtke, Department of Sociology, University of Toronto
*!
*! version 0.9 - 4/25/2019 - GT Wodtke
*!

program define rwrmedbs, rclass
	
	version 13	

	syntax varname(numeric), avar(varname numeric) mvar(varname numeric) ///
			[cvars(varlist numeric)] [lvars(varlist numeric)] ///
			a(real) astar(real) m(real) mreg(string) ///
			[NOINTERaction] ///
			
	local yvar `varlist'

	marksample touse	
		
	local mregtypes linear logistic poisson
	local nmreg : list posof "`mreg'" in mregtypes
	if !`nmreg' {
		display as error "Error: mreg must be chosen from: `mregtypes'."
		error 198		
		}
	else {
		local mreg : word `nmreg' of `mregtypes'
		}
	
	local cvar `cvars'
	local lvar `lvars'
	
	local nc : word count `cvar'

	if "`nointeraction'" == "" {
		local interaction true	
		local inter_var_names "_`avar'_X_`mvar' _`mvar'_X_`avar' _`avar_X_`mvar'_001 _`avar'_X_`mvar'_010 _`avar'_X_`mvar'_100 _`mvar'_X_`avar'_001 _`mvar'_X_`avar'_010 _`mvar'_X_`avar'_100"
		foreach name of local inter_var_names {
			capture confirm new variable `name'
			if !_rc {
				local inter `name'
				continue, break
				}
			}
		if _rc {
			display as error "{p 0 0 5 0}The command needs to create an interaction variable "
			display as error "with one of the following names: `inter_var_names', "
			display as error "but these variables have already been defined.{p_end}"
			error 110
			}
		gen `inter' = `avar'*`mvar'
		}
	else {
		local interaction false
		}

	***************************************
	///block 1: for mreg=linear
	if ("`mreg'"=="linear") {
		
		if (("`cvar'"!="") & ("`lvar'"!="")) {
		
			foreach c in `cvar' {
				capture confirm new variable `c'_r001
				if _rc {
					display as error "{p 0 0 5 0}The command needs to create a residualized variable "
					display as error "with the following name: `c'_r001, "
					display as error "but this variable has already been defined.{p_end}"
					error 110
					}
				}
				
			local cvar_r ""
			foreach c in `cvar' {
				regress `c' 
				predict `c'_r001, resid
				local cvar_r `cvar_r' `c'_r001
				}

			foreach l in `lvar' {
				capture confirm new variable `l'_r001
				if _rc {
					display as error "{p 0 0 5 0}The command needs to create a residualized variable "
					display as error "with the following name: `l'_r001, "
					display as error "but this variable has already been defined.{p_end}"
					error 110
					}
				}

			local lvar_r ""
			foreach l in `lvar' {
				regress `l' `avar' `cvar'
				predict `l'_r001, resid
				local lvar_r `lvar_r' `l'_r001
				}
				
			if ("`interaction'"=="false") {
				regress `yvar' `avar' `mvar' `cvar_r' `lvar_r'
				matrix beta = e(b)
				regress `mvar' `avar' `cvar_r'
				matrix theta = e(b)
				return scalar cde = beta[1,1]*(`astar'-`a')
				return scalar rnde = beta[1,1]*(`astar'-`a')
				return scalar rnie = theta[1,1]*beta[1,2]*(`astar'-`a')
				return scalar rate = (beta[1,1]*(`astar'-`a'))+(theta[1,1]*beta[1,2]*(`astar'-`a'))
				}
		
			if ("`interaction'"=="true") {
				regress `yvar' `avar' `mvar' `inter' `cvar_r' `lvar_r'
				matrix beta = e(b)
				regress `mvar' `avar' `cvar_r'
				matrix theta = e(b)
				return scalar cde = (beta[1,1]+beta[1,3]*`m')*(`astar'-`a')
				return scalar rnde = (beta[1,1]+beta[1,3]*(theta[1,`nc'+2]+theta[1,1]*`a'))*(`astar'-`a')
				return scalar rnie = theta[1,1]*(beta[1,2]+beta[1,3]*`astar')*(`astar'-`a')
				return scalar rate = ((beta[1,1]+beta[1,3]*(theta[1,`nc'+2]+theta[1,1]*`a'))*(`astar'-`a'))+ ///
									 (theta[1,1]*(beta[1,2]+beta[1,3]*`astar')*(`astar'-`a'))
				}
			}

		if (("`cvar'"!="") & ("`lvar'"=="")) {
			
			foreach c in `cvar' {
				capture confirm new variable `c'_r001
				if _rc {
					display as error "{p 0 0 5 0}The command needs to create a residualized variable "
					display as error "with the following name: `c'_r001, "
					display as error "but this variable has already been defined.{p_end}"
					error 110
					}
				}
				
			local cvar_r ""
			foreach c in `cvar' {
				regress `c' 
				predict `c'_r001, resid
				local cvar_r `cvar_r' `c'_r001
				}
					
			if ("`interaction'"=="false") {
				regress `yvar' `avar' `mvar' `cvar_r'
				matrix beta = e(b)
				regress `mvar' `avar' `cvar_r'
				matrix theta = e(b)
				return scalar cde = beta[1,1]*(`astar'-`a')
				return scalar nde = beta[1,1]*(`astar'-`a')
				return scalar nie = theta[1,1]*beta[1,2]*(`astar'-`a')
				return scalar ate = (beta[1,1]*(`astar'-`a'))+(theta[1,1]*beta[1,2]*(`astar'-`a'))
				}
		
			if ("`interaction'"=="true") {
				regress `yvar' `avar' `mvar' `inter' `cvar_r'
				matrix beta = e(b)
				regress `mvar' `avar' `cvar_r'
				matrix theta = e(b)
				return scalar cde = (beta[1,1]+beta[1,3]*`m')*(`astar'-`a')
				return scalar nde = (beta[1,1]+beta[1,3]*(theta[1,`nc'+2]+theta[1,1]*`a'))*(`astar'-`a')
				return scalar nie = theta[1,1]*(beta[1,2]+beta[1,3]*`astar')*(`astar'-`a')
				return scalar ate = ((beta[1,1]+beta[1,3]*(theta[1,`nc'+2]+theta[1,1]*`a'))*(`astar'-`a'))+(theta[1,1]*(beta[1,2]+beta[1,3]*`astar')*(`astar'-`a'))
				}
			}

		if (("`cvar'"=="") & ("`lvar'"!="")) {
			
			foreach l in `lvar' {
				capture confirm new variable `l'_r001
				if _rc {
					display as error "{p 0 0 5 0}The command needs to create a residualized variable "
					display as error "with the following name: `l'_r001, "
					display as error "but this variable has already been defined.{p_end}"
					error 110
					}
				}

			local lvar_r ""
			foreach l in `lvar' {
				regress `l' `avar'
				predict `l'_r001, resid
				local lvar_r `lvar_r' `l'_r001
				}
				
			if ("`interaction'"=="false") {
				regress `yvar' `avar' `mvar' `lvar_r'
				matrix beta = e(b)
				regress `mvar' `avar'
				matrix theta = e(b)
				return scalar cde = beta[1,1]*(`astar'-`a')
				return scalar rnde = beta[1,1]*(`astar'-`a')
				return scalar rnie = theta[1,1]*beta[1,2]*(`astar'-`a')
				return scalar rate = (beta[1,1]*(`astar'-`a'))+(theta[1,1]*beta[1,2]*(`astar'-`a'))
				}
		
			if ("`interaction'"=="true") {
				regress `yvar' `avar' `mvar' `inter' `lvar_r'
				matrix beta = e(b)
				regress `mvar' `avar' 
				matrix theta = e(b)
				return scalar cde = (beta[1,1]+beta[1,3]*`m')*(`astar'-`a')
				return scalar rnde = (beta[1,1]+beta[1,3]*(theta[1,2]+theta[1,1]*`a'))*(`astar'-`a')
				return scalar rnie = theta[1,1]*(beta[1,2]+beta[1,3]*`astar')*(`astar'-`a')
				return scalar rate = ((beta[1,1]+beta[1,3]*(theta[1,2]+theta[1,1]*`a'))*(`astar'-`a'))+ ///
									 (theta[1,1]*(beta[1,2]+beta[1,3]*`astar')*(`astar'-`a'))
				}
			}

		if (("`cvar'"=="") & ("`lvar'"=="")) {
			
			if ("`interaction'"=="false") {
				regress `yvar' `avar' `mvar'
				matrix beta = e(b)
				regress `mvar' `avar'
				matrix theta = e(b)
				return scalar cde = beta[1,1]*(`astar'-`a')
				return scalar nde = beta[1,1]*(`astar'-`a')
				return scalar nie = theta[1,1]*beta[1,2]*(`astar'-`a')
				return scalar ate = (beta[1,1]*(`astar'-`a'))+(theta[1,1]*beta[1,2]*(`astar'-`a'))
				}
		
			if ("`interaction'"=="true") {
				regress `yvar' `avar' `mvar' `inter'
				matrix beta = e(b)
				regress `mvar' `avar' 
				matrix theta = e(b)
				return scalar cde = (beta[1,1]+beta[1,3]*`m')*(`astar'-`a')
				return scalar nde = (beta[1,1]+beta[1,3]*(theta[1,2]+theta[1,1]*`a'))*(`astar'-`a')
				return scalar nie = theta[1,1]*(beta[1,2]+beta[1,3]*`astar')*(`astar'-`a')
				return scalar ate = ((beta[1,1]+beta[1,3]*(theta[1,2]+theta[1,1]*`a'))*(`astar'-`a'))+ ///
									 (theta[1,1]*(beta[1,2]+beta[1,3]*`astar')*(`astar'-`a'))
				}
			}
		}
	///end block 1
	***************************************
	
	***************************************
	///block 2: for mreg=logistic
	if ("`mreg'"=="logistic") {
		
		if (("`cvar'"!="") & ("`lvar'"!="")) {
		
			foreach c in `cvar' {
				capture confirm new variable `c'_r001
				if _rc {
					display as error "{p 0 0 5 0}The command needs to create a residualized variable "
					display as error "with the following name: `c'_r001, "
					display as error "but this variable has already been defined.{p_end}"
					error 110
					}
				}
				
			local cvar_r ""
			foreach c in `cvar' {
				regress `c' 
				predict `c'_r001, resid
				local cvar_r `cvar_r' `c'_r001
				}

			foreach l in `lvar' {
				capture confirm new variable `l'_r001
				if _rc {
					display as error "{p 0 0 5 0}The command needs to create a residualized variable "
					display as error "with the following name: `l'_r001, "
					display as error "but this variable has already been defined.{p_end}"
					error 110
					}
				}

			local lvar_r ""
			foreach l in `lvar' {
				regress `l' `avar' `cvar'
				predict `l'_r001, resid
				local lvar_r `lvar_r' `l'_r001
				}
				
			if ("`interaction'"=="false") {
				regress `yvar' `avar' `mvar' `cvar_r' `lvar_r'
				matrix beta = e(b)
				logit `mvar' `avar' `cvar_r'
				matrix theta = e(b)
				return scalar cde = beta[1,1]*(`astar'-`a')
				return scalar rnde_c = beta[1,1]*(`astar'-`a')
				return scalar rnie_c = beta[1,2]* ///
									   ((exp(theta[1,`nc'+2]+theta[1,1]*`astar')/(1+exp(theta[1,`nc'+2]+theta[1,1]*`astar')))- ///
									   (exp(theta[1,`nc'+2]+theta[1,1]*`a')/(1+exp(theta[1,`nc'+2]+theta[1,1]*`a'))))
				return scalar rate_c = beta[1,1]*(`astar'-`a')+ ///
									   beta[1,2]* ///
									   ((exp(theta[1,`nc'+2]+theta[1,1]*`astar')/(1+exp(theta[1,`nc'+2]+theta[1,1]*`astar')))- ///
									   (exp(theta[1,`nc'+2]+theta[1,1]*`a')/(1+exp(theta[1,`nc'+2]+theta[1,1]*`a'))))
				}
		
			if ("`interaction'"=="true") {
				regress `yvar' `avar' `mvar' `inter' `cvar_r' `lvar_r'
				matrix beta = e(b)
				logit `mvar' `avar' `cvar_r'
				matrix theta = e(b)
				return scalar cde = (beta[1,1]+beta[1,3]*`m')*(`astar'-`a')
				return scalar rnde_c = (beta[1,1]+beta[1,3]*(exp(theta[1,`nc'+2]+theta[1,1]*`a')/(1+exp(theta[1,`nc'+2]+theta[1,1]*`a'))))*(`astar'-`a')
				return scalar rnie_c = (beta[1,2]+beta[1,3]*`astar')* ///
									 ((exp(theta[1,`nc'+2]+theta[1,1]*`astar')/(1+exp(theta[1,`nc'+2]+theta[1,1]*`astar')))- ///
									 (exp(theta[1,`nc'+2]+theta[1,1]*`a')/(1+exp(theta[1,`nc'+2]+theta[1,1]*`a'))))
				return scalar rate_c = (beta[1,1]+beta[1,3]*(exp(theta[1,`nc'+2]+theta[1,1]*`a')/(1+exp(theta[1,`nc'+2]+theta[1,1]*`a'))))*(`astar'-`a')+ ///
									 (beta[1,2]+beta[1,3]*`astar')* ///
									 ((exp(theta[1,`nc'+2]+theta[1,1]*`astar')/(1+exp(theta[1,`nc'+2]+theta[1,1]*`astar')))- ///
									 (exp(theta[1,`nc'+2]+theta[1,1]*`a')/(1+exp(theta[1,`nc'+2]+theta[1,1]*`a'))))
				}
			}

		if (("`cvar'"!="") & ("`lvar'"=="")) {
			
			foreach c in `cvar' {
				capture confirm new variable `c'_r001
				if _rc {
					display as error "{p 0 0 5 0}The command needs to create a residualized variable "
					display as error "with the following name: `c'_r001, "
					display as error "but this variable has already been defined.{p_end}"
					error 110
					}
				}
				
			local cvar_r ""
			foreach c in `cvar' {
				regress `c' 
				predict `c'_r001, resid
				local cvar_r `cvar_r' `c'_r001
				}
					
			if ("`interaction'"=="false") {
				regress `yvar' `avar' `mvar' `cvar_r'
				matrix beta = e(b)
				logit `mvar' `avar' `cvar_r'
				matrix theta = e(b)
				return scalar cde = beta[1,1]*(`astar'-`a')
				return scalar nde_c = beta[1,1]*(`astar'-`a')
				return scalar nie_c = beta[1,2]* ///
									   ((exp(theta[1,`nc'+2]+theta[1,1]*`astar')/(1+exp(theta[1,`nc'+2]+theta[1,1]*`astar')))- ///
									   (exp(theta[1,`nc'+2]+theta[1,1]*`a')/(1+exp(theta[1,`nc'+2]+theta[1,1]*`a'))))
				return scalar ate_c = beta[1,1]*(`astar'-`a')+ ///
									   beta[1,2]* ///
									   ((exp(theta[1,`nc'+2]+theta[1,1]*`astar')/(1+exp(theta[1,`nc'+2]+theta[1,1]*`astar')))- ///
									   (exp(theta[1,`nc'+2]+theta[1,1]*`a')/(1+exp(theta[1,`nc'+2]+theta[1,1]*`a'))))
				}
		
			if ("`interaction'"=="true") {
				regress `yvar' `avar' `mvar' `inter' `cvar_r'
				matrix beta = e(b)
				logit `mvar' `avar' `cvar_r'
				matrix theta = e(b)
				return scalar cde = (beta[1,1]+beta[1,3]*`m')*(`astar'-`a')
				return scalar nde_c = (beta[1,1]+beta[1,3]*(exp(theta[1,`nc'+2]+theta[1,1]*`a')/(1+exp(theta[1,`nc'+2]+theta[1,1]*`a'))))*(`astar'-`a')
				return scalar nie_c = (beta[1,2]+beta[1,3]*`astar')* ///
									 ((exp(theta[1,`nc'+2]+theta[1,1]*`astar')/(1+exp(theta[1,`nc'+2]+theta[1,1]*`astar')))- ///
									 (exp(theta[1,`nc'+2]+theta[1,1]*`a')/(1+exp(theta[1,`nc'+2]+theta[1,1]*`a'))))
				return scalar ate_c = (beta[1,1]+beta[1,3]*(exp(theta[1,`nc'+2]+theta[1,1]*`a')/(1+exp(theta[1,`nc'+2]+theta[1,1]*`a'))))*(`astar'-`a')+ ///
									 (beta[1,2]+beta[1,3]*`astar')* ///
									 ((exp(theta[1,`nc'+2]+theta[1,1]*`astar')/(1+exp(theta[1,`nc'+2]+theta[1,1]*`astar')))- ///
									 (exp(theta[1,`nc'+2]+theta[1,1]*`a')/(1+exp(theta[1,`nc'+2]+theta[1,1]*`a'))))
				}
			}

		if (("`cvar'"=="") & ("`lvar'"!="")) {
			
			foreach l in `lvar' {
				capture confirm new variable `l'_r001
				if _rc {
					display as error "{p 0 0 5 0}The command needs to create a residualized variable "
					display as error "with the following name: `l'_r001, "
					display as error "but this variable has already been defined.{p_end}"
					error 110
					}
				}

			local lvar_r ""
			foreach l in `lvar' {
				regress `l' `avar'
				predict `l'_r001, resid
				local lvar_r `lvar_r' `l'_r001
				}
				
			if ("`interaction'"=="false") {
				regress `yvar' `avar' `mvar' `lvar_r'
				matrix beta = e(b)
				logit `mvar' `avar'
				matrix theta = e(b)
				return scalar cde = beta[1,1]*(`astar'-`a')
				return scalar rnde = beta[1,1]*(`astar'-`a')
				return scalar rnie = beta[1,2]* ///
									((exp(theta[1,2]+theta[1,1]*`astar')/(1+exp(theta[1,2]+theta[1,1]*`astar')))- ///
									(exp(theta[1,2]+theta[1,1]*`a')/(1+exp(theta[1,2]+theta[1,1]*`a'))))
				return scalar rate = beta[1,1]*(`astar'-`a')+ ///
									 beta[1,2]* ///
									 ((exp(theta[1,2]+theta[1,1]*`astar')/(1+exp(theta[1,2]+theta[1,1]*`astar')))- ///
									 (exp(theta[1,2]+theta[1,1]*`a')/(1+exp(theta[1,2]+theta[1,1]*`a'))))
									   
				}

			if ("`interaction'"=="true") {
				regress `yvar' `avar' `mvar' `inter' `lvar_r'
				matrix beta = e(b)
				logit `mvar' `avar' 
				matrix theta = e(b)
				return scalar cde = (beta[1,1]+beta[1,3]*`m')*(`astar'-`a')
				return scalar rnde = (beta[1,1]+beta[1,3]*(exp(theta[1,2]+theta[1,1]*`a')/(1+exp(theta[1,2]+theta[1,1]*`a'))))*(`astar'-`a')
				return scalar rnie = (beta[1,2]+beta[1,3]*`astar')* ///
									 ((exp(theta[1,2]+theta[1,1]*`astar')/(1+exp(theta[1,2]+theta[1,1]*`astar')))- ///
									 (exp(theta[1,2]+theta[1,1]*`a')/(1+exp(theta[1,2]+theta[1,1]*`a'))))
				return scalar rate = (beta[1,1]+beta[1,3]*(exp(theta[1,2]+theta[1,1]*`a')/(1+exp(theta[1,2]+theta[1,1]*`a'))))*(`astar'-`a')+ ///
									 (beta[1,2]+beta[1,3]*`astar')* ///
									 ((exp(theta[1,2]+theta[1,1]*`astar')/(1+exp(theta[1,2]+theta[1,1]*`astar')))- ///
									 (exp(theta[1,2]+theta[1,1]*`a')/(1+exp(theta[1,2]+theta[1,1]*`a'))))
				}
			}

		if (("`cvar'"=="") & ("`lvar'"=="")) {
			
			if ("`interaction'"=="false") {
				regress `yvar' `avar' `mvar'
				matrix beta = e(b)
				logit `mvar' `avar'
				matrix theta = e(b)
				return scalar cde = beta[1,1]*(`astar'-`a')
				return scalar nde = beta[1,1]*(`astar'-`a')
				return scalar nie = beta[1,2]* ///
									((exp(theta[1,2]+theta[1,1]*`astar')/(1+exp(theta[1,2]+theta[1,1]*`astar')))- ///
									(exp(theta[1,2]+theta[1,1]*`a')/(1+exp(theta[1,2]+theta[1,1]*`a'))))
				return scalar ate = beta[1,1]*(`astar'-`a')+ ///
									 beta[1,2]* ///
									 ((exp(theta[1,2]+theta[1,1]*`astar')/(1+exp(theta[1,2]+theta[1,1]*`astar')))- ///
									 (exp(theta[1,2]+theta[1,1]*`a')/(1+exp(theta[1,2]+theta[1,1]*`a'))))
				}
		
			if ("`interaction'"=="true") {
				regress `yvar' `avar' `mvar' `inter'
				matrix beta = e(b)
				logit `mvar' `avar' 
				matrix theta = e(b)
				return scalar cde = (beta[1,1]+beta[1,3]*`m')*(`astar'-`a')
				return scalar nde = (beta[1,1]+beta[1,3]*(exp(theta[1,2]+theta[1,1]*`a')/(1+exp(theta[1,2]+theta[1,1]*`a'))))*(`astar'-`a')
				return scalar nie = (beta[1,2]+beta[1,3]*`astar')* ///
									 ((exp(theta[1,2]+theta[1,1]*`astar')/(1+exp(theta[1,2]+theta[1,1]*`astar')))- ///
									 (exp(theta[1,2]+theta[1,1]*`a')/(1+exp(theta[1,2]+theta[1,1]*`a'))))
				return scalar ate = (beta[1,1]+beta[1,3]*(exp(theta[1,2]+theta[1,1]*`a')/(1+exp(theta[1,2]+theta[1,1]*`a'))))*(`astar'-`a')+ ///
									 (beta[1,2]+beta[1,3]*`astar')* ///
									 ((exp(theta[1,2]+theta[1,1]*`astar')/(1+exp(theta[1,2]+theta[1,1]*`astar')))- ///
									 (exp(theta[1,2]+theta[1,1]*`a')/(1+exp(theta[1,2]+theta[1,1]*`a'))))
				}
			}
		}
	///end block 2
	***************************************

	***************************************
	///block 3: for mreg=poisson
	if ("`mreg'"=="poisson") {
		
		if (("`cvar'"!="") & ("`lvar'"!="")) {
		
			foreach c in `cvar' {
				capture confirm new variable `c'_r001
				if _rc {
					display as error "{p 0 0 5 0}The command needs to create a residualized variable "
					display as error "with the following name: `c'_r001, "
					display as error "but this variable has already been defined.{p_end}"
					error 110
					}
				}
				
			local cvar_r ""
			foreach c in `cvar' {
				regress `c' 
				predict `c'_r001, resid
				local cvar_r `cvar_r' `c'_r001
				}

			foreach l in `lvar' {
				capture confirm new variable `l'_r001
				if _rc {
					display as error "{p 0 0 5 0}The command needs to create a residualized variable "
					display as error "with the following name: `l'_r001, "
					display as error "but this variable has already been defined.{p_end}"
					error 110
					}
				}

			local lvar_r ""
			foreach l in `lvar' {
				regress `l' `avar' `cvar'
				predict `l'_r001, resid
				local lvar_r `lvar_r' `l'_r001
				}

			if ("`interaction'"=="false") {
				regress `yvar' `avar' `mvar' `cvar_r' `lvar_r'
				matrix beta = e(b)
				poisson `mvar' `avar' `cvar_r'
				matrix theta = e(b)
				return scalar cde = beta[1,1]*(`astar'-`a')
				return scalar rnde_c = beta[1,1]*(`astar'-`a')
				return scalar rnie_c = beta[1,2]*((exp(theta[1,`nc'+2]+theta[1,1]*`astar'))-(exp(theta[1,`nc'+2]+theta[1,1]*`a')))
				return scalar rate_c = beta[1,1]*(`astar'-`a')+ ///
									   beta[1,2]*((exp(theta[1,`nc'+2]+theta[1,1]*`astar'))-(exp(theta[1,`nc'+2]+theta[1,1]*`a')))
				}
		
			if ("`interaction'"=="true") {
				regress `yvar' `avar' `mvar' `inter' `cvar_r' `lvar_r'
				matrix beta = e(b)
				poisson `mvar' `avar' `cvar_r'
				matrix theta = e(b)
				return scalar cde = (beta[1,1]+beta[1,3]*`m')*(`astar'-`a')
				return scalar rnde_c = (beta[1,1]+beta[1,3]*(exp(theta[1,`nc'+2]+theta[1,1]*`a')))*(`astar'-`a')
				return scalar rnie_c = (beta[1,2]+beta[1,3]*`astar')*((exp(theta[1,`nc'+2]+theta[1,1]*`astar'))-(exp(theta[1,`nc'+2]+theta[1,1]*`a')))
				return scalar rate_c = (beta[1,1]+beta[1,3]*(exp(theta[1,`nc'+2]+theta[1,1]*`a')))*(`astar'-`a')+ ///
									   (beta[1,2]+beta[1,3]*`astar')*((exp(theta[1,`nc'+2]+theta[1,1]*`astar'))-(exp(theta[1,`nc'+2]+theta[1,1]*`a')))
				}
			}

		if (("`cvar'"!="") & ("`lvar'"=="")) {
			
			foreach c in `cvar' {
				capture confirm new variable `c'_r001
				if _rc {
					display as error "{p 0 0 5 0}The command needs to create a residualized variable "
					display as error "with the following name: `c'_r001, "
					display as error "but this variable has already been defined.{p_end}"
					error 110
					}
				}
				
			local cvar_r ""
			foreach c in `cvar' {
				regress `c' 
				predict `c'_r001, resid
				local cvar_r `cvar_r' `c'_r001
				}
					
			if ("`interaction'"=="false") {
				regress `yvar' `avar' `mvar' `cvar_r'
				matrix beta = e(b)
				poisson `mvar' `avar' `cvar_r'
				matrix theta = e(b)
				return scalar cde = beta[1,1]*(`astar'-`a')
				return scalar nde_c = beta[1,1]*(`astar'-`a')
				return scalar nie_c = beta[1,2]*((exp(theta[1,`nc'+2]+theta[1,1]*`astar'))-(exp(theta[1,`nc'+2]+theta[1,1]*`a')))
				return scalar ate_c = beta[1,1]*(`astar'-`a')+ ///
									  beta[1,2]*((exp(theta[1,`nc'+2]+theta[1,1]*`astar'))-(exp(theta[1,`nc'+2]+theta[1,1]*`a')))									   
				}

			if ("`interaction'"=="true") {
				regress `yvar' `avar' `mvar' `inter' `cvar_r'
				matrix beta = e(b)
				poisson `mvar' `avar' `cvar_r'
				matrix theta = e(b)
				return scalar cde = (beta[1,1]+beta[1,3]*`m')*(`astar'-`a')
				return scalar nde_c = (beta[1,1]+beta[1,3]*(exp(theta[1,`nc'+2]+theta[1,1]*`a')))*(`astar'-`a')
				return scalar nie_c = (beta[1,2]+beta[1,3]*`astar')*((exp(theta[1,`nc'+2]+theta[1,1]*`astar'))-(exp(theta[1,`nc'+2]+theta[1,1]*`a')))
				return scalar ate_c = (beta[1,1]+beta[1,3]*(exp(theta[1,`nc'+2]+theta[1,1]*`a')))*(`astar'-`a')+ ///
									  (beta[1,2]+beta[1,3]*`astar')*((exp(theta[1,`nc'+2]+theta[1,1]*`astar'))-(exp(theta[1,`nc'+2]+theta[1,1]*`a')))									 
				}
			}

		if (("`cvar'"=="") & ("`lvar'"!="")) {
			
			foreach l in `lvar' {
				capture confirm new variable `l'_r001
				if _rc {
					display as error "{p 0 0 5 0}The command needs to create a residualized variable "
					display as error "with the following name: `l'_r001, "
					display as error "but this variable has already been defined.{p_end}"
					error 110
					}
				}

			local lvar_r ""
			foreach l in `lvar' {
				regress `l' `avar'
				predict `l'_r001, resid
				local lvar_r `lvar_r' `l'_r001
				}

			if ("`interaction'"=="false") {
				regress `yvar' `avar' `mvar' `lvar_r'
				matrix beta = e(b)
				poisson `mvar' `avar'
				matrix theta = e(b)
				return scalar cde = beta[1,1]*(`astar'-`a')
				return scalar rnde = beta[1,1]*(`astar'-`a')
				return scalar rnie = beta[1,2]*((exp(theta[1,2]+theta[1,1]*`astar'))-(exp(theta[1,2]+theta[1,1]*`a')))
				return scalar rate = beta[1,1]*(`astar'-`a')+ ///
									 beta[1,2]*((exp(theta[1,2]+theta[1,1]*`astar'))-(exp(theta[1,2]+theta[1,1]*`a')))
				}

			if ("`interaction'"=="true") {
				regress `yvar' `avar' `mvar' `inter' `lvar_r'
				matrix beta = e(b)
				poisson `mvar' `avar' 
				matrix theta = e(b)
				return scalar cde = (beta[1,1]+beta[1,3]*`m')*(`astar'-`a')
				return scalar rnde = (beta[1,1]+beta[1,3]*(exp(theta[1,2]+theta[1,1]*`a')))*(`astar'-`a')
				return scalar rnie = (beta[1,2]+beta[1,3]*`astar')*((exp(theta[1,2]+theta[1,1]*`astar'))-(exp(theta[1,2]+theta[1,1]*`a')))
				return scalar rate = (beta[1,1]+beta[1,3]*(exp(theta[1,2]+theta[1,1]*`a')))*(`astar'-`a')+ ///
									 (beta[1,2]+beta[1,3]*`astar')*((exp(theta[1,2]+theta[1,1]*`astar'))-(exp(theta[1,2]+theta[1,1]*`a')))
				}
			}

		if (("`cvar'"=="") & ("`lvar'"=="")) {
			
			if ("`interaction'"=="false") {
				regress `yvar' `avar' `mvar'
				matrix beta = e(b)
				poisson `mvar' `avar'
				matrix theta = e(b)
				return scalar cde = beta[1,1]*(`astar'-`a')
				return scalar nde = beta[1,1]*(`astar'-`a')
				return scalar nie = beta[1,2]*((exp(theta[1,2]+theta[1,1]*`astar'))-(exp(theta[1,2]+theta[1,1]*`a')))
				return scalar ate = beta[1,1]*(`astar'-`a')+ ///
									beta[1,2]*((exp(theta[1,2]+theta[1,1]*`astar'))-(exp(theta[1,2]+theta[1,1]*`a')))									 
				}
		
			if ("`interaction'"=="true") {
				regress `yvar' `avar' `mvar' `inter'
				matrix beta = e(b)
				poisson `mvar' `avar' 
				matrix theta = e(b)
				return scalar cde = (beta[1,1]+beta[1,3]*`m')*(`astar'-`a')
				return scalar nde = (beta[1,1]+beta[1,3]*(exp(theta[1,2]+theta[1,1]*`a')))*(`astar'-`a')
				return scalar nie = (beta[1,2]+beta[1,3]*`astar')*((exp(theta[1,2]+theta[1,1]*`astar'))-(exp(theta[1,2]+theta[1,1]*`a')))
				return scalar ate = (beta[1,1]+beta[1,3]*(exp(theta[1,2]+theta[1,1]*`a')))*(`astar'-`a')+ ///
									(beta[1,2]+beta[1,3]*`astar')*((exp(theta[1,2]+theta[1,1]*`astar'))-(exp(theta[1,2]+theta[1,1]*`a')))
				}
			}
		}
	///end block 3
	***************************************
	
	if ("`interaction'"=="true") {
		drop `inter'
		}
	
	if ("`cvar'"!="") {
		drop `cvar_r'
		}
		
	if ("`lvar'"!="") {
		drop `lvar_r'
		}

end paramedbs
