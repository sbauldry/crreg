*! v1.0.0, S Bauldry, 15nov2016

capture program drop crreg
program crreg, properties(swml svyb svyj svyr mi)
	version 14.1
	if replay() {
		if (`"`e(cmd)'"' != "crreg") error 301
		Replay `0'
	}
	else Estimate `0'
end


capture program drop Estimate
program Estimate, eclass sortpreserve
	syntax varlist(numeric) [if] [in] [pweight fweight aweight iweight] [, /// 	       
		   PRop(varlist)     /// vars with proportionality constraint
		   FRee(varlist)     /// vars with no constraint
		   LINK(string)      /// link function (logit, probit, or cloglog)
		   Level(cilevel)    /// display option for confidence intervals
		   vce(string)]      /// robust and cluster robust standard errors              
		   
	marksample touse
	
	* identify DV and IVs
	gettoken Y X : varlist
	_fv_check_depvar `Y'
	
	* set globals for # categories
	qui tab `Y'
	global nCat   = r(r)
	global nCatm1 = r(r) - 1
	
		* too few or two many categories
		if ( $nCat < 3 ) {
			dis ""
			dis as error "{yellow}`Y'{red} has $nCat categories - a minimum" ///
			             " of 3 is required"
			exit 148
		}
		else if ( $nCat > 10 ) {
			dis ""
			dis as error "{yellow}`Y'{red} has $nCat categories - a maximum" ///
			             " of 10 is allowed"
			exit 149
		}
		
	* check for link function (default to logit)
	if ( "`link'" == "logit" | "`link'" == "l" | "`link'" == "" ) {
		global Link       "logit"
		local  link_title "Ordered Logit Estimates"
	}
	else if ( "`link'" == "probit" | "`link'" == "p" ) {
		global Link       "probit"
		local  link_title "Ordered Probit Estimates"
	} 
	else if ( "`link'" == "cloglog" | "`link'" == "c" ) {
		global Link       "cloglog"
		local  link_title "Ordered Complementary Log-Log Estimates"
	} 
	else {
		dis ""
		dis as error "{yellow}`link'{red} is not a supported link function"
		exit 198
	}
	
	* parse VCE statement
	if ( `"`vce'"' != "" ) {
		my_vce_parse , vce(`vce')
		local vcetype    "robust"
		local clustervar "`r(clustervar)'"
		if "`clustervar'" != "" {
			markout `touse' `clustervar'
		}
	}
	
	* parse weight statement
	if ( "`weight'" != "" ) {
		local wgt "[`weight'`exp']"
	}
		
	* prepare IVs
	if ( "`X'" != "" ) {
		fvexpand `X'
		local cnIV `r(varlist)'
	}
	
	* prepare no constraint variables
	if ( "`free'" != "" ) {
		fvexpand `free'
		local frIV `r(varlist)'
		
		* verify subset of IVs
		local frchk : list local(frIV) - local(cnIV)
		if ( "`frchk'" != "" ) {
			dis ""
			dis as error "free{yellow}(`frchk'){red} is not included in" ///
			             " the list of independent variables: {yellow}`IV'"
			exit 198
		}
		
		* remove from list of IVs
		local cnIV : list local(cnIV) - local(frIV)
	}
	
	* prepare proportionality constraint variables
	if ( "`prop'" != "" ) {
		fvexpand `prop'
		local prIV `r(varlist)'
		
		* verify subset of IVs
		local prchk : list local(prIV) - local(cnIV)
		if ( "`prchk'" != "" ) {
			dis ""
			dis as error "prop{yellow}(`prchk'){red} is not included in" ///
			             " the list of independent variables: {yellow}`IV'"
			exit 198
		}
		
		* remove from list of IVs
		local cnIV : list local(cnIV) - local(prIV)
	}
	

	* create ML model statements
	
	* case 1: all variables with parallel assumption
	if ( "`free'" == "" & "`prop'" == "" & "`cnIV'" != "" ) {
		dis "case1"
		local model "(cns: `Y' = `cnIV', nocons)"
	
		forval i = 1/$nCatm1 {
			local model "`model' /tau`i'"
		}
	
		* obtain ML estimates
		ml model lf cr_lf `model' `wgt' if `touse', title(`link_title') ///
		   vce(`vcetype') maximize
			
		* replace current b, V, and eqnames matrices
		tempname b v
		mat `b' = e(b)
		mat `v' = e(V)
	}
	
	* case 2: all variables with non-parallel assumption
	if ( "`free'" != "" & "`prop'" == "" & "`cnIV'" == "" ) {
		dis "case2"
		local model "(eq1: `Y' = `free')"
		
		forval i = 2/$nCatm1 {
			local model "`model' (eq`i': `free')"
		}
		
		* obtain ML estimates
		ml model lf cr_np_lf `model' `wgt' if `touse', title(`link_title') ///
		   vce(`vcetype') maximize
		
		* replace current b, V, and eqnames matrices
		tempname b v
		mat `b' = e(b)
		mat `v' = e(V)
	}
	
	* case 3: subset of variables with non-parallel assumption
	if ( "`free'" != "" & "`prop'" == "" & "`cnIV'" != "" )  {
		dis "case3"
		local model "(cns: `Y' = `cnIV', nocons)"
		
		forval i = 1/$nCatm1 {
			local model "`model' (eq`i': `free')"
		}
		
		* obtain ML estimates
		ml model lf cr_pnp_lf `model' `wgt' if `touse', title(`link_title') ///
		   vce(`vcetype') maximize
		
		* replace current b, V, and eqnames matrices
		tempname b v
		mat `b' = e(b)
		mat `v' = e(V)
	}
	
	* case 4: subset of variables with proportionality assumption
	if ( "`free'" == "" & "`prop'" != "" & "`cnIV'" != "" )  {
		dis "case4"
		local model "(cns: `Y' = `cnIV', nocons) (prp: `prIV')"
		
		forval i = 2/$nCatm1 {
			local model "`model' /phi`i'"
		}	
		
		* obtain ML estimates
		ml model lf cr_pc_lf `model' `wgt' if `touse', title(`link_title') ///
		   vce(`vcetype') maximize
		
		* replace current b, V, and eqnames matrices
		tempname b v
		mat `b' = e(b)
		mat `v' = e(V)
	}
	
	* case 5: subset of variables with non-parallel assumption and 
	*         proportionality assumption
	if ( "`free'" != "" & "`prop'" != "" & "`cnIV'" != "" )  {
		dis "case5"
		local model "(cns: `Y' = `cnIV', nocons) (prp: `prIV', nocons)"
		
		forval i = 1/$nCatm1 {
			local model "`model' (eq`i': `free')"
		}
		
		forval i = 2/$nCatm1 {
			local model "`model' /phi`i'"
		}
		
		* obtain ML estimates
		ml model lf cr_ppc_lf `model' `wgt' if `touse', title(`link_title') ///
		   vce(`vcetype') maximize
		
		* replace current b, V, and eqnames matrices
		tempname b v
		mat `b' = e(b)
		mat `v' = e(V)
	}	

	* return and display results
	ereturn local cmd crreg
	ereturn repost b = `b' V = `v', rename
	Replay, level(`level')
end



capture program drop Replay
program Replay
	syntax [, Level(cilevel) or]
	ml display, level(`level')
end



capture program drop my_vce_parse
program my_vce_parse, rclass
	syntax [, vce(string) ]
	
	local case : word count `vce'
	
	if ( `case' > 2 ) {
		dis `"{red}{bf:vce(`vce')} invalid"'
		exit 498
	}
	
	local 0 `", `vce'"'
	syntax [, Robust Cluster * ]
	
	if ( `case' == 2 ) {
		if "`robust'" == "robust" | "`cluster'" == "" {
			dis `"{red}{bf:vce(`vce')} invalid"'
			exit 498
		}
		
		capture confirm numeric variable `options'
		if _rc {
			dis `"{red}{bf:vce(`vce')} invalid"'
			exit 498
		}
		local clustervar "`options'"
	}
	else {
		if ( "`robust'" == "" ) {
			dis `"{red}{bf:vce(`vce')} invalid"'
			exit 498
		}
	}
	
	return clear
	return local clustervar "`clustervar'"
end

/* History
1.0.0  11.15.16  initial program for arbitrary number of categories






