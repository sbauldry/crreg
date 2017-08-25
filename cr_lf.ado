*! v1.1.0, S Bauldry, 25aug2017

capture program drop cr_lf
program cr_lf
  version 14
  
  * creating arguments based on number of categories of Y stored in globals
  local arguments "lnf xb"
  
  forval i = 1/$nCatm1 {
    local arguments "`arguments' f`i'"
  }
	
  args `arguments'
	
  * tempvars for thresholds
  forval j = 1/$nCatm1 {
    tempvar tau`j'
    qui gen double `tau`j'' = `f`j''
  }
	
  * likelihood function for logit link
  if ( "$Link" == "logit" ) {	
		
    * equation for first value of Y
    qui replace `lnf' = ln(invlogit(`tau1' - `xb')) if $ML_y == 1
		
    * build equations for middle values of Y
    forval k = 2/$nCatm1 {
      local meqn_b `" ln(invlogit(`tau`k'' - `xb')) "'
    
	  local meqn_a ""
      forval m = 1/`k' {
        local meqn_a `" `eqn_a' ln(1 - invlogit(`tau`m'' - `xb')) + "'
      }
	
      local meqn `" `meqn_a' `meqn_b' "'
      qui replace `lnf' = `meqn' if $ML_y == `k'
    }
	
	* build equation for last value of Y
	local eqn `" ln(1 - invlogit(`tau1' - `xb')) "'
	forval n = 2/$nCatm1 {
      local eqn `" `eqn' + ln(1 - invlogit(`tau`n'' - `xb')) "'
    }
	qui replace `lnf' = `eqn' if $ML_y == $nCat
  }
	
	
	* likelihood function for probit link
	if ( "$Link" == "probit" ) {	
		qui replace `lnf' = ln(normal(`tau1' - `xb')) if $ML_y == 1
		
		qui replace `lnf' = ln(1 - normal(`tau1' - `xb')) + ///
		                    ln(normal(`tau2' - `xb')) if $ML_y == 2
		
		if( $nCat > 3 ) {
			qui replace `lnf' = ln(1 - normal(`tau1' - `xb')) + ///
			                    ln(1 - normal(`tau2' - `xb')) + ///
		                        ln(normal(`tau3' - `xb')) if $ML_y == 3
		}
		
		if( $nCat > 4 ) {
			qui replace `lnf' = ln(1 - normal(`tau1' - `xb')) + ///
			                    ln(1 - normal(`tau2' - `xb')) + ///
								ln(1 - normal(`tau3' - `xb')) + ///
		                        ln(normal(`tau4' - `xb')) if $ML_y == 4
		}
		
		if( $nCat > 5 ) {
			qui replace `lnf' = ln(1 - normal(`tau1' - `xb')) + ///
			                    ln(1 - normal(`tau2' - `xb')) + ///
								ln(1 - normal(`tau3' - `xb')) + ///
								ln(1 - normal(`tau4' - `xb')) + ///
		                        ln(normal(`tau5' - `xb')) if $ML_y == 5
		}
		
		if( $nCat > 6 ) {
			qui replace `lnf' = ln(1 - normal(`tau1' - `xb')) + ///
			                    ln(1 - normal(`tau2' - `xb')) + ///
								ln(1 - normal(`tau3' - `xb')) + ///
								ln(1 - normal(`tau4' - `xb')) + ///
								ln(1 - normal(`tau5' - `xb')) + ///
		                        ln(normal(`tau6' - `xb')) if $ML_y == 6
		}
		
		if( $nCat > 7 ) {
			qui replace `lnf' = ln(1 - normal(`tau1' - `xb')) + ///
			                    ln(1 - normal(`tau2' - `xb')) + ///
								ln(1 - normal(`tau3' - `xb')) + ///
								ln(1 - normal(`tau4' - `xb')) + ///
								ln(1 - normal(`tau5' - `xb')) + ///
								ln(1 - normal(`tau6' - `xb')) + ///
		                        ln(normal(`tau7' - `xb')) if $ML_y == 7
		}
		
		if( $nCat > 8 ) {
			qui replace `lnf' = ln(1 - normal(`tau1' - `xb')) + ///
			                    ln(1 - normal(`tau2' - `xb')) + ///
								ln(1 - normal(`tau3' - `xb')) + ///
								ln(1 - normal(`tau4' - `xb')) + ///
								ln(1 - normal(`tau5' - `xb')) + ///
								ln(1 - normal(`tau6' - `xb')) + ///
								ln(1 - normal(`tau7' - `xb')) + ///
		                        ln(normal(`tau8' - `xb')) if $ML_y == 8
		}
		
		if( $nCat > 9 ) {
			qui replace `lnf' = ln(1 - normal(`tau1' - `xb')) + ///
			                    ln(1 - normal(`tau2' - `xb')) + ///
								ln(1 - normal(`tau3' - `xb')) + ///
								ln(1 - normal(`tau4' - `xb')) + ///
								ln(1 - normal(`tau5' - `xb')) + ///
								ln(1 - normal(`tau6' - `xb')) + ///
								ln(1 - normal(`tau7' - `xb')) + ///
								ln(1 - normal(`tau8' - `xb')) + ///
		                        ln(normal(`tau9' - `xb')) if $ML_y == 9
		}
		
		* build equation for last value of Y
		local eqn `" ln(1 - normal(`tau1' - `xb')) "'
		forval i = 2/$nCatm1 {
			local eqn `" `eqn' + ln(1 - normal(`tau`i'' - `xb')) "'
		}
		qui replace `lnf' = `eqn' if $ML_y == $nCat
	}
	
	
	* likelihood function for complementary log-log link
	if ( "$Link" == "cloglog" ) {	
		qui replace `lnf' = ln(1 - exp(-exp(`tau1' - `xb'))) if $ML_y == 1
		
		qui replace `lnf' = ln(exp(-exp(`tau1' - `xb'))) + ///
		                    ln(1 - exp(-exp(`tau2' - `xb'))) if $ML_y == 2
		
		if( $nCat > 3 ) {
			qui replace `lnf' = ln(exp(-exp(`tau1' - `xb'))) + ///
			                    ln(exp(-exp(`tau2' - `xb'))) + ///
		                        ln(1 - exp(-exp(`tau3' - `xb'))) if $ML_y == 3
		}
		
		if( $nCat > 4 ) {
			qui replace `lnf' = ln(exp(-exp(`tau1' - `xb'))) + ///
			                    ln(exp(-exp(`tau2' - `xb'))) + ///
								ln(exp(-exp(`tau3' - `xb'))) + ///
		                        ln(1 - exp(-exp(`tau4' - `xb'))) if $ML_y == 4
		}
		
		if( $nCat > 5 ) {
			qui replace `lnf' = ln(exp(-exp(`tau1' - `xb'))) + ///
			                    ln(exp(-exp(`tau2' - `xb'))) + ///
								ln(exp(-exp(`tau3' - `xb'))) + ///
								ln(exp(-exp(`tau4' - `xb'))) + ///
		                        ln(1 - exp(-exp(`tau5' - `xb'))) if $ML_y == 5
		}
		
		if( $nCat > 6 ) {
			qui replace `lnf' = ln(exp(-exp(`tau1' - `xb'))) + ///
			                    ln(exp(-exp(`tau2' - `xb'))) + ///
								ln(exp(-exp(`tau3' - `xb'))) + ///
								ln(exp(-exp(`tau4' - `xb'))) + ///
								ln(exp(-exp(`tau5' - `xb'))) + ///
		                        ln(1 - exp(-exp(`tau6' - `xb'))) if $ML_y == 6
		}
		
		if( $nCat > 7 ) {
			qui replace `lnf' = ln(exp(-exp(`tau1' - `xb'))) + ///
			                    ln(exp(-exp(`tau2' - `xb'))) + ///
								ln(exp(-exp(`tau3' - `xb'))) + ///
								ln(exp(-exp(`tau4' - `xb'))) + ///
								ln(exp(-exp(`tau5' - `xb'))) + ///
								ln(exp(-exp(`tau6' - `xb'))) + ///
		                        ln(1 - exp(-exp(`tau7' - `xb'))) if $ML_y == 7
		}
		
		if( $nCat > 8 ) {
			qui replace `lnf' = ln(exp(-exp(`tau1' - `xb'))) + ///
			                    ln(exp(-exp(`tau2' - `xb'))) + ///
								ln(exp(-exp(`tau3' - `xb'))) + ///
								ln(exp(-exp(`tau4' - `xb'))) + ///
								ln(exp(-exp(`tau5' - `xb'))) + ///
								ln(exp(-exp(`tau6' - `xb'))) + ///
								ln(exp(-exp(`tau7' - `xb'))) + ///
		                        ln(1 - exp(-exp(`tau8' - `xb'))) if $ML_y == 8
		}
		
		if( $nCat > 9 ) {
			qui replace `lnf' = ln(exp(-exp(`tau1' - `xb'))) + ///
			                    ln(exp(-exp(`tau2' - `xb'))) + ///
								ln(exp(-exp(`tau3' - `xb'))) + ///
								ln(exp(-exp(`tau4' - `xb'))) + ///
								ln(exp(-exp(`tau5' - `xb'))) + ///
								ln(exp(-exp(`tau6' - `xb'))) + ///
								ln(exp(-exp(`tau7' - `xb'))) + ///
								ln(exp(-exp(`tau8' - `xb'))) + ///
		                        ln(1 - exp(-exp(`tau9' - `xb'))) if $ML_y == 9
		}
		
		* build equation for last value of Y
		local eqn `" ln(exp(-exp(`tau1' - `xb'))) "'
		forval i = 2/$nCatm1 {
			local eqn `" `eqn' + ln(exp(-exp(`tau`i'' - `xb'))) "'
		}
		qui replace `lnf' = `eqn' if $ML_y == $nCat
	}
	
end


/* History
1.0.0  11.15.16  initial likelihood program for arbitrary number of categories
1.1.0  08.25.17  generalized program for unlimited number of categories
*/
