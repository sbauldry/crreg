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
	
  *** likelihood function for logit link
  if ( "$Link" == "logit" ) {	
		
    * equation for first value of Y
    qui replace `lnf' = ln(invlogit(`tau1' - `xb')) if $ML_y == 1
		
    * build equations for middle values of Y
    forval k = 2/$nCatm1 {
      local meqn_b `" ln(invlogit(`tau`k'' - `xb')) "'
    
	  local meqn_a ""
	  local m = `k' - 1
      forval n = 1/`m' {
        local meqn_a `" `meqn_a' ln(1 - invlogit(`tau`n'' - `xb')) + "'
      }
	
      local meqn `" `meqn_a' `meqn_b' "'
      qui replace `lnf' = `meqn' if $ML_y == `k'
    }
	
	* build equation for last value of Y
	local eqn `" ln(1 - invlogit(`tau1' - `xb')) "'
	forval o = 2/$nCatm1 {
      local eqn `" `eqn' + ln(1 - invlogit(`tau`o'' - `xb')) "'
    }
	qui replace `lnf' = `eqn' if $ML_y == $nCat
  }
  
  *** likelihood function for probit link
  if ( "$Link" == "probit" ) {	
		
    * equation for first value of Y
    qui replace `lnf' = ln(normal(`tau1' - `xb')) if $ML_y == 1
		
    * build equations for middle values of Y
    forval k = 2/$nCatm1 {
      local meqn_b `" ln(normal(`tau`k'' - `xb')) "'
    
	  local meqn_a ""
	  local m = `k' - 1
      forval n = 1/`m' {
        local meqn_a `" `meqn_a' ln(1 - normal(`tau`n'' - `xb')) + "'
      }
	
      local meqn `" `meqn_a' `meqn_b' "'
      qui replace `lnf' = `meqn' if $ML_y == `k'
    }
	
	* build equation for last value of Y
	local eqn `" ln(1 - normal(`tau1' - `xb')) "'
	forval o = 2/$nCatm1 {
      local eqn `" `eqn' + ln(1 - normal(`tau`o'' - `xb')) "'
    }
	qui replace `lnf' = `eqn' if $ML_y == $nCat
  }
  
  *** likelihood function for complementary log-log link
  if ( "$Link" == "cloglog" ) {	
		
    * equation for first value of Y
    qui replace `lnf' = ln(1 - exp(-exp(`tau1' - `xb')) if $ML_y == 1
		
    * build equations for middle values of Y
    forval k = 2/$nCatm1 {
      local meqn_b `" ln(1 - exp(-exp(`tau`k'' - `xb')) "'
    
	  local meqn_a ""
	  local m = `k' - 1
      forval n = 1/`m' {
        local meqn_a `" `meqn_a' ln(exp(-exp(`tau`n'' - `xb')) + "'
      }
	
      local meqn `" `meqn_a' `meqn_b' "'
      qui replace `lnf' = `meqn' if $ML_y == `k'
    }
	
	* build equation for last value of Y
	local eqn `" ln(exp(-exp(`tau1' - `xb')) "'
	forval o = 2/$nCatm1 {
      local eqn `" `eqn' + ln(exp(-exp(`tau`o'' - `xb')) "'
    }
	qui replace `lnf' = `eqn' if $ML_y == $nCat
  }
	
end


/* History
1.0.0  11.15.16  initial likelihood program for arbitrary number of categories
1.1.0  08.25.17  generalized program for unlimited number of categories
*/
