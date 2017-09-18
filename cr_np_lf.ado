*! v1.2.0, S Bauldry, 18sep2017

capture program drop cr_np_lf
program cr_np_lf
  version 14
		
  * creating arguments based on number of categories of Y stored in globals
  local arguments "lnf"
		
  forval i = 1/$nCatm1 {
    local arguments "`arguments' xb`i'"
  }
	
  args `arguments'
  
  * setting values for y
  forval j = 1/$nCat {
    local y_`j' ${y_`j'}
  }
  local M $nCat
	
  *** likelihood function for logit link
  if ( "$Link" == "logit" ) {	
    
	* equation for first value of Y
	qui replace `lnf' = ln(invlogit(-`xb1')) if $ML_y == `y_1'
		
	* build equations for middle value of Y
	forval k = 2/$nCatm1 {
      local meqn_b `" ln(invlogit(- `xb`k'')) "'
    
	  local meqn_a ""
	  local m = `k' - 1
      forval n = 1/`m' {
        local meqn_a `" `meqn_a' ln(1 - invlogit(- `xb`n'')) + "'
      }
	
      local meqn `" `meqn_a' `meqn_b' "'
      qui replace `lnf' = `meqn' if $ML_y == `y_`k''
    }
	
	* build equation for last value of Y
	local eqn `" ln(1 - invlogit(-`xb1')) "'
	forval o = 2/$nCatm1 {
      local eqn `" `eqn' + ln(1 - invlogit(-`xb`o'')) "'
	}
	qui replace `lnf' = `eqn' if $ML_y == `y_`M''
  }
  
  *** likelihood function for probit link
  if ( "$Link" == "probit" ) {	
    
	* equation for first value of Y
	qui replace `lnf' = ln(normal(-`xb1')) if $ML_y == `y_1'
		
	* build equations for middle value of Y
	forval k = 2/$nCatm1 {
      local meqn_b `" ln(normal(- `xb`k'')) "'
    
	  local meqn_a ""
	  local m = `k' - 1
      forval n = 1/`m' {
        local meqn_a `" `meqn_a' ln(1 - normal(- `xb`n'')) + "'
      }
	
      local meqn `" `meqn_a' `meqn_b' "'
      qui replace `lnf' = `meqn' if $ML_y == `y_`k''
    }
	
	* build equation for last value of Y
	local eqn `" ln(1 - normal(-`xb1')) "'
	forval o = 2/$nCatm1 {
      local eqn `" `eqn' + ln(1 - normal(-`xb`o'')) "'
	}
	qui replace `lnf' = `eqn' if $ML_y == `y_`M''
  }
	
  *** likelihood function for complementary log-log link
  if ( "$Link" == "cloglog" ) {	
    
	* equation for first value of Y
	qui replace `lnf' = ln(1 - exp(-exp(-`xb1'))) if $ML_y == `y_1'
		
	* build equations for middle value of Y
	forval k = 2/$nCatm1 {
      local meqn_b `" ln(1 - exp(-exp(- `xb`k''))) "'
    
	  local meqn_a ""
	  local m = `k' - 1
      forval n = 1/`m' {
        local meqn_a `" `meqn_a' ln(exp(-exp(- `xb`n''))) + "'
      }
	
      local meqn `" `meqn_a' `meqn_b' "'
      qui replace `lnf' = `meqn' if $ML_y == `y_`k''
    }
	
	* build equation for last value of Y
	local eqn `" ln(exp(-exp(-`xb1'))) "'
	forval o = 2/$nCatm1 {
      local eqn `" `eqn' + ln(exp(-exp(-`xb`o''))) "'
	}
	qui replace `lnf' = `eqn' if $ML_y == `y_`M''
  }
	
end


/* History
1.0.0  11.16.16  initial likelihood program for arbitrary number of categories
1.1.0  08.25.17  generalized program for unlimited number of categories
1.2.0  09.18.17  fixed bug with non-standard values for Y
*/
