*! v1.1.0, S Bauldry, 25aug2017

capture program drop cr_pnp_lf
program cr_pnp_lf
  version 14
		
  * creating arguments based on number of categories of Y stored in globals
  local arguments "lnf xb_c"
		
  forval i = 1/$nCatm1 {
    local arguments "`arguments' xb_f`i'"
  }
	
  args `arguments'
	
  *** likelihood function for logit link
  if ( "$Link" == "logit" ) {	

    * equation for first value of Y
    qui replace `lnf' = ln(invlogit(-`xb_c' - `xb_f1')) if $ML_y == 1
		
	* build equations for middle values of Y
    forval k = 2/$nCatm1 {
      local meqn_b `" ln(invlogit(-`xb_c' - `xb_f`k'')) "'
    
	  local meqn_a ""
	  local m = `k' - 1
      forval n = 1/`m' {
        local meqn_a `" `meqn_a' ln(1 - invlogit(-`xb_c' - `xb_f`n'')) + "'
      }
	
      local meqn `" `meqn_a' `meqn_b' "'
      qui replace `lnf' = `meqn' if $ML_y == `k'
    }	
	
	* build equation for last value of Y
	local eqn `" ln(1 - invlogit(-`xb_c' - `xb_f1')) "'
	forval o = 2/$nCatm1 {
	  local eqn `" `eqn' + ln(1 - invlogit(-`xb_c' - `xb_f`o'')) "'
	}
	qui replace `lnf' = `eqn' if $ML_y == $nCat
  }

  *** likelihood function for probit link
  if ( "$Link" == "probit" ) {	

    * equation for first value of Y
    qui replace `lnf' = ln(normal(-`xb_c' - `xb_f1')) if $ML_y == 1
		
	* build equations for middle values of Y
    forval k = 2/$nCatm1 {
      local meqn_b `" ln(normal(-`xb_c' - `xb_f`k'')) "'
    
	  local meqn_a ""
	  local m = `k' - 1
      forval n = 1/`m' {
        local meqn_a `" `meqn_a' ln(1 - normal(-`xb_c' - `xb_f`n'')) + "'
      }
	
      local meqn `" `meqn_a' `meqn_b' "'
      qui replace `lnf' = `meqn' if $ML_y == `k'
    }	
	
	* build equation for last value of Y
	local eqn `" ln(1 - normal(-`xb_c' - `xb_f1')) "'
	forval o = 2/$nCatm1 {
	  local eqn `" `eqn' + ln(1 - normal(-`xb_c' - `xb_f`o'')) "'
	}
	qui replace `lnf' = `eqn' if $ML_y == $nCat
  }
  
  
  *** likelihood function for complementary log-log link
  if ( "$Link" == "cloglog" ) {	

    * equation for first value of Y
    qui replace `lnf' = ln(1 - exp(-exp(-`xb_c' - `xb_f1'))) if $ML_y == 1
		
	* build equations for middle values of Y
    forval k = 2/$nCatm1 {
      local meqn_b `" ln(1 - exp(-exp(-`xb_c' - `xb_f`k''))) "'
    
	  local meqn_a ""
	  local m = `k' - 1
      forval n = 1/`m' {
        local meqn_a `" `meqn_a' ln(exp(-exp(-`xb_c' - `xb_f`n''))) + "'
      }
	
      local meqn `" `meqn_a' `meqn_b' "'
      qui replace `lnf' = `meqn' if $ML_y == `k'
    }	
	
	* build equation for last value of Y
	local eqn `" ln(exp(-exp(-`xb_c' - `xb_f1'))) "'
	forval o = 2/$nCatm1 {
	  local eqn `" `eqn' + ln(exp(-exp(-`xb_c' - `xb_f`o''))) "'
	}
	qui replace `lnf' = `eqn' if $ML_y == $nCat
  }
  	
end


/* History
1.0.0  11.17.16  initial likelihood program for arbitrary number of categories
1.1.0  08.25.17  generalized program for unlimited number of categories
*/
