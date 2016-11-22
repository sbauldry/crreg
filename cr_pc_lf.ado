*! v1.0.0, S Bauldry, 22nov2016

capture program drop cr_pc_lf
program cr_pc_lf
	version 14

	* creating arguments based on number of categories of Y stored in globals
	local arguments "lnf xb_c xb_p"
		
	forval i = 2/$nCatm1 {
		local arguments "`arguments' f`i'"
	}
	
	args `arguments'
	
	forval i = 2/$nCatm1 {
		tempvar phi`i'
		qui gen double `phi`i'' = `f`i''
	}
	
	* likelihood function for logit link
	if ( "$Link" == "logit" ) {	
		qui replace `lnf' = ln(invlogit(-`xb_c' - `xb_p')) if $ML_y == 1
		
		qui replace `lnf' = ln(1 - invlogit(-`xb_c' - `xb_p')) +  ///
		                    ln(invlogit(-`xb_c' - `xb_p'*`phi2')) if $ML_y == 2
		
		if( $nCat > 3 ) {
			qui replace `lnf' = ln(1 - invlogit(-`xb_c' - `xb_p')) + ///
			                    ln(1 - invlogit(-`xb_c' - `xb_p'*`phi2')) + ///
		                        ln(invlogit(-`xb_c' - `xb_p'*`phi3')) if $ML_y == 3
		}
		
		if( $nCat > 4 ) {
			qui replace `lnf' = ln(1 - invlogit(-`xb_c' - `xb_p')) + ///
			                    ln(1 - invlogit(-`xb_c' - `xb_p'*`phi2')) + ///
								ln(1 - invlogit(-`xb_c' - `xb_p'*`phi3')) + ///
		                        ln(invlogit(-`xb_c' - `xb_p'*`phi4')) if $ML_y == 4
		}
		
		if( $nCat > 5 ) {
			qui replace `lnf' = ln(1 - invlogit(-`xb_c' - `xb_p')) + ///
			                    ln(1 - invlogit(-`xb_c' - `xb_p'*`phi2')) + ///
								ln(1 - invlogit(-`xb_c' - `xb_p'*`phi3')) + ///
								ln(1 - invlogit(-`xb_c' - `xb_p'*`phi4')) + ///
		                        ln(invlogit(-`xb_c' - `xb_p'*`phi5')) if $ML_y == 5
		}
		
		if( $nCat > 6 ) {
			qui replace `lnf' = ln(1 - invlogit(-`xb_c' - `xb_p')) + ///
			                    ln(1 - invlogit(-`xb_c' - `xb_p'*`phi2')) + ///
								ln(1 - invlogit(-`xb_c' - `xb_p'*`phi3')) + ///
								ln(1 - invlogit(-`xb_c' - `xb_p'*`phi4')) + ///
								ln(1 - invlogit(-`xb_c' - `xb_p'*`phi5')) + ///
		                        ln(invlogit(-`xb_c' - `xb_p'*`phi6')) if $ML_y == 6
		}
		
		if( $nCat > 7 ) {
			qui replace `lnf' = ln(1 - invlogit(-`xb_c' - `xb_p')) + ///
			                    ln(1 - invlogit(-`xb_c' - `xb_p'*`phi2')) + ///
								ln(1 - invlogit(-`xb_c' - `xb_p'*`phi3')) + ///
								ln(1 - invlogit(-`xb_c' - `xb_p'*`phi4')) + ///
								ln(1 - invlogit(-`xb_c' - `xb_p'*`phi5')) + ///
								ln(1 - invlogit(-`xb_c' - `xb_p'*`phi6')) + ///
		                        ln(invlogit(-`xb_c' - `xb_p'*`phi7')) if $ML_y == 7
		}
		
		if( $nCat > 8 ) {
			qui replace `lnf' = ln(1 - invlogit(-`xb_c' - `xb_p')) + ///
			                    ln(1 - invlogit(-`xb_c' - `xb_p'*`phi2')) + ///
								ln(1 - invlogit(-`xb_c' - `xb_p'*`phi3')) + ///
								ln(1 - invlogit(-`xb_c' - `xb_p'*`phi4')) + ///
								ln(1 - invlogit(-`xb_c' - `xb_p'*`phi5')) + ///
								ln(1 - invlogit(-`xb_c' - `xb_p'*`phi6')) + ///
								ln(1 - invlogit(-`xb_c' - `xb_p'*`phi7')) + ///
		                        ln(invlogit(-`xb_c' - `xb_p'*`phi8')) if $ML_y == 8
		}
		
		if( $nCat > 9 ) {
			qui replace `lnf' = ln(1 - invlogit(-`xb_c' - `xb_p')) + ///
			                    ln(1 - invlogit(-`xb_c' - `xb_p'*`phi2')) + ///
								ln(1 - invlogit(-`xb_c' - `xb_p'*`phi3')) + ///
								ln(1 - invlogit(-`xb_c' - `xb_p'*`phi4')) + ///
								ln(1 - invlogit(-`xb_c' - `xb_p'*`phi5')) + ///
								ln(1 - invlogit(-`xb_c' - `xb_p'*`phi6')) + ///
								ln(1 - invlogit(-`xb_c' - `xb_p'*`phi7')) + ///
								ln(1 - invlogit(-`xb_c' - `xb_p'*`phi8')) + ///
		                        ln(invlogit(-`xb_c' - `xb_p'*`phi9')) if $ML_y == 9
		}
		
		* build equation for last value of Y
		local eqn `" ln(1 - invlogit(-`xb_c' - `xb_p')) "'
		forval i = 2/$nCatm1 {
			local eqn `" `eqn' + ln(1 - invlogit(-`xb_c' - `xb_p'*`phi`i'')) "'
		}
		qui replace `lnf' = `eqn' if $ML_y == $nCat
	}
	
	
	* likelihood function for probit link
	if ( "$Link" == "probit" ) {	
		qui replace `lnf' = ln(normal(-`xb_c' - `xb_p')) if $ML_y == 1
		
		qui replace `lnf' = ln(1 - normal(-`xb_c' - `xb_p')) +  ///
		                    ln(normal(-`xb_c' - `xb_p'*`phi2')) if $ML_y == 2
		
		if( $nCat > 3 ) {
			qui replace `lnf' = ln(1 - normal(-`xb_c' - `xb_p')) + ///
			                    ln(1 - normal(-`xb_c' - `xb_p'*`phi2')) + ///
		                        ln(normal(-`xb_c' - `xb_p'*`phi3')) if $ML_y == 3
		}
		
		if( $nCat > 4 ) {
			qui replace `lnf' = ln(1 - normal(-`xb_c' - `xb_p')) + ///
			                    ln(1 - normal(-`xb_c' - `xb_p'*`phi2')) + ///
								ln(1 - normal(-`xb_c' - `xb_p'*`phi3')) + ///
		                        ln(normal(-`xb_c' - `xb_p'*`phi4')) if $ML_y == 4
		}
		
		if( $nCat > 5 ) {
			qui replace `lnf' = ln(1 - normal(-`xb_c' - `xb_p')) + ///
			                    ln(1 - normal(-`xb_c' - `xb_p'*`phi2')) + ///
								ln(1 - normal(-`xb_c' - `xb_p'*`phi3')) + ///
								ln(1 - normal(-`xb_c' - `xb_p'*`phi4')) + ///
		                        ln(normal(-`xb_c' - `xb_p'*`phi5')) if $ML_y == 5
		}
		
		if( $nCat > 6 ) {
			qui replace `lnf' = ln(1 - normal(-`xb_c' - `xb_p')) + ///
			                    ln(1 - normal(-`xb_c' - `xb_p'*`phi2')) + ///
								ln(1 - normal(-`xb_c' - `xb_p'*`phi3')) + ///
								ln(1 - normal(-`xb_c' - `xb_p'*`phi4')) + ///
								ln(1 - normal(-`xb_c' - `xb_p'*`phi5')) + ///
		                        ln(normal(-`xb_c' - `xb_p'*`phi6')) if $ML_y == 6
		}
		
		if( $nCat > 7 ) {
			qui replace `lnf' = ln(1 - normal(-`xb_c' - `xb_p')) + ///
			                    ln(1 - normal(-`xb_c' - `xb_p'*`phi2')) + ///
								ln(1 - normal(-`xb_c' - `xb_p'*`phi3')) + ///
								ln(1 - normal(-`xb_c' - `xb_p'*`phi4')) + ///
								ln(1 - normal(-`xb_c' - `xb_p'*`phi5')) + ///
								ln(1 - normal(-`xb_c' - `xb_p'*`phi6')) + ///
		                        ln(normal(-`xb_c' - `xb_p'*`phi7')) if $ML_y == 7
		}
		
		if( $nCat > 8 ) {
			qui replace `lnf' = ln(1 - normal(-`xb_c' - `xb_p')) + ///
			                    ln(1 - normal(-`xb_c' - `xb_p'*`phi2')) + ///
								ln(1 - normal(-`xb_c' - `xb_p'*`phi3')) + ///
								ln(1 - normal(-`xb_c' - `xb_p'*`phi4')) + ///
								ln(1 - normal(-`xb_c' - `xb_p'*`phi5')) + ///
								ln(1 - normal(-`xb_c' - `xb_p'*`phi6')) + ///
								ln(1 - normal(-`xb_c' - `xb_p'*`phi7')) + ///
		                        ln(normal(-`xb_c' - `xb_p'*`phi8')) if $ML_y == 8
		}
		
		if( $nCat > 9 ) {
			qui replace `lnf' = ln(1 - normal(-`xb_c' - `xb_p')) + ///
			                    ln(1 - normal(-`xb_c' - `xb_p'*`phi2')) + ///
								ln(1 - normal(-`xb_c' - `xb_p'*`phi3')) + ///
								ln(1 - normal(-`xb_c' - `xb_p'*`phi4')) + ///
								ln(1 - normal(-`xb_c' - `xb_p'*`phi5')) + ///
								ln(1 - normal(-`xb_c' - `xb_p'*`phi6')) + ///
								ln(1 - normal(-`xb_c' - `xb_p'*`phi7')) + ///
								ln(1 - normal(-`xb_c' - `xb_p'*`phi8')) + ///
		                        ln(normal(-`xb_c' - `xb_p'*`phi9')) if $ML_y == 9
		}
		
		* build equation for last value of Y
		local eqn `" ln(1 - normal(-`xb_c' - `xb_p')) "'
		forval i = 2/$nCatm1 {
			local eqn `" `eqn' + ln(1 - normal(-`xb_c' - `xb_p'*`phi`i'')) "'
		}
		qui replace `lnf' = `eqn' if $ML_y == $nCat
	}
	
	
	* likelihood function for complementary log-log link
	if ( "$Link" == "cloglog" ) {	
		qui replace `lnf' = ln(exp(-exp(-`xb_c' - `xb_p'))) if $ML_y == 1
		
		qui replace `lnf' = ln(1 - exp(-exp(-`xb_c' - `xb_p'))) +  ///
		                    ln(exp(-exp(-`xb_c' - `xb_p'*`phi2'))) if $ML_y == 2
		
		if( $nCat > 3 ) {
			qui replace `lnf' = ln(1 - exp(-exp(-`xb_c' - `xb_p'))) + ///
			                    ln(1 - exp(-exp(-`xb_c' - `xb_p'*`phi2'))) + ///
		                        ln(exp(-exp(-`xb_c' - `xb_p'*`phi3'))) if $ML_y == 3
		}
		
		if( $nCat > 4 ) {
			qui replace `lnf' = ln(1 - exp(-exp(-`xb_c' - `xb_p'))) + ///
			                    ln(1 - exp(-exp(-`xb_c' - `xb_p'*`phi2'))) + ///
								ln(1 - exp(-exp(-`xb_c' - `xb_p'*`phi3'))) + ///
		                        ln(exp(-exp(-`xb_c' - `xb_p'*`phi4'))) if $ML_y == 4
		}
		
		if( $nCat > 5 ) {
			qui replace `lnf' = ln(1 - exp(-exp(-`xb_c' - `xb_p'))) + ///
			                    ln(1 - exp(-exp(-`xb_c' - `xb_p'*`phi2'))) + ///
								ln(1 - exp(-exp(-`xb_c' - `xb_p'*`phi3'))) + ///
								ln(1 - exp(-exp(-`xb_c' - `xb_p'*`phi4'))) + ///
		                        ln(exp(-exp(-`xb_c' - `xb_p'*`phi5'))) if $ML_y == 5
		}
		
		if( $nCat > 6 ) {
			qui replace `lnf' = ln(1 - exp(-exp(-`xb_c' - `xb_p'))) + ///
			                    ln(1 - exp(-exp(-`xb_c' - `xb_p'*`phi2'))) + ///
								ln(1 - exp(-exp(-`xb_c' - `xb_p'*`phi3'))) + ///
								ln(1 - exp(-exp(-`xb_c' - `xb_p'*`phi4'))) + ///
								ln(1 - exp(-exp(-`xb_c' - `xb_p'*`phi5'))) + ///
		                        ln(exp(-exp(-`xb_c' - `xb_p'*`phi6'))) if $ML_y == 6
		}
		
		if( $nCat > 7 ) {
			qui replace `lnf' = ln(1 - exp(-exp(-`xb_c' - `xb_p'))) + ///
			                    ln(1 - exp(-exp(-`xb_c' - `xb_p'*`phi2'))) + ///
								ln(1 - exp(-exp(-`xb_c' - `xb_p'*`phi3'))) + ///
								ln(1 - exp(-exp(-`xb_c' - `xb_p'*`phi4'))) + ///
								ln(1 - exp(-exp(-`xb_c' - `xb_p'*`phi5'))) + ///
								ln(1 - exp(-exp(-`xb_c' - `xb_p'*`phi6'))) + ///
		                        ln(exp(-exp(-`xb_c' - `xb_p'*`phi7'))) if $ML_y == 7
		}
		
		if( $nCat > 8 ) {
			qui replace `lnf' = ln(1 - exp(-exp(-`xb_c' - `xb_p'))) + ///
			                    ln(1 - exp(-exp(-`xb_c' - `xb_p'*`phi2'))) + ///
								ln(1 - exp(-exp(-`xb_c' - `xb_p'*`phi3'))) + ///
								ln(1 - exp(-exp(-`xb_c' - `xb_p'*`phi4'))) + ///
								ln(1 - exp(-exp(-`xb_c' - `xb_p'*`phi5'))) + ///
								ln(1 - exp(-exp(-`xb_c' - `xb_p'*`phi6'))) + ///
								ln(1 - exp(-exp(-`xb_c' - `xb_p'*`phi7'))) + ///
		                        ln(exp(-exp(-`xb_c' - `xb_p'*`phi8'))) if $ML_y == 8
		}
		
		if( $nCat > 9 ) {
			qui replace `lnf' = ln(1 - exp(-exp(-`xb_c' - `xb_p'))) + ///
			                    ln(1 - exp(-exp(-`xb_c' - `xb_p'*`phi2'))) + ///
								ln(1 - exp(-exp(-`xb_c' - `xb_p'*`phi3'))) + ///
								ln(1 - exp(-exp(-`xb_c' - `xb_p'*`phi4'))) + ///
								ln(1 - exp(-exp(-`xb_c' - `xb_p'*`phi5'))) + ///
								ln(1 - exp(-exp(-`xb_c' - `xb_p'*`phi6'))) + ///
								ln(1 - exp(-exp(-`xb_c' - `xb_p'*`phi7'))) + ///
								ln(1 - exp(-exp(-`xb_c' - `xb_p'*`phi8'))) + ///
		                        ln(exp(-exp(-`xb_c' - `xb_p'*`phi9'))) if $ML_y == 9
		}
		
		* build equation for last value of Y
		local eqn `" ln(1 - exp(-exp(-`xb_c' - `xb_p'))) "'
		forval i = 2/$nCatm1 {
			local eqn `" `eqn' + ln(1 - exp(-exp(-`xb_c' - `xb_p'*`phi`i''))) "'
		}
		qui replace `lnf' = `eqn' if $ML_y == $nCat
	}
	
end


/* History
1.0.0  11.22.16  initial likelihood program for arbitrary number of categories

