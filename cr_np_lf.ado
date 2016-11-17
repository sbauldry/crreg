*! v1.0.0, S Bauldry, 15nov2016

capture program drop cr_np_lf
program cr_np_lf
	version 14
		
	* creating arguments based on number of categories of Y stored in globals
	local arguments "lnf"
		
	forval i = 1/$nCatm1 {
		local arguments "`arguments' xb`i'"
	}
	
	args `arguments'
	
	* likelihood function for logit link
	if ( "$Link" == "logit" ) {	
		qui replace `lnf' = ln(invlogit(-`xb1')) if $ML_y == 1
		
		qui replace `lnf' = ln(1 - invlogit(-`xb1')) + ///
		                    ln(invlogit(-`xb2')) if $ML_y == 2
		
		if( $nCat > 3 ) {
			qui replace `lnf' = ln(1 - invlogit(-`xb1')) + ///
			                    ln(1 - invlogit(-`xb2')) + ///
		                        ln(invlogit(-`xb3')) if $ML_y == 3
		}
		
		if( $nCat > 4 ) {
			qui replace `lnf' = ln(1 - invlogit(-`xb1')) + ///
			                    ln(1 - invlogit(-`xb2')) + ///
								ln(1 - invlogit(-`xb3')) + ///
		                        ln(invlogit(-`xb4')) if $ML_y == 4
		}
		
		if( $nCat > 5 ) {
			qui replace `lnf' = ln(1 - invlogit(-`xb1')) + ///
			                    ln(1 - invlogit(-`xb2')) + ///
								ln(1 - invlogit(-`xb3')) + ///
								ln(1 - invlogit(-`xb4')) + ///
		                        ln(invlogit(-`xb5')) if $ML_y == 5
		}
		
		if( $nCat > 6 ) {
			qui replace `lnf' = ln(1 - invlogit(-`xb1')) + ///
			                    ln(1 - invlogit(-`xb2')) + ///
								ln(1 - invlogit(-`xb3')) + ///
								ln(1 - invlogit(-`xb4')) + ///
								ln(1 - invlogit(-`xb5')) + ///
		                        ln(invlogit(-`xb6')) if $ML_y == 6
		}
		
		if( $nCat > 7 ) {
			qui replace `lnf' = ln(1 - invlogit(-`xb1')) + ///
			                    ln(1 - invlogit(-`xb2')) + ///
								ln(1 - invlogit(-`xb3')) + ///
								ln(1 - invlogit(-`xb4')) + ///
								ln(1 - invlogit(-`xb5')) + ///
								ln(1 - invlogit(-`xb6')) + ///
		                        ln(invlogit(-`xb7')) if $ML_y == 7
		}
		
		if( $nCat > 8 ) {
			qui replace `lnf' = ln(1 - invlogit(-`xb1')) + ///
			                    ln(1 - invlogit(-`xb2')) + ///
								ln(1 - invlogit(-`xb3')) + ///
								ln(1 - invlogit(-`xb4')) + ///
								ln(1 - invlogit(-`xb5')) + ///
								ln(1 - invlogit(-`xb6')) + ///
								ln(1 - invlogit(-`xb7')) + ///
		                        ln(invlogit(-`xb8')) if $ML_y == 8
		}
		
		if( $nCat > 9 ) {
			qui replace `lnf' = ln(1 - invlogit(-`xb1')) + ///
			                    ln(1 - invlogit(-`xb2')) + ///
								ln(1 - invlogit(-`xb3')) + ///
								ln(1 - invlogit(-`xb4')) + ///
								ln(1 - invlogit(-`xb5')) + ///
								ln(1 - invlogit(-`xb6')) + ///
								ln(1 - invlogit(-`xb7')) + ///
								ln(1 - invlogit(-`xb8')) + ///
		                        ln(invlogit(-`xb9')) if $ML_y == 9
		}
		
		* build equation for last value of Y
		local eqn `" ln(1 - invlogit(-`xb1')) "'
		forval i = 2/$nCatm1 {
			local eqn `" `eqn' + ln(1 - invlogit(-`xb`i'')) "'
		}
		qui replace `lnf' = `eqn' if $ML_y == $nCat
	}
	
	
	* likelihood function for probit link
	if ( "$Link" == "probit" ) {	
		qui replace `lnf' = ln(normal(-`xb1')) if $ML_y == 1
		
		qui replace `lnf' = ln(1 - normal(-`xb1')) + ///
		                    ln(normal(-`xb2')) if $ML_y == 2
		
		if( $nCat > 3 ) {
			qui replace `lnf' = ln(1 - normal(-`xb1')) + ///
			                    ln(1 - normal(-`xb2')) + ///
		                        ln(normal(-`xb3')) if $ML_y == 3
		}
		
		if( $nCat > 4 ) {
			qui replace `lnf' = ln(1 - normal(-`xb1')) + ///
			                    ln(1 - normal(-`xb2')) + ///
								ln(1 - normal(-`xb3')) + ///
		                        ln(normal(-`xb4')) if $ML_y == 4
		}
		
		if( $nCat > 5 ) {
			qui replace `lnf' = ln(1 - normal(-`xb1')) + ///
			                    ln(1 - normal(-`xb2')) + ///
								ln(1 - normal(-`xb3')) + ///
								ln(1 - normal(-`xb4')) + ///
		                        ln(normal(-`xb5')) if $ML_y == 5
		}
		
		if( $nCat > 6 ) {
			qui replace `lnf' = ln(1 - normal(-`xb1')) + ///
			                    ln(1 - normal(-`xb2')) + ///
								ln(1 - normal(-`xb3')) + ///
								ln(1 - normal(-`xb4')) + ///
								ln(1 - normal(-`xb5')) + ///
		                        ln(normal(-`xb6')) if $ML_y == 6
		}
		
		if( $nCat > 7 ) {
			qui replace `lnf' = ln(1 - normal(-`xb1')) + ///
			                    ln(1 - normal(-`xb2')) + ///
								ln(1 - normal(-`xb3')) + ///
								ln(1 - normal(-`xb4')) + ///
								ln(1 - normal(-`xb5')) + ///
								ln(1 - normal(-`xb6')) + ///
		                        ln(normal(-`xb7')) if $ML_y == 7
		}
		
		if( $nCat > 8 ) {
			qui replace `lnf' = ln(1 - normal(-`xb1')) + ///
			                    ln(1 - normal(-`xb2')) + ///
								ln(1 - normal(-`xb3')) + ///
								ln(1 - normal(-`xb4')) + ///
								ln(1 - normal(-`xb5')) + ///
								ln(1 - normal(-`xb6')) + ///
								ln(1 - normal(-`xb7')) + ///
		                        ln(normal(-`xb8')) if $ML_y == 8
		}
		
		if( $nCat > 9 ) {
			qui replace `lnf' = ln(1 - normal(-`xb1')) + ///
			                    ln(1 - normal(-`xb2')) + ///
								ln(1 - normal(-`xb3')) + ///
								ln(1 - normal(-`xb4')) + ///
								ln(1 - normal(-`xb5')) + ///
								ln(1 - normal(-`xb6')) + ///
								ln(1 - normal(-`xb7')) + ///
								ln(1 - normal(-`xb8')) + ///
		                        ln(normal(-`xb9')) if $ML_y == 9
		}
		
		* build equation for last value of Y
		local eqn `" ln(1 - normal(-`xb1')) "'
		forval i = 2/$nCatm1 {
			local eqn `" `eqn' + ln(1 - normal(-`xb`i'')) "'
		}
		qui replace `lnf' = `eqn' if $ML_y == $nCat
	}
	
	
	* likelihood function for complementary log-log link
	if ( "$Link" == "cloglog" ) {	
		qui replace `lnf' = ln(exp(-exp(-`xb1')) if $ML_y == 1
		
		qui replace `lnf' = ln(1 - exp(-exp(-`xb1')) + ///
		                    ln(exp(-exp(-`xb2')) if $ML_y == 2
		
		if( $nCat > 3 ) {
			qui replace `lnf' = ln(1 - exp(-exp(-`xb1')) + ///
			                    ln(1 - exp(-exp(-`xb2')) + ///
		                        ln(exp(-exp(-`xb3')) if $ML_y == 3
		}
		
		if( $nCat > 4 ) {
			qui replace `lnf' = ln(1 - exp(-exp(-`xb1')) + ///
			                    ln(1 - exp(-exp(-`xb2')) + ///
								ln(1 - exp(-exp(-`xb3')) + ///
		                        ln(exp(-exp(-`xb4')) if $ML_y == 4
		}
		
		if( $nCat > 5 ) {
			qui replace `lnf' = ln(1 - exp(-exp(-`xb1')) + ///
			                    ln(1 - exp(-exp(-`xb2')) + ///
								ln(1 - exp(-exp(-`xb3')) + ///
								ln(1 - exp(-exp(-`xb4')) + ///
		                        ln(exp(-exp(-`xb5')) if $ML_y == 5
		}
		
		if( $nCat > 6 ) {
			qui replace `lnf' = ln(1 - exp(-exp(-`xb1')) + ///
			                    ln(1 - exp(-exp(-`xb2')) + ///
								ln(1 - exp(-exp(-`xb3')) + ///
								ln(1 - exp(-exp(-`xb4')) + ///
								ln(1 - exp(-exp(-`xb5')) + ///
		                        ln(exp(-exp(-`xb6')) if $ML_y == 6
		}
		
		if( $nCat > 7 ) {
			qui replace `lnf' = ln(1 - exp(-exp(-`xb1')) + ///
			                    ln(1 - exp(-exp(-`xb2')) + ///
								ln(1 - exp(-exp(-`xb3')) + ///
								ln(1 - exp(-exp(-`xb4')) + ///
								ln(1 - exp(-exp(-`xb5')) + ///
								ln(1 - exp(-exp(-`xb6')) + ///
		                        ln(exp(-exp(-`xb7')) if $ML_y == 7
		}
		
		if( $nCat > 8 ) {
			qui replace `lnf' = ln(1 - exp(-exp(-`xb1')) + ///
			                    ln(1 - exp(-exp(-`xb2')) + ///
								ln(1 - exp(-exp(-`xb3')) + ///
								ln(1 - exp(-exp(-`xb4')) + ///
								ln(1 - exp(-exp(-`xb5')) + ///
								ln(1 - exp(-exp(-`xb6')) + ///
								ln(1 - exp(-exp(-`xb7')) + ///
		                        ln(exp(-exp(-`xb8')) if $ML_y == 8
		}
		
		if( $nCat > 9 ) {
			qui replace `lnf' = ln(1 - exp(-exp(-`xb1')) + ///
			                    ln(1 - exp(-exp(-`xb2')) + ///
								ln(1 - exp(-exp(-`xb3')) + ///
								ln(1 - exp(-exp(-`xb4')) + ///
								ln(1 - exp(-exp(-`xb5')) + ///
								ln(1 - exp(-exp(-`xb6')) + ///
								ln(1 - exp(-exp(-`xb7')) + ///
								ln(1 - exp(-exp(-`xb8')) + ///
		                        ln(exp(-exp(-`xb9')) if $ML_y == 9
		}
		
		* build equation for last value of Y
		local eqn `" ln(1 - exp(-exp(-`xb1')) "'
		forval i = 2/$nCatm1 {
			local eqn `" `eqn' + ln(1 - exp(-exp(-`xb`i'')) "'
		}
		qui replace `lnf' = `eqn' if $ML_y == $nCat
	}
	
end


/* History
1.0.0  11.16.16  initial likelihood program for arbitrary number of categories

