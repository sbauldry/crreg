*! v1.0.0, S Bauldry, 15nov2016

capture program drop cr_pnp_lf
program cr_pnp_lf
	version 14
		
	* creating arguments based on number of categories of Y stored in globals
	local arguments "lnf xb_c"
		
	forval i = 1/$nCatm1 {
		local arguments "`arguments' xb_f`i'"
	}
	
	args `arguments'
	
	* likelihood function for logit link
	if ( "$Link" == "logit" ) {	
		qui replace `lnf' = ln(invlogit(-`xb_c' - `xb_f1')) if $ML_y == 1
		
		qui replace `lnf' = ln(1 - invlogit(-`xb_c' - `xb_f1')) + ///
		                    ln(invlogit(-`xb_c' - `xb_f2')) if $ML_y == 2
		
		if( $nCat > 3 ) {
			qui replace `lnf' = ln(1 - invlogit(-`xb_c' - `xb_f1')) + ///
			                    ln(1 - invlogit(-`xb_c' - `xb_f2')) + ///
		                        ln(invlogit(-`xb_c' - `xb_f3')) if $ML_y == 3
		}
		
		if( $nCat > 4 ) {
			qui replace `lnf' = ln(1 - invlogit(-`xb_c' - `xb_f1')) + ///
			                    ln(1 - invlogit(-`xb_c' - `xb_f2')) + ///
								ln(1 - invlogit(-`xb_c' - `xb_f3')) + ///
		                        ln(invlogit(-`xb_c' - `xb_f4')) if $ML_y == 4
		}
		
		if( $nCat > 5 ) {
			qui replace `lnf' = ln(1 - invlogit(-`xb_c' - `xb_f1')) + ///
			                    ln(1 - invlogit(-`xb_c' - `xb_f2')) + ///
								ln(1 - invlogit(-`xb_c' - `xb_f3')) + ///
								ln(1 - invlogit(-`xb_c' - `xb_f4')) + ///
		                        ln(invlogit(-`xb_c' - `xb_f5')) if $ML_y == 5
		}
		
		if( $nCat > 6 ) {
			qui replace `lnf' = ln(1 - invlogit(-`xb_c' - `xb_f1')) + ///
			                    ln(1 - invlogit(-`xb_c' - `xb_f2')) + ///
								ln(1 - invlogit(-`xb_c' - `xb_f3')) + ///
								ln(1 - invlogit(-`xb_c' - `xb_f4')) + ///
								ln(1 - invlogit(-`xb_c' - `xb_f5')) + ///
		                        ln(invlogit(-`xb_c' - `xb_f6')) if $ML_y == 6
		}
		
		if( $nCat > 7 ) {
			qui replace `lnf' = ln(1 - invlogit(-`xb_c' - `xb_f1')) + ///
			                    ln(1 - invlogit(-`xb_c' - `xb_f2')) + ///
								ln(1 - invlogit(-`xb_c' - `xb_f3')) + ///
								ln(1 - invlogit(-`xb_c' - `xb_f4')) + ///
								ln(1 - invlogit(-`xb_c' - `xb_f5')) + ///
								ln(1 - invlogit(-`xb_c' - `xb_f6')) + ///
		                        ln(invlogit(-`xb_c' - `xb_f7')) if $ML_y == 7
		}
		
		if( $nCat > 8 ) {
			qui replace `lnf' = ln(1 - invlogit(-`xb_c' - `xb_f1')) + ///
			                    ln(1 - invlogit(-`xb_c' - `xb_f2')) + ///
								ln(1 - invlogit(-`xb_c' - `xb_f3')) + ///
								ln(1 - invlogit(-`xb_c' - `xb_f4')) + ///
								ln(1 - invlogit(-`xb_c' - `xb_f5')) + ///
								ln(1 - invlogit(-`xb_c' - `xb_f6')) + ///
								ln(1 - invlogit(-`xb_c' - `xb_f7')) + ///
		                        ln(invlogit(-`xb_c' - `xb_f8')) if $ML_y == 8
		}
		
		if( $nCat > 9 ) {
			qui replace `lnf' = ln(1 - invlogit(-`xb_c' - `xb_f1')) + ///
			                    ln(1 - invlogit(-`xb_c' - `xb_f2')) + ///
								ln(1 - invlogit(-`xb_c' - `xb_f3')) + ///
								ln(1 - invlogit(-`xb_c' - `xb_f4')) + ///
								ln(1 - invlogit(-`xb_c' - `xb_f5')) + ///
								ln(1 - invlogit(-`xb_c' - `xb_f6')) + ///
								ln(1 - invlogit(-`xb_c' - `xb_f7')) + ///
								ln(1 - invlogit(-`xb_c' - `xb_f8')) + ///
		                        ln(invlogit(-`xb_c' - `xb_f9')) if $ML_y == 9
		}
		
		* build equation for last value of Y
		local eqn `" ln(1 - invlogit(-`xb_c' - `xb_f1')) "'
		forval i = 2/$nCatm1 {
			local eqn `" `eqn' + ln(1 - invlogit(-`xb_c' - `xb_f`i'')) "'
		}
		qui replace `lnf' = `eqn' if $ML_y == $nCat
	}
	
	
	* likelihood function for probit link
	if ( "$Link" == "probit" ) {	
		qui replace `lnf' = ln(normal(-`xb_c' - `xb_f1')) if $ML_y == 1
		
		qui replace `lnf' = ln(1 - normal(-`xb_c' - `xb_f1')) + ///
		                    ln(normal(-`xb_c' - `xb_f2')) if $ML_y == 2
		
		if( $nCat > 3 ) {
			qui replace `lnf' = ln(1 - normal(-`xb_c' - `xb_f1')) + ///
			                    ln(1 - normal(-`xb_c' - `xb_f2')) + ///
		                        ln(normal(-`xb_c' - `xb_f3')) if $ML_y == 3
		}
		
		if( $nCat > 4 ) {
			qui replace `lnf' = ln(1 - normal(-`xb_c' - `xb_f1')) + ///
			                    ln(1 - normal(-`xb_c' - `xb_f2')) + ///
								ln(1 - normal(-`xb_c' - `xb_f3')) + ///
		                        ln(normal(-`xb_c' - `xb_f4')) if $ML_y == 4
		}
		
		if( $nCat > 5 ) {
			qui replace `lnf' = ln(1 - normal(-`xb_c' - `xb_f1')) + ///
			                    ln(1 - normal(-`xb_c' - `xb_f2')) + ///
								ln(1 - normal(-`xb_c' - `xb_f3')) + ///
								ln(1 - normal(-`xb_c' - `xb_f4')) + ///
		                        ln(normal(-`xb_c' - `xb_f5')) if $ML_y == 5
		}
		
		if( $nCat > 6 ) {
			qui replace `lnf' = ln(1 - normal(-`xb_c' - `xb_f1')) + ///
			                    ln(1 - normal(-`xb_c' - `xb_f2')) + ///
								ln(1 - normal(-`xb_c' - `xb_f3')) + ///
								ln(1 - normal(-`xb_c' - `xb_f4')) + ///
								ln(1 - normal(-`xb_c' - `xb_f5')) + ///
		                        ln(normal(-`xb_c' - `xb_f6')) if $ML_y == 6
		}
		
		if( $nCat > 7 ) {
			qui replace `lnf' = ln(1 - normal(-`xb_c' - `xb_f1')) + ///
			                    ln(1 - normal(-`xb_c' - `xb_f2')) + ///
								ln(1 - normal(-`xb_c' - `xb_f3')) + ///
								ln(1 - normal(-`xb_c' - `xb_f4')) + ///
								ln(1 - normal(-`xb_c' - `xb_f5')) + ///
								ln(1 - normal(-`xb_c' - `xb_f6')) + ///
		                        ln(normal(-`xb_c' - `xb_f7')) if $ML_y == 7
		}
		
		if( $nCat > 8 ) {
			qui replace `lnf' = ln(1 - normal(-`xb_c' - `xb_f1')) + ///
			                    ln(1 - normal(-`xb_c' - `xb_f2')) + ///
								ln(1 - normal(-`xb_c' - `xb_f3')) + ///
								ln(1 - normal(-`xb_c' - `xb_f4')) + ///
								ln(1 - normal(-`xb_c' - `xb_f5')) + ///
								ln(1 - normal(-`xb_c' - `xb_f6')) + ///
								ln(1 - normal(-`xb_c' - `xb_f7')) + ///
		                        ln(normal(-`xb_c' - `xb_f8')) if $ML_y == 8
		}
		
		if( $nCat > 9 ) {
			qui replace `lnf' = ln(1 - normal(-`xb_c' - `xb_f1')) + ///
			                    ln(1 - normal(-`xb_c' - `xb_f2')) + ///
								ln(1 - normal(-`xb_c' - `xb_f3')) + ///
								ln(1 - normal(-`xb_c' - `xb_f4')) + ///
								ln(1 - normal(-`xb_c' - `xb_f5')) + ///
								ln(1 - normal(-`xb_c' - `xb_f6')) + ///
								ln(1 - normal(-`xb_c' - `xb_f7')) + ///
								ln(1 - normal(-`xb_c' - `xb_f8')) + ///
		                        ln(normal(-`xb_c' - `xb_f9')) if $ML_y == 9
		}
		
		* build equation for last value of Y
		local eqn `" ln(1 - normal(-`xb_c' - `xb_f1')) "'
		forval i = 2/$nCatm1 {
			local eqn `" `eqn' + ln(1 - normal(-`xb_c' - `xb_f`i'')) "'
		}
		qui replace `lnf' = `eqn' if $ML_y == $nCat
	}
	
	
	* likelihood function for complementary log-log link
	if ( "$Link" == "cloglog" ) {	
		qui replace `lnf' = ln(exp(-exp(-`xb_c' - `xb_f1')) if $ML_y == 1
		
		qui replace `lnf' = ln(1 - exp(-exp(-`xb_c' - `xb_f1')) + ///
		                    ln(exp(-exp(-`xb_c' - `xb_f2')) if $ML_y == 2
		
		if( $nCat > 3 ) {
			qui replace `lnf' = ln(1 - exp(-exp(-`xb_c' - `xb_f1')) + ///
			                    ln(1 - exp(-exp(-`xb_c' - `xb_f2')) + ///
		                        ln(exp(-exp(-`xb_c' - `xb_f3')) if $ML_y == 3
		}
		
		if( $nCat > 4 ) {
			qui replace `lnf' = ln(1 - exp(-exp(-`xb_c' - `xb_f1')) + ///
			                    ln(1 - exp(-exp(-`xb_c' - `xb_f2')) + ///
								ln(1 - exp(-exp(-`xb_c' - `xb_f3')) + ///
		                        ln(exp(-exp(-`xb_c' - `xb_f4')) if $ML_y == 4
		}
		
		if( $nCat > 5 ) {
			qui replace `lnf' = ln(1 - exp(-exp(-`xb_c' - `xb_f1')) + ///
			                    ln(1 - exp(-exp(-`xb_c' - `xb_f2')) + ///
								ln(1 - exp(-exp(-`xb_c' - `xb_f3')) + ///
								ln(1 - exp(-exp(-`xb_c' - `xb_f4')) + ///
		                        ln(exp(-exp(-`xb_c' - `xb_f5')) if $ML_y == 5
		}
		
		if( $nCat > 6 ) {
			qui replace `lnf' = ln(1 - exp(-exp(-`xb_c' - `xb_f1')) + ///
			                    ln(1 - exp(-exp(-`xb_c' - `xb_f2')) + ///
								ln(1 - exp(-exp(-`xb_c' - `xb_f3')) + ///
								ln(1 - exp(-exp(-`xb_c' - `xb_f4')) + ///
								ln(1 - exp(-exp(-`xb_c' - `xb_f5')) + ///
		                        ln(exp(-exp(-`xb_c' - `xb_f6')) if $ML_y == 6
		}
		
		if( $nCat > 7 ) {
			qui replace `lnf' = ln(1 - exp(-exp(-`xb_c' - `xb_f1')) + ///
			                    ln(1 - exp(-exp(-`xb_c' - `xb_f2')) + ///
								ln(1 - exp(-exp(-`xb_c' - `xb_f3')) + ///
								ln(1 - exp(-exp(-`xb_c' - `xb_f4')) + ///
								ln(1 - exp(-exp(-`xb_c' - `xb_f5')) + ///
								ln(1 - exp(-exp(-`xb_c' - `xb_f6')) + ///
		                        ln(exp(-exp(-`xb_c' - `xb_f7')) if $ML_y == 7
		}
		
		if( $nCat > 8 ) {
			qui replace `lnf' = ln(1 - exp(-exp(-`xb_c' - `xb_f1')) + ///
			                    ln(1 - exp(-exp(-`xb_c' - `xb_f2')) + ///
								ln(1 - exp(-exp(-`xb_c' - `xb_f3')) + ///
								ln(1 - exp(-exp(-`xb_c' - `xb_f4')) + ///
								ln(1 - exp(-exp(-`xb_c' - `xb_f5')) + ///
								ln(1 - exp(-exp(-`xb_c' - `xb_f6')) + ///
								ln(1 - exp(-exp(-`xb_c' - `xb_f7')) + ///
		                        ln(exp(-exp(-`xb_c' - `xb_f8')) if $ML_y == 8
		}
		
		if( $nCat > 9 ) {
			qui replace `lnf' = ln(1 - exp(-exp(-`xb_c' - `xb_f1')) + ///
			                    ln(1 - exp(-exp(-`xb_c' - `xb_f2')) + ///
								ln(1 - exp(-exp(-`xb_c' - `xb_f3')) + ///
								ln(1 - exp(-exp(-`xb_c' - `xb_f4')) + ///
								ln(1 - exp(-exp(-`xb_c' - `xb_f5')) + ///
								ln(1 - exp(-exp(-`xb_c' - `xb_f6')) + ///
								ln(1 - exp(-exp(-`xb_c' - `xb_f7')) + ///
								ln(1 - exp(-exp(-`xb_c' - `xb_f8')) + ///
		                        ln(exp(-exp(-`xb_c' - `xb_f9')) if $ML_y == 9
		}
		
		* build equation for last value of Y
		local eqn `" ln(1 - exp(-exp(-`xb_c' - `xb_f1')) "'
		forval i = 2/$nCatm1 {
			local eqn `" `eqn' + ln(1 - exp(-exp(-`xb_c' - `xb_f`i'')) "'
		}
		qui replace `lnf' = `eqn' if $ML_y == $nCat
	}
	
end


/* History
1.0.0  11.17.16  initial likelihood program for arbitrary number of categories

