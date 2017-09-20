{smcl}
{* Last revised 11/27/16}{…}
{title:Title}

{p2colset 8 18 24 2}{...}
{p2col:{cmd:cratreg}} Continuation Ratio Regression Models {p_end}
{p2colreset}{...}

{marker syn}
{title:Syntax}

{p 8 18 2}
{cmd:cratreg} {it:depvar} [{it:indepvars}] 
{ifin} {weight} 
[{cmd:,}
{cmdab:pr:op(}{it:varlist}{cmd:)}
{cmdab:fr:ee(}{it:varlist}{cmd:)}
{cmdab:link(}{it:string}{cmd:)}
{cmdab:vce(}{it:vcetype}{cmd:)}
{cmdab:or} {cmdab:ef:orm}
{it:display_options} {it:maximize_options}
{p_end}


{title:Description}

{p 4 4 2} {cmd:cratreg} is a user-written command that fits generalized continuation ratio models (Fullerton and Xu 2016). This class of models permits three types of covariates: (1) covariates whose coefficients are constrained to be equal across cutpoint equations, (2) covariates whose coefficients vary across cutpoint equations by a common factor (listed using the prop option), and (3) covariates whose coefficients freely vary across cutpoint equations (listed using the free option). 

{p 4 4 2} {cmd:cratreg} supports factor variables, the {cmdab:svy} prefix, and the {cmdab:mi estimate} prefix.


{title:Options}

{p 4 8 2} {cmdab:pr:op(}{it:varlist}{cmd:)} specifies which, if any, of the independent variables have coefficients that vary across cutpoint equations by a common factor.

{p 4 8 2} {cmdab:fr:ee(}{it:varlist}{cmd:)} specifies which, if any, of the independent variables have coefficients that freely vary across cutpoint equations.

{p 4 8 2} {cmdab:link(}{it:string}{cmd:)} specifies the link function to be used. The default link function is the logit link. Users may also specify a probit or a cloglog link.

{p 4 8 2} {cmdab:vce(}{it:vcetype}{cmd:)} specifies the type of standard error reported. Users may specify {cmdab:r:obust} or {cmdab:cl:uster} {it:clustvar} standard errors. 

{p 4 8 2} {cmdab:or} report odds ratios

{p 4 8 2} {cmdab:ef:orm} report exponentiated coefficients



{marker aut}
{title:Author}

{p 5 5}
Shawn Bauldry {break}
Purdue University {break}
Department of Sociology {break}
sbauldry@purdue.edu {break}
{browse "https://github.com/sbauldry/crreg"}



{title:References}

{p 4 8 2}Fullerton, A and Xu, J. 2016. {it:Ordered Regression Models: Parallel, Partial, and Non-Parallel Alternatives}. New York: CRC Press.
