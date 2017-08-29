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
{it:display_options}
{p_end}


{title:Description}

{p 4 4 2} {cmd:cratreg} is a user-written command that fits continuation ratio models (Fullerton and Xu 2016). This class of models permits three types of covariates: (1) covariates whose effects are constrained to be equal across cutpoint equations, (2) covariates whose effects vary across cutpoint equations by a common factor (listed using the prop option), and (3) covariates whose effects freely vary across cutpoint equations (listed using the free option). 

{p 4 4 2} {cmd:cratreg} supports factor variables, the {cmdab:svy} prefix, and the {cmdab:mi estimate} prefix.


{title:Options}

{p 4 8 2} {cmdab:pr:op(}{it:varlist}{cmd:)} specifies which, if any, of the independent variables have effects that vary across cutpoint equations by a common factor.

{p 4 8 2} {cmdab:fr:ee(}{it:varlist}{cmd:)} specifies which, if any, of the independent variables have effects that freely vary across cutpoint equations.

{p 4 8 2} {cmdab:link(}{it:string}{cmd:)} specifies the link function to be used. The default link function is the logit link. Users may also specify a probit or a cloglog link.

{p 4 8 2} {cmdab:vce(}{it:vcetype}{cmd:)} specifies the type of standard error reported. Users may specify {cmdab:r:obust} or {cmdab:cl:uster} {it:clustvar} standard errors. 



{marker exa}
{dlgtab: Examples}



{marker aut}
{title:Author}

{p 5 5}
Shawn Bauldry {break}
Purdue University {break}
Department of Sociology {break}
sbauldry@purdue.edu {break}
{browse "https://github.com/sbauldry/cratreg"}


{title:Acknowledgements}



{title:References}
