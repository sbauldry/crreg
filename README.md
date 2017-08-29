# crreg
This repository contains a set of commands for fitting continuation ratio models. The command allows for effects to be freely estimated across thresholds, to have a proportionality constraint across thresholds, or to be constrained to be equal across thresholds. Type "net install cratreg, from(https://github.com/sbauldry/crreg/raw/master) replace" without the quotes in Stata's command window to install the commands.

#Files
1. cratreg.ado:     primary command file
2. cr_lf.ado:     likelihood function for constrained model
3. cr_np_lf.ado:  likelihood function for free model
4. cr_pc_lf.ado:  likelihood function for proportional constraint model
5. cr_pnp_lf.ado: likelihood function for partial free model
6. cr_ppc_lf.ado: likelihood function for partial ppc model
7. cratreg.sthlp:   help file for command
8. cratreg.pkg:     Stata package file
9. stata.toc:     Stata toc file
