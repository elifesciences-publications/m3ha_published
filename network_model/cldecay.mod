TITLE Fast mechanism for submembranal Cl- concentration (cli)
:
: Takes into account:
:
:	- increase of cli due to chloride currents
:	- extrusion of chloride with a simple first order equation
:
: Parameters:
:
:	- depth: depth of the shell just beneath the membrane (in um)
:	- clinf: equilibrium concentration of chloride (10 mM)
:	- tauKCC2: time constant of chloride extrusion (must be fast)
:
: 2017-02-07 Adapted from cadecay.mod by Alain Destexhe, Salk Institute, 1995
: 2017-02-22 This is not used; cldecay2.mod is used instead
: 2017-03-07 Initialize cli from hoc
: 2017-03-12 This is now an option
: 2017-03-15 Made drive_channel & drive_extrusion RANGE variables
: 2017-03-15 Forced drive_extrusion to be negative instead
: 2017-03-31 Changed units of tauKCC2 to seconds
:

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX cld
	USEION cl READ icl, cli WRITE cli VALENCE -1
	RANGE depth, tauKCC2, clinf
	RANGE drive_channel, drive_extrusion
}

UNITS {
	(molar) = (1/liter)		: moles do not appear in units
	(mM)	= (millimolar)
	(um)	= (micron)
	(mA)	= (milliamp)
	(msM)	= (ms mM)
	FARADAY = (faraday) (coulombs)
}

PARAMETER {
	: RANGE variables whose values are specified in hoc
	depth	= .1	(um)		: depth of shell for Cl- (um)
	tauKCC2	= 30	(s)		: Cl- removal time constant (s), Peter's value (Jedlicka et al 2011 used 3 s)
	clinf	= 8	(mM)		: steady state intracellular [Cl-] (mM)
}

ASSIGNED {
	: Variables that are assigned outside the mod file
	icl		(mA/cm2)

	: RANGE variables that are assigned in the DERIVATIVE block
	drive_channel	(mM/ms)		: driving Cl- flux (mM/ms) due to channel opening
	drive_extrusion	(mM/ms)		: driving Cl- flux (mM/ms) due to leak/extrusion
}
	
STATE {
	cli		(mM)		: intracellular [Cl-] (mM)
}

BREAKPOINT {
	SOLVE state METHOD derivimplicit
}

DERIVATIVE state { 

	: Calculate driving flux due to channel opening
	drive_channel =  - (10000) * icl / ((-1) * FARADAY * depth)

	: Calculate driving flux due to leak/extrusion
	drive_extrusion = (clinf-cli)/(tauKCC2*(1000))

	: Calculate change in chloride concentration
	cli' = drive_channel + (clinf-cli)/(tauKCC2*(1000))
}

COMMENT
OLD CODE:
	FARADAY = 96489		(coul)		: moles do not appear in units
CONSTANT {
	FARADAY = 96485.3329	(coul)	: moles do not appear in units
}

INITIAL {
	cli = clinf
}

	: if (drive_channel <= 0.) { drive_channel = 0. }	: since the pump is outward, do not let cli decrease from icl
	if (drive_extrusion > 0.) { drive_extrusion = 0. }	: Cl- pump is outward
	tauKCC2	= 30000	(ms)		: Cl- removal time constant (ms), Peter's value (Jedlicka et al 2011 used 3000 ms)

ENDCOMMENT
