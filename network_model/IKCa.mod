TITLE Linear calcium-dependent potassium current
: This is a simple model in which the potassium conductance is just
: proportional to the intracellular calcium concentration
: gbar = 2e-3 reduces spikes/burst from 22 to 16 when P_Ts = .002
:
: 2017-02-07 Parameters modified to be consistent with Sohal & Huguenard 2003
: 2017-02-10 Removed celsius from the PARAMETER block
: 2017-03-04 Moved cai from PARAMETER to ASSIGNED
: 2017-03-05 Moved some variables from PARAMETER to ASSIGNED

NEURON {
	SUFFIX kca
	USEION ca READ cai
	USEION k READ ek WRITE ik
	RANGE gbar
	RANGE ik
}

UNITS {
	(mV)	=	(millivolt)
	(mA)	=	(milliamp)
	(molar)	=	(/liter)
	(mM)	=	(millimolar)
}

PARAMETER {
	: RANGE variables whose values are specified in hoc
	ek			(mV)
	gbar	= 1.5e-3	(mho/cm2 mM)	: [Ca++]i-dependent K+ conductance, Sohal & Huguenard 2003
}

ASSIGNED {
	: Variables that are assigned outside the mod file
	v			(mV)
	cai			(mM)		: intracellular [Ca++] (mM)

	: RANGE variables that are assigned in the BREAKPOINT block
	ik			(mA/cm2)
}

BREAKPOINT {
	ik = gbar * cai * (v - ek)
}

COMMENT
OLD CODE:

	gbar	= 2e-3	(mho/cm2 mM)	: cai-dependent conductance, Dextexhe et al 1996
	celsius = 37		(degC)	: Dextexhe et al 1996
	celsius = 34		(degC)	: Sohal & Huguenard 2003; however, this is overriden by NEURON

INITIAL {

}

ENDCOMMENT

