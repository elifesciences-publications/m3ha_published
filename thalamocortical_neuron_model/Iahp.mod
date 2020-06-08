TITLE Potassium AHP type current Iahp

COMMENT

    From RD Traub, J Neurophysiol 89:909-921, 2003
    Implemented by Maciej Lazarewicz 2003 (mlazarew@seas.upenn.edu)
    Modified by Yimy Amarillo, 2014 (Amarillo et al., J Neurophysiol, 2014)

    2017-07-25 AL - Added units and comments
    2017-07-25 AL - Now initializes m to minf instead of 0
    2017-07-27 AL - Renamed suffix kahp -> Iahp
    2017-07-27 AL - Renamed gkmax -> gkbar

ENDCOMMENT

NEURON { 
	SUFFIX Iahp
	USEION k READ ek WRITE ik
	USEION ca READ cai
	RANGE gkbar
	RANGE alpha, beta, ik
}

UNITS { 
	(mV)    = (millivolt) 
	(mA)    = (milliamp)
	(mM)    = (milli/liter)
    (S)     = (siemens)
} 

PARAMETER { 
	: RANGE variables whose values are specified in hoc
	gkbar = 1.5e-5  (S/cm2)	    : default maximum conductance of kahp
}
 
ASSIGNED { 
    : Variables that are assigned outside the mod file
    v               (mV)
	ek 		        (mV)    : K+ reversal potential (mV)
	cai		        (mM)

	: RANGE variables that are assigned in the INITIAL & DERIVATIVE blocks
	alpha           (/ms)   : rate of activation
	beta	        (/ms)   : rate of inactivation

	: RANGE variables that are assigned in the BREAKPOINT block
	ik              (mA/cm2): potassium current (mA/cm2) generated
}
 
STATE {
	m               (1)     : steady state value of activation gating variable
}

INITIAL { 
	: Initialize state variables to steady state values
	rates(cai)
	m = alpha / ( alpha + beta )
}
 
BREAKPOINT {
	: Update gating variable m from voltage v
	SOLVE states METHOD cnexp
	
    : Calculate current ik from m, v
	ik = gkbar * m * ( v - ek )
}
 
DERIVATIVE states { 
    : Update alpha, beta from cai
	rates(cai)

    : Update m from alpha, beta
	m' = alpha * ( 1 - m ) - beta * m 
}

PROCEDURE rates(cai (mM)) { 
	if( cai < 0.0001 (mM) ) {
		alpha = cai / ( 10000 (mM ms) )
	}else{
		alpha = 0.01 (/ms)
	}
	beta = 0.01 (/ms)
}

