TITLE Potassium C type current IKCa

COMMENT

    from RD Traub, J Neurophysiol 89:909-921, 2003
    implemented by Maciej Lazarewicz 2003 (mlazarew@seas.upenn.edu)
    Modified by Yimy Amarillo, 2014 (Amarillo et al., J Neurophysiol, 2014)
 
    2017-07-25 AL - Added units and comments
    2017-07-25 AL - Now initializes m to minf instead of 0
    2017-07-27 AL - Renamed gkmax -> gkbar

ENDCOMMENT

INDEPENDENT { t FROM 0 TO 1 WITH 1 (ms) }

NEURON { 
	SUFFIX IKCa
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
	gkbar = 1.0e-4  (S/cm2) : default maximum conductance of IKCa
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
	settables(v) 
	m = alpha / ( alpha + beta )
}
 
BREAKPOINT { 
	: Update gating variable m from voltage v
	SOLVE states METHOD cnexp
	
    : Calculate current ik from m, cai, v
	if ( 0.004 (/mM) * cai < 0.000001 ) {
		ik = gkbar * m * ( 0.004 (/mM) ) * cai * ( v - ek ) 
	} else {
		ik = gkbar * m * ( v - ek ) 
	}
}
 
DERIVATIVE states { 
    : Update alpha, beta from v
	settables(v) 

    : Update m from alpha, beta
	m' = alpha * ( 1 - m ) - beta * m 
}

PROCEDURE settables(v (mV)) { 
	TABLE alpha, beta FROM -120 TO 40 WITH 641

	if ( v < -10.0 (mV) ) {
		alpha = 2 (/ms) / 37.95 * ( exp( ( v + 50 (mV) ) / 11 (mV) 
		                            - ( v + 53.5 (mV) ) / 27 (mV) ) )
		
		: Note that there is a typo in the paper:
		:   missing minus sign in the front of 'v'
		beta  = 2 (/ms) * exp( ( v + 53.5 (mV) ) / -27 (mV) ) - alpha
	} else {
		alpha = 2 (/ms) * exp( ( v + 53.5 (mV) ) / -27 (mV) )
		beta  = 0 (/ms)
	}
}

