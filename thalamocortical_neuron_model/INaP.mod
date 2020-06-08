TITLE Persistent sodium current INaP

COMMENT

    History: Modified from INaP.mod of Amarillo et al., J Neurophysiol, 2014.
                Based on the model by Wu et al, 2005 on
                mesencephalic trigeminal sensory neurons

    Current-voltage relationship: Described by Ohm’s Law.

    Gating: Uses 1 activation gate and 1 inactivation gate (mh).
            The activation is instantaneous whereas 
            the inactivation is slow and time-dependent.
            Voltage dependence and kinetics of activation/inactivation 
            at 23 degC from voltage-clamp data (whole cell patch clamp) of 
            Wu et al., 2005.

    2017-07-25 AL - Added units and comments
    2017-07-25 AL - Removed ina & ena from range variables
    2017-07-25 AL - Made qh a global variable
    2017-07-27 AL - Renamed gnabar -> gnabar
    2017-08-01 AL - Renamed tadj -> phih

ENDCOMMENT

NEURON {
	SUFFIX INaP
	USEION na READ ena WRITE ina
    GLOBAL qh
    GLOBAL phih
	RANGE gnabar
	RANGE minf, hinf, tauh, ina
}

UNITS {
	(mV)    = (millivolt)
	(mA)    = (milliamp)
	(S)     = (siemens)
}

PARAMETER {
    : GLOBAL variables whose values are fixed
    qh      = 3     (1)     : Q10 for inactivation, assumed by Amarillo et al

	: RANGE variables whose values are specified in hoc
	gnabar  = 5.5e-6(S/cm2)	: default maximum conductance of INaP
}

ASSIGNED { 
    : Variables that are assigned outside the mod file
    v               (mV)
    celsius         (degC)
	ena             (mV)    : Na+ reversal potential

	: GLOBAL variables that are assigned in the INITIAL block
    phih            (1)     : temperature adjustment to taum

	: RANGE variables that are assigned in the INITIAL & DERIVATIVE blocks
    minf            (1)     : steady state value of activation gating variable
    hinf            (1)     : steady state value of inactivation gating variable
    tauh            (ms)    : time constant for inactivation

	: RANGE variables that are assigned in the BREAKPOINT block
	ina	            (mA/cm2): sodium current (mA/cm2) generated
}
	
STATE {
	h		        (1)		: inactivation gating variable
}

INITIAL {
	: Calculate constants from temperature
    phih = qh ^ ( (celsius - 23 (degC) ) / 10 (degC) )

	: Initialize state variables
    setvalues(v)			: calculate steady state values of state variables
	h = hinf
}

BREAKPOINT {
	: Update gating variable h from voltage v
	SOLVE states METHOD cnexp

    : Calculate current ih from minf, h, ena, v
	ina = gnabar * minf * h * (v - ena)
}

DERIVATIVE states {
    : Update minf, hinf, tauh from v
    setvalues(v)

    : Update h from hinf, tauh
	h' = (hinf - h) / tauh
}

PROCEDURE setvalues(v (mV)) { 
	: Update minf, hinf, taum from v
	minf = 1 / ( 1 + exp( ( v + 57.9 (mV) ) / -6.4 (mV) ) )
	hinf = 1 / ( 1 + exp( ( v + 58.7 (mV) ) / 14.2 (mV) ) )
	tauh = (1 / phih) * (1000 (ms) + 
	            (10000 (ms) / (1 + exp( ( v + 60 (mV) ) / 10 (mV) ) ) ) )
}

