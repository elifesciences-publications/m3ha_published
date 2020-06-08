TITLE Fast transient potassium current IA

COMMENT

    History: Modified from IA.mod of Amarillo et al., J Neurophysiol, 2014
                based on the model of Huguenard & McCormick, 
                J Neurophysiol 68: 1373-1383, 1992.

    Current-voltage relationship: Described by Ohm’s Law.

    Gating: Uses 4 activation gates and 1 inactivation gate (m4h). 
            There are two types of activation gates, each paired with
            a type of inactivation gates. The ratio of contribution is 3:2.
            Voltage dependence and kinetics of activation/inactivation 
            at 23 degC from voltage-clamp data (whole cell patch clamp) of 
            Huguenard & Prince, J. Neurosci. 12: 3804-3817, 1992.

    2017-07-25 AL - Added units and comments
    2017-07-27 AL - Renamed gkmax -> gkbar
    2017-08-01 AL - Renamed tadj -> phi
    2017-08-01 AL - Changed placement of minus signs

ENDCOMMENT

NEURON {
    SUFFIX IA
    USEION k READ ek WRITE ik
    GLOBAL q10
    GLOBAL phi
    RANGE gkbar
    RANGE m1inf, m2inf, hinf, taum, tauh1, tauh2, ik
}

UNITS {
    (mV)    = (millivolt)
    (mA)    = (milliamp)
    (S)     = (siemens)
}

PARAMETER {
    : GLOBAL variables whose values are fixed
    q10     = 2.8   (1)     : Q10 for both activation and inactivation,
                            :   from Huguenard et al., 1991
    
	: RANGE variables whose values are specified in hoc
    gkbar   = 5.5e-3(S/cm2) : default maximum conductance of IA
}

ASSIGNED { 
    : Variables that are assigned outside the mod file
    v               (mV)
    celsius         (degC)
    ek              (mV)    : K+ reversal potential

	: GLOBAL variables that are assigned in the INITIAL block
    phi             (1)     : temperature adjustment to taum & tauh

	: RANGE variables that are assigned in the INITIAL & DERIVATIVE blocks
    m1inf   (1)     : steady state value of activation gating variable #1
    m2inf   (1)     : steady state value of activation gating variable #2
    hinf    (1)     : steady state value of both inactivation gating variables
    taum    (ms)    : time constant for both activation gating variables
    tauh1   (ms)    : time constant for inactivation gating variable #1
    tauh2   (ms)    : time constant for inactivation gating variable #1

	: RANGE variables that are assigned in the BREAKPOINT block
	ik              (mA/cm2): potassium current (mA/cm2) generated
}

STATE {
	m1		        (1)		: activation gating variable #1
	m2		        (1)		: activation gating variable #2
	h1		        (1)		: inactivation gating variable #1
	h2		        (1)		: inactivation gating variable #2
}

INITIAL {
	: Calculate constants from temperature
    phi = q10 ^ ( ( celsius - 23 (degC) ) / 10 (degC) )

	: Initialize state variables to steady state values
    settables(v)			: calculate steady state values of state variables
    m1 = m1inf
    m2 = m2inf
    h1 = hinf
    h2 = hinf
}

BREAKPOINT {
	: Update gating variables m1, m2, h1, h2 from voltage v
    SOLVE states METHOD cnexp

    : Calculate current ik from  m1, m2, h1, h2, ek, v
    ik = gkbar * ( 0.6 * m1^4 * h1 + 0.4 * m2^4 * h2 ) * ( v - ek )
}

DERIVATIVE states {  
    : Update m1inf, m2inf, hinf, taum, tauh1, tauh2 from v
    settables(v)

    : Update m1, m2, h1, h2 from m1inf, m2inf, hinf, taum, tauh1, tauh2
    m1' = (m1inf - m1) / taum
    m2' = (m2inf - m2) / taum
    h1' = (hinf - h1) / tauh1
    h2' = (hinf - h2) / tauh2
}

PROCEDURE settables(v (mV)) { LOCAL tauh
    : Update m1inf, m2inf, hinf, taum, tauh1, tauh2 from v
    m1inf = 1 / ( 1 + exp( ( v + 60 (mV) ) / -8.5 (mV) ) )
    m2inf = 1 / ( 1 + exp( ( v + 36 (mV) ) / -20 (mV) ) )
    hinf  = 1 / ( 1 + exp( ( v + 78 (mV) ) / 6.0 (mV) ) )

    : Update taum from v
    taum  = (1 / phi) * 
            ( 0.37 (ms) + 
                1.0 (ms) / (exp( ( v + 35.8 (mV) ) / 19.7 (mV) ) + 
                            exp( ( v + 79.7 (mV) ) / -12.7 (mV) ) ) ) 

    : Update tauh1, tauh2 from v
    tauh  = (1 / phi) *
            ( 1.0 (ms) / ( exp( ( v + 46 (mV) ) / 5.0 (mV) ) + 
                            exp( ( v + 238 (mV) ) / -37.5 (mV) ) ) )
    if ( v < -63 (mV) ) {
        tauh1 = tauh
    } else {
        tauh1 = (19 (ms) / phi)
    }
    if ( v < -73 (mV) ) {
        tauh2 = tauh
    } else {
        tauh2 = (60 (ms) / phi)
    }
}

