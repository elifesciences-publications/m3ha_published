TITLE High threshold calcium current IL

COMMENT

    Based on the model by Huguenard and & McCormick, J Neurophysiol, 1992
    including errata corrections
    Written by Yimy Amarillo, 2014 (Amarillo et al., J Neurophysiol, 2014)

    2017-07-25 AL - Added units and comments
    2017-07-25 AL - Now initializes mL to mLinf instead of 0
    2017-07-27 AL - Renamed pLmax -> pcabar
    2017-08-01 AL - Renamed tadj -> phimL

ENDCOMMENT

NEURON {
    SUFFIX IL
    USEION ca READ cai, cao WRITE ica
    GLOBAL q10
    GLOBAL phimL
    RANGE pcabar
    RANGE mLinf, taumL, ica
}

UNITS {
    (molar) = (1/liter)
    (mV)    = (millivolt)
    (mA)    = (milliamp)
    (mM)    = (millimolar)

    FARADAY = (faraday) (coulombs)
    R       = (k-mole) (joule/degC)
}

PARAMETER {
    : GLOBAL variables whose values are fixed
    q10     = 3     (1)     : Q10 for activation

	: RANGE variables whose values are specified in hoc
    pcabar   = 1.0e-4(cm/s)	: default maximum permeability of IL
}

ASSIGNED { 
    : Variables that are assigned outside the mod file
    v               (mV)
    celsius         (degC)
	cai		        (mM)    : intracellular [Ca++] (mM)
	cao		        (mM)    : extracellular [Ca++] (mM)

	: GLOBAL variables that are assigned in the INITIAL block
    phimL            (1)     : temperature adjustment to taumL

	: RANGE variables that are assigned in the INITIAL & DERIVATIVE blocks
    mLinf           (1)     : steady state value of activation gating variable
    taumL           (ms)    : time constant for activation

	: RANGE variables that are assigned in the BREAKPOINT block
    ica             (mA/cm2): calcium current (mA/cm2) generated
}

STATE {
    mL		        (1)		: activation gating variable
}

INITIAL {
	: Calculate constants from temperature
    phimL = q10 ^ ( ( celsius - 23.5 (degC) ) / 10 (degC) )

	: Initialize state variables to steady state values
    settables(v)			: calculate steady state values of state variables
    mL = mLinf
}

BREAKPOINT {
	: Update gating variable mL from voltage v
    SOLVE states METHOD cnexp
    
    : Calculate current ica from mL, v
    ica = pcabar * mL^2 * ghk(v, cai, cao)
}

DERIVATIVE states {  
    : Update mLinf, taumL from v
    settables(v)
    
    : Update mL from mLinf, taumL
    mL' = (mLinf - mL) / taumL
}

PROCEDURE settables(v (mV)) { LOCAL mLalpha, mLbeta
    : Update mLinf from v
    mLinf = 1 / ( 1 + exp( ( v + 10 (mV) ) / -10 (mV) ) )

    : Update taumL from v
    mLalpha = 1.6 (/ms) / ( 1 + exp( -0.072 (/mV) * ( v - 5 (mV) ) ) )
    mLbeta = 1.0 (/ms) * 0.02 (/mV) * ( v - 1.31 (mV) ) 
                / ( exp ( ( v - 1.31 (mV) ) / 5.36 (mV) ) - 1 )
    taumL = ( 1 / (mLalpha + mLbeta) ) / phimL
}

FUNCTION ghk(v (mV), ci (mM), co (mM)) (.001 coul/cm3) {
    LOCAL z, eci, eco
    
    z = (1e-3) * 2 * FARADAY * v / (R * (celsius + 273.15 (degC) ) )
                            : this is ZFV/RT, which is dimensionless 
                            : after applying conversion factor 1e-3 V/mV
    eco = co * efun(z)      : this has units of [mM] = [umol/cm3]
    eci = ci * efun(-z)     : this has units of [mM] = [umol/cm3]
    ghk = (1e-3) * 2 * FARADAY * (eci - eco)
                            : this has units of [mC/cm3]
                            : after applying conversion factor 1e-3 mmol/umol
}

FUNCTION efun(z) {
	if (fabs(z) < 1e-4) {
		efun = 1 - z/2      : 1st order Taylor approximation of z/(exp(z) - 1)
	} else {
		efun = z / (exp(z) - 1)
	}
}

FUNCTION nongat(v (mV), cai (mM), cao (mM)) (mA/cm2) {    : non gating variabled current
	nongat = pcabar * ghk(v, cai, cao)
}

