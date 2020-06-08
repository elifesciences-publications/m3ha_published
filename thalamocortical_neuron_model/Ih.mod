TITLE Hyperpolarization-activated nonspecific cationic current Ih

COMMENT

    History: Modified from Ih.mod of Amarillo et al., J Neurophysiol, 2014. 
                Based on the model of Huguenard & McCormick, 
                J Neurophysiol 68: 1373-1383, 1992, with updated
                kinetics from Santoro et al., 2000 & Amarillo et al., 2014.

    Current-voltage relationship: Described by Ohm’s Law.

    Gating: Uses 1 activation gate (m). 
            Voltage dependence and kinetics of activation at 34 degC
            from Amarillo et al., 2014.
            Huguenard & McCormick originally had V_1/2 = -75 mV and k = 5.5 mV
            Santoro et al. had V_1/2 = -82 mV

    Permeability ratio: K:Na is about 3:1~4:1, Santoro et al., 1999.

    Approximate eh:
        Based on [Na+]out = 127.25 mM, [Na+]in = 4.5 mM
        [K+]out = 2.5 mM, [K+]in = 113 mM & celsius = 33 degC
        the GHK voltage equation yields -24 ~ -32 mV
        Santoro et al., 1999 had -35 mV
        Amarillo et al., 2014 used -43 mV

    Identity: HCN channels (Hyperpolarization-activated cyclic-nucleotide 
                dependent cation-nonspecific channels).
                mHCN2 & mHCN4 found in thalamocortical relay neurons.
                See Santoro et al., 2000.

    2017-07-23 AL - Reorganized code
    2017-07-23 AL - Made shift a depolarizing shift
    2017-07-23 AL - Added qm
    2017-07-25 AL - Removed ih and added minf, taum from RANGE
    2017-07-27 AL - Renamed suffix iar -> Ih
    2017-07-28 AL - Renamed shift -> shiftm
    2017-08-01 AL - Renamed tadj -> phim
    2017-08-01 AL - Changed placement of minus signs
    2017-08-01 AL - Changed k of minf 5.49 -> 5.5

ENDCOMMENT

NEURON {
    SUFFIX Ih
    NONSPECIFIC_CURRENT ih
    GLOBAL qm
    GLOBAL phim
    RANGE ghbar, eh, shiftm
    RANGE minf, taum
}

UNITS {
    (mV)    = (millivolt)
    (mA)    = (milliamp)
    (S)     = (siemens)
}

PARAMETER {
    : GLOBAL variables whose values are fixed
    qm      = 4.0   (1)     : Q10 for activation, Santoro & Tibbs, 1999, 
                            :   based on values of 3.13 (sheep Purkinje fibers), 
                            :   4.5 (rat CA1 pyramidal neurons), 
                            :   5 (guinea pig CA1 pyramidal neurons)

    : RANGE variables whose values are specified in hoc
    ghbar   = 2.2e-5(S/cm2) : default maximum conductance of Ih
    eh      = -43   (mV)    : reversal potential of Ih
    shiftm  = 0     (mV)    : depolarizing shift of activation curve
}

ASSIGNED { 
    : Variables that are assigned outside the mod file
    v               (mV)
    celsius         (degC)

    : GLOBAL variables that are assigned in the INITIAL block
    phim            (1)     : temperature adjustment to taum

    : RANGE variables that are assigned in the INITIAL & DERIVATIVE blocks
    minf            (1)     : steady state value of activation gating variable
    taum            (ms)    : time constant for activation

    : RANGE variables that are assigned in the BREAKPOINT block
    ih              (mA/cm2): current (mA/cm2) generated
}

STATE {
    m               (1)     : activation gating variable
}

INITIAL {
    : Calculate constants from temperature
    phim = qm ^ ( ( celsius - 34 (degC) ) / 10 (degC) )

	: Initialize state variables
    settables(v)            : calculate steady state values of state variables
    m = minf
}

BREAKPOINT {
    : Update gating variable m from voltage v
    SOLVE states METHOD cnexp

    : Calculate current ih from m, eh, v
    ih = ghbar * m * (v - eh)
}

DERIVATIVE states {  
    : Update minf, taum from v
    settables(v)
    
    : Update m from minf, taum
    m' = (minf - m) / taum
}

PROCEDURE settables(v (mV)) { 
    : Update minf, taum from v
    minf  = 1 / (1 + exp( (v + 82 (mV) - shiftm ) / 5.5 (mV) ) )           
    taum  = ( 1.0 (ms) / 
                ( ( 8e-4 + 3.5e-6 * exp(-0.05787 (/mV) * (v - shiftm) ) ) 
                    + exp(-1.87 + 0.0701 (/mV) * (v - shiftm) ) ) ) / phim      : TODO: Where did these numbers come from?
}

