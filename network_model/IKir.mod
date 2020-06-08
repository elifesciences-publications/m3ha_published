TITLE Potassium strong inward rectifier current IKir

COMMENT

    History: Modified from IA.mod of Amarillo et al., J Neurophysiol, 2014

    Current-voltage relationship: Described by Ohm’s Law.

    Gating: Uses an instantaneous activation gate. 
            Voltage dependence from Amarillo et al., J Neurophysiol, 2014
    
    2017-07-25 AL - Added units and comments
    2017-07-25 AL - Removed ik & ek from range variables
    2017-07-27 AL - Renamed suffix Ikir -> IKir
    2017-07-27 AL - Renamed gkmax -> gkbar
    2019-10-13 AL - Added minf in the NEURON block

ENDCOMMENT

NEURON {
    SUFFIX IKir
    USEION k READ ek WRITE ik
    RANGE gkbar
    RANGE minf, ik
}

UNITS {
    (mV)    = (millivolt)
    (mA)    = (milliamp)
    (S)     = (siemens)
}

PARAMETER {
	: RANGE variables whose values are specified in hoc
    gkbar   = 2.0e-5(S/cm2) : default maximum conductance of Ikir
}

ASSIGNED { 
    : Variables that are assigned outside the mod file
    v               (mV)
    ek              (mV)    : K+ reversal potential

	: RANGE variables that are assigned in the BREAKPOINT block
    minf            (1)     : steady state value of activation gating variable
    ik              (mA/cm2)
}

BREAKPOINT {    LOCAL bir
	: Calculate (instantaneous) activation gate
    minf = 1 / (1 + exp( (v + 97.9 (mV)) / 9.7 (mV) ) )

	: Calculate current ik from v, ek
    ik  = gkbar * minf * (v - ek)
}
