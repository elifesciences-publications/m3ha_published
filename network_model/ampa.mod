TITLE simple AMPA receptors

COMMENT

    Adapted from Sohal et al 2000 with
        updated kinetics reflecting the measurements of Deleuze & Huguenard 2016

    2017-11-01 Changed default Ninputs to 1

ENDCOMMENT

NEURON {
    POINT_PROCESS ampa
    NONSPECIFIC_CURRENT i
    GLOBAL q10, Trise, Tfall, Erev
    RANGE gmax, Ninputs
    GLOBAL phi, factor
    RANGE act, total, g, Ron, Roff
}

UNITS {
    (nA)    = (nanoamp)
    (mV)    = (millivolt)
    (uS)    = (micromho)
    (mM)    = (milli/liter)
}

PARAMETER {
    : GLOBAL variables whose values are specified in hoc
    q10     = 1     (1)         : Q10 for all phases
    Trise   = 1.1   (ms)        : rise time constant
    Tfall   = 6.0   (ms)        : fall time constant
    Erev    = 0     (mV)        : reversal potential

	: RANGE variables whose values are specified in hoc
    gmax    = 0     (uS)        : maximum conductance
    Ninputs = 1     (1)         : number of input streams
}

ASSIGNED {
    : Variables that are assigned outside the mod file
    v               (mV)        : postsynaptic membrane potential
    celsius         (degC)      : temperature

	: GLOBAL variables that are assigned in the INITIAL block
    phi             (1)         : temperature correction for rates
    factor          (1)         : normalization factor for the conductance

    : RANGE variables that are assigned in the NET_RECEIVE block
    act             (1)         : amount of new activation
    total           (1)         : total activation received

    : RANGE variables that are assigned in the BREAKPOINT block
    g               (uS)        : total conductance (uS) generated
    i               (nA)        : total current (nA) generated
}

STATE {
	Ron             (1)         : rise variable
	                            :   with maximum value 1 for an isolated IPSC
	Roff            (1)	        : decay variable 
	                            :   with maximum value 1 for an isolated IPSC
                                : Ron - Roff yields dual exp timecourse
}

INITIAL { LOCAL tp, mv
	: Calculate constants from temperature
	phi = q10 ^ ( ( celsius - 24 (degC) ) / 10 (degC))

    : Calculate normalization factor from phi, Trise & Tfall
    tp = (Trise * Tfall) / ((Tfall - Trise) * phi) * log(Tfall/Trise)    
                                : t at which exp(-t/(Tfall/phi))-exp(-t/(Trise/phi)) reaches maximum
    mv = exp(-phi * tp / Tfall) - exp(-phi * tp / Trise)    
                                : maximum value of exp(-t/(Tfall/phi))-exp(-t/(Trise/phi))
    factor = 1/mv               : normalization factor for exp(-t/(Tfall/phi))-exp(-t/(Trise/phi))

    : Initialize state variables
    Ron = 0
    Roff = 0
}

BREAKPOINT {
    : Calculate rising and falling phases of the conductance
    SOLVE release METHOD cnexp

    : Calculate the conductance
    g = (Roff - Ron) * factor * gmax

    : Calculate the total current through the receptor
    i = g * (v - Erev)
}

DERIVATIVE release {
    : Update Ron, Roff
    Ron' = -Ron / (Trise / phi)
    Roff' = -Roff / (Tfall / phi)
}

NET_RECEIVE(weight (1)) {
    : weight    - amount of activation added by this synaptic event,
    :               relative to the total activation by the simultaneous
    :               activation of all synapses at this site

    : Compute new amount of activation
    act = weight / Ninputs

    : Add activation to current state variables
    Ron = Ron + act             : add to the rise variable
    Roff = Roff + act           : add to the decay variable

    : Increment the total activation received by this synapse
    total = total + act
}

COMMENT
OLD CODE:

    i = gmax*(Roff + Ron)*(v-Erev)

    NET_RECEIVE(weight, Rlast, Tlast (ms)) { LOCAL act
    :weight represents the fraction of available channels activated by a spike
        Rlast = Rlast * norm * (Exp((Tlast-t)/Tfall) - Exp((Tlast-t)/Trise))
        Tlast = t
        act = weight * (1 - Rlast) / Ninputs
        state_discontinuity(Roff, Roff + act)
        state_discontinuity(Ron, Ron - act)
    }

    PROCEDURE Exp(x) { LOCAL y
        if (x > -100 && x < 100) {
            y = exp(x)
            VERBATIM
            return _ly;
            ENDVERBATIM
        } else {
            VERBATIM
            return 0;
            ENDVERBATIM
        }
    }

ENDCOMMENT
