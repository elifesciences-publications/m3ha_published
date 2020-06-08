TITLE Simple GABA-B receptor

COMMENT

    The kinetics reflect the measurements of Mark Beenhakker.

    Used by:    
        cd/singleneuron4compgabab.hoc

    2009-07-22 CKL  Edited to have two decay components (fast and slow) - 
    ????-??-?? CKL  Edited to make conductance (g) a public variable
    2016-10-07 AL   Changed Erev to -115 mV because that's what Christine
                    used in the dclamp experiments, effectively
    2016-10-21 AL   Changed name of file to gabab_m3ha.mod
    2017-07-25 AL - Reorganized code
    2017-07-25 AL   Removed Rdummy
    2017-07-25 AL   Changed Ninputs to 1
    2017-07-25 AL   Made Erev GLOBAL
    2017-07-25 AL   Changed units of amp to umho and made Rs dimensionless
    2017-08-02 AL   Added q10 and changed value 2.2 -> 2.1 (Otis et al, 1993)
    2017-08-02 AL   Added phi in the NET_RECEIVE block as well
    2017-08-03 AL   Removed Rlast
    2017-11-01 AL   Made Erev, Trise, TfallFast, TfallSlow, w GLOBAL
    2017-11-05 AL   Replaced Rlast -> no change
    2017-11-06 AL   Fixed nonlinearity: 
                    Changed state variables to RFast0, RFast1, etc.

ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
    POINT_PROCESS gabab
    NONSPECIFIC_CURRENT i
    GLOBAL q10, Erev, Trise, TfallFast, TfallSlow, w
    RANGE amp, Ninputs
    GLOBAL phi
    RANGE act, total, g
    RANGE RFast0, RFast1, RFast2, RFast3, RFast4
    RANGE RFast5, RFast6, RFast7, RFast8
    RANGE RSlow0, RSlow1, RSlow2, RSlow3, RSlow4
    RANGE RSlow5, RSlow6, RSlow7, RSlow8
}

UNITS {
    (nA)    = (nanoamp)
    (mV)    = (millivolt)
    (uS)    = (micromho)
    (mM)    = (milli/liter)
}

PARAMETER {
    : GLOBAL variables whose values are fixed
    q10         = 2.1   (1)     : Q10 for all phases, from Otis et al, 1993
    Erev        = -115  (mV)    : reversal potential; this is probably what
                                :     was used in the dynamic clamp experiments
    Trise       = 52    (ms)    : rise time constant
    TfallFast   = 90.1  (ms)    : fast decay time constant
    TfallSlow   = 1073.2(ms)    : slow decay time constant
    w           = 0.952 (1)     : weight of fast decay (out of 1; 
                                :    weight of slow decay is 1 - w)

    : RANGE variables whose values are specified in hoc
    amp         = 0.016 (uS)    : maximum amplitude of conductance
    Ninputs     = 1     (1)     : number of input streams
}

ASSIGNED {
    : Variables that are assigned outside the mod file
    v                   (mV)    : postsynaptic membrane potential
    celsius             (degC)  : temperature

    : GLOBAL variables that are assigned in the INITIAL block
    phi                 (1)     : temperature adjustment for rates

    : RANGE variables that are assigned in the NET_RECEIVE block
    act                 (1)     : amount of new activation
    total               (1)     : total activation received

    : RANGE variables that are assigned in the BREAKPOINT block
    g                   (uS)    : total conductance (uS) generated
    i                   (nA)    : total current (nA) generated
}

STATE {
    : These variables have a maximum value 1 for an isolated IPSC
    RFast0              (1)
    RFast1              (1)
    RFast2              (1)
    RFast3              (1)
    RFast4              (1)
    RFast5              (1)
    RFast6              (1)
    RFast7              (1)
    RFast8              (1)
    RSlow0              (1)
    RSlow1              (1)
    RSlow2              (1)
    RSlow3              (1)
    RSlow4              (1)
    RSlow5              (1)
    RSlow6              (1)
    RSlow7              (1)
    RSlow8              (1)
}

INITIAL { LOCAL tpeak
    : Calculate constants from temperature
    phi = q10 ^ ( ( celsius - 33 (degC) ) / 10 (degC))

    : Initialize state variables
    RFast0 = 0
    RFast1 = 0
    RFast2 = 0
    RFast3 = 0
    RFast4 = 0
    RFast5 = 0
    RFast6 = 0
    RFast7 = 0
    RFast8 = 0
    RSlow0 = 0
    RSlow1 = 0
    RSlow2 = 0
    RSlow3 = 0
    RSlow4 = 0
    RSlow5 = 0
    RSlow6 = 0
    RSlow7 = 0
    RSlow8 = 0
    total = 0
}

BREAKPOINT {
    : Update state variables from TfallSlow, TfallFast, Trise
    SOLVE release METHOD cnexp
    
    : Compute g from state variables Ron, RoffFast, RoffSlow
    g = amp * (w * ( RFast0 - 8 * RFast1 + 28 * RFast2
            - 56 * RFast3 + 70 * RFast4 - 56 * RFast5
            + 28 * RFast6 - 8 * RFast7 + RFast8 )
            + (1-w) * ( RSlow0 - 8 * RSlow1 + 28 * RSlow2
            - 56 * RSlow3 + 70 * RSlow4 - 56 * RSlow5
            + 28 * RSlow6 - 8 * RSlow7 + RSlow8 ) )

    : Compute i from g, v
    i = g * (v - Erev)
}

DERIVATIVE release {
    : Update state variables
    RFast0' = -RFast0 * (0/Trise + 1/TfallFast) * phi
    RFast1' = -RFast1 * (1/Trise + 1/TfallFast) * phi
    RFast2' = -RFast2 * (2/Trise + 1/TfallFast) * phi
    RFast3' = -RFast3 * (3/Trise + 1/TfallFast) * phi
    RFast4' = -RFast4 * (4/Trise + 1/TfallFast) * phi
    RFast5' = -RFast5 * (5/Trise + 1/TfallFast) * phi
    RFast6' = -RFast6 * (6/Trise + 1/TfallFast) * phi
    RFast7' = -RFast7 * (7/Trise + 1/TfallFast) * phi
    RFast8' = -RFast8 * (8/Trise + 1/TfallFast) * phi
    RSlow0' = -RSlow0 * (0/Trise + 1/TfallSlow) * phi
    RSlow1' = -RSlow1 * (1/Trise + 1/TfallSlow) * phi
    RSlow2' = -RSlow2 * (2/Trise + 1/TfallSlow) * phi
    RSlow3' = -RSlow3 * (3/Trise + 1/TfallSlow) * phi
    RSlow4' = -RSlow4 * (4/Trise + 1/TfallSlow) * phi
    RSlow5' = -RSlow5 * (5/Trise + 1/TfallSlow) * phi
    RSlow6' = -RSlow6 * (6/Trise + 1/TfallSlow) * phi
    RSlow7' = -RSlow7 * (7/Trise + 1/TfallSlow) * phi
    RSlow8' = -RSlow8 * (8/Trise + 1/TfallSlow) * phi
}

NET_RECEIVE(weight (1)) {
    : weight    - amount of activation added by this synaptic event,
    :               relative to the total activation by the simultaneous
    :               activation of all synapses at this site

    : Compute new amount of activation
    act = weight / Ninputs

    : Add activation to current state variables
    RFast0 = RFast0 + act
    RFast1 = RFast1 + act
    RFast2 = RFast2 + act
    RFast3 = RFast3 + act
    RFast4 = RFast4 + act
    RFast5 = RFast5 + act
    RFast6 = RFast6 + act
    RFast7 = RFast7 + act
    RFast8 = RFast8 + act
    RSlow0 = RSlow0 + act
    RSlow1 = RSlow1 + act
    RSlow2 = RSlow2 + act
    RSlow3 = RSlow3 + act
    RSlow4 = RSlow4 + act
    RSlow5 = RSlow5 + act
    RSlow6 = RSlow6 + act
    RSlow7 = RSlow7 + act
    RSlow8 = RSlow8 + act

    : Increment the total activation received by this synapse
    total = total + act
}

COMMENT

    : Add activation to current state variables
    RoffFast = RoffFast + act
    RoffSlow = RoffSlow + act
    Ron = Ron + act

    p           = 8     (1)     : power of rising phase
    g = amp * (1 - Ron)^8 * ( w * RoffFast + (1-w) * RoffSlow )

    RoffSlow            (1)     : slow decay variable 
                                :   with maximum value 1 for an isolated IPSC
    RoffFast            (1)     : fast decay variable
                                :   with maximum value 1 for an isolated IPSC
    Ron                 (1)     : rise variable
                                :   with maximum value 1 for an isolated IPSC   
    : Update Ron, RoffFast, RoffSlow
    Ron' = -Ron / (Trise / phi)
    RoffFast' = -RoffFast / (TfallFast / phi)
    RoffSlow' = -RoffSlow / (TfallSlow / phi)

    Rlast               (1)     : amount of activation remaining from before
    Tlast               (ms)    : last time an event was received

ENDCOMMENT
