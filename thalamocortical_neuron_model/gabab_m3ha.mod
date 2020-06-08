TITLE Simple GABA-B receptor

COMMENT

	The kinetics reflect the measurements of Mark Beenhakker.

	Used by:	
		cd/singleneuron4compgabab.hoc

	2009-07-22 CKL	Edited to have two decay components (fast and slow) - 
	????-??-?? CKL	Edited to make conductance (g) a public variable
	2016-10-07 AL	Changed Erev to -115 mV because that's what Christine
			        used in the dclamp experiments, effectively
	2016-10-21 AL	Changed name of file to gabab_m3ha.mod
    2017-07-25 AL - Reorganized code
    2017-07-25 AL   Removed Rdummy
    2017-07-25 AL   Changed Ninputs to 1
    2017-07-25 AL   Made Erev GLOBAL
    2017-07-25 AL   Changed units of amp to umho and made Rs dimensionless
    2017-08-02 AL   Added q10 and changed value 2.2 -> 2.1 (Otis et al, 1993)
    2017-08-02 AL   Added phi in the NET_RECEIVE block as well
    2017-08-03 AL   Removed Rlast

ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS gabab
	NONSPECIFIC_CURRENT i
	GLOBAL p
	GLOBAL phi
	RANGE Erev, amp, Trise, TfallFast, TfallSlow, w, Ninputs
	RANGE Ron, RoffSlow, RoffFast, g
}

UNITS {
	(nA)    = (nanoamp)
	(mV)    = (millivolt)
	(umho)  = (micromho)
	(mM)    = (milli/liter)
}

PARAMETER {
    : GLOBAL variables whose values are fixed
	p           = 8		(1)		: power of rising phase
    q10         = 2.1   (1)     : Q10 for all phases, from Otis et al, 1993

	: RANGE variables whose values are specified in hoc
	Erev        = -115  (mV)	: reversal potential; this is probably what
					            : 	was used in the dynamic clamp experiments
	amp         = 15.92	(umho)  : maximum amplitude of gabab conductance
	Trise       = 52    (ms)	: rise time constant
	TfallFast   = 140.02(ms)	: fast decay time constant
	TfallSlow   = 1073  (ms)	: slow decay time constant
	w           = 0.952	(1)		: weight of fast decay (out of 1; 
					            :	weight of slow decay is 1 - w)
	Ninputs	    = 1		(1)	    : number of input streams
}

ASSIGNED {
    : Variables that are assigned outside the mod file
    v                   (mV)    : postsynaptic membrane potential
    celsius             (degC)  : temperature

	: GLOBAL variables that are assigned in the INITIAL block
	phi                 (1)     : temperature adjustment for rates

	: RANGE variables that are assigned in the NET_RECEIVE block
    act                 (1)     : amount of new activation

	: RANGE variables that are assigned in the BREAKPOINT block
	g                   (umho)  : conductance generated
	i 		            (nA)    : current generated = g * (v - Erev)
}

STATE {
	RoffSlow            (1)	    : slow decay variable 
	                            :   with maximum value 1 for an isolated IPSC
	RoffFast            (1)	    : fast decay variable
	                            :   with maximum value 1 for an isolated IPSC
	Ron                 (1)     : rise variable
	                            :   with maximum value 1 for an isolated IPSC
}

INITIAL { LOCAL tpeak
	: Calculate constants from temperature
	phi = q10 ^ ( ( celsius - 33 (degC) ) / 10 (degC))
	         		            : Christine did the experiments at 33 degC

	: Initialize state variables
	Ron = 0
	RoffSlow = 0
	RoffFast = 0
}

BREAKPOINT {
	: Update state variables from TfallSlow, TfallFast, Trise
	SOLVE release METHOD cnexp
	
	: Compute g from state variables Ron, RoffFast, RoffSlow
	g = amp * (1 - Ron)^p * ( w * RoffFast + (1-w) * RoffSlow )

	: Compute i from g, v
	i = g * (v - Erev)
}

DERIVATIVE release {
	: Update Ron, RoffFast, RoffSlow
	Ron' = -Ron / (Trise / phi)
	RoffFast' = -RoffFast / (TfallFast / phi)
	RoffSlow' = -RoffSlow / (TfallSlow / phi)
}

NET_RECEIVE(weight (1)) { 
    : weight    - amount of activation added by this synaptic event,
    :               relative to the total activation by the simultaneous
    :               activation of all GABAB synapses at this site

    : Compute new amount of activation (presynaptic potentiation?)
	act = weight / Ninputs

    : Add activation to current state variables
	RoffFast = RoffFast + act
	RoffSlow = RoffSlow + act
	Ron = Ron + act
}

