TITLE T-type calcium current responsible for low-threshold spikes (LTS) IT

COMMENT

    History: Modified from ITGHK.mod of the Destexhe et al 1998a model, 
                based on the model of Huguenard & McCormick, 
                J Neurophysiol 68: 1373-1383, 1992.

    Current-voltage relationship: Described by Goldman-Hodgkin-Katz equations.

    Gating: Uses 2 activation gates and 1 inactivation gate (m²h). 
            Voltage dependence and kinetics of activation/inactivation 
            at 23 degC from voltage-clamp data (whole cell patch clamp) of 
            Huguenard & Prince, J. Neurosci. 12: 3804-3817, 1992. 
            Updated to reflect values in Destexhe et al, 1998.

            The activation and inactivation functions can be 
            empirically corrected to account for the contamination 
            of inactivation, to compensate for screening charge, etc. 
            The correction terms are denoted shiftm and shifth and 
            cause depolarizing (rightward) shifts.

            The steepness of the activation and inactivation functions 
            can be varied with the parameters slopem and slopeh, respectively.

    This model is described in detail in:
    Destexhe A, Neubig M, Ulrich D and Huguenard JR.  
    Dendritic low-threshold calcium currents in thalamic relay cells.  
    Journal of Neuroscience 18: 3574-3588, 1998.
    (a postscript version of this paper, including figures, is available on
    the Internet at http://cns.fmed.ulaval.ca)

    - shift parameter for screening charge
    - empirical correction for contamination by inactivation (Huguenard)
    - GHK equations

    Written by Alain Destexhe, Laval University, 1995

    Modifications by Christine:
    1) Decoupled shift terms for ease of parameter fitting
       i.e. shift (screening charge parameter) and actshift are now changed to
            shifth and shiftm.  
            shiftm = shift+actshift, shifth = shift
       both shiftm > 0 and shifth > 0 cause depolarizing (rightward) shifts
            suffix is still itGHK
            [CKL 2011-03-19]

    2) Made the activation and inactivation curve "steepness" terms into 
       range variables (slopeh, slopem) where steepness is a scaling factor for 
       the original values (m: -6.2, h: 4.0)

    2016-10-21 AL - Removed UNITSOFF and UNITSON
    2017-07-19 AL - Removed cai & cao from parameters and make them assigned
    2017-07-22 AL - Changed numerical method from euler to cnexp
    2017-07-22 AL - Reorganized code
    2017-07-27 AL - Renamed suffix itGHK -> IT
    2017-07-31 AL - Fixed default Dextexhe shiftm from 2 -> -2
    2017-08-01 AL - Changed placement of minus signs
    2017-08-01 AL - Changed q10 base point from 24 -> 23 degC
    2017-08-06 AL - Changed default values for shiftm 2 -> 1
    2017-08-06 AL - Changed default values for shiftm 0 -> 1
    2020-01-02 AL - Added tauhmode

ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
    SUFFIX IT
    USEION ca READ cai, cao WRITE ica
    GLOBAL qm, qh
    GLOBAL phim, phih
    RANGE pcabar, shiftm, shifth, slopem, slopeh, tauhmode
    RANGE minf, hinf, taum, tauh, ica
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
    qm      = 3.6   (1)     : Q10 for activation, based on Coulter et al 1989
    qh      = 2.8   (1)     : Q10 for inactivation, based on Coulter et al 1989

	: RANGE variables whose values are specified in hoc
    pcabar  = .2e-3 (cm/s)  : default maximum Ca++ permeability
    shiftm  = 1     (mV)    : depolarizing shift of activation curve
    shifth  = 1     (mV)    : depolarizing shift of inactivation curve
    slopem  = 1     (1)     : scaling factor for slope of activation curve
    slopeh  = 1     (1)     : scaling factor for slope of inactivation curve
    tauhmode = 0    (1)     : tauh mode:
                            :   0 - original curve
                            :   1 - the same as taum
                            :   2 - 10 times smaller amplitude
                            :   3 - 10 times larger amplitude
                            :   4 - 5 times smaller amplitude
                            :   5 - 5 times larger amplitude
}

ASSIGNED {
    : Variables that are assigned outside the mod file
    v               (mV)
    celsius         (degC)
	cai		        (mM)    : intracellular [Ca++] (mM)
	cao		        (mM)    : extracellular [Ca++] (mM)

	: GLOBAL variables that are assigned in the INITIAL block
    phim            (1)     : temperature adjustment to taum
    phih            (1)     : temperature adjustment to tauh

	: RANGE variables that are assigned in the INITIAL & DERIVATIVE blocks
    minf            (1)     : steady state value of activation gating variable
    hinf            (1)     : steady state value of inactivation gating variable
    taum            (ms)    : time constant for activation
    tauh            (ms)    : time constant for inactivation

	: RANGE variables that are assigned in the BREAKPOINT block
    ica             (mA/cm2): calcium current (mA/cm2) generated
}

STATE {
	m		        (1)		: activation gating variable
	h		        (1)		: inactivation gating variable
}

INITIAL {
    : Equations for taum & tauh were taken from Huguenard & McCormick, 1992 
    :   and were at 23 degC
    : Transformation to celsius using Q10 (qm & qh)

	: Calculate constants from temperature
    phim = qm ^ ( ( celsius - 23 (degC) ) / 10 (degC) )
    phih = qh ^ ( ( celsius - 23 (degC) ) / 10 (degC) )

	: Initialize state variables to steady state values
	evaluate_fct(v)			: calculate steady state values of state variables
    m = minf
    h = hinf
}

BREAKPOINT {
	: Update gating variables m & h from voltage v
    SOLVE castate METHOD cnexp

	: Calculate calcium current ica from m, h, cai, cao, v
    ica = pcabar * m*m*h * ghk(v, cai, cao)
}

DERIVATIVE castate {
	: Update minf, hinf, taum, tauh from v
    evaluate_fct(v)

	: Update m & h from minf, hinf, taum, tauh
    m' = (minf - m) / taum
    h' = (hinf - h) / tauh
}

PROCEDURE evaluate_fct(v (mV)) {
    LOCAL ampFactor

	: Update minf, hinf from v
    minf = 1.0 / ( 1 + exp( (v + 57 (mV) - shiftm) / (-6.2 (mV) * slopem) ) )
    hinf = 1.0 / ( 1 + exp( (v + 81 (mV) - shifth) / (4.0 (mV) * slopeh) ) )

	: Update taum from v
    taum = ( 0.612 (ms) + 1.0 (ms) / 
                                ( exp( (v + 132 (mV) - shiftm) / -16.7 (mV)) 
                                + exp( (v + 16.8 (mV) - shiftm) / 18.2 (mV)) ) ) 
            / (phim * slopem)

    : Update tauh from v
    if (tauhmode == 1) {
        tauh = taum
    } else {
        if (tauhmode == 0) {
            ampFactor = 1
        } else {
            if (tauhmode == 2) {
                ampFactor = 0.1
            }
            if (tauhmode == 3) {
                ampFactor = 10
            }
            if (tauhmode == 4) {
                ampFactor = 2/3
            }
            if (tauhmode == 5) {
                ampFactor = 1.5
            }
            if (tauhmode == 6) {
                ampFactor = 1/2
            }
            if (tauhmode == 7) {
                ampFactor = 2
            }
        }

        if ( v < -80 + shifth) {
            tauh = ampFactor * 
                    1.0 (ms) * exp((v + 467 (mV) - shifth) / 66.6 (mV))
                    / (phih * slopeh)
        } else {
            tauh = ampFactor * 
            ( 28 (ms) + 1.0 (ms) * exp( (v + 22 (mV) - shifth) / -10.5 (mV)) ) 
                    / (phih * slopeh)
        }        
    }
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

FUNCTION nongat(v (mV), cai (mM), cao (mM)) (mA/cm2) {    : non gated current
    nongat = pcabar * ghk(v, cai, cao)
}

