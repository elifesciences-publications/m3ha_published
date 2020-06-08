TITLE Fast mechanism for submembranal Ca++ concentration (cai) cad

COMMENT

    Takes into account:
    
        - increase of cai due to calcium currents
        - extrusion of calcium with a simple first order equation
    
    This mechanism is compatible with the calcium pump "cad" and has the 
    same name and parameters; however the parameters specific to the pump
    are dummy here.
    
    Parameters:
        - depth: depth of the shell just beneath the membrane (in um)
        - cainf: equilibrium concentration of calcium (2e-4 mM)
        - taur: time constant of calcium extrusion (must be fast)
    
    Written by Alain Destexhe, Salk Institute, 1995

    Update kinetics from Sohal & Huguenard 2003:
    "These kinetics reproduced the burst AHP observed in 
    several studies of RE cells (Avanzini et al., 1989; 
    Bal and McCormick, 1993; Contreras et al., 1993)."

    
    2017-02-07 Modified the Faraday constant
    2017-02-07 Parameters modified to be consistent with Sohal & Huguenard 2003
    2017-03-07 Initialize cai from hoc
    2017-07-22 Copied from RTCl
    2017-07-22 Changed tabs to spaces
    2017-07-25 AL - Added units and comments

ENDCOMMENT

NEURON {
    SUFFIX cad
    USEION ca READ ica, cai WRITE cai
    RANGE depth, taur, cainf
    RANGE drive_channel
}

UNITS {
    (molar) = (1/liter)     : moles do not appear in units
    (mM)    = (millimolar)
    (um)    = (micron)
    (mA)    = (milliamp)
    (msM)   = (ms mM)
    FARADAY = (faraday) (coulombs)
}

PARAMETER {
    : RANGE variables whose values are specified in hoc
    depth   = .1    (um)    : depth of shell just beneath the membrane (um)
    taur    = 24    (ms)    : Ca++ removal time constant (ms), 
                            :   Sohal & Huguenard 2003 (@ 34 degC)
                            :   Note: Destexhe used 5 ms, Amarillo used 1 ms
    cainf   = 2.4e-4(mM)    : steady state intracellular [Ca++] (mM)

}

ASSIGNED {
    : Variables that are assigned outside the mod file
    ica             (mA/cm2): calcium current

    : RANGE variables that are assigned in the DERIVATIVE block
    drive_channel   (mM/ms) : calcium flux due to ica
}
    
STATE {
    cai             (mM)    : submembranal [Ca++] (mM)
}

BREAKPOINT {
    : Calculate cai from ica
    SOLVE state METHOD derivimplicit
}

DERIVATIVE state { 
    : Calculate flux due to ica
    drive_channel =  - (10000) * ica / (2 * FARADAY * depth)

    : Since the pump is outward, do not let cai decrease from ica
    if (drive_channel <= 0.) { 
        drive_channel = 0. 
    }

    : Calculate change in calcium concentration
    cai' = drive_channel + (cainf - cai) / taur
}

