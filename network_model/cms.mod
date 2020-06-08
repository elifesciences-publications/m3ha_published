TITLE A mechanism for checking for maximum slope and firing events

COMMENT

    Helper TC neuron for determining when to send spikes to connected neurons

    2018-03-28 Adapted from CheckMaxSlope.tem and HelperTC.mod
    2018-04-04 Fixed bugs
    2018-04-04 Forced firstDelayTime and intervalToFire to be nonnegative
    2018-04-26 Forced intervalToFire to be at least 1 ms

ENDCOMMENT

NEURON {
    POINT_PROCESS cms
    GLOBAL a1, b1, c1, a2, b2, c2, a3, b3, c3
    RANGE nSpikesToFire, firstDelayTime, intervalToFire
    RANGE inRefractory, inSlopeWatching, periodRefractory, t0
    RANGE lastIPSCtime, delay, slope, slopeLast
    RANGE vLast2, vLast, dtLast
}

UNITS {
    (mV)    = (millivolt)
}

PARAMETER {
    : GLOBAL variables whose values are fixed
    : Linear model coefficients taken from data_dclamp/linear_models.pdf
    a1 = -0.6809448     (1)         : spikes per peak intercept
    b1 = 1.3957288      (s/V)       : spikes per peak per maxslope increase
    c1 = 8.3791814e-4   (/ms)       : spikes per peak per delay increase
    a2 = 2.8054103      (ms)        : spike latency intercept
    b2 = -0.5685347     (ms s/V)    : spike latency per maxslope increase
    c2 = -2.0618266e-4  (1)         : spike latency per delay increase
    a3 = 333.8204802    (/s)        : spike frequency intercept
    b3 = 8.7678243      (/V)        : spike frequency per maxslope increase
    c3 = -0.0819026     (/s /ms)    : spike frequency per delay increase
}

ASSIGNED {
    : Variables that are assigned outside the mod file
    v                   (mV)        : postsynaptic membrane potential
    dt                  (ms)        : interval of the current time step

    : RANGE variables that are assigned in the INITIAL & NET_RECEIVE blocks
    nSpikesToFire       (1)         : number of spikes yet to be fired
    firstDelayTime      (ms)        : delay from max slope to first spike
    intervalToFire      (ms)        : interval between successive spikes
    inRefractory        (1)         : whether in a "refractory period"
    inSlopeWatching     (1)         : whether in a "slope watching period"
    periodRefractory    (ms)        : refractory period length
    t0                  (ms)        : time of latest spike sent
    lastIPSCtime        (ms)        : the time of the last IPSC
    delay               (ms)        : the time since the last IPSC
    slope               (V/s)       : the slope of the current time step
    slopeLast           (V/s)       : the slope of the last time step

    : RANGE variables that are assigned in the INITIAL & BREAKPOINT blocks
    vLast2              (mV)        : the last last voltage
    vLast               (mV)        : the last voltage
    dtLast              (ms)        : interval of the last time step
}

INITIAL {
    : Initialize variables
    nSpikesToFire    = 0
    firstDelayTime   = 0
    intervalToFire   = 0          
    inRefractory     = 0
    inSlopeWatching  = 0
    periodRefractory = 0
    t0               = 0
    lastIPSCtime     = 0
    slope            = 0
    slopeLast        = 0
    vLast2           = 0
    vLast            = 0
    dtLast           = 0
}

BREAKPOINT {
    : Compute the slope (dv/dt) over the current time step
    :   Note: dt is the magnitude of the current time step
    :   vLast is the voltage value saved from the last time step
    slope = (v - vLast) / dt

    : Compute the slope (dv/dt) over the last time step
    slopeLast = (vLast - vLast2) / dtLast

    : Store the voltage, dt and slope values for future use
    vLast2 = vLast
    vLast = v
    dtLast = dt
}

NET_RECEIVE (weight (1)) {
    : Note: flag is an internal variable in NET_RECEIVE that is modified by
    :       net_send(). By default (external events), flag == 0
    : Weight is not used here

    : Deal with external events
    if (flag == 0 && inRefractory == 0) {   : An external event arrives while
                                            :   not in a "refractory period"
        : Start the "refractory period"
        inRefractory = 1

        : Update the last IPSC time (ms)
        lastIPSCtime = t

        : Send a self event after 250 ms with flag == 1 to start slope watching
        :   Note: The smallest maximum slope time in the dynamic clamp data
        :           is 288.2 ms
        net_send(250, 1)
    } else if (flag == 0 && inSlopeWatching == 0) {
        : Update the last IPSC time (ms)
        lastIPSCtime = t
    } else {
        : Do nothing
    }

    : Deal with self events
    if (flag == 1) {                : start slope watching
        : Start the "slope watching period"
        inSlopeWatching = 1

        : Send a self event with flag == 2 when the slope reaches a maximum
        :   This happens when the slope decreases
        WATCH (slope < slopeLast) 2
    } else if (flag == 2) {         : slope reaches a maximum
        : Terminate the "slope watching period"
        inSlopeWatching = 0

        : Compute the delay time (ms) since the last IPSC
        delay = t - lastIPSCtime

        : Compute the number of spikes yet to be fired
        nSpikesToFire = spikesperpeak(slopeLast, delay)

        : Compute the delay time for the first spike
        if (nSpikesToFire > 0) {
            firstDelayTime = spikelatency(slopeLast, delay)
            : Force the delay time to be nonnegative
            if (firstDelayTime < 0) {
               firstDelayTime = 0
            }
        }

        : Compute the interval between spikes
        if (nSpikesToFire > 1) {
            intervalToFire = 1000 (ms/s) / 
                                spikefrequency(slopeLast, delay)
            : Force the interval to fire to be positive
            if (intervalToFire < 0) {
                : At least 1 ms
                intervalToFire = 1
            }
        }

        : Compute the length of the refractory period
        periodRefractory = firstDelayTime + 
                            (nSpikesToFire - 1) * intervalToFire

        : Fire an action potentials according to the computed delay times
        : Note: nSpikesToFire will be 0 after this block
        if (nSpikesToFire > 0) {
            : Fire the first spike
            net_event(t + firstDelayTime)

            : Update the time of latest spike sent
            t0 = t + firstDelayTime

            : Decrement the number of spikes to fire
            nSpikesToFire = nSpikesToFire - 1

            : Fire subsequent spikes if any
            while (nSpikesToFire > 0) {
                : Fire the next spike
                net_event(t0 + intervalToFire)

                : Update the time of latest spike sent
                t0 = t0 + intervalToFire

                : Decrement the number of spikes to fire
                nSpikesToFire = nSpikesToFire - 1                
            }
        }

        : Send a self event with flag == 3 and a delay of periodRefractory, 
        :   which could be 0
        net_send(periodRefractory, 3)
    } else if (flag == 3) {                 : refractory period ends
        : Terminate the "refractory period"
        inRefractory = 0

        : Reset number of spikes to fire and interval to fire
        nSpikesToFire = 0
        firstDelayTime = 0
        intervalToFire = 0
        periodRefractory = 0
    } else {
        : Do nothing
    }
}

FUNCTION spikesperpeak(maxslope (V/s), delay (ms)) (1) {
    spikesperpeak = a1 + b1 * maxslope + c1 * delay
}

FUNCTION spikelatency(maxslope (V/s), delay (ms)) (ms) {
    spikelatency = a2 + b2 * maxslope + c2 * delay
}

FUNCTION spikefrequency(maxslope (V/s), delay (ms)) (/s) {
    spikefrequency = a3 + b3 * maxslope + c3 * delay
}

COMMENT

OLD CODE:
    : Send a self event to initiate the WATCH statement 
    :   See https://www.neuron.yale.edu/phpBB/viewtopic.php?t=1044
    net_send(0, 1)

    t                   (ms)        : current time
    : Moved here because of this error:
    `Multiple declaration of t at line 200 in file cms.mod`

    : Fire a spike now (for DEBUG)
    net_event(t)

    : For cmstest
    nSpikesToFire = 2
    firstDelayTime = 300 (ms)
    intervalToFire = 1000 (ms/s) / 300 (/s)

    : Force the interval to fire to be nonnegative
    if (intervalToFire < 0) {
       intervalToFire = 0
    }

ENDCOMMENT