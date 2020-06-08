TITLE Helper TC neuron

COMMENT

    Helper TC neuron for determining when to send spikes to connected neurons

    2018-03-24 Created by Adam Lu

ENDCOMMENT

NEURON {
    ARTIFICIAL_CELL HelperTC
    GLOBAL a1, b1, c1, a2, b2, c2, a3, b3, c3
    RANGE nSpikesToFire, firstDelayTime, intervalToFire
    RANGE inRefractory, periodRefractory, t0
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
    : RANGE variables that are assigned in the INITIAL & NET_RECEIVE blocks
    nSpikesToFire       (1)         : number of spikes yet to be fired
    firstDelayTime      (ms)        : delay from max slope to first spike (ms)
    intervalToFire      (ms)        : interval between successive spikes (ms)
    inRefractory        (1)         : whether in a "refractory period"
    periodRefractory    (ms)        : refractory period length (ms)
    t0                  (ms)        : time of latest spike sent
}

INITIAL {
    : Initialize variables
    nSpikesToFire    = 0
    firstDelayTime   = 0
    intervalToFire   = 0          
    inRefractory     = 0
    periodRefractory = 0
    t0               = 0
}

NET_RECEIVE (maxslope (V/s), delay (ms)) {
    : maxslope  - value of maximum slope (V/s)
    : delay     - time since last IPSC (ms)
    : Note: flag is an internal variable in NET_RECEIVE that is modified by
    :       net_send(). By default (external events), flag == 0

    if (flag == 0 && inRefractory == 0) {   : An external event received 
                                            :   while not in a refractory period
        : Compute the number of spike to fire
        nSpikesToFire = spikesperpeak(maxslope, delay)

        : Compute the delay time for the first spike
        if (nSpikesToFire > 0) {
            firstDelayTime = spikelatency(maxslope, delay)
        }

        : Compute the interval between spikes
        if (nSpikesToFire > 1) {
            intervalToFire = 1000 (ms/s) / 
                                spikefrequency(maxslope, delay)
        }

        : Fire an action potentials according to the computed delay times
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

        : If spiking, put the cell in a refractory period
        :   and send a self event when it is supposed to end
        if (nSpikesToFire > 0) {
            : Compute the length of the refractory period
            periodRefractory = firstDelayTime + 
                                (nSpikesToFire - 1) * intervalToFire

            : Send a self event with flag == 1 and a delay of periodRefractory
            net_send(periodRefractory, 1)        
        }
    } else if (flag == 1) {                 : A self event is received
        : Terminate the refractory period
        inRefractory = 0

        : Reset number of spikes to fire and interval to fire
        nSpikesToFire = 0
        firstDelayTime = 0
        intervalToFire = 0
        periodRefractory = 0
    } else {                                : In a "refractory period""
        : Ignore the external event
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
