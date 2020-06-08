: dummy.mod
: Makes functions written in C or NMODL available to hoc
: 
: 2017-02-22 Combined with asin.mod, which was adapted from Jedlicka et al 2011

VERBATIM
#include <math.h>
#include <time.h>
ENDVERBATIM

NEURON {
   SUFFIX nothing
}

: Rounds up z
: Call it gceil to avoid naming conflicts
FUNCTION gceil(z) {
   gceil = ceil(z)
}

: Rounds down z
FUNCTION gfloor(z) {
   gfloor = floor(z)
}

: Computes the inverse sine function
FUNCTION my_asin(arg) {
	VERBATIM
	double ret;
	ret=asin(*getarg(1));
	ENDVERBATIM
	my_asin=ret
}

: Computes the sine function
FUNCTION my_sin(arg) {
	VERBATIM
	double ret;
	ret=sin(*getarg(1));
	ENDVERBATIM
	my_sin=ret
}

: %%% TO EXAMINE
FUNCTION exp_i(arg) {
	VERBATIM
	double ret=1.0;
	double euler=exp(1.0);
	int n=0;
	for (n=0;n<(*getarg(1));n++) {
		ret*=euler;
	}
	ENDVERBATIM
	exp_i=ret
}

: %%% TO EXAMINE
FUNCTION mytime() {
	VERBATIM
	double ret = 0.0;
	ret = (double)time(0)/3600.0;
	ENDVERBATIM
	mytime = ret
}

