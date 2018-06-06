#pragma once
#include <iomanip>

namespace lbfgs {
	using namespace std;
	#include "lbfgsb.c.h"
	#include "optimizer.types.h"
}

#include <vector>

// Class to wrap the lbfgsb routine
class optimizer {
public:
	typedef lbfgs::vecD vecD;
	struct limit { double *vmin, *vmax; };
	typedef lbfgs::vecL  vecL;
	typedef unsigned long ulong;
	typedef lbfgs::options options;			// user settable options for the optimizer
	struct ResultProvonance {
		double	objective,objectivePrev;
		vecD	pVals,pDerivs;
		ulong	iters,numEvals;
		double	pgtol, factr;				// tolerances used to find extrema factr is in units of Double machine precision (usually 1E2 to 1E8), pgtol is a minimum gradient.
		char	lastmsg[32];
		ResultProvonance() { clear(); };
		void	clear() { iters = numEvals = 0; pgtol = factr = objective = objectivePrev = 0.0; pVals.clear(); pDerivs.clear(); lastmsg[0] = '\0'; };
	};
private:
	typedef lbfgs::wrkBuffer wrkBuffer;		// internal working buffer
	bool optimize_iterate(wrkBuffer &w, const options &opts);
public:
	// Useful utilities....
	inline static int sstrcmp(const char *A, const char *B) { return(strncmp(A, B, strlen(B))); }
	template< class T >	inline static void vec2array(const std::vector<T> &V, T *A) { for (uint32_t i = 0; i < V.size(); i++) A[i] = V[i]; }
	template< class T >	inline static void array2vec(const T *A, std::vector<T> &V) { for (uint32_t i = 0; i < V.size(); i++) V[i] = A[i]; }

	// at least calcValue must be implimented... others can be derived from it, if needed
	virtual bool calcValue(const void *data, const vecD &parameters, double &Value) { return false; };

	// override with an analytic deravitive to improve accuracy and speed calculation.
	virtual bool calcDeriv(const void *data, const vecD &parameters, vecD &Deravitives); 

	// Sometimes it is most efficient to calculate deriv & value at same time. If so, override, if not, deligate...
	virtual bool calcValueAndDeriv(const void *data, const vecD &parameters,  double &Value, vecD &Deravitives ) {
		bool ok = calcValue(data, parameters, Value);
		if (ok) ok = calcDeriv(data, parameters, Deravitives);
		return(ok);
	}

	static inline double mylog(double p) { return(p > .5 ? log1p(p - 1) : log(p)); };
	bool optimize( const void *data, const vecD &paramIn, const vecL &limits, const options &opts, vecD &paramOut, ResultProvonance *prov = NULL );
	//
};

// override with an analytic deravitive to improve accuracy and speed calculation.
bool optimizer::calcDeriv(const void *data, const vecD &parameters, vecD &Deravitives) {
	// TODO: a parabolic approximation requiring 3 evaluations might be MUCH better here....
	bool ok=true;
	double delta = 1.0e-8, p, p1, p2, f1, f2;
	vecD	tmpP = parameters;
	for (ulong i = 0; i < parameters.size(); i++) {
		p = tmpP[i]; p1 = p - p*delta; p2 = p + p*delta;
		if (p1 == p2) { p1 = -delta; p2 = delta; };
		tmpP[i] = p1; ok = calcValue(data, tmpP, f1); if (!ok) return(ok);
		tmpP[i] = p2; ok = calcValue(data, tmpP, f2); if (!ok) return(ok);
		tmpP[i] = p; // restore original value...
		Deravitives[i] = (f2 - f1) / (p2 - p1);
	}
	return(ok);
}

bool optimizer::optimize_iterate(optimizer::wrkBuffer &w, const optimizer::options &opts) {
	bool done = true;
	//int r;
	double factr = opts.factr, pgtol = opts.pgtol;	// needed to get around "const" compatibility issue with fortran

	// return value is always 0, it is s ok to ignore it.
	(void) lbfgs::setulb_(w.n, w.m, w.parameters, w.boundl, w.boundu, w.boundf, &w.f, w.gradient, &factr, &pgtol, w.wa, w.wai, w.task, w.csave, w.lsave, w.isave, w.dsave);
	if (w.isTask("FG") || w.isTask("NEW_X")) done = false;
	return(done);
}

bool optimizer::optimize(const void *data, const vecD &paramIn, const vecL &limits, const options &opts, vecD &paramOut, ResultProvonance *prov ) {
	bool ok, done;
	wrkBuffer wrk;
	int	icount = 0;
	vecD best_param = paramIn, best_grad=paramIn;
	double best_obj = 0;


	// temporary values
	std::vector<double> wrkparam = paramIn;
	std::vector<double> wrkgrad = paramIn;
	int N= (int)wrkparam.size();

	ok = wrk.init(paramIn, limits, opts.m); 
	if (!ok) return(false);	// TODO handle error

	// set initial task, value and gradient,.,
	ok = calcValueAndDeriv(data, wrkparam, wrk.f, wrkgrad);
	if (!ok) return(false);	// TODO handle error
	best_obj = wrk.f; best_param = wrkparam; best_grad = wrkgrad;
	vec2array<double>(wrkgrad, wrk.gradient);

	wrk.setTask("START"); done = false;
	while (!done) {
		done = optimize_iterate(wrk, opts);
		icount++;
		// most common results, 
		if (wrk.isTask("FG") || wrk.isTask("NEW_X") ) {
			// Optimizer wants new F(X) and dF/dx
			if (opts.verbose) { std::cerr << "FG-befor " << icount << "\t" << wrk.f << "\t"; 
				for (int i = 0; i < N; i++)  std::cerr << std::setprecision(16) << wrkparam[i] << "\t";
				for (int i = 0; i < N; i++)  std::cerr << wrkgrad[i] << "\t";
				std::cerr << wrk.task << std::endl; };
			if (wrk.f < best_obj) { best_obj = wrk.f; best_param = wrkparam; best_grad = wrkgrad; };
			array2vec<double>(wrk.parameters, wrkparam);
			ok = calcValueAndDeriv(data, wrkparam, wrk.f, wrkgrad);
			// std::cerr << "FG " << icount << "\t" << wrk.f << "\t";
			if (!ok) return(ok);
			if (wrk.f < best_obj) { best_obj = wrk.f; best_param = wrkparam; best_grad = wrkgrad; };
			vec2array<double>(wrkgrad, wrk.gradient);
			if (opts.verbose) { std::cerr << "FG-after " << icount << "\t" << wrk.f << "\t"; 
				for (int i = 0; i < N; i++)  std::cerr << wrkparam[i] << "\t";
				for (int i = 0; i < N; i++)  std::cerr << wrkgrad[i] << "\t";
				std::cerr << std::endl; };
			ok = true; done = false;
		} else if (wrk.isTask("NEW_X")) {
			// An iteration has concluded, user can stop or continue...
			// TODO: Add verbose statement
			if (opts.verbose) { std::cerr << "NEW_X " << icount << "\t" << wrk.f << "\n"; };
			if (opts.verbose > 2) {}
		} else if (wrk.isTask("CONV")) {
			// CONV Results converged, we are done
			if (opts.verbose) { std::cerr << "CONV " << icount << "\t" << wrk.f << "\n"; };
			ok = true; done = true;
		}
		else if (wrk.isTask("ABNO")) {
			// ABNO - abnormal termination
			if (opts.verbose) { std::cerr << "ABNO  " << icount << "\t" << wrk.f << "\n"; };
			ok = false; done = true;
		} else if (wrk.isTask("ERROR")) {
			// ERROR - error, generally in input parameters...
			if (opts.verbose) { std::cerr << "ERR " << icount << "\t" << wrk.f << "\n"; };
			ok = false; done = true;
		} else {
			// Unknown response, this is an error
			std::cerr << "butils::optimizer::optimize::optimize_iterate() returns: Unknown response (UNKO) - this is an error. " << icount << "\t" << wrk.f << "\n";
			ok = false; done = true;
		}

		// Alternatively, estimator has not reached termination with existing gradient, continue estimation... Never seen it take more than 100, so kill it at 1000.
		if (icount > 1000) { 
			std::cerr << "butils::optimizer::optimize More than 1000 iterations - this is an error. " << icount << "\t" << wrk.f << "\n";
			ok = false; done = true; }
	}

	if (ok) {
		array2vec<double>(wrk.parameters, paramOut);
	} else {
		// if we technically fail to converge, return the best solution that we found...
		wrk.f= best_obj; 
		paramOut = wrkparam = best_param; 
		wrkgrad = best_grad;
	};

	if (prov != NULL) {
		prov->clear();
		strncpy(prov->lastmsg, wrk.task, 32); prov->lastmsg[31] = '\0';
		prov->objective = wrk.f; prov->objectivePrev = wrk.dsave[2];
		prov->iters = wrk.isave[30]; prov->numEvals = wrk.isave[34];
		prov->pVals = wrkparam; prov->pDerivs = wrkgrad;
	};
	return(ok);
}

