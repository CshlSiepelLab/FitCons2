#pragma once

#include <vector>
// Derive a class from optimizer for inferring Betas...
#include "butils/optimizer.h"
#include "BlockSet.h"					// data type definitions
#include "modelInsight.h"				// base model class
#include "modelInsightAccelerator.h"	// accelerated likelihood calculations

#include <iostream>	// debugging

class optBeta : public butils::optimizer {
public:
	typedef butils::optimizer optimizer;
	typedef optimizer::vecD vecD;
	typedef optimizer::vecL vecL;
	typedef inscomp::modelInsight modelInsight;
	typedef inscomp::modelInsightAcceleratorBeta modelInsightAcceleratorBeta;

	struct  options : optimizer::options { };
private:
	modelInsightAcceleratorBeta mAccel;
	void clear() { mAccel.clear(); };
public:
	optBeta( const modelInsight &M ) : mAccel(M) {	}
	virtual ~optBeta() { clear(); };

	virtual bool calcValue(const void *data, const vecD &parameters, double &Value) {
		if (!mAccel.initialized) return false;
		modelInsightAcceleratorBeta::parameters p; p = mAccel.getparameters(); p.beta.b3 = parameters[0]; p.beta.b1 = 1.0 - p.beta.b2 - p.beta.b3; mAccel.setparameters(p);
		Value = mAccel.calcNLL(); 
		return true;
	};

	// override with an analytic deravitive to improve accuracy and speed calculation.
	virtual bool calcDeriv(const void *data, const vecD &parameters, vecD &Deravitives) {
		double tmp_val;
		bool ok = calcValueAndDeriv(data, parameters, tmp_val, Deravitives);
		return( ok );
	};

	// Sometimes it is most efficient to calculate deriv & value at same time. If so, override, if not, delegate...
	virtual bool calcValueAndDeriv(const void *data, const vecD &parameters, double &Value, vecD &Deravitives) {
		if (!mAccel.initialized) return false;
		modelInsightAcceleratorBeta::parameters p; p=mAccel.getparameters(); p.beta.b3 = parameters[0]; p.beta.b1=1.0- p.beta.b2- p.beta.b3; mAccel.setparameters( p );
		Value = mAccel.calcNLL( &p ); Deravitives[0] = p.beta.b3;
		return true;
	};

	double estimateBeta2(const BlockSet_t &bs, double PriorWeight = 0, const modelInsight *modelPrior = NULL) const {
		double beta2 = -1, hf=0, lf=0;	// an impossible value.

		// calculate beta2 directly...
		{	uint32_t h = 0, l = 0; bs.polyCountsSites(h, l); hf = (double) h; lf = (double) l;	}

		double polycounts = hf + lf;
		if ( (PriorWeight > 0) && (modelPrior != NULL) ) {
			hf			+= PriorWeight * modelPrior->getparameters().beta.b2;
			polycounts	+= PriorWeight;
		};
		if (polycounts > 0)	beta2 = hf / polycounts; // usually around .20 (.15 to .25)
		return beta2;
	}

	bool optimize(const BlockSet_t &bs, modelInsight &model, double PriorWeight = 0, const modelInsight *modelPrior = NULL, void *provonance = NULL) {
		typedef butils::optimizer::ResultProvonance ResultProvonance;
		typedef butils::lbfgs::vecD	vecD;
		clear();

		// extract parameters from source model (init) and destimation model (refine)
		modelInsight::parameters params = model.getparameters();

		// estimate beta2
		params.beta.b2 = estimateBeta2( bs, PriorWeight, modelPrior );
		if ((params.beta.b2<=0) || (params.beta.b2 >=1)) { model.setparameters(params); return(false); };
		// renormalize refiend values....
		params.beta.b1 = (1- params.beta.b2) * params.beta.b1 / (params.beta.b1 + params.beta.b3);
		params.beta.b3 = (1- params.beta.b2) * params.beta.b3 / (params.beta.b1 + params.beta.b3);

		// Initialize high performance likelihood / AP calculation to find beta3
		mAccel.setparameters( params );
		mAccel.initialize( bs, PriorWeight, modelPrior );

		double vmin=0, vmax = 1.0 - params.beta.b2;									// range for beta3
		vecL limits(1); limits[0].vmin = &vmin; limits[0].vmax = &vmax;	

		vecD paramIn(1), paramOut(1); paramIn[0] = params.beta.b3; paramOut[0]=0;	// initial values
		ResultProvonance prov_tmp, *prov = (provonance == NULL ? &prov_tmp : (ResultProvonance *)provonance);
		options op;	op.factr = 1000.0; op.pgtol = 1e-8;								// factr, tolerance relative to machine precision, pgtol, terminate on low gradient.

		// Call optimizer to find ML / MAP value for beta3
		bool ok = optimizer::optimize((const void *) &mAccel, paramIn, limits, op, paramOut, prov);
		// we may fail a technical test for convergence, but if the solution looks good, use it.
		if (!ok && ((prov->objectivePrev / prov->objective)< 1.0e-6) ) ok = true;
		if (ok) {
			// save results in modelRefine
			params.beta.b3 = paramOut[0]; params.beta.b1 = 1.0 - params.beta.b2 - params.beta.b3;
			model.setparameters( params );
		};
		
		clear(); return(ok);						// release AUX values
	}
};





#ifdef OLD_CODE
#ifdef OLD_CODE
// precomputed cofactors for each term
struct beta_aux_t {
#ifndef V0
	double A;		// cofactor of beta3
	double B;		// constant term
	double Nobs;	// Number of times this was oberved in current block.
	inline void load(double lambda, double zxmaj, double zxmin, double nobs) {
		double asb1, asb3;
		// calculate cofactor for RELATIVE betas, that is beta1 & beta3, where (beta1+beta3=1), renorm later to account for beta2
		// asb1 = zxmaj * (1.0 - lambda) + (1.0 - zxmaj)*lambda / 3.0;	// cofactor for beta1
		// asb3 = zxmin * (1.0 - lambda) + (1.0 - zxmin)*lambda / 3.0;	// cofactor for beta3 (beta3=1-beta1)
		asb1 = zxmaj * (1.0 - lambda) + (zxmin)*lambda / 3.0;	// cofactor for beta1
		asb3 = zxmin * (1.0 - lambda) + (zxmaj)*lambda / 3.0;	// cofactor for beta3 (beta3=1-beta1)
																// usually beta3 is on the order of .04, (.02-.07), so looking for beta3 is more stable than looking for beta1
		A = (asb3 - asb1); B = asb1; Nobs = nobs;
		// A = asb1; B = asb3; Nobs = nobs;
	}
#else
	double lam;		// cofactor of beta3
	double NzeqXmaj;		// constant term
	double NzeqXmin;	// Number of times this was oberved in current block.
	inline void load(double lambda, double zxmaj, double zxmin, double nobs) {
		lam = lambda;
		NzeqXmaj = zxmaj * nobs; NzeqXmin = zxmin * nobs;
	}
#endif
};
// we dont need a double vector... we abastract away all block info
//	so one long vector will do!
std::vector<std::vector<beta_aux_t>> auxen;
#endif

bool optBeta::initialize(const void *data, const vecD &paramIn, int V) {
	using butils::ulong;
	const BlockSet_t &bsn = *(BlockSet_t *)data;	// type casting
	ulong nblocks = (ulong)bsn.blocks.size();		// 64-> 32 bit

	if (initialized) { clear(); };					// clears initialized flag...
	auxen.resize(nblocks);							// Deeply appreciating vectors....

	BlockSet_t::AncPriInd_t		imaj, imin;
	BlockSet_t::AncPriCount_t	nobs;
	double zxmaj, zxmin, lambda;
	double count = 0, wcount = 0, ab1s = 0, ab3s = 0;
	double sxma = 0, sxmi = 0, slam = 0, sobs = 0;
	for (ulong iblock = 0; iblock < nblocks; iblock++) {
		auto &bl = bsn.blocks[iblock];								// get current block from database
		auto &ab = auxen[iblock]; ab.resize(bl.polyLNumPat());		// resize aux block to same number of positions
		bl.getHeader(NULL, NULL, NULL, &lambda);					// get lambda from block... can provide first args if debugging is desired...
		for (ulong ipoly = 0; ipoly < ab.size(); ipoly++) {
			bl.polyLVal(ipoly, imaj, imin, nobs);
			zxmaj = *(bsn.ancPri[imaj]); zxmin = *(bsn.ancPri[imin]);
			ab[ipoly].load(lambda, zxmaj, zxmin, nobs);
			sxma += nobs * zxmaj; sxmi += nobs * zxmin; slam += lambda * nobs; sobs += nobs;
		}
	}

	sxma /= sobs; sxmi /= sobs; slam /= sobs;
	std::cout << "expected: xma " << sxma << "  smi " << sxmi << " slam " << slam << " nobs = " << sobs << std::endl;

	// add pseudocounts to account for prior.
	if (op.priWeight > 0.0) {
		std::vector<beta_aux_t> tmp_pri; tmp_pri.clear();
		beta_aux_t				tmp_pos;
		double b, b1, b2, b3, zxmaj, zxmin, zxoth, lambdaT = op.priLambdaT, thetaA = op.priThetaA;	// really, we should use expected lambda on NEUTRAL poly blocks.
		b = op.priBeta1 + op.priBeta2 + op.priBeta3; b1 = op.priBeta1 / b; b2 = op.priBeta2 / b; b3 = op.priBeta3 / b;
		// calculate distributions based on prior-values, the ML of parameters value for these observations should BE the prior parameters themselves.
		// optsinight:Likelihiid, class = false
		zxmaj = 3.0*((1 - lambdaT)*b1 + lambdaT*b3 / 3.0)*thetaA / 3.0; zxmin = 3.0*((1 - lambdaT)*b3 + lambdaT*b1 / 3.0)*thetaA / 3.0; zxoth = 6.0*(b1 + b3)*(lambdaT / 3.0)*(thetaA / 3.0);
		// optsinight:Likelihiid, class = true
		// double zpost=zxmaj*3.0+zxmin*3.0+zxoth*6.0, zscale = (zxmaj*3.0 + zxmin*3.0)/zpost;
		double z = zxmaj + zxmin + zxoth; zxmaj /= z; zxmin /= z; zxoth /= z;
		double pcounts = op.priWeight *(b1 + b3);

		// proWeight is hte weight of the prior, in observations
		//	(b1+b3) is the fraction of obs that are Low freq (over High freq)
		//	zscale is the fraction of obs that are NOT Z=Xoth (that is Z=Xmaj or Z=Xmin)
		//tmp_pos.load( lambdaT, zxmaj,zxmin, op.priWeight *(b1+b3)*zscale); tmp_pri.push_back( tmp_pos );
		tmp_pos.load(lambdaT, 1.0, 0.0, zxmaj * pcounts); tmp_pri.push_back(tmp_pos);
		tmp_pos.load(lambdaT, 0.0, 1.0, zxmin * pcounts); tmp_pri.push_back(tmp_pos);
		auxen.push_back(tmp_pri);
	};

	initialized = true;
	return(true);
}


// parameters beta3: Really, this is beta3 / (beta1 + beta3)

// auxStats[site].beta1 = pZeqXmaj*(1 - lambda[block]) + (1 - pZeqXmaj)*lambda[block] / 3;
// auxStats[site].beta3 = pZeqXmin*(1 - lambda[block]) + (1 - pZeqXmin)*lambda[block] / 3;
// qBeta1 = beta1       * auxStats[site].beta1;
// qBeta3 = (1 - beta1) * auxStats[site].beta3;
// qSum = qBeta1 + qBeta3;
// We want to rewrite this in qSum = A*beta3 +  B form for eae of optimization
bool optBeta::calcValue(const void *data, const vecD &parameters, double &Value) {
	// calc val
	double beta3 = parameters[0], cum = 0;
	ulong nblocks = (ulong)auxen.size();	// 64-> 32 bit...
	for (ulong iblock = 0; iblock < nblocks; iblock++) {
		auto &ab = auxen[iblock];
		for (ulong ipoly = 0; ipoly < ab.size(); ipoly++) {
			const beta_aux_t &v = ab[ipoly];
			// can speed this up with prob hack, multiple untill resutls would denormalize 
			// then log add.... mult is much faster... and may get a hundred mult per log...
#ifndef V0
			cum += mylog(beta3 * ab[ipoly].A + ab[ipoly].B) * ab[ipoly].Nobs;
#else
			cum += (v.NzeqXmaj * mylog((1 - v.lam) + beta3*(v.lam * 4 / 3 - 1)) + v.NzeqXmin*mylog(v.lam / 3 - (4 * v.lam / 3 - 1)*beta3));
#endif
		}
	}

	Value = -cum;
	return(true);
}


bool optBeta::calcDeriv(const void *data, const vecD &parameters, vecD &Derivs) {
	// d val / d b; where val = sum_i log( A * b + B ) -> sum_i ([ 1/(A*b+B)]*(A))
	double beta3 = parameters[0], cum = 0;
	ulong nblocks = (ulong)auxen.size(); // 64 -> 32 bit
	double tmp = 0;
	for (ulong iblock = 0; iblock < nblocks; iblock++) {
		auto &ab = auxen[iblock];
		for (ulong ipoly = 0; ipoly < ab.size(); ipoly++) {
			const beta_aux_t &v = ab[ipoly];
#ifndef V0
			cum += (ab[ipoly].A / (beta3 * ab[ipoly].A + ab[ipoly].B))*ab[ipoly].Nobs;
#else
			tmp = 4 * v.lam / 3 - 1;
			cum += (v.NzeqXmaj * tmp / ((1 - v.lam) + beta3*tmp) + v.NzeqXmin * tmp / (v.lam / 3 - beta3*tmp));
#endif
		}
	}

	Derivs[0] = -cum;
	return(true);
};

bool optBeta::calcValueAndDeriv(const void *data, const vecD &parameters, double &Value, vecD &Derivs) {
	// calc val
	double beta3 = parameters[0], cumv = 0.0, cumd = 0.0;
	ulong nblocks = (ulong)auxen.size();	// 64 -> 32 bit
	for (ulong iblock = 0; iblock < nblocks; iblock++) {
		auto &ab = auxen[iblock];
		for (ulong ipoly = 0; ipoly < ab.size(); ipoly++) {
			// can speed this up with prob hack, multiple untill resutls would denormalize 
			// then log add.... mult is much faster... and may get a hundred mult per log...
#ifndef V0
			double t;
			t = beta3 * ab[ipoly].A + ab[ipoly].B;
			cumv += mylog(t) * ab[ipoly].Nobs;
			cumd += (ab[ipoly].A / (t))*ab[ipoly].Nobs;
#else
			const beta_aux_t &v = ab[ipoly];
			double tmp = 4 * v.lam / 3 - 1;
			cumv += (v.NzeqXmaj * mylog((1 - v.lam) + beta3*tmp) + v.NzeqXmin*mylog(v.lam / 3 - tmp*beta3));
			cumd += (v.NzeqXmaj * tmp / ((1 - v.lam) + beta3*tmp) + v.NzeqXmin * tmp / (v.lam / 3 - beta3*tmp));
#endif
		}
	}

	Value = -cumv;
	Derivs[0] = -cumd;
	return true;
};
#endif