#pragma once
#include <cassert>
#include <vector>
#include <algorithm>	//std::sort
#include <iostream>		// Used for slow / debugging  IO

// Derive a class from optimizer for inferring Betas...
#include "butils/butils.h"		// includes optimizer...
#include "inscomp/inscomp.h"
#include "BlockSet.h"

#include <iostream>	// debugging
#include <cmath>	// sqrt...

class optInsight : public butils::optimizer {
public:
	typedef butils::ezmatrix ezmatrix;
	typedef butils::optimizer optimizer;
	typedef optimizer::vecD vecD;
	typedef optimizer::vecL vecL;
	typedef inscomp::BlockSet_t BlockSet_t;
	typedef inscomp::modelInsight modelInsight;
	typedef modelInsight::parameters parameters;
	struct insightProvonance : optimizer::ResultProvonance {
		double	insNll;
		vecD	insDeriv;
		insightProvonance() { clear(); };
		void clear() { ResultProvonance::clear(); insNll = 0; insDeriv.clear(); }
	};
	
	struct  options : optimizer::options {	// Add INSIGHT specific options to generic optimizer options 
		bool		fixRho, fixEta , fixGam ;	// do not refine these values during optimization.
		double		minEta;						// experimental!
		int			v;							// verbosity, 0 = silent, higher for runtime debugging and monitoring.
		options() { clear(); }
		void clear() { fixRho = fixEta = fixGam = false; minEta = 0.0; optimizer::options::clear(); v = 0;  return; };
	};

private:
	// precomputed cofactors for each term. These could likely be converted to floats
	//typedef float rstore;		// real storage type

	// the block information gets incorproated into the aux variables as sufficient statistics that are observation count weighted!
	// once we have it we no longer need the original block structure and can free that memory.
	options					opt;
	parameters				mFix;
	void clear(bool ClearAccel = true) { opt.clear(); mFix.clear(); if (ClearAccel) mAccel.clear(); };

public:
	inscomp::modelInsightAccelerator	mAccel;

	optInsight( const modelInsight &M ) : mAccel(M) {	clear( false); }
	virtual ~optInsight() { clear(); };

	// Calculate the LL for a given parameter set, this is an optimized routine for high performance...
	virtual bool calcValue(const void *data, const vecD &params, double &Value) {
		if (!mAccel.initialized) return false;
		parameters p=mAccel.getparameters(); p.rho=params[0]; p.eta = params[1]; p.gam = params[2];
		if (opt.fixRho) p.rho = mFix.rho;
		if (opt.fixEta) p.eta = mFix.eta;
		if (opt.fixGam) p.gam = mFix.gam;
		mAccel.setparameters(p); Value = mAccel.calcNLL();
		return true;
	};

	// override first difference method with an analytic deravitive to improve accuracy and speed calculation.
	virtual bool calcDeriv(const void *data, const vecD &params, vecD &Deravitives) {
		double tmp_val = 0;
		bool ok = calcValueAndDeriv(data, params, tmp_val, Deravitives);	// this is no slower than just calculating derivs...
		return ok;
	};

	// Sometimes it is most efficient to calculate deriv & value at same time. If so, override, if not, delegate...
	virtual bool calcValueAndDeriv(const void *data, const vecD &params, double &Value, vecD &Deravitives) {
		if (!mAccel.initialized) return false;
		parameters p = mAccel.getparameters(); p.rho = params[0]; p.eta = params[1]; p.gam = params[2];
		if (opt.fixRho) p.rho = mFix.rho;
		if (opt.fixEta) p.eta = mFix.eta;
		if (opt.fixGam) p.gam = mFix.gam;
		mAccel.setparameters(p); Value = mAccel.calcNLL( &p ); // p becomes deravitives...
		if (opt.fixRho) p.rho = 0;
		if (opt.fixEta) p.eta = 0;
		if (opt.fixGam) p.gam = 0;
		Deravitives[0] = p.rho; Deravitives[1] = p.eta; Deravitives[2] = p.gam;
		return true;
	}

	// central Maximum likelihood finder, find the params that minimize NLL
	bool optimize(const BlockSet_t &bs, modelInsight &model, const options &opts, double PriorWeight = 0, const modelInsight *modelPrior = NULL, void *provonance = NULL) {
		typedef butils::optimizer::ResultProvonance ResultProvonance;
		typedef butils::lbfgs::vecD	vecD;

		bool ok = false;
        //# very large MAX values destabalize likelihood esitmates, but soimetimes low rho & high weak/pos selectio ndrives large eta / gamma for a given dw , PA
        //So allow adaptive relaxation of  these limits
        bool   maxFactor_ok=false;
        double maxFactor_eta = 1.0, maxFactor_gam = 1.0, maxFactor_max=600.0;
		while (! maxFactor_ok ) {
  			

			ok = false;
			clear(); opt = opts; opt.verbose = opt.v; mFix = model.getparameters();
			mAccel.setparameters( model.getparameters() );
			mAccel.initialize( bs, PriorWeight, modelPrior );

			vecL limits(3); 
			// double eps = 1e-6; // epsilon for when we need to make sure parameter values are not on boundaries...
			// hardcoded , from manual curationexpereince....
    		double rhoMin = .001,			rhoMax	= .999;
    		double etaMin = opts.minEta,	etaMax	= 2.0 * maxFactor_eta; // 20; Very large values can drive infinite deravitive... and destablize likelihood estimates...
    		double gamMin = 0,				gamMax	= 5.0 * maxFactor_gam; // 20;
    		limits[0].vmin = &rhoMin; limits[0].vmax = &rhoMax;
    		limits[1].vmin = &etaMin; limits[1].vmax = &etaMax;
    		limits[2].vmin = &gamMin; limits[2].vmax = &gamMax;
    
    		vecD paramIn(3), paramOut(3);
    		int iter=0;	// start with tight tolerance, if it fails, relax tolerance
    		while (!ok && iter <= 4 ) {
    			switch (iter) {
    				case 0: break;											// use default, this is high precision 1e2x machine precision, no slope limit this si about 13 sig digits
    				case 1: opt.factr = 1.0e4 ;                    break;	// 1e4 * machine precision, no slope limit, this is about 11 sig digits still excellent
    				case 2: opt.factr = 1.0e6 ; opt.pgtol = 1e-12; break;	// 1e6 * machine presion, weaker slope limit, 9 sig digits, still very good
    				case 3: opt.factr = 1.0e8 ; opt.pgtol = 1e-10; break;	// 1e8 * this gives us only ~7 significant digits... perfectly servicable
    				case 4: opt.factr = 1.0e10; opt.pgtol = 1e-08; break;	// 1e10 * this gives us only ~5 significant digits... adequate but minimally acceptable. anything worse is a real problem
    			}
    			paramIn[0] = mFix.rho; paramIn[1] = mFix.eta; paramIn[2] = mFix.gam; paramOut = paramIn;
    			ResultProvonance prov_tmp, *prov = (provonance == NULL ? &prov_tmp : (ResultProvonance *)provonance);
    
    			// opt.factr = 100.0; opt.pgtol = 1e-9;								// factr, tolerance relative to machine precision, pgtol, terminate on low gradient.
    			ok = optimizer::optimize( (const void *) &bs, paramIn, limits, opt, paramOut, prov);
    			prov->pgtol = opt.pgtol; prov->factr = opt.factr;
    			iter++;
    		}
    		mFix.rho = paramOut[0]; mFix.eta = paramOut[1]; mFix.gam = paramOut[2]; model.setparameters(mFix);
    		if (opts.v > 0 || opts.verbose > 0 || iter > 2 ) {
    			std::cerr << "optInsight:optimize() found solution on iteration " << iter << " factr = " << opt.factr << "   pgtol = " << opt.pgtol << std::endl;
    		}
			

			maxFactor_ok = true;
    		// lift limits if solution is at bound.
            if ( paramOut[1]>(.9*etaMax) ) { maxFactor_eta *= 3.0; maxFactor_ok = false; };
            if ( paramOut[2]>(.9*gamMax) ) { maxFactor_gam *= 3.0; maxFactor_ok = false; };
			// prevent infinite loostening
			if ( (maxFactor_eta>maxFactor_max) || (maxFactor_gam>maxFactor_max) ) {
				// this is an ERROR condition!
    			std::cerr << "ERROR: optInsight:optimize() solution found at parameter boundry, but adaptive boundary at max. Solution may not be optimal. eta " << paramOut[1] << "  eta bound " << etaMax << "  gamma " << paramOut[2] << "  gamma bound " << gamMax << "  ." << std::endl;
				maxFactor_ok = true;
			}
    		if ( (opts.v > 0 || opts.verbose > 0) && !maxFactor_ok ) {
    			std::cerr << "optInsight:optimize() solution found at parameter boundry, loosening boundary. eta " << paramOut[1] << "  eta bound " << etaMax << "  gamma " << paramOut[2] << "  gamma bound " << gamMax << "  ." << std::endl;
    		}
		}

		//clear();									// release AUX values? ok.. for now, but they are expensive.. needed for err & Supp stat calculation, keep them around...
		return ok;
	};

};

#ifdef OLD_CODE
struct posteriordist {	// Admissible posteriors over hidden variables (S,A) as S and allowable relations between (Z, X & A).
	double cSelNodiv;	// 0) S == 1, A == Z | Y=M,  (inf sites -> Z==Xma==A)
	double cSelDiv;		// 1) S == 1, A != Z | Y=M,  (inf sites -> Z!=Xma==A)
	double cSelLow;		// 2) S == 1, A == Xma | Y=L, Z==Xma==A
	double cNeutNodiv;  // 3) S == 0, A == Z | Y=M,  (inf sites -> Z==Xma==A) - See Supp table S1.
	double cNeutDiv;	// 4) S == 0, A != Z | Y=M,  (inf sites -> Z!=Xma==A)
	double cNeutLow;	// 5) S == 0, A == Xma | Y=L, (Z==Xma==A)+(Z!=Xma==A) Low  Derived allele Freq 
	double cNeutHigh;	// 6) S == 0, A == Xmi | Y=L, (Z==Xmi==A)+(Z!=Xmi==A) High Derived allele Freq 				
	double cHigh;		// 7) S == 0, A == Xma or Xmi  | Y=H, any Z (Y=H -> S=0), by definition Neut...
	double cCnt;		// number of counts, shoud be equal to sum of all counts...
	double probSel() { double c = count(); return((c <= 0 ? 0 : (cSelNodiv + cSelDiv + cSelLow) / c)); };	// S==1
	double probWeak() { double c = count(); return((c <= 0 ? 0 : cSelLow / c)); };							// Y==L, S==1
	double probAdap() { double c = count(); return((c <= 0 ? 0 : cSelDiv / c)); };							// Y==M, Z!=A,S==1
	double count() { return(cSelNodiv + cSelDiv + cSelLow + cNeutNodiv + cNeutDiv + cNeutLow + cNeutHigh + cHigh); };
	void clear() { cSelNodiv = cSelDiv = cSelLow = cNeutNodiv = cNeutDiv = cNeutLow = cNeutHigh = cHigh = cCnt = 0.0; return; }
	void renorm() { double c = count(); cSelNodiv /= c; cSelDiv /= c; cSelLow /= c; cNeutNodiv /= c; cNeutDiv /= c; cNeutLow /= c; cNeutHigh /= c; cHigh = 0; cHigh = 1.0 - count(); cCnt = 1.0; };
	posteriordist() { clear(); }
	std::string header() { return("SelNoD\tSelDiv\tSelLow\tNeutNoD\tNeutDiv\tNeutLow\tNeutHi\tHi"); };
	std::string toStr(const char *f = "%7.5lf") { return(std::string("") + fmt(cSelNodiv, f) + "\t" + fmt(cSelDiv, f) + "\t" + fmt(cSelLow, f) + "\t" + fmt(cNeutNodiv, f) + "\t" + fmt(cNeutDiv, f) + "\t" + fmt(cNeutLow, f) + "\t" + fmt(cNeutHigh, f) + "\t" + fmt(cHigh, f)); };
	inline char *fmt(double a, const char *f = "%7.5lf") { static char buf[32];  sprintf(buf, (a < 1e-4 ? "%7le" : f), a); return(buf); };
};

struct parameters {
	typedef std::string string;
	struct derived_t { double alpha, tau, dp, pw; };	// values derived from rho, eta,gamma over all positions in all blocks
	struct block_t { double lambdaT, theta; };			// values typically derivbed from a single polymorphism block
	struct betas_t { double b1, b2, b3; };				// values estimated from neutral poly statistics over collection of blocks
	double rho, eta, gamma;								// essential model parameters, derived from all positions in all blocks
	block_t block;										// values typically derivbed from a single polymorphism block
	betas_t beta;										// values estimated from neutral poly statistics over all blocks
	derived_t sup;										// values derived directly from other parameters
	parameters() { clear(); };
	void clear() { rho = eta = gamma = -1; beta.b1 = beta.b1 = beta.b3 = -1; sup.alpha = sup.tau = sup.dp = sup.pw = -1; block.lambdaT = block.theta = -1.0; };
	string fmt(double v, const char *f = "%12.6lf") { static char b[32]; sprintf(b, (v >= 1e-4 ? f : "%.6le"), v); return(b); };
	string toStr() {
		string s = "Rho:\t" + fmt(rho) + "\nEta:\t" + fmt(eta) + "\nGamma:\t" + fmt(gamma) + "\n";
		s += string("Lambda:\t") + fmt(block.lambdaT) + "\n"; s += "Theta:\t" + fmt(block.theta) + "\n";
		s += "Beta1:\t" + fmt(beta.b1) + "\nBeta2:\t" + fmt(beta.b2) + "\nBeta3:\t" + fmt(beta.b3) + "\n";
		s += "Alpha:\t" + fmt(sup.alpha) + "\nTau:\t" + fmt(sup.tau) + "\tDp:\t" + fmt(sup.dp) + "\nPw:\t" + fmt(sup.pw) + "\n";
		return s;
	};
};

struct postElem {
	double rho, eta, gam;			// parameter values
	double drho, deta, dgam;		// dimensions of cell (differential)
	double lik;						// data negative log likelihood, base e
	double priorDens;				// the point samples prior, actualy a prior DENSITY, NLL base e
	double posterior, posteriorDens;	// data negative log likelihood, base e
	};
struct postSet {
	std::vector<postElem> elements;
	vecD gridRho, gridEta, gridGam;		// sampling grids, finer near ML/MAP values
	double cumPost;						// cumulative unnormalized posterior, thisis the normalization factor
	double expRho, expEta, expGam;		// expected values
	double expLam, expThet;				// expected valeus for Theta and LambdaT
	double beta1, beta2, beta3;			// point estiamtes for beta1,2&3
	ezmatrix priorObs;					// population distribution of input states, so P(param) \porp [P(X|param)]^N, for N pseudo counts.

};

enum class SelectionClass { first = 0, Mono = first, PolyL, PolyH, last = PolyH, size = (last + 1) };	// Polymorphism class, Monomorphic, low frequence Polymorphisn, or High Frequency Polymotphism
enum class AncestralAllele { first = 0, Xmaj = first, Xmin, Xoth, last = Xoth, size = (last + 1) };	// Ancestral State (Z) matches population Major Allele, Minor Allele or neither (oth).

static double wattersonsA(uint32_t N) { double a = 0.0; for (int i = 2; i <= (int)N; i++) a += (1.0 / (double)(i - 1)); return(a); }
{ from ioptimize

bool optimize(const BlockSet_t &data, const vecD &paramIn, const vecL &limits, const options &opts, vecD &paramOut, insightProvonance *provonance) {
	this->opts = &opts;							// save a pointer to the parameter optimization options, this is referenced in getVals / getDerivs callbacks...
	bool ok = initialize(data, paramIn, opts);	// Calculate AUX values, store options.
	{
		parameters pi = opts.modelCurrent;
		for (double r = .001; r < 0.5; r += .01) {
			pi.rho = r;
			pi.eta = 0.001; pi.gamma = .54;
			const inscomp::BlockSet_t &bs = *((const inscomp::BlockSet_t *)data);
			printf("%.4lf\t%.4lf\t%.4lf\t%.6lf\n", pi.rho, pi.eta, pi.gamma, dataNLL(pi, bs));
		}
	};

	ok = true;
}

double dataNLL(const parameters &Param, const inscomp::BlockSet_t &bs);

// calculate supplimentry stats, including Dp, Pw alpha & Tau
bool suppStats(const void *data, const vecD &params, const vecD &pErrs, const options &opts, vecD &deriv, vecD &dErrs);

// standard error value for parameters using LL curvature method..
bool estimateErr(const void *data, const vecD &paramIn, const options &opts, vecD &errs, double delta = 1.0e-6);

// Calculate prior probability of Model as expected likelihood of Nobs pseudocounts drawn from Population (under Model)
double Prior(const parameters &Pop, const parameters &Proposed, ulong Nallele, double Nobs = 1.0, ezmatrix *PriorCache = NULL);
bool calcDerivRaw(const parameters &Model, ulong Nallele, SelectionClass Icls, double ZeqXmaj, double ZeqXmin, vecD &Derivs);


// Calculate Likelihood, that is Class probability based in parameters, P(Y,X|Model).
//	Lambda is genrally multiplied by T before storage in the database... Nallele is used to genreate theta*A
//	X is provided as an enunerated constant representing one sufficient condition: Xmaj=Z, Xmin=Z or Z!=Xmaj,Xmin
// ByClass =false generates valeus for a single allele observation (1 position). ByClass=true generates weighted 
//	values by input class: Y, Xmaj=Z, Xmin=Z or Z!=Xmaj,Xmin.
static double likelihood(const parameters &Param, ulong Nallele, SelectionClass Y, AncestralAllele Z, bool ByClass = false, double *CondS0 = NULL, double *CondS1 = NULL);
static void likelihoodMatrix(const parameters &Param, ulong Nallele, bool ByClass, ezmatrix &Joint);
static void likelihoodMatrixLog(const parameters &Param, ulong Nallele, bool ByClass, ezmatrix &Joint, bool Proper = false);
static double likelihoodExpected(const parameters &Param, ulong Nallele, SelectionClass Y, double pZeqXmaj, double pZeqXmin, bool ByClass = true) {
	return(pZeqXmaj * likelihood(Param, Nallele, Y, AncestralAllele::Xmaj, ByClass) + pZeqXmin * likelihood(Param, Nallele, Y, AncestralAllele::Xmin, ByClass) + (1.0 - pZeqXmaj - pZeqXmaj) * likelihood(Param, Nallele, Y, AncestralAllele::Xoth, ByClass));
}

// Expected probabilities of hidden variables S & A (A is provided as a sufficent set of relationships between A, X and Z. Expectation over Z.)
static bool posteriorExpected(const parameters &Param, const void *data, posteriordist &Post, ulong V);							// Cumulative over whole data set, fast.
static bool posteriorIndividual(const parameters &Param, ulong Nallele, SelectionClass Y, double ZeqXmaj, double ZeqXmin, posteriordist &Post);	// One site, slow

																																				// Public interface for high performance likelihood calculations
bool posteriorInit(const void *data, const options &Opt) { vecD unused_params;	return(initialize(data, unused_params, Opt)); };
bool parameterPosterior(const void *data, const parameters &model, const options &opts, postSet &Posterior);


// Prior derived from intuition provided by pseudocounts, given rigor by Dirchlet priors. Start with an initial uniform prior on model parameters M.
// The likelihood of an observation d is thus P(d|M)=P(M|d)P(d)/P(M). As P(d) is a constant, the uniform prior provides that P(M|d) \porp P(d|M).
// For an expected d of D, I use a M generated from the whole genome (population) and generate an expected observation (D) then set the prior
// to P(M') = P(D|M'). To apply weight consitent with N pseudocounts, simple take P(D|M')^N. This prior, related to the KL distance between M and M'
// biases models towars those that represent the global population, inducing statistical shrinkage.
double optInsight::Prior(const parameters &Pop, const parameters &Proposed, butils::ulong Nallele, double Nobs, ezmatrix *PriorCache) {
	ezmatrix *pc = PriorCache; double prob = 0.0;
	if (Nobs <= 0.0) return(1.0);
	if (PriorCache == NULL) {	// gerate distribution of input states, according to global (population) model. This is the expected distribution over states.
		pc = new ezmatrix(); likelihoodMatrix(Pop, Nallele, true, *pc);
	}
	// Calculate expected probbaility of global observation, under proposed model. This estimates KL distance...
	for (int i = (int)SelectionClass::first; i <= (int)SelectionClass::last; i++) {
		for (int j = (int)AncestralAllele::first; j <= (int)AncestralAllele::last; j++) {
			prob += pc->v(i, j) * likelihood(Proposed, Nallele, (SelectionClass)i, (AncestralAllele)j);
		};
	};
	// exponentiate to generate relative probabiltiy of N observations (pseudocounts) drawn from Population.
	if (Nobs != 1.0) prob = pow(prob, Nobs);
	if (PriorCache == NULL) { delete pc; pc = NULL; }
	return prob;
}

double optInsight::likelihood(const parameters &Param, ulong NAllele, SelectionClass Y, AncestralAllele Z, bool ByClass, double *CondS0, double *CondS1) {
	// P(d|M) directly From Table3, doi:10.1093/molbev/mst019 . To get values conditional on S=0, set Rho=0.0. To get values conditional on S=1, set Rho=1.0. 
	// CAREFULL HERE: a data observation d, consists if an allele{a,c,g,t} Observation for Z,Xmaj,Xmin and classY.
	//   expected observations become distributions over the 16 allele combinations an d 3 classes of Y... that sum to 1.
	//	The sufficient properties to calculate the likelihood of that data are Y,Z=Xmaj,Z=Xmin... while a specific d, determines the 
	//	values for these sufficient properties (SP(d)\in SP={Y,ZeqX), so we can calcute P(d|model) as P(SP(d)|Model), this probability is NOT
	//	P(sp \in SP| Model), because diffeing d map to the same SP. Compare to 6 sided die model M with 1-4 colored blue and 5-6 colored red.
	//    and for a roll outrcome d P(d)=1/8 iff blue, 1/4 if red P(d=3) can be represented as p(d|M) = p(d is a specific blue number |M) = 1/8
	//    but p(any blue number) = 1/2, not 1/8. We are dealing with exactlu the same thing. ZeqX(d) is a sufficient property to calculate P(d) but
	//		this probability is NOT P(ZeqX).
	// The provide probabilties of a particular observation of observed variables (considered class Y, and alleles X_maj, X_min). 
	// While this probbaility distribution is convienetly parameterized in terms of Y, X_maj==Z and X_min==Z, this is NOT the probability
	//	of Y, X_min==Z or X_maj==Z... for that we need to account for multiple alleles (priors). Careful! See likelihoodByClass().
	// Also, I know that Z is not actually observed, but we estiamte it in a seperate step and treat it as observed for the purposes of INSIGHT.
	double LambdaT = Param.block.lambdaT, ThetaA = Param.block.theta * wattersonsA(NAllele);
	double Rho = Param.rho, Eta = Param.eta, Gamma = Param.gamma;
	double Beta1 = Param.beta.b1, Beta2 = Param.beta.b2, Beta3 = Param.beta.b3;
	double ps0 = 0.0, ps1 = 0.0, p = 0.0, w = 1.0;		// Combinations have probability 0 unless identified below.....
	switch (Y) {
	case SelectionClass::Mono:
		switch (Z) {
		case AncestralAllele::Xmaj:	ps0 = (1.0 - LambdaT)*(1 - ThetaA);		ps1 = (1.0 - Eta*LambdaT)*(1.0 - Gamma*ThetaA);		w = 1.0;	break;
		case AncestralAllele::Xmin:	ps0 = 0.0;								ps1 = 0.0;											w = 0.0;	break;	// M means there is NO minor X allele.
		case AncestralAllele::Xoth:	ps0 = (1 - ThetaA)*LambdaT / 3.0;		ps1 = Eta*LambdaT / 3.0;							w = 3.0;	break;
		}; break;
	case SelectionClass::PolyL:
		switch (Z) {
		case AncestralAllele::Xmaj:	ps0 = ((1.0 - LambdaT)*Beta1 + LambdaT*Beta3 / 3.0)*ThetaA / 3.0;	ps1 = (1.0 - Eta*LambdaT)*Gamma*ThetaA / 3.0;	w = 3.0; break;
		case AncestralAllele::Xmin:	ps0 = ((1.0 - LambdaT)*Beta3 + LambdaT*Beta1 / 3.0)*ThetaA / 3.0;	ps1 = 0.0;										w = 3.0; break;
		case AncestralAllele::Xoth:	ps0 = (LambdaT / 3.0)*(Beta1 + Beta3)*ThetaA / 3.0;					ps1 = 0.0;										w = 6.0; break;
		}; break;
	case SelectionClass::PolyH:
		ps1 = 0.0;
		switch (Z) {	// logic for weights is nuanced here... see Normalization table.
		case AncestralAllele::Xmaj:	ps0 = (1.0 - 2.0*LambdaT / 3.0)*Beta2*ThetaA / 3.0;	ps1 = 0.0;	w = 1.5; break;
		case AncestralAllele::Xmin:	ps0 = (1.0 - 2.0*LambdaT / 3.0)*Beta2*ThetaA / 3.0;	ps1 = 0.0;	w = 1.5; break;
		case AncestralAllele::Xoth:	ps0 = (2.0*LambdaT*Beta2 / 3.0)*ThetaA / 3.0;		ps1 = 0.0;	w = 3.0; break;
		}; break;
	};
	if (ByClass) { ps0 *= w; ps1 *= w; };	// Upweight by class, to get normalized probability, by class (Selec / Ancestral Allele)
	p = (Rho * ps1) + (1.0 - Rho)*ps0;
	if (CondS0 != NULL) *CondS0 = ps0;
	if (CondS1 != NULL) *CondS1 = ps1;
	return p;
};

double optInsight::dataNLL(const parameters &Param, const inscomp::BlockSet_t &bs) {
	double zxmaj, zxmin, ll = 0.0;
	ulong nallele = bs.numAlleles;
	ezmatrix like;							// data log likelihood for current block
	parameters tmp_param = Param;			// parameters for current block
	inscomp::AncPriIdx_t idxMa = 0, idxMi = 0;	// index into probility table
	inscomp::PatternCounts_t nobs;			// count of observations of this pattern

	ll = 0.0; tmp_param = Param;
	for (ulong iblock = 0; iblock < bs.numBlocks(); iblock++) {
		const inscomp::Block_t &bl = bs.blocks[iblock];
		bl.getHeader(NULL, NULL, NULL, &tmp_param.block.lambdaT, &tmp_param.block.theta);
		if (true) {
			double l;
			// This is log( expected [likelihood])
			likelihoodMatrix(tmp_param, nallele, false, like); //like.logElem();
			for (ulong i = 0; i < bl.monoNumPat(); i++) {
				bl.monoVal(i, idxMa, nobs); zxmaj = *(bs.ancPri[idxMa]);
				l = like.v((int)SelectionClass::Mono, (int)AncestralAllele::Xmaj) * (zxmaj)+like.v((int)SelectionClass::Mono, (int)AncestralAllele::Xoth) * ((1.0 - zxmaj));
				if (l>0) ll += nobs * mylog(l);
			};
			for (ulong i = 0; i < bl.polyLNumPat(); i++) {
				bl.polyLVal(i, idxMa, idxMi, nobs); zxmaj = *(bs.ancPri[idxMa]); zxmin = *(bs.ancPri[idxMi]);
				l = like.v((int)SelectionClass::PolyL, (int)AncestralAllele::Xmaj) * (zxmaj)+
					like.v((int)SelectionClass::PolyL, (int)AncestralAllele::Xmin) * (zxmin)+
					like.v((int)SelectionClass::PolyL, (int)AncestralAllele::Xoth) * ((1.0 - zxmaj - zxmin));
				if (l>0) ll += nobs * mylog(l);
			};
			for (ulong i = 0; i < bl.polyHNumPat(); i++) {
				bl.polyHVal(i, idxMa, idxMi, nobs); zxmaj = *(bs.ancPri[idxMa]); zxmin = *(bs.ancPri[idxMi]);
				l = like.v((int)SelectionClass::PolyH, (int)AncestralAllele::Xmaj) * (zxmaj)+
					like.v((int)SelectionClass::PolyH, (int)AncestralAllele::Xmin) * (zxmin)+
					like.v((int)SelectionClass::PolyH, (int)AncestralAllele::Xoth) * ((1.0 - zxmaj - zxmin));
				if (l>0) ll += nobs * mylog(l);
			};
		} else {

			// this is expected [ log(likelihood) ]
			likelihoodMatrix(tmp_param, nallele, false, like); like.logElem();
			//likelihoodMatrixLog(tmp_param, nallele, false, like,false); // get actual data log-likelihoods for this block...
			for (ulong i = 0; i < bl.monoNumPat(); i++) {
				bl.monoVal(i, idxMa, nobs); zxmaj = *(bs.ancPri[idxMa]);
				ll += like.v((int)SelectionClass::Mono, (int)AncestralAllele::Xmaj) * (nobs * zxmaj);			assert(ll < -1e-100);
				ll += like.v((int)SelectionClass::Mono, (int)AncestralAllele::Xoth) * (nobs * (1.0 - zxmaj));	assert(ll < -1e-100);
			};
			for (ulong i = 0; i < bl.polyLNumPat(); i++) {
				bl.polyLVal(i, idxMa, idxMi, nobs); zxmaj = *(bs.ancPri[idxMa]); zxmin = *(bs.ancPri[idxMi]);
				ll += like.v((int)SelectionClass::PolyL, (int)AncestralAllele::Xmaj) * (nobs * zxmaj);			assert(ll < -1e-100);
				ll += like.v((int)SelectionClass::PolyL, (int)AncestralAllele::Xmin) * (nobs * zxmin);			assert(ll < -1e-100);
				ll += like.v((int)SelectionClass::PolyL, (int)AncestralAllele::Xoth) * (nobs * (1.0 - zxmaj - zxmin));	assert(ll < -1e-100);
			};
			for (ulong i = 0; i < bl.polyHNumPat(); i++) {
				bl.polyHVal(i, idxMa, idxMi, nobs); zxmaj = *(bs.ancPri[idxMa]); zxmin = *(bs.ancPri[idxMi]);	assert(ll < -1e-100);
				ll += like.v((int)SelectionClass::PolyH, (int)AncestralAllele::Xmaj) * (nobs * zxmaj);			assert(ll < -1e-100);
				ll += like.v((int)SelectionClass::PolyH, (int)AncestralAllele::Xmin) * (nobs * zxmin);			assert(ll < -1e-100);
				ll += like.v((int)SelectionClass::PolyH, (int)AncestralAllele::Xoth) * (nobs * (1.0 - zxmaj - zxmin));	assert(ll < -1e-100);
			};
		}
	};
	if (false) {
		ezmatrix t_lik, t_prob;
		likelihoodMatrix(opts->modelPrior, nallele, true, t_prob);	// calculate data distribution
		parameters t_param = opts->modelPrior; t_param.rho = Param.rho;
		likelihoodMatrix(t_param, nallele, false, t_lik); t_lik.logElem();	// get actual data log-likelihoods using expected block values.

		if (opts && opts->priorWeight > 0) {
			// pseudocount bias the parameters toward those of the generating distribution, in this case population-wide values.
			tmp_param = Param;	// use expected model values to calculate data likelihood for pseudocounts
								// likelihoodMatrix(tmp_param, nallele, false, like); like.logElem();	// get actual data log-likelihoods using expected block values.
			double tll = 0, pseudocounts = 0.0;
			for (int sc = (int)SelectionClass::first; sc <= (int)SelectionClass::last; sc++) {
				for (int zxc = (int)AncestralAllele::first; zxc <= (int)AncestralAllele::last; zxc++) {
					tll = t_lik.v((int)sc, (int)zxc); pseudocounts = opts->priorWeight * t_prob.v(sc, zxc);
					// tll = like.v((int)sc, (int)zxc); pseudocounts = opts->priorWeight * popExpInp.v(sc, zxc);
					if (pseudocounts > 0) ll += tll * pseudocounts;
				};
			};
		};
	};
	return (-ll);
};

// Calculate log expected likelihood (improper), or expected log likelihood (proper)
void optInsight::likelihoodMatrixLog(const parameters &Param, ulong Nallele, bool ByClass, ezmatrix &Joint, bool Proper) {
	//	bool ret;
	if (!Proper) {
		likelihoodMatrix(Param, Nallele, ByClass, Joint); Joint.logElem();
	} else {
		parameters p = Param;
		double rho = Param.rho;
		ezmatrix j0, j1;
		p.rho = 0;  likelihoodMatrix(p, Nallele, ByClass, j0); j0.logElem();
		p.rho = 1;  likelihoodMatrix(p, Nallele, ByClass, j1); j1.logElem();
		j0.mult(1 - rho); j1.mult(rho);
		Joint = j0; Joint.add(j1);
	}
	return;
}

void optInsight::likelihoodMatrix(const parameters &Param, ulong Nallele, bool ByClass, ezmatrix &Joint) {
	Joint.resize((int)SelectionClass::size, (int)AncestralAllele::size);
	for (int i = (int)SelectionClass::first; i <= (int)SelectionClass::last; i++) {
		for (int j = (int)AncestralAllele::first; j <= (int)AncestralAllele::last; j++) {
			Joint.v(i, j) = likelihood(Param, Nallele, (SelectionClass)i, (AncestralAllele)j, ByClass);
		}
	}
}

bool optInsight::posteriorIndividual(const parameters &Param, ulong Nallele, SelectionClass Y, double ZeqXmaj, double ZeqXmin, posteriordist &Post) {
	// The provide probabilties of UNOBSERVED (latent) variables (considered S & A). Yes, I know that the allele distribution over A is not explicitely 
	//	treated, and this is a little confusing (legacy). For efficency, we calculate the distribution over allowable S,A combinations for input
	//  observations. Input observation consists of the sufficient characterization of Selection Class Y {M,L,H} and the prob that Z=Xmaj & Z=Xmin.
	Post.clear();
	// informative alias for parameters
	double rho = Param.rho, eta = Param.eta, gamma = Param.gamma;
	double lambda = Param.block.lambdaT, theta = Param.block.theta * wattersonsA(Nallele);
	double beta1 = Param.beta.b1, beta3 = Param.beta.b3;
	double pzxmaj = ZeqXmaj, pzxmin = ZeqXmin;
	switch (Y) {
	case SelectionClass::Mono:
		Post.cNeutNodiv = (1 - rho)*(1 - theta)*pzxmaj*(1 - lambda);			// 3) S == 0, A == Z | Y=M,  (inf sites -> Z==Xma==A) - See Supp table S1.
		Post.cNeutDiv = (1 - rho)*(1 - theta)*(1 - pzxmaj)*lambda / 3;		// 4) S == 0, A != Z | Y=M,  (inf sites -> Z!=Xma==A)
		Post.cSelNodiv = rho*(1 - gamma*theta)*pzxmaj*(1 - eta*lambda);		// 0) S == 1, A == Z | Y=M,  (inf sites -> Z==Xma==A)
		Post.cSelDiv = rho*(1 - pzxmaj)*eta*lambda / 3;						// 1) S == 1, A != Z | Y=M,  (inf sites -> Z!=Xma==A)
		break;
	case SelectionClass::PolyL:
		Post.cSelLow = rho * gamma*theta / 3 * pzxmaj * (1 - eta*lambda);								// 2) S == 1, A == Xma | Y=L, Z==Xma==A
		Post.cNeutLow = (1 - rho) * beta1 * theta / 3 * ((1 - pzxmaj)*lambda / 3 + pzxmaj*(1 - lambda));	// 5) S == 0, A == Xma | Y=L, (Z==Xma==A)+(Z!=Xma==A) Low  Derived allele Freq 
		Post.cNeutHigh = (1 - rho) * beta3 * theta / 3 * ((1 - pzxmin)*lambda / 3 + pzxmin*(1 - lambda));	// 6) S == 0, A == Xmi | Y=L, (Z==Xmi==A)+(Z!=Xmi==A) High Derived allele Freq 				ptot      = pSelLow + pNeutLow + pNeutHigh;
		break;
	case SelectionClass::PolyH:
		Post.cHigh = 1.0;	break;										// 7) S=0, A=Xma or A=Xmi | Y = H (all HF pollys are Neut, by def)
	default: return(false);
	};
	Post.cCnt = 1.0;
	return true;
}
// Assess the posterior probabilities of the hidden variables (conditioned on a model and observed data).
// Formally, the hidden variables are A, Z & S, but we treat the deep ancestor (Z) as observed, and take the expectation over it, or 
// more precisely, over a distribution over admissible combinations of Xmaj, Xmin and Z. This leaves us with S and A. Rather than explicitely
// assessing the posterior of A over {a,c,g,t}, we look at the implication of the distibution over A for 
// over allowable combinations of Z, A & S. This gives us a compact representation over the 8 admissable states, all other combinations
// have 0 probability under the INSIGHT model (or due to the infinite sites assumption applied for the population derived from A).
// 
bool optInsight::posteriorExpected(const parameters &paramIn, const void *data, posteriordist &pos, ulong V) {
	using butils::ulong;
	using inscomp::BlockSet_t;
	const BlockSet_t &bs = *(BlockSet_t *)data;		// type casting
	ulong nblocks = (ulong)bs.blocks.size();		// 64-> 32 bit
	double pNeutNodiv, pNeutDiv, pSelNodiv, pSelDiv, pSelLow, pNeutLow, pNeutHigh;
	double cNeutNodiv = 0, cNeutDiv = 0, cSelNodiv = 0, cSelDiv = 0, cSelLow = 0, cNeutLow = 0, cNeutHigh = 0, cHigh = 0, cCnt = 0;
	double lambda = 0, theta = 0, a = 0, ptot;

	// informative alias for parameters
	double rho = paramIn.rho, eta = paramIn.eta, gamma = paramIn.gamma, beta1 = paramIn.beta.b1, beta2 = paramIn.beta.b2, beta3 = paramIn.beta.b3;

	inscomp::PatternCounts_t	nobs = 0;
	inscomp::AncPriIdx_t		i_ancpriMaj = 0, i_ancpriMin = 0;
	a = optInsight::wattersonsA(bs.numAlleles);
	for (ulong iblock = 0; iblock < nblocks; iblock++) {
		auto &bl = bs.blocks[iblock];					// get current block from database
		bl.getHeader(NULL, NULL, NULL, &lambda, &theta); theta *= a;
		for (ulong isite = 0; isite < bl.monoNumPat(); isite++) {
			bl.monoVal(isite, i_ancpriMaj, nobs);
			if (i_ancpriMaj > 0) {
				auto pzxmaj = *(bs.ancPri[i_ancpriMaj]);							// floating point type
				pNeutNodiv = (1 - rho)*(1 - theta)*pzxmaj*(1 - lambda);			// 3) S == 0, A == Z | Y=M,  (inf sites -> Z==Xma==A) - See Supp table S1.
				pNeutDiv = (1 - rho)*(1 - theta)*(1 - pzxmaj)*lambda / 3;		// 4) S == 0, A != Z | Y=M,  (inf sites -> Z!=Xma==A)
				pSelNodiv = rho*(1 - gamma*theta)*pzxmaj*(1 - eta*lambda);		// 0) S == 1, A == Z | Y=M,  (inf sites -> Z==Xma==A)
				pSelDiv = rho*(1 - pzxmaj)*eta*lambda / 3;						// 1) S == 1, A != Z | Y=M,  (inf sites -> Z!=Xma==A)
				ptot = pNeutNodiv + pNeutDiv + pSelNodiv + pSelDiv;
				if (ptot > 0) {
					ptot = (1.0 / ptot) * nobs;
					cNeutNodiv += pNeutNodiv * ptot;
					cNeutDiv += pNeutDiv	  * ptot;
					cSelNodiv += pSelNodiv  * ptot;
					cSelDiv += pSelDiv	  * ptot;
					cCnt += nobs;
				};
			};
		}
		for (ulong isite = 0; isite < bl.polyLNumPat(); isite++) {
			bl.polyLVal(isite, i_ancpriMaj, i_ancpriMin, nobs);
			if (i_ancpriMaj > 0 && i_ancpriMin > 0) {
				auto pzxmaj = *(bs.ancPri[i_ancpriMaj]);		// floating point type
				auto pzxmin = *(bs.ancPri[i_ancpriMin]);		// floating point type
				pSelLow = rho * gamma*theta / 3 * pzxmaj * (1 - eta*lambda);										// 2) S == 1, A == Xma | Y=L, Z==Xma==A
				pNeutLow = (1 - rho) * beta1 * theta / 3 * ((1 - pzxmaj)*lambda / 3 + pzxmaj*(1 - lambda));		// 5) S == 0, A == Xma | Y=L, (Z==Xma==A)+(Z!=Xma==A) Low  Derived allele Freq 
				pNeutHigh = (1 - rho) * beta3 * theta / 3 * ((1 - pzxmin)*lambda / 3 + pzxmin*(1 - lambda));		// 6) S == 0, A == Xmi | Y=L, (Z==Xmi==A)+(Z!=Xmi==A) High Derived allele Freq 				ptot      = pSelLow + pNeutLow + pNeutHigh;
				if (ptot > 0) {
					ptot = (1.0 / ptot) * nobs;
					cSelLow += pSelLow   * ptot;
					cNeutLow += pNeutLow  * ptot;
					cNeutHigh += pNeutHigh * ptot;
					cCnt += nobs;
				};
			};
		};

		for (ulong isite = 0; isite < bl.polyHNumPat(); isite++) {
			bl.polyHVal(isite, i_ancpriMaj, i_ancpriMin, nobs);
			if (i_ancpriMaj > 0 && i_ancpriMin > 0) { cHigh += nobs; cCnt += nobs; }								// 7) S=0, A == Xma or Xmi  | Y=H, any Z (Y=H -> S=0)
		};
	};  // loop or blocks...

		// Recapitulate Illans stats, but add an item for total counts...
	pos.cSelNodiv = cSelNodiv;	// 0 
	pos.cSelDiv = cSelDiv;		// 1
	pos.cSelLow = cSelLow;		// 2
	pos.cNeutNodiv = cNeutNodiv;	// 3
	pos.cNeutDiv = cNeutDiv;		// 4 
	pos.cNeutLow = cNeutLow;		// 5
	pos.cNeutHigh = cNeutHigh;	// 6
	pos.cHigh = cHigh;		// 7
	pos.cCnt = cCnt;			// 8
	return true;
};

// Get paramater error estimates using the hessian...
bool optInsight::estimateErr(const void *data, const vecD &paramIn, const options &opts, vecD &errs, double delta) {
	ezmatrix a, ainv;
	bool ok = hessian(data, paramIn, opts, a, delta);
	if (!ok) return(ok);
	a.inv(ainv); if (errs.size() < paramIn.size()) errs.resize(paramIn.size());
	if (false) { // debugging
		ezmatrix m = a;
		std::cerr << "Inverse Hessian is  " << std::endl;
		std::cerr << ainv.toStr() << std::endl;
		a.mult(ainv, m);
		std::cerr << "Is Normal?: " << std::endl;
		std::cerr << m.toStr() << std::endl;
	};
	for (int i = 0; i < paramIn.size(); i++) errs[i] = sqrt(ainv.v(i, i));
	return ok;
};

// calculate the second deravitive matrix using finite differences of analytic first deravitive...
bool optInsight::hessian(const void *data, const vecD &paramIn, const options &opts, ezmatrix &Hess, double delta) {
	bool ok = true;
	vecD d1 = paramIn, d2 = paramIn, ptmp;
	if (!initialized) ok = initialize(data, paramIn, opts);
	if (!ok) return(ok);
	Hess.resize((int)paramIn.size());				// square matrix
	for (int i = 0; i < Hess.nrow(); i++) {		// use finite difference of first deravitives to get second deravitives
		ptmp = paramIn; ptmp[i] -= delta; calcDeriv(data, ptmp, d1);
		ptmp = paramIn; ptmp[i] += delta; calcDeriv(data, ptmp, d2);
		for (int j = 0; j < Hess.ncol(); j++) Hess.v(i, j) = (d2[j] - d1[j]) / (2 * delta);
	}
	if (false) { // debugging
		std::cerr << "Estimating Hessian: " << std::endl;
		std::cerr << Hess.toStr() << std::endl;
	};
	return ok;
}


#ifdef V0
bool optInsight::initialize(const void *data, const vecD &paramIn, const options &Opts) {
	using butils::ulong;
	using inscomp::BlockSet_t;

	const BlockSet_t			&bs = *(BlockSet_t *)data;			// type casting
	ulong						nblocks = (ulong)bs.blocks.size();	// 64-> 32 bit, many positions, but fewer blocks.
	BlockSet_t::AncPriInd_t		xz_maji = 0, xz_mini = 0;			// indices
	BlockSet_t::AncPriVal_t		xz_maj = 0.0, xz_min = 0.0;		// values
	BlockSet_t::AncPriCount_t	nobs = 0;						// counts = # of M/L/H observations in data set.
	parameters					tmp_param = Opts.fixVals;			// currently developed model parameters...
	double						wa = wattersonsA(bs.numAlleles);	// Assumed constant (full data)

	aux_h_t	tmp_h; aux_l_t	tmp_l; aux_m_t	tmp_m;					// space optimized data structures for holding sufficient statistics and weights.

	clear(); aux_m.reserve(bs.numMonoPat());	aux_l.reserve(bs.numPolyLPat());	aux_h.reserve(bs.numPolyHPat()); // clear and reallocate memory...

																													 // loop over blocks, for each block handle mono, PolyL and polyH in subloops.
	for (ulong iblock = 0; iblock < nblocks; iblock++) {
		auto &bl = bs.blocks[iblock];									// get current block from database
		bl.getHeader(NULL, NULL, NULL, &tmp_param.block.lambdaT, &tmp_param.block.theta); // get lambda from block... can provide first args if debugging is desired...
		tmp_param.block.theta *= wa;	// This is the only time we inclide wa in a param::theta value.
		for (ulong i = 0; i < bl.polyHNumPat(); i++) {
			// get index into prob ability table for z=Maj/Min human allele
			bl.polyHVal(i, xz_maji, xz_mini, nobs);
			// lookup indices and convert nobs into double
			xz_maj = *(bs.ancPri[xz_maji]); xz_min = *(bs.ancPri[xz_mini]);
			// calculate summary statistics, yes this can be optimized, but it is not really a bottleneck. Long form flows directly from paper.
			tmp_h.load(tmp_param, xz_maj, xz_min, nobs);
			aux_h.push_back(tmp_h);
		};
		for (ulong i = 0; i < bl.polyLNumPat(); i++) {
			// get index into prob ability table for z=Maj/Min human allele
			bl.polyLVal(i, xz_maji, xz_mini, nobs);
			// lookup indices and convert nobs into double
			xz_maj = *(bs.ancPri[xz_maji]); xz_min = *(bs.ancPri[xz_mini]);
			// calculate summary statistics, yes this can be optimized, but it is not really a bottleneck. Long form flows directly from paper.
			tmp_l.load(tmp_param, xz_maj, xz_min, nobs);
			aux_l.push_back(tmp_l);
		};
		for (ulong i = 0; i < bl.monoNumPat(); i++) {
			// get index into prob ability table for z=Maj human allele
			bl.monoVal(i, xz_maji, nobs);
			// lookup indices and convert nobs into double
			xz_maj = *(bs.ancPri[xz_maji]);
			// calculate cofactors of parameter rho,eta,gamma, rho*eta, rho*gamma, rho*eta*gamma
			tmp_m.load(tmp_param, xz_maj, nobs);
			aux_m.push_back(tmp_m);
		};
	};

	if (Opts.priorWeight > 0) {
		// use priors to generate posterior (normalized) data likelihood for each observation class 
		tmp_param.block.theta = Opts.fixVals.block.theta * wa; tmp_param.block.lambdaT = Opts.fixVals.block.lambdaT;
		butils::ezmatrix counts;

		parameters tp2 = Opts.priorVal;
		likelihoodMatrix(tp2, bs.numAlleles, true, counts);	// get relative counts

		SelectionClass sc;
		double xma = 0, xmi = 0, xot = 0, x = 0;	// P(Y) = x = xma+xmi+xot. Probs of Z==Xma, Z==Xmi, Z==Xother
		sc = SelectionClass::Mono;  xma = counts.v((int)sc, (int)AncestralAllele::Xmaj); xmi = counts.v((int)sc, (int)AncestralAllele::Xmin); xot = counts.v((int)sc, (int)AncestralAllele::Xoth); x = xma + xmi + xot;
		if (x>0) {
			tmp_m.load(tmp_param, 1.0, xma*Opts.priorWeight); aux_m.push_back(tmp_m);
			tmp_m.load(tmp_param, 0.0, xmi*Opts.priorWeight); aux_m.push_back(tmp_m);
			//tmp_m.load(tmp_param, xma / x, x*Opts.priorWeight); aux_m.push_back(tmp_m);
		}

		sc = SelectionClass::PolyL; xma = counts.v((int)sc, (int)AncestralAllele::Xmaj); xmi = counts.v((int)sc, (int)AncestralAllele::Xmin); xot = counts.v((int)sc, (int)AncestralAllele::Xoth); x = xma + xmi + xot;
		if (x>0) {
			tmp_l.load(tmp_param, 1.0, 0.0, xma * Opts.priorWeight); aux_l.push_back(tmp_l);
			tmp_l.load(tmp_param, 0.0, 1.0, xmi * Opts.priorWeight); aux_l.push_back(tmp_l);
			tmp_l.load(tmp_param, 0.0, 0.0, xot * Opts.priorWeight); aux_l.push_back(tmp_l);
			//tmp_l.load(tmp_param, xma / x, xmi / x, x * Opts.priorWeight); aux_l.push_back(tmp_l);
		}

		sc = SelectionClass::PolyH; xma = counts.v((int)sc, (int)AncestralAllele::Xmaj); xmi = counts.v((int)sc, (int)AncestralAllele::Xmin); xot = counts.v((int)sc, (int)AncestralAllele::Xoth); x = xma + xmi + xot;
		if (x>0) {
			tmp_h.load(tmp_param, 1.0, 0.0, xma * Opts.priorWeight); aux_h.push_back(tmp_h);
			tmp_h.load(tmp_param, 0.0, 1.0, xmi * Opts.priorWeight); aux_h.push_back(tmp_h);
			tmp_h.load(tmp_param, 0.0, 0.0, xot * Opts.priorWeight); aux_h.push_back(tmp_h);
			//tmp_h.load(tmp_param, xma / x, xmi / x, x * Opts.priorWeight); aux_h.push_back(tmp_h); 
		}
	}

	opts = &Opts;  initialized = true;
	return(true);
}
#else
bool optInsight::initialize(const void *data, const vecD &paramIn, const options &Opts) {
	using butils::ulong;
	using inscomp::BlockSet_t;
	const BlockSet_t			&bs = *(BlockSet_t *)data;			// recover type, we may no longer need this....

	if (Opts.priorWeight > 0) {
		// use priors to generate posterior (normalized) data likelihood for each observation class. Weight*popExpInp = Pseudocounts.
		likelihoodMatrix(Opts.modelPrior, bs.numAlleles, true, popExpInp);	// get normalized frequency of counts in population.
	}

	opts = &Opts;  initialized = true;
	return(true);
}
#endif

#ifdef V0
// The sufficient statistics transform the likelihood estimates into a simple form that utilises
//	A + B * Rho + C * Rho*Eta + D * Rho*Gamma + E *  Rho * Eta * Gamma. For each of H, L, M dat
//  may be implicitely 0.
bool optInsight::calcValue(const void *data, const vecD &parameters, double &Value) {
	// calc val
	double rho = parameters[0], eta = parameters[1], gamma = parameters[2];
	if (opts != NULL) {
		if (opts->fixRho) rho = opts->fixVals.rho;
		if (opts->fixEta) eta = opts->fixVals.eta;
		if (opts->fixGam) gamma = opts->fixVals.gamma;
	}
	double rhoeta = rho*eta, rhogamma = rho*gamma, rhoetagamma = rho*eta*gamma;
	// double c_con = 0, c_rho = 0, c_rhoeta = 0, c_rhogamma = 0, c_rhoetagamma = 0, c_tot = 0;	// keep counts for debugging
	double prob = 1.0, cum = 0;
	ulong ns = 0;
	// Process Mono
	ns = (ulong)aux_m.size();	// 64-> 32 bit...
	for (ulong isite = 0; isite < ns; isite++) {
		aux_m_t &site = aux_m[isite];
		prob = site.fCon + site.fRho * rho + site.fRhoEta*rhoeta + site.fRhoGamma * rhogamma + site.fRhoEtaGamma * rhoetagamma;
		//c_con			+= site.Nobs * site.fCon;
		if (prob > 0) {
			cum += mylog(prob)*site.Nobs; // We can likely speed this up by multiplying probs ratehr than adding logs...
		};
	}
	// Process PolyH
	ns = (ulong)aux_h.size();	// 64-> 32 bit...
	for (ulong isite = 0; isite < ns; isite++) {
		aux_h_t &site = aux_h[isite];
		prob = site.fCon + site.fRho * rho;
		// c_rho += site.Nobs * site.fRho;
		if (prob > 0) {
			cum += mylog(prob)*site.Nobs; // We can likely speed this up by multiplying probs ratehr than adding logs...
		};
	}
	// Process PolyL
	ns = (ulong)aux_l.size();	// 64-> 32 bit...
	for (ulong isite = 0; isite < ns; isite++) {
		aux_l_t &site = aux_l[isite];
		prob = site.fCon + site.fRho * rho + site.fRhoGamma * rhogamma + site.fRhoEtaGamma * rhoetagamma;
		if (prob > 0) {
			cum += mylog(prob)*site.Nobs; // We can likely speed this up by multiplying probs ratehr than adding logs...
		};
	}

	Value = -cum;
	return(true);
}


bool optInsight::calcDeriv(const void *data, const vecD &parameters, vecD &Derivs) {
	// calc val
	double rho = parameters[0], eta = parameters[1], gamma = parameters[2];
	if (opts != NULL) {
		if (opts->fixRho) rho = opts->fixVals.rho;
		if (opts->fixEta) eta = opts->fixVals.eta;
		if (opts->fixGam) gamma = opts->fixVals.gamma;
	}
	double rhoeta = rho*eta, rhogamma = rho*gamma, rhoetagamma = rho*eta*gamma, etagamma = eta * gamma;
	double prob = 1.0, pinv = 1, cum = 0;
	Derivs.resize(parameters.size());
	for (ulong i = 0; i < Derivs.size(); i++) Derivs[i] = 0;

	// Process Mono
	ulong ns = (ulong)aux_m.size();	// 64-> 32 bit...
	for (ulong isite = 0; isite < ns; isite++) {
		aux_m_t &site = aux_m[isite];
		// prob = site.fCon + site.fRho * rho + site.fRhoEta*rhoeta + site.fRhoGamma * rhogamma + site.fRhoEtaGamma*rhoetagamma;	// Implicit cofactor of 1 for rhoetagamma
		prob = site.fCon + site.fRho * rho + site.fRhoEta*rhoeta + site.fRhoGamma * rhogamma + site.fRhoEtaGamma * rhoetagamma;	// Implicit cofactor of 1 for rhoetagamma
		if (prob > 0) {
			pinv = 1.0 / prob;
			Derivs[0] += site.Nobs * (-pinv * (site.fRho + site.fRhoEta * eta + site.fRhoGamma * gamma + site.fRhoEtaGamma * etagamma));	// d(log(prob))/drho
			Derivs[1] += site.Nobs * (-pinv * (site.fRhoEta * rho + site.fRhoEtaGamma * rhogamma));	// d(log(prob))/deta
			Derivs[2] += site.Nobs * (-pinv * (site.fRhoGamma * rho + site.fRhoEtaGamma * rhoeta));		// d(log(prob))/dgamma
																										//cum += site.Nobs * log(prob); // We can likely speed this up by multiplying probs ratehr than adding logs...
		};
	}
	// Process PolyH
	ns = (ulong)aux_h.size();	// 64-> 32 bit...
	for (ulong isite = 0; isite < ns; isite++) {
		aux_h_t &site = aux_h[isite];
		prob = site.fCon + site.fRho * rho;
		if (prob > 0) {
			pinv = 1.0 / prob;
			Derivs[0] += site.Nobs * (-pinv * (site.fRho));	// d(log(prob))/drho
															//Derivs[1] += site.Nobs * (pinv + (rho * site.fRhoEta + rhogamma));	// d(log(prob))/deta
															//Derivs[2] += site.Nobs * (pinv + (rho * site.fRhoGamma + rhoeta));		// d(log(prob))/dgamma
															//cum += site.Nobs * log(prob); // We can likely speed this up by multiplying probs ratehr than adding logs...
		};
	}
	// Process PolyL
	ns = (ulong)aux_l.size();	// 64-> 32 bit...
	for (ulong isite = 0; isite < ns; isite++) {
		aux_l_t &site = aux_l[isite];
		prob = site.fCon + site.fRho * rho + site.fRhoGamma * rhogamma + site.fRhoEtaGamma * rhoetagamma;
		if (prob > 0) {
			pinv = 1.0 / prob;
			Derivs[0] += site.Nobs * (-pinv * (site.fRho + site.fRhoGamma * gamma + site.fRhoEtaGamma * etagamma));	// d(log(prob))/drho
			Derivs[1] += site.Nobs * (-pinv * (site.fRhoEtaGamma * rhogamma));	// d(log(prob))/deta
			Derivs[2] += site.Nobs * (-pinv * (site.fRhoGamma * rho + site.fRhoEtaGamma * rhoeta));	// d(log(prob))/dgamma
																									//cum += site.Nobs * log(prob); // We can likely speed this up by multiplying probs ratehr than adding logs...
		};
	}

	if (opts != NULL) {
		if (opts->fixRho) Derivs[0] = 0.0;
		if (opts->fixEta) Derivs[1] = 0.0;
		if (opts->fixGam) Derivs[2] = 0.0;
	}
	return(true);
};

bool optInsight::calcValueAndDeriv(const void *data, const vecD &parameters, double &Value, vecD &Derivs) {
	// calc val
	double rho = parameters[0], eta = parameters[1], gamma = parameters[2];
	if (opts != NULL) {
		if (opts->fixRho) rho = opts->fixVals.rho;
		if (opts->fixEta) eta = opts->fixVals.eta;
		if (opts->fixGam) gamma = opts->fixVals.gamma;
	}
	double rhoeta = rho*eta, rhogamma = rho*gamma, rhoetagamma = rho*eta*gamma, etagamma = eta * gamma;
	double prob = 1.0, pinv = 1, cum = 0;
	ulong ns = 0;
	Derivs.resize(parameters.size());
	for (ulong i = 0; i < Derivs.size(); i++) Derivs[i] = 0;

	// Process Mono
	ns = (ulong)aux_m.size();	// 64-> 32 bit...
	for (ulong isite = 0; isite < ns; isite++) {
		aux_m_t &site = aux_m[isite];
		// prob = site.fCon + site.fRho * rho + site.fRhoEta*rhoeta + site.fRhoGamma * rhogamma + site.fRhoEtaGamma*rhoetagamma;	// Implicit cofactor of 1 for rhoetagamma
		prob = site.fCon + site.fRho * rho + site.fRhoEta*rhoeta + site.fRhoGamma * rhogamma + site.fRhoEtaGamma * rhoetagamma;	// Implicit cofactor of 1 for rhoetagamma
		if (prob > 0) {
			pinv = 1.0 / prob;
			Derivs[0] += site.Nobs * (-pinv * (site.fRho + site.fRhoEta * eta + site.fRhoGamma * gamma + site.fRhoEtaGamma * etagamma));	// d(log(prob))/drho
			Derivs[1] += site.Nobs * (-pinv * (site.fRhoEta * rho + site.fRhoEtaGamma * rhogamma));	// d(log(prob))/deta
			Derivs[2] += site.Nobs * (-pinv * (site.fRhoGamma * rho + site.fRhoEtaGamma * rhoeta));	// d(log(prob))/dgamma
			cum += site.Nobs * mylog(prob); // We can likely speed this up by multiplying probs ratehr than adding logs...
		}
	}
	// Process PolyH
	ns = (ulong)aux_h.size();	// 64-> 32 bit...
	for (ulong isite = 0; isite < ns; isite++) {
		aux_h_t &site = aux_h[isite];
		prob = site.fCon + site.fRho * rho;
		if (prob > 0) {
			pinv = 1.0 / prob;
			Derivs[0] += site.Nobs * (-pinv * (site.fRho));	// d(log(prob))/drho
															//Derivs[1] += site.Nobs * (pinv + (rho * site.fRhoEta + rhogamma));	// d(log(prob))/deta
															//Derivs[2] += site.Nobs * (pinv + (rho * site.fRhoGamma + rhoeta));		// d(log(prob))/dgamma
			cum += site.Nobs * mylog(prob); // We can likely speed this up by multiplying probs ratehr than adding logs...
		};
	}
	// Process PolyL
	ns = (ulong)aux_l.size();	// 64-> 32 bit...
	for (ulong isite = 0; isite < ns; isite++) {
		aux_l_t &site = aux_l[isite];
		prob = site.fCon + site.fRho * rho + site.fRhoGamma * rhogamma + site.fRhoEtaGamma * rhoetagamma;
		if (prob > 0) {
			pinv = 1.0 / prob;
			Derivs[0] += site.Nobs * (-pinv * (site.fRho + site.fRhoGamma * gamma + site.fRhoEtaGamma * etagamma));	// d(log(prob))/drho
			Derivs[1] += site.Nobs * (-pinv * (site.fRhoEtaGamma * rhogamma));	// d(log(prob))/deta
			Derivs[2] += site.Nobs * (-pinv * (site.fRhoGamma * rho + site.fRhoEtaGamma * rhoeta));	// d(log(prob))/dgamma
			cum += site.Nobs * mylog(prob); // We can likely speed this up by multiplying probs ratehr than adding logs...
		};
	}

	Value = -cum;

	if (opts != NULL) {
		if (opts->fixRho) Derivs[0] = 0.0;
		if (opts->fixEta) Derivs[1] = 0.0;
		if (opts->fixGam) Derivs[2] = 0.0;
	}
	return true;
};
#else
bool optInsight::calcValue(const void *data, const vecD &parameters, double &Value) {
	const inscomp::BlockSet_t &bs = *((const inscomp::BlockSet_t *) data);
	optInsight::parameters tmp_param = opts->modelCurrent;	// fill in betas etc..
	tmp_param.rho = parameters[0];	tmp_param.eta = parameters[1];	tmp_param.gamma = parameters[2];
	if (opts->fixRho) tmp_param.rho = opts->modelFix.rho;
	if (opts->fixEta) tmp_param.eta = opts->modelFix.eta;
	if (opts->fixGam) tmp_param.gamma = opts->modelFix.gamma;
	Value = dataNLL(tmp_param, bs);
	return true;
}
#endif

// To get raw derives, call this 3 times for each SelectionClass, with (ZeqXmaj,ZeqXmin) = {1,0), {0,1}, {0,0}...
bool optInsight::calcDerivRaw(const parameters &Model, ulong Nallele, SelectionClass Icls, double ZeqXmaj, double ZeqXmin, vecD &Derivs) {
	// Code is drawn from original EM code, so it is a little messy... apologies... (see Initialzie() and CalcDeriv()-- Brad 
	double beta1 = Model.beta.b1, beta2 = Model.beta.b2, beta3 = Model.beta.b3;
	double lambda = Model.block.lambdaT, theta = Model.block.theta * wattersonsA(Nallele);
	double xz_maj = ZeqXmaj, xz_min = ZeqXmin;
	double rho = Model.rho, eta = Model.eta, gamma = Model.gamma;
	double rhoeta = rho*eta, rhogamma = rho*gamma, rhoetagamma = rho*eta*gamma, etagamma = eta * gamma;
	double p = 0, pinv = 0;
	Derivs.resize(3); for (auto i = 0; i < Derivs.size(); i++) Derivs[i] = 0.0;
	switch (Icls) {
	case SelectionClass::Mono: { aux_m_t c;
		c.fCon = xz_maj         * (1.0 - lambda)*(1 - theta);
		c.fCon += (1.0 - xz_maj) * (1.0 - theta)*lambda / 3.0;
		c.fRho = (xz_maj)-c.fCon;
		c.fRhoEta = (1.0 - xz_maj) * lambda / 3.0;
		c.fRhoEta += (xz_maj)*-lambda;
		c.fRhoGamma = xz_maj*-theta;
		c.fRhoEtaGamma = xz_maj*lambda*theta;
		p = c.fCon + c.fRho * rho + c.fRhoEta*rhoeta + c.fRhoGamma * rhogamma + c.fRhoEtaGamma * rhoetagamma;	// Implicit cofactor of 1 for rhoetagamma
		if (p <= 0) return false;
		pinv = 1.0 / p;
		Derivs[0] += -pinv * (c.fRho + c.fRhoEta * eta + c.fRhoGamma * gamma + c.fRhoEtaGamma * etagamma);	// d(log(prob))/drho
		Derivs[1] += -pinv * (c.fRhoEta * rho + c.fRhoEtaGamma * rhogamma);	// d(log(prob))/deta
		Derivs[2] += -pinv * (c.fRhoGamma * rho + c.fRhoEtaGamma * rhoeta);	// d(log(prob))/dgamma
	}; break;
	case SelectionClass::PolyL: { aux_l_t c;
		c.fCon = xz_maj                 *((1.0 - lambda)*beta1 + lambda*beta3 / 3.0)*theta / 3.0;
		c.fCon += xz_min                 *((1.0 - lambda)*beta3 + lambda*beta1 / 3.0)*theta / 3.0;
		c.fCon += (1.0 - (xz_min + xz_maj))*((lambda / 3.0)*(beta1 + beta3)*(theta / 3.0));
		c.fRho = -c.fCon;
		c.fRhoGamma = xz_maj*theta / 3.0;
		c.fRhoEtaGamma = -xz_maj*lambda*theta / 3.0;
		p = c.fCon + c.fRho * rho + c.fRhoGamma * rhogamma + c.fRhoEtaGamma * rhoetagamma;
		if (p <= 0) return false;
		pinv = 1.0 / p;
		Derivs[0] += -pinv * (c.fRho + c.fRhoGamma * gamma + c.fRhoEtaGamma * etagamma);	// d(log(prob))/drho
		Derivs[1] += -pinv * (c.fRhoEtaGamma * rhogamma);	// d(log(prob))/deta
		Derivs[2] += -pinv * (c.fRhoGamma * rho + c.fRhoEtaGamma * rhoeta);	// d(log(prob))/dgamma
	}; break;
	case SelectionClass::PolyH: { aux_h_t c;
		c.fCon = (xz_maj + xz_min)         * (1.0 - 2.0*lambda / 3.0)*beta2*theta / 3.0;
		c.fCon += (1.0 - (xz_maj + xz_min)) * (2.0*lambda / 3.0)*beta2*theta / 3.0;
		c.fRho = -c.fCon;
		p = c.fCon + c.fRho * rho;
		if (p <= 0) return false;
		pinv = 1.0 / p;
		Derivs[0] += -pinv * (c.fRho);	// d(log(prob))/drho
	}; break;
	};
	return true;
}

// Calculate supplimentary stats
bool optInsight::suppStats(const void *data, const vecD &params, const vecD &pErrs, const options &opts, vecD &deriv, vecD &dErrs) {
	using butils::ulong;
	using inscomp::BlockSet_t;
	const BlockSet_t &bs = *(BlockSet_t *)data;	// type casting

												// get data-weighted estimates of lambda / theta across blocks.
	double lambda = 0, theta = 0, lambdatheta = 0, nsites = 0;
	double cum_count = 0, cum_lam = 0, cum_thet = 0, cum_lamthet = 0;
	for (ulong iblock = 0; iblock < bs.blocks.size(); iblock++) {
		const inscomp::Block_t &b = bs.blocks[iblock];
		b.getHeader(NULL, NULL, NULL, &lambda, &theta);  nsites = b.NumSites();
		cum_lam += lambda * nsites; cum_thet += theta * nsites; cum_lamthet += theta * lambda * nsites; cum_count += nsites;
	}
	{
		double a = wattersonsA(bs.numAlleles);
		lambda = cum_lam / cum_count; theta = a* cum_thet / cum_count; lambdatheta = a* cum_lamthet / cum_count;
	};

	// relable useful values....
	double rho = params[0], eta = params[1], gamma = params[2];
	double a, b, c; // temp values...

					// Calculate primary uncertainty in parameters (variance is ainv)
	ezmatrix ah, ainv;
	bool ok = hessian(data, params, opts, ah, 1e-6);	if (!ok) return(ok);
	ah.inv(ainv);
	if (false)  std::cerr << "\nCovariance Matrix:\n" << ainv.toStr() << std::endl;

	// insure elements of ainv are positive and symmetric...
	for (int i = 0; i < ainv.nrow(); i++) {
		for (int j = i; j < ainv.ncol(); j++) {
			double v = (ainv.v(i, j) + ainv.v(j, i)) / 2.0;
			ainv.v(i, j) = ainv.v(j, i) = v;
		}
	}

	// Dp (sites under positive selection)
	double eDp, Dp = rho * eta * lambda;
	if (rho == 0 || eta == 0) { Dp = 0; eDp = -1.0; } else {
		eDp = ainv.v(0, 0) / (rho*rho) + ainv.v(1, 1) / (eta*eta) + 2 * ainv.v(0, 1) / (eta*rho);
		eDp = rho*eta*lambda * sqrt(eDp);
	}
	// Pw (sites udner weak selection)
	double ePw, Pw = rho * gamma * (theta - eta * lambdatheta);
	if (rho == 0 || gamma == 0 || (theta - eta * lambdatheta) == 0) { Pw = 0; ePw = -1.0; } else {
		a = 1.0 / rho;  c = 1.0 / gamma;  b = -lambdatheta / (theta - eta * lambdatheta);
		ePw = a*a*ainv.v(0, 0) + b*b*ainv.v(1, 1) + c * c * ainv.v(2, 2);
		ePw += 2.0 * a * b *ainv.v(0, 1) + 2.0*a*ainv.v(0, 2) + 2.0 * b * c * ainv.v(1, 2);
		ePw = Pw * sqrt(ePw);
	}
	// Alpha
	double ealpha, alpha;
	if (rho*eta == 0) { alpha = 0.0; ealpha = -1.0; } else {
		alpha = 1.0 / (1.0 + (1.0 - rho) / (rho*eta));
		a = 1.0 / (rho * rho * eta); b = (1.0 - rho) / (rho * eta * eta);
		ealpha = a*a*ainv.v(0, 0) + b*b*ainv.v(1, 1) + 2.0 * a*b*ainv.v(0, 1);
		ealpha = alpha * alpha * sqrt(ealpha);
	}
	// Tau
	double etau, tau;
	if (rho*eta == 0) { tau = 0.0; etau = -1.0; } else {
		tau = 1.0 / (1.0 + (1.0 - rho) / (rho*gamma));
		a = 1.0 / (rho * rho * gamma); b = (1.0 - rho) / (rho * gamma * gamma);
		etau = a * a * ainv.v(0, 0) + b * b * ainv.v(2, 2) + 2.0 * a * b * ainv.v(0, 2);
		etau = tau * tau * sqrt(etau);
	}
	// return values
	deriv.resize(4); dErrs.resize(4);
	deriv[0] = Dp * 1000;  deriv[1] = Pw * 1000;  deriv[2] = alpha;  deriv[3] = tau;
	dErrs[0] = eDp * 1000; dErrs[1] = ePw * 1000; dErrs[2] = ealpha; dErrs[3] = etau;
	return true;
};

#endif
