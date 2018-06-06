#pragma once

#include "modelInsight.h"

// Likelihood calculation for a fixed data set, over varying models, is frequently used both in
//	MAP / ML calculation of optimal parameter values, but also in posterior probability calculations.
//	This code greatly accelerates the calculation of likelihoods by computing sufficient statistics for
//	each block which form a function that takes rho, eta and gamma to a likelihood value.
//
// The likelihood relationship is nonlinear and requires rho*eta, rho*gamma and rho*eta*gamma as well as rho, eta and gamma.
//	Analytic calculation of the gradient of log likelihood curve also requires eta*gamma, rounding out all combinations
//	of one, two or three unuqie factors.
//

struct modelInsightAccelerator : modelInsight {
	typedef double rstore;		// try float if these structures become too big....

	struct aux_factors { 
		rstore rho, eta, gamma, etagamma, rhoeta, rhogamma, rhoetagamma; 
		aux_factors() { clear(); }; 
		void load( double R, double E, double G ) {	rho=R; eta=E; gamma=G; etagamma=E*G; rhoeta=R*E; rhogamma=R*G; rhoetagamma=R*E*G; }
		void clear () { load( 0,0,0); }; 
	};

	struct aux_m_t {			// cofcators from Mono Sites
		rstore fCon;			// constant
		rstore fRho;			// cofactor of rho
		rstore fRhoEta;			// cofactor of rho*eta
		rstore fRhoGamma;		// cofactor of rho*gamma
		rstore fRhoEtaGamma;	// cofactor of rho*eta*gamma
		rstore Nobs;			// Number of times this was oberved in current block. might be ok as be uint16
		inline void load(const parameters &model, double ZeqXmaj, double nobs) {
			double lambda = model.block.lambdaT, theta = model.block.theta;
			fCon = ZeqXmaj         * (1.0 - lambda)*(1 - theta);
			fCon += (1.0 - ZeqXmaj) * (1.0 - theta)*lambda / 3.0;
			fRho = (ZeqXmaj)-fCon;
			fRhoEta = (1.0 - ZeqXmaj) * lambda / 3.0;
			fRhoEta += (ZeqXmaj)*-lambda;
			fRhoGamma = ZeqXmaj*-theta;
			fRhoEtaGamma = ZeqXmaj*lambda*theta;
			Nobs = nobs;
		};
		inline double calcL1( const aux_factors &x) const { double d=fCon + fRho * x.rho + fRhoEta * x.rhoeta + fRhoGamma * x.rhogamma + fRhoEtaGamma * x.rhoetagamma; return(d<0||d>1?DBL_MIN:d); };
		inline double calcL(  const aux_factors &x ) const { return( pow( calcL1(x),Nobs) ); };
		inline double calcLL( const aux_factors &x ) const { return( Nobs * mylog( calcL1(x) )); };
		inline double calcLLD(const aux_factors &x, parameters &d ) const { 
			double p = calcL1(x); d.clear();
			if (p<=0) {return( 0 ); };
			double pi = 1.0 / p;
			d.rho = Nobs * (pi * (fRho + fRhoEta * x.eta + fRhoGamma * x.gamma + fRhoEtaGamma * x.etagamma));	// d(log(prob))/drho
			d.eta = Nobs * (pi * (fRhoEta * x.rho + fRhoEtaGamma * x.rhogamma));								// d(log(prob))/deta
			d.gam = Nobs * (pi * (fRhoGamma * x.rho + fRhoEtaGamma * x.rhoeta));								// d(log(prob))/dgamma
			return(Nobs * mylog(p));
		};
	};

	struct aux_h_t {			// cofcators from HF polys
		rstore fCon;			// constant
		rstore fRho;			// cofactor of rho
		rstore Nobs;			// Number of times this was oberved in current block.
		inline void load(const parameters &model, double ZeqXmaj, double ZeqXmin, double nobs) {
			double beta2 = model.beta.b2, lambda = model.block.lambdaT, theta = model.block.theta;
			fCon = (ZeqXmaj + ZeqXmin)          * (1.0 - 2.0*lambda / 3.0)*beta2*theta / 3.0;
			fCon += (1.0 - (ZeqXmaj + ZeqXmin)) * (2.0*lambda / 3.0)*beta2*theta / 3.0;
			fRho = -fCon;
			Nobs = nobs;
		};
		inline double calcL1(const aux_factors &x) const { double d= fCon + fRho * x.rho; return (d<0||d>1?DBL_MIN:d); };		// likelihood of 1 observation
		inline double calcL(const aux_factors &x) const { return(pow(calcL1(x),Nobs)); };		// likelihood of N independent observations
		inline double calcLL(const aux_factors &x) const { return(Nobs * mylog(calcL1(x))); };	// Log likelihood of N independent observations
		inline double calcLLD(const aux_factors &x, parameters &d) const {						// Vanue and Deravitive of Log likelihood of N independent observations
			double p=calcL1(x);	d.clear();
			if (p <= 0) return( 0.0 );
			double pinv = 1.0 / p;
			d.rho = Nobs * (pinv * (fRho));// d(log(prob))/drho	
			return Nobs * mylog( p );
		}
	};

	struct aux_l_t {			// cofactors for low frequency (LF) polys...
		rstore fCon;			// constant
		rstore fRho;			// cofactor of rho
		rstore fRhoGamma;		// cofactor of rho*gamma
		rstore fRhoEtaGamma;	// cofactor of rho*eta*gamma
		rstore Nobs;			// Number of times this was oberved in current block. might be ok as be uint16
		inline void load(const parameters &model, double ZeqXmaj, double ZeqXmin, double nobs) {
			double beta1 = model.beta.b1, beta3 = model.beta.b3, lambda = model.block.lambdaT, theta = model.block.theta;
			fCon = ZeqXmaj                    *((1.0 - lambda)*beta1 + lambda*beta3 / 3.0)*theta / 3.0;
			fCon += ZeqXmin                    *((1.0 - lambda)*beta3 + lambda*beta1 / 3.0)*theta / 3.0;
			fCon += (1.0 - (ZeqXmin + ZeqXmaj))*((lambda / 3.0)*(beta1 + beta3)*(theta / 3.0));
			fRho = -fCon;
			fRhoGamma = ZeqXmaj*theta / 3.0;
			fRhoEtaGamma = -ZeqXmaj*lambda*theta / 3.0;
			Nobs = nobs;
		};
		inline double calcL1(const aux_factors &x) const { double d= fCon + fRho * x.rho + fRhoGamma * x.rhogamma + fRhoEtaGamma * x.rhoetagamma; return(d<0||d>1? DBL_MIN : d); };
		inline double calcL(const aux_factors &x) const { return(pow(calcL1(x), Nobs)); };
		inline double calcLL(const aux_factors &x) const { return(Nobs * mylog(calcL1(x))); };
		inline double calcLLD(const aux_factors &x, parameters &d) const {
			double p = calcL1(x); d.clear();
			if (p <= 0) { return(0); };
			double pi = 1.0 / p;
			d.rho = Nobs * (pi * (fRho + fRhoGamma * x.gamma + fRhoEtaGamma * x.etagamma));	// d(log(prob))/drho
			d.eta = Nobs * (pi * (fRhoEtaGamma * x.rhogamma));		// d(log(prob))/dgamma
			d.gam = Nobs * (pi * (fRhoGamma * x.rho + fRhoEtaGamma * x.rhoeta));	// d(log(prob))/deta
			return(Nobs * mylog(p));
		};

	};
	bool initialized;
	std::vector<aux_m_t>	aux_m;			// precomputed data for monomorphic sites
	std::vector<aux_l_t>	aux_l;			// precomptued data for low  frequencey polymorphisms
	std::vector<aux_h_t>	aux_h;			// precomputed data for high frequency polymorphisms.

	// Model is only used to grab global paremeters like beta.
	bool	initialize( const BlockSet_t &Data, double PriorWeight = 0.0, const modelInsight *Prior=NULL );
	double	calcNLL(parameters *Derivs = NULL , const parameters *OverrideParams=NULL ) const;
	void	clear( ) { initialized = false; aux_m.clear(); aux_l.clear(); aux_h.clear();  };
	modelInsightAccelerator(const modelInsight &M) : modelInsight(M) { clear( ); };
};

struct modelInsightAcceleratorBeta : modelInsight {
	// accelerated relative likelihood estimator for P(d | S=0, Y=L, Z\in{Xmin,Xmaj}, this is used to find MAp / ML values for
	//	beta1 and beta3, generally over blocks containing all neutral polymorphisms.... only data with Y=L is sensative to 
	//	beta1 and beta3 values. As beta1 + beta2 + beta3 =1 and beta2 is calculated seperately, beta1+beta3 is constant.
	// The only classes sensative to the porportionate valeus of beta1 and beta2 have ZeqXmaj or ZeqXmin, so we condition
	//	observations on THOSE clases (entries 3 & 4 in table 3, doi:10.1093/molbev/mst019).
	typedef modelInsightAccelerator::rstore rstore;		// generally double. Try float if these structures become too big....

	struct aux_beta {
		rstore A;		// cofactor of beta3 
		rstore B;		// constant
		double Nobs;	// Number of times this was oberved in current block.
		inline void load( const parameters &model, double zxmaj, double zxmin, double nobs) {
			double asb1, asb3, lambda = model.block.lambdaT, beta2 = model.beta.b2;
			// Sum of weightes entries for P(d|S=0,Y=L) 3, 4 & 5 in table 3. P(Z=Xmin)*L(d|model,Z)+P(Z=Xmax)*L(d|model,Z).
			//	For efficency remove constants, then any shared factors (theta) and rearragne in terms of Beta1*abs1 + Bta3*abs3 + const.
			asb1 = zxmaj * (1.0 - lambda) + zxmin*lambda / 3.0;	// cofactor for beta1
			asb3 = zxmin * (1.0 - lambda) + zxmaj*lambda / 3.0;	// cofactor for beta3 
			// b1 + b2 + b3 = 1. Likelihood porp asb1 * b1 + asb3 * b3 = asb1 * (1-b2-b3) + asb3 * b3 = b3 * (asb3 - asb1) + asb1 * (1.0-b2).
			// this is no longer true likelihood(beta3), but it is porportionate to it.
			A = (asb3 - asb1); B = asb1 * (1.0-beta2); Nobs = nobs;				// beta3 ~ .04, (.02-.07), so looking for beta3 is more numerically stable than looking for beta1
		}
		inline double calcL1( double beta3 ) const { return( A * beta3 + B ); };
		inline double calcL(double beta3) const { return(pow(calcL1(beta3), Nobs)); };
		inline double calcLL(double beta3) const { return(Nobs * mylog(calcL1(beta3))); };
		inline double calcLLD(double beta3, double &d ) const {
			double p = calcL1(beta3); d=0.0;
			if (p <= 0) { return(0); };
			double pi = 1.0 / p;
			d = Nobs * (pi * A );	// d(log(prob))/dbeta3
			return(Nobs * mylog(p));
		};
	};

	std::vector<aux_beta>	auxen;			// precomputed data for monomorphic sites
	bool initialized;

	bool	initialize(const BlockSet_t &data, double PriorWeight = 0.0, const modelInsight *Prior = NULL);
	double	calcNLL(parameters *Derivs = NULL) const;
	void	clear( ) { initialized = false; auxen.clear(); };
	modelInsightAcceleratorBeta(const modelInsight &M) : modelInsight(M) { clear(); };
};

bool modelInsightAccelerator::initialize(const BlockSet_t &Data, double PriorWeight , const modelInsight *Prior) {
	// set Priorweight == 0, or Prior = NULL (the defaults) to disable use of priors and just calculate ML values rather than MAp
	const BlockSet_t			&bs = Data;
	ulong						nblocks = (ulong) bs.blocks.size();	// 64-> 32 bit, many positions, but fewer blocks.
	BlockSet_t::AncPriInd_t		xz_maji = 0,   xz_mini = 0;			// indices
	BlockSet_t::AncPriVal_t		xz_maj  = 0.0, xz_min  = 0.0;		// values
	BlockSet_t::AncPriCount_t	nobs = 0;							// counts = # of M/L/H observations in data set.
	double						wa = wattersonsA(bs.numAlleles);	// Assumed constant (full data)
	double						obsweight = bs.meanWeight;			// scale factor representing 1.0/(max nobs per genomic position)

	aux_h_t	tmp_h; aux_l_t	tmp_l; aux_m_t	tmp_m;					// space optimized data structures for holding sufficient statistics and weights.

	clear(); aux_m.reserve(bs.numMonoPat()+3);	aux_l.reserve(bs.numPolyLPat()+3);	aux_h.reserve(bs.numPolyHPat()+3); // clear and reallocate memory...

																														// loop over blocks, for each block handle mono, PolyL and polyH in subloops.
	for (ulong iblock = 0; iblock < nblocks; iblock++) {
		auto &bl = bs.blocks[iblock];									// get current block from database
		bl.getHeader(NULL, NULL, NULL, &model.block.lambdaT, &model.block.theta); // get lambda from block... can provide first args if debugging is desired...
		model.block.theta *= wa;	// This is the only time we inclide wa in a param::theta value.
		for (ulong i = 0; i < bl.polyHNumPat(); i++) {
			bl.polyHVal(i, xz_maji, xz_mini, nobs);								// get index into prob ability table for z=Maj/Min human allele
			xz_maj = *(bs.ancPri[xz_maji]); xz_min = *(bs.ancPri[xz_mini]);		// lookup indices and convert nobs into double
			tmp_h.load(model, xz_maj, xz_min, nobs*obsweight); aux_h.push_back(tmp_h);	// Calculate succifient statistics and save...
		};
		for (ulong i = 0; i < bl.polyLNumPat(); i++) {
			bl.polyLVal(i, xz_maji, xz_mini, nobs);
			xz_maj = *(bs.ancPri[xz_maji]); xz_min = *(bs.ancPri[xz_mini]);
			tmp_l.load(model, xz_maj, xz_min, nobs*obsweight); aux_l.push_back(tmp_l);
		};
		for (ulong i = 0; i < bl.monoNumPat(); i++) {
			bl.monoVal(i, xz_maji, nobs);
			xz_maj = *(bs.ancPri[xz_maji]);
			tmp_m.load(model, xz_maj, nobs*obsweight); aux_m.push_back(tmp_m);
		};
	};

	// Dont weigh priors. Priors are demoninaetd in number of physiocal observations, not fractional observations.
	if (PriorWeight > 0) {
		// use priors to generate posterior (normalized) data likelihood for each observation class 
		model.block.theta = Prior->getparameters().block.theta * wa; model.block.lambdaT = Prior->getparameters().block.lambdaT;
		butils::ezmatrix counts;

		Prior->likelihoodMatrix(bs.numAlleles, true, counts);	// get prior distribution over input states. This is normalized!

		// Push back psudocounts, as expected observations of input data states (combinatiosn of sufficient properties)
		SelectionClass sc;
		double xma = 0, xmi = 0, xot = 0;	// P(Y) = x = xma+xmi+xot. Probs of Z==Xma, Z==Xmi, Z==Xother
		sc = SelectionClass::Mono;  xma = counts.v((int)sc, (int)AncestralAllele::Xmaj); xmi = counts.v((int)sc, (int)AncestralAllele::Xmin); xot = counts.v((int)sc, (int)AncestralAllele::Xoth);
		if (xma>0) { tmp_m.load(model, 1.0, xma*PriorWeight); aux_m.push_back(tmp_m); };
		//if (xmi>0) { tmp_m.load(model, ??, xma*PriorWeight); aux_m.push_back(tmp_m); }; // class = Monomorphic... there IS no minor allele state, only Major and Other. xmi should always be 0!
		if (xot>0) { tmp_m.load(model, 0.0, xot*PriorWeight); aux_m.push_back(tmp_m); };	

		sc = SelectionClass::PolyL; xma = counts.v((int)sc, (int)AncestralAllele::Xmaj); xmi = counts.v((int)sc, (int)AncestralAllele::Xmin); xot = counts.v((int)sc, (int)AncestralAllele::Xoth); 
		if (xma>0) { tmp_l.load(model, 1.0, 0.0, xma*PriorWeight); aux_l.push_back(tmp_l); };
		if (xmi>0) { tmp_l.load(model, 0.0, 1.0, xmi*PriorWeight); aux_l.push_back(tmp_l); };
		if (xot>0) { tmp_l.load(model, 0.0, 0.0, xot*PriorWeight); aux_l.push_back(tmp_l); };

		sc = SelectionClass::PolyH; xma = counts.v((int)sc, (int)AncestralAllele::Xmaj); xmi = counts.v((int)sc, (int)AncestralAllele::Xmin); xot = counts.v((int)sc, (int)AncestralAllele::Xoth);
		if (xma>0) { tmp_h.load(model, 1.0, 0.0, xma*PriorWeight); aux_h.push_back(tmp_h); };
		if (xmi>0) { tmp_h.load(model, 0.0, 1.0, xmi*PriorWeight); aux_h.push_back(tmp_h); };
		if (xot>0) { tmp_h.load(model, 0.0, 0.0, xot*PriorWeight); aux_h.push_back(tmp_h); };
	}

	initialized = ((aux_m.size()+aux_l.size()+aux_h.size()) > 1);
	return(initialized);
}
double modelInsightAccelerator::calcNLL(parameters *Derivs, const parameters *NewParams ) const {	
	// Pass NULL in to skip deravitive calculation (faster).
	//	Pass in override parameters to perserve exsiting model...
	// The sufficient statistics transform the likelihood estimates into a simple form that utilises
	//	A + B * Rho + C * Rho*Eta + D * Rho*Gamma + E *  Rho * Eta * Gamma. For each of H, L, M dat
	//  may be implicitely 0.
	if (!initialized) return false ;

	aux_factors cofac; 
	if (NewParams != NULL)	cofac.load(NewParams->rho, NewParams->eta, NewParams->gam );
	else					cofac.load(model.rho, model.eta, model.gam);
	double ll = 0.0; ulong ns;
	if (Derivs == NULL ) {
		ns = (ulong)aux_m.size();	// Process Mono
		for (ulong isite = 0; isite < ns; isite++) ll += aux_m[isite].calcLL(cofac);

		ns = (ulong)aux_l.size();	// Process PolyL
		for (ulong isite = 0; isite < ns; isite++) ll += aux_l[isite].calcLL(cofac);

		ns = (ulong)aux_h.size();	// Process PolyH
		for (ulong isite = 0; isite < ns; isite++) ll += aux_h[isite].calcLL(cofac);
	} else {
		parameters tp, &d=(*Derivs); tp.clear(); d.clear();
		ns = (ulong)aux_m.size();	// Process Mono
		for (ulong isite = 0; isite < ns; isite++) { 
			ll += aux_m[isite].calcLLD(cofac,tp);
			d.rho += tp.rho; d.eta += tp.eta; d.gam += tp.gam; 
			// TODO we can get tp.rho = -inf if an egregious value of eta or gamma is provided... p->DBLMIN, (1/p)*NumDeriv->inf... limit search space a bit.
			//	This should be handled better. --Brad
		};

		ns = (ulong)aux_l.size();	// Process PolyL
		for (ulong isite = 0; isite < ns; isite++) { ll += aux_l[isite].calcLLD(cofac,tp); d.rho += tp.rho; d.eta += tp.eta; d.gam += tp.gam; };

		ns = (ulong)aux_h.size();	// Process PolyH
		for (ulong isite = 0; isite < ns; isite++) { ll += aux_h[isite].calcLLD(cofac,tp); d.rho += tp.rho; d.eta += tp.eta; d.gam += tp.gam; };

		d.rho = -d.rho; d.eta = -d.eta; d.gam = -d.gam;	// We want deravitives of NLL not LL.
	};
	return(-ll);
}

bool modelInsightAcceleratorBeta::initialize(const BlockSet_t &Data, double PriorWeight, const modelInsight *Prior) {
	// This is for neutral poly blocks! NOT for primary blocksets. That is how beta is defined.
	// set Priorweight == 0, or Prior = NULL (the defaults) to disable use of priors and just calculate ML values rather than MAp
	const BlockSet_t			&bs = Data;
	ulong						nblocks = (ulong)bs.blocks.size();	// 64-> 32 bit, many positions, but fewer blocks.
	BlockSet_t::AncPriInd_t		xz_maji = 0, xz_mini = 0;			// indices
	BlockSet_t::AncPriVal_t		xz_maj = 0.0, xz_min = 0.0;			// values
	BlockSet_t::AncPriCount_t	nobs = 0;							// counts = # of M/L/H observations in data set.
	//double						wa = wattersonsA(bs.numAlleles);	// Assumed constant (full data)

	aux_beta tmp_b;							// space optimized data structures for holding sufficient statistics and weights.

	clear(); auxen.reserve(bs.numPolyLPat() + 3);	// clear and reallocate memory...
																							 // loop over blocks, for each block handle mono, PolyL and polyH in subloops.
	for (ulong iblock = 0; iblock < nblocks; iblock++) {
		auto &bl = bs.blocks[iblock];									// get current block from database
		bl.getHeader(NULL, NULL, NULL, &model.block.lambdaT, NULL);		// get lambda from block... can provide first args if debugging is desired...
		// model.block.theta *= wa;	// This is the only time we inclide wa in a param::theta value. Beta calculations dont use thetaA, so we don't have to calculate it..
		for (ulong i = 0; i < bl.polyLNumPat(); i++) {
			bl.polyLVal(i, xz_maji, xz_mini, nobs);
			xz_maj = *(bs.ancPri[xz_maji]); xz_min = *(bs.ancPri[xz_mini]);
			tmp_b.load(model, xz_maj, xz_min, nobs); auxen.push_back(tmp_b);
		};
	};

	if (PriorWeight > 0) {
		// use priors to generate posterior (normalized) data likelihood for each observation class 
		model.block.lambdaT = Prior->getparameters().block.lambdaT;  // model.block.theta = Prior->getparameters().block.theta * wa;
		butils::ezmatrix counts;

		Prior->likelihoodMatrix(bs.numAlleles, true, counts);	// get prior distribution over input states. This is normalized!

		// eliminate monomorphics sites, and renormalize to polymorphic sites. Prior Model beta2 will redistribute some counts to PolyH. Thats Good.
		for (int j = 0; j < counts.ncol(); j++) counts.v((int) SelectionClass::Mono,j)=0;
		counts.norm1();

		// Push back psudocounts, as expected observations of input data states (combinatiosn of sufficient properties)
		SelectionClass sc;
		double xma = 0, xmi = 0;	// P(Y) = x = xma+xmi+xot. Probs of Z==Xma, Z==Xmi, Z==Xother
		sc = SelectionClass::PolyL; xma = counts.v((int)sc, (int)AncestralAllele::Xmaj); xmi = counts.v((int)sc, (int)AncestralAllele::Xmin); 
		// xot = counts.v((int)sc, (int)AncestralAllele::Xoth); // not used, so dont set...
		if (xma>0) { tmp_b.load(Prior->getparameters(), 1.0, 0.0, xma*PriorWeight); auxen.push_back(tmp_b); };
		if (xmi>0) { tmp_b.load(Prior->getparameters(), 0.0, 1.0, xmi*PriorWeight); auxen.push_back(tmp_b); };
		// tmp_b.load(model, 0.0, 0.0, xot*PriorWeight); auxen.push_back(tmp_b); // constant term ignored in likelihood maximization...
	};
	initialized = (auxen.size() > 1);	// ML can fail if there is only one data point...
	return(initialized);
}

double modelInsightAcceleratorBeta::calcNLL(parameters *Derivs) const {	
	double beta3 = model.beta.b3;
	double ll = 0.0; ulong ns;
	if (Derivs == NULL) {
		ns = (ulong)auxen.size();	// Process PolyL
		for (ulong isite = 0; isite < ns; isite++) ll += auxen[isite].calcLL(beta3);
	} else {
		double tp, &d = Derivs->beta.b3; tp = d = 0.0;
		ns = (ulong)auxen.size();	// Process PolyL
		for (ulong isite = 0; isite < ns; isite++) { ll += auxen[isite].calcLLD(beta3, tp); d += tp; };
		d = -d;	// we want deravitive of NLL, not LL
	}
	return(-ll);	// return NLL not LL.
}
