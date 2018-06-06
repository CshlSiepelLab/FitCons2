#pragma once
#pragma once
#include <cassert>
#include <vector>
#include <algorithm>	//std::sort
#include <iostream>		// Used for slow / debugging  IO
#include <string>

// Derive a class from optimizer for inferring Betas...
#include "butils/butils.h"		// includes optimizer...
#include "inscomp/inscomp.h"
#include "BlockSet.h"

#include <iostream>	// debugging
#include <cmath>	// sqrt...

class modelInsight  {
public:
	typedef butils::ezmatrix	ezmatrix;
	typedef butils::ulong		ulong;		// central place to define long int type. for 1 cell type 32 bits is ok, otherwise got with 64.
	typedef std::string			string;
	typedef inscomp::BlockSet_t	BlockSet_t;
	typedef std::vector<double>	vecD;

	// Polymorphism class, Monomorphic, low frequence Polymorphisn, or High Frequency Polymotphism
	enum class SelectionClass { first = 0, Mono = first, PolyL, PolyH, last = PolyH, size = (last + 1) };
	// Ancestral State (Z) matches population Major Allele, Minor Allele or neither (oth).
	enum class AncestralAllele { first = 0, Xmaj = first, Xmin, Xoth, last = Xoth, size = (last + 1) };

	// Key model parameters.
	struct parameters {
		// Data Types..
		struct block_t { double lambdaT, theta; void clear() {lambdaT=theta=-1; }; };	// values typically derivbed from a single polymorphism block
		struct betas_t { double b1, b2, b3; void clear() {b1=b2=b3=-1;}; string toStr() const { return( std::to_string(b1) + "\t" + std::to_string(b2) + "\t" + std::to_string(b3) ); }; };				// values estimated from neutral poly statistics over collection of blocks

		// Data Elements
		double rho, eta, gam;								// essential model parameters, derived from all positions in all blocks
		block_t block;										// values typically derivbed from a single polymorphism block
		betas_t beta;										// values estimated from neutral poly statistics over all blocks

		// Functions
		parameters() { clear(); };
		void clear() { rho = eta = gam = 0.0; beta.clear(); block.clear(); };
		string fmt(double v, const char *f = "%25.20lf") const { static char b[32]; sprintf(b, ((v==0) || v >= 1e-6 ? f : "%25.15le"), v); return(b); };
		string toStr() const {
			string s = "Rho:\t" + fmt(rho) + "\nEta:\t" + fmt(eta) + "\nGamma:\t" + fmt(gam) + "\n";
			s += string("Lambda:\t") + fmt(block.lambdaT) + "\n"; s += "Theta:\t" + fmt(block.theta) + "\n";
			s += "Beta1:\t" + fmt(beta.b1) + "\nBeta2:\t" + fmt(beta.b2) + "\nBeta3:\t" + fmt(beta.b3) + "\n";
			// s += "Alpha:\t" + fmt(sup.alpha) + "\nTau:\t" + fmt(sup.tau) + "\tDp:\t" + fmt(sup.dp) + "\nPw:\t" + fmt(sup.pw) + "\n";
			return s;
		};
	};
protected:
	parameters				model;			// as we refine model, gather aggregate statistics here, including E[lambda] E[theta] etc.....
	inline static double mylog( double d ) { return (butils::mathplus::mylog( d )); };
	inline	void	likelihoodRaw(ulong Nallele, SelectionClass Y, AncestralAllele Z, double &PS0, double &PS1, double &Weight) const;
public:
	void clear() { model.clear(); };
	inline const parameters &getparameters() const { return model; };
	inline void setparameters( const parameters &p ) { model = p; return;};
	string toStr() const { return( model.toStr() ); };

	static double wattersonsA(uint32_t N) { 
		static std::vector<double>	w;
		if (N<2) return( 0.0);
		if (N >= w.size() ) {	// generally this only happens once....
			w.resize(N+1); w[0]=0.0; w[1]=0.0;
			for (int i = 2; i <= (int)N; i++) w[i] = w[i-1] + (1.0 / (double)(i - 1)); 
		};
		return( w[N] );
	};
	
	// Calculate prior probability of propsoed as expected likelihood of Nobs pseudocounts drawn from *this (Population) Model.
	double Prior(const modelInsight &PriorModel, ulong Nallele, double Nobs = 1.0, const ezmatrix *PriorDist = NULL) const;
	double PriorLog(const modelInsight &PriorModel, ulong Nallele, double Nobs = 1.0, const ezmatrix *PriorDist = NULL) const;

	// Calculate Likelihood, that is Class probability based in parameters, P(Y,X|Model).
	//	Lambda is genrally multiplied by T before storage in the database... Nallele is used to genreate theta*A
	//	X is provided as an enunerated constant representing one sufficient condition: Xmaj=Z, Xmin=Z or Z!=Xmaj,Xmin
	// ByClass =false generates valeus for a single allele observation (1 position). ByClass=true generates weighted 
	//	values by input class: Y, Xmaj=Z, Xmin=Z or Z!=Xmaj,Xmin.
	inline	double	likelihood(ulong Nallele, SelectionClass Y, AncestralAllele Z, bool ByClass = false, parameters *Deravitive= NULL ) const;
	inline	void	likelihoodMatrix(ulong Nallele, bool ByClass, ezmatrix &Joint) const;
	inline	double	likelihoodExpected(ulong Nallele, SelectionClass Y, double pZeqXmaj, double pZeqXmin, bool ByClass = false) const  {
		return(pZeqXmaj * likelihood(Nallele, Y, AncestralAllele::Xmaj, ByClass) + pZeqXmin * likelihood(Nallele, Y, AncestralAllele::Xmin, ByClass) + (1.0 - pZeqXmaj - pZeqXmaj) * likelihood(Nallele, Y, AncestralAllele::Xoth, ByClass));
	}
	double	likelihoodDataLog(const BlockSet_t &bs, double PriorWeight = 0.0, const modelInsight *Prior=NULL ) const;
	void posteriorByClass(ulong NAllele, SelectionClass Y, ezmatrix &pd_S0, ezmatrix &pd_S1 ) const ;
	void posteriorByObs(ulong NAllele, SelectionClass Y, ezmatrix &pd_S0, ezmatrix &pd_S1) const;
	void posteriorMatrix(ulong NAllele, ezmatrix &pd_S0, ezmatrix &pd_S1, bool ByClass = true) const ;
	void posteriorJointAlleleAndS(ulong NAllele, SelectionClass Y, ezmatrix &dist) const;

	// Beta functions are CONDITIONED on neutral (S=0) polymorphic sites (PolyH or PolyL). Generally these are only applied to 
	//	data sets consisting of neutral polymorphisms, used to infer beta values for the model
	inline	void	likelihoodBetaMatrix(ulong Nallele, bool ByClass, ezmatrix &Joint) const;
	inline	double	likelihoodBeta(ulong Nallele, SelectionClass Y, AncestralAllele Z, bool ByClass = false, parameters *Deravitive = NULL) const;
			double	likelihoodBetaDataLog(const BlockSet_t &bs, double PriorWeight, const modelInsight *Prior) const;
};

// Proportionate to the prior probability density of the current model, under the PriorModel.
double	modelInsight::PriorLog(const modelInsight  &PriorModel, ulong Nallele, double Nobs, const ezmatrix *PriorDist) const {
	// RELATIVE Prior density (unnoralized) derived from intuition provided by pseudocounts, given rigor by Dirchlet priors.
	//	Start with an initial uniform prior on model parameters in M.
	//	The likelihood of an observation d is thus P(d|M)=P(M|d)P(d)/P(M). As P(d) is a constant, the uniform prior provides that P(M|d) \porp P(d|M).
	//	For an expected d of D, I use a M generated from the whole genome (population) and generate expected distribution over data d (or at least,
	//	sufficient properties of D) then set the prior to P(M') = P(D|M'). To apply weight consitent with N pseudocounts, simply
	//	take P(D|M')^N. This prior, related to the KL distance between M and M' biases models towars those that represent the global
	//	population, inducing statistical shrinkage.
	const ezmatrix *pc = PriorDist; double ll=0.0;
	if (Nobs <= 0.0) return(1.0);
	if (pc == NULL) {	// gerate distribution of input states, according to global (population) model. This is the expected distribution over states.
		ezmatrix *pc_tmp = new ezmatrix();
		PriorModel.likelihoodMatrix(Nallele, true, *pc_tmp); pc = pc_tmp;	// if doign a of this, calculate ED once at parent level and pass in...
	}

	// Here is where we get into trouble. Do we want expected log likelihood (A), or log expected likelihood? (B). E[F(X)] >= F(E[X])
	// A) sum over [ Nobs * pc->v(i, j) * log( Proposed->v(i,j) ) ]
	// B) Nobs * log( [ sum pc->v(i, j) * Proposed->v(i,j)  ] )...
	// Tecall: over classes i, N=sum_Ni, P(D|P) = Prod_i[Pi^Ni]*(Combinatorics in Ni and N). The latter term is constant if expected Distrib doesnt change.
	//		This turns into something porportionate to A), sum Ni Log(Pi) = sum N*pi log(Pi), where N*pi are pseudocounts....
	//	Go with the former..... (ie pseudocounts)
	// Calculate expected probbaility of global observation, under proposed model. This estimates KL distance...
	for (int i = (int)SelectionClass::first; i <= (int)SelectionClass::last; i++) {
		for (int j = (int)AncestralAllele::first; j <= (int)AncestralAllele::last; j++) {
			double pseudocounts = Nobs * pc->v(i, j);
			ll += pseudocounts * mylog( likelihood(Nallele, (SelectionClass)i, (AncestralAllele)j ) );	
		};
	};
	if (PriorDist == NULL) { delete pc; pc = NULL; }
	return ll;
}

double	modelInsight::Prior(const modelInsight  &PriorModel, ulong Nallele, double Nobs, const ezmatrix *PriorDist ) const {
	// while prob of 1 observation is generally > 1e-9, prob^Nobs might underflow in unlogged format, so primal form is Logspace...
	return exp( PriorLog( PriorModel, Nallele, Nobs, PriorDist ));
}

// Private! Encapsulates table 3 from paper...
void	modelInsight::likelihoodRaw(ulong NAllele, SelectionClass Y, AncestralAllele Z, double &PS0, double &PS1, double &Weight) const {
	// P(d|M) directly From Table3, doi:10.1093/molbev/mst019 . 
	double LambdaT = model.block.lambdaT, ThetaA = model.block.theta * wattersonsA(NAllele);
	double Eta = model.eta, Gamma = model.gam;
	double Beta1 = model.beta.b1, Beta2 = model.beta.b2, Beta3 = model.beta.b3;
	PS0 = PS1 = Weight = 0.0;
	switch (Y) {
		case SelectionClass::Mono:
			switch (Z) {
				case AncestralAllele::Xmaj:	PS0 = (1.0 - LambdaT)*(1 - ThetaA);		PS1 = (1.0 - Eta*LambdaT)*(1.0 - Gamma*ThetaA);		Weight = 1.0;	break;
				case AncestralAllele::Xmin:	PS0 = 0.0;								PS1 = 0.0;											Weight = 0.0;	break;	// M means there is NO minor X allele.
				case AncestralAllele::Xoth:	PS0 = (1 - ThetaA)*LambdaT / 3.0;		PS1 = Eta*LambdaT / 3.0;							Weight = 3.0;	break;
				default: break;
			}; break;
		case SelectionClass::PolyL:
			switch (Z) {
				case AncestralAllele::Xmaj:	PS0 = ((1.0 - LambdaT)*Beta1 + LambdaT*Beta3 / 3.0)*ThetaA / 3.0;	PS1 = (1.0 - Eta*LambdaT)*Gamma*ThetaA / 3.0;	Weight = 3.0; break;
				case AncestralAllele::Xmin:	PS0 = ((1.0 - LambdaT)*Beta3 + LambdaT*Beta1 / 3.0)*ThetaA / 3.0;	PS1 = 0.0;										Weight = 3.0; break;
				case AncestralAllele::Xoth:	PS0 = (LambdaT / 3.0)*(Beta1 + Beta3)*ThetaA / 3.0;					PS1 = 0.0;										Weight = 6.0; break;
				default: break;
			}; break;
		case SelectionClass::PolyH:
			PS1 = 0.0;
			switch (Z) {	// logic for weights is nuanced here... see Normalization table.
				case AncestralAllele::Xmaj:	PS0 = (1.0 - 2.0*LambdaT / 3.0)*Beta2*ThetaA / 3.0;	PS1 = 0.0;	Weight = 1.5; break;
				case AncestralAllele::Xmin:	PS0 = (1.0 - 2.0*LambdaT / 3.0)*Beta2*ThetaA / 3.0;	PS1 = 0.0;	Weight = 1.5; break;
				case AncestralAllele::Xoth:	PS0 = (2.0*LambdaT*Beta2 / 3.0)*ThetaA / 3.0;		PS1 = 0.0;	Weight = 3.0; break;
				default: break;
			}; break;
		default: break;
	};
	return;
}

// get one row of matrix P(X,Y|M,S) for a single Selection Class Y.
void modelInsight::posteriorJointAlleleAndS(ulong NAllele, SelectionClass Y, ezmatrix &dist) const {
	double w=0;
	ezmatrix s0(1,(int)AncestralAllele::size), s1(1,(int)AncestralAllele::size);
	dist.resize(2, (int) AncestralAllele::size);	// joint posterior distribution S\in{0,1} & Allele
	//posteriorByClass(NAllele, Y, s0, s1);			// This is P(SufProp(X),Y,S,Z|Model), for one Y..
	posteriorByObs(NAllele, Y, s0, s1);				// This is simple P(X,Y,S|Z,Model), for a particular X for one Y....
	for (int i = (int)AncestralAllele::first; i <= (int)AncestralAllele::last; ++i) { 
		// There is always 1 Xmaj, but M has 0 Xmin (1 otherwise) and M has 3 Xoth, while L,H have 2 Xoth... distribution must sum to 1....
		w=1.0;  if (i == (int)AncestralAllele::Xoth) w = (Y == SelectionClass::Mono ? 3.0 : 2.0);
		dist.v(0,i)=s0.v(0,i)*w; dist.v(1, i) = s1.v(0, i)*w; 
	}
	dist.norm1();									// Renormalize to get posterior P(X,S|Y,M) for one site...
}

// get one row of matrix P(X,Y,S|M) for a single Selection Class Y.
void modelInsight::posteriorByObs(ulong NAllele, SelectionClass Y, ezmatrix &pd_S0, ezmatrix &pd_S1) const {
	double s0 = 0, s1 = 0, w = 0;
	pd_S0.resize(1, (int)AncestralAllele::size); pd_S1.resize(1, (int)AncestralAllele::size);
	for (int i = (int)AncestralAllele::first; i <= (int)AncestralAllele::last; ++i) {
		likelihoodRaw(NAllele, Y, (AncestralAllele)i, s0, s1, w);
		pd_S0.v(0, i) = s0*(1-model.rho); pd_S1.v(0, i) = s1 * model.rho;
	}
}

// get one row of matrix P(X,Y,S|M) for a single Selection Class Y.
void modelInsight::posteriorByClass(ulong NAllele, SelectionClass Y, ezmatrix &pd_S0, ezmatrix &pd_S1) const {
	double s0=0, s1=0, w=0;
	pd_S0.resize(1,(int)AncestralAllele::size); pd_S1.resize(1,(int) AncestralAllele::size);
	for (int i = (int)AncestralAllele::first; i <= (int)AncestralAllele::last; ++i) {
		likelihoodRaw( NAllele, Y, (AncestralAllele) i, s0, s1, w); s0*=w;s1*=w;
		pd_S0.v(0,i)=s0*(1 - model.rho); pd_S1.v(0,i)=s1 * model.rho;
	}
}

// get full matrix of P(X,Y|M,S)
void modelInsight::posteriorMatrix(ulong NAllele, ezmatrix &pd_S0, ezmatrix &pd_S1, bool ByClass ) const {
	ezmatrix s0, s1;
	pd_S0.resize((int) SelectionClass::size, (int) AncestralAllele::size); pd_S1.resize((int)SelectionClass::size, (int)AncestralAllele::size);
	for (int i = (int)SelectionClass::first; i <= (int)SelectionClass::last; ++i) {
		if (ByClass ) {		posteriorByClass(NAllele, (SelectionClass) i , s0, s1);	
		} else {			posteriorByObs(NAllele, (SelectionClass) i , s0, s1);	};
		for ( int j= (int)AncestralAllele::first; j <= (int)AncestralAllele::last; ++j) { pd_S0.v(i,j) = s0.v(0,j); pd_S1.v(i, j) = s1.v(0,j); };
	}
}


// Get data likelihood for an observation with the given sufficient properties N, Y, Z.
double	modelInsight::likelihood(ulong NAllele, SelectionClass Y, AncestralAllele Z, bool ByClass, parameters *Deravitive) const  {
	// To get values conditional on S=0, set Rho=0.0. To get values conditional on S=1, set Rho=1.0. 
	// A data observation d, consists if an allele{a,c,g,t} Observation of alleles for Z,Xmaj,Xmin and a polymorphism class Y.
	//   expected observations become distributions over the 16 allele combinations an d 3 classes of Y... that sum to 1.
	// The sufficient properties (SP) to calculate the likelihood of that data are Y,Z=Xmaj,Z=Xmin,Z=Xoth. For a specific d, these properties are 
	//	sufficient to determines for P(d|model). However this probability is NOT P(SP(d)|Model) because diffeing d map to the same SP.
	//	Consider a 6 sided die M with 1-4 colored blue and 5-6 colored red. Allow sufficient properties P(d)=1/8 iff blue, 1/4 if red.
	//	Thus, p(d=3) = 1/8 but this is NOT P(SP(d)|M) = P(blue|M) = 1/2.
	//	Similarly, ZeqX(d) is a sufficient property to calculate P(d) but this probability is NOT P(ZeqX).
	// The following table provides P(d|M) via the SP Y, X_maj==Z and X_min==Z, when Class == false, and P(sp(d)|M) when Class = true.
	// Also, I know that Z is not actually observed, but we estiamte it in a seperate step and treat it as observed for the purposes of INSIGHT.
	double ps0 = 0.0, ps1 = 0.0, p = 0.0, w = 1.0;		// Combinations have probability 0 unless identified below.....
	likelihoodRaw( NAllele, Y, Z, ps0, ps1,w);
	if (ByClass) { ps0 *= w; ps1 *= w; };	// Upweight by class, to get normalized probability, by class (Selec / Ancestral Allele)
	p = (model.rho * ps1) + (1.0 - model.rho)*ps0;

	if (Deravitive != NULL) {
		// Very SLOW, and numerical, if we use it much, make this an analytic deravitive.
		modelInsight mt = *this; double *pr=NULL, delta=1e-6; parameters &d=(*Deravitive);
		// take a numerical deravitive......
		auto numd = [&] ( modelInsight &m, double *par, bool ubound ) {
			double s=*par,x1=s-delta,x2=s+delta,y1,y2; if (x1<0) x1=0; if (ubound&&(x2>1))x2=1;
			*par=x1; y1=m.likelihood( NAllele,Y,Z,ByClass); *par = x2; y2 = m.likelihood(NAllele, Y, Z, ByClass); *par=s;
			return( (y2-y1)/(x2-x1) ); };
		pr = &mt.model.rho;		d.rho = numd(mt, pr, true);
		pr = &mt.model.eta;		d.eta = numd(mt, pr, false);
		pr = &mt.model.gam;		d.gam = numd(mt, pr, false);
		pr = &mt.model.beta.b1; d.beta.b1 = numd(mt, pr, true);
		pr = &mt.model.beta.b2; d.beta.b2 = numd(mt, pr, true);
		pr = &mt.model.beta.b3; d.beta.b3 = numd(mt, pr, true);
	}
	return p;
}

// More efficent routine for getting Posterior Likelihood of data. Set PriorWeight=0 to get simple data likeihood.
double	modelInsight::likelihoodDataLog(const BlockSet_t &bs, double PriorWeight, const modelInsight *Prior ) const {
	double zxmaj, zxmin, l, ll = 0.0;
	ulong nallele = bs.numAlleles;
	ezmatrix like;								// data log likelihood for current block
	parameters tmp_param;						// parameters for current block
	inscomp::AncPriIdx_t idxMa = 0, idxMi = 0;	// index into probility table
	inscomp::PatternCounts_t nobs;				// count of observations of this pattern

	ll = 0.0; tmp_param = model;
	for (ulong iblock = 0; iblock < bs.numBlocks(); iblock++) {
		const inscomp::Block_t &bl = bs.blocks[iblock];
		bl.getHeader(NULL, NULL, NULL, &tmp_param.block.lambdaT, &tmp_param.block.theta);
		// This is log( expected [likelihood])
		likelihoodMatrix(nallele, false, like); //like.logElem();
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
	};
	if ( PriorWeight >0 && (Prior != NULL) ) {
		ll += Prior->PriorLog( *this, nallele, PriorWeight);
	};
	return ll;
}

// Data likelihood (ByClass=false), or sufficient class probability (ByClass=true)
void	modelInsight::likelihoodMatrix(ulong Nallele, bool ByClass, ezmatrix &Joint) const  {
	Joint.resize((int)SelectionClass::size, (int)AncestralAllele::size);
	for (int i = (int)SelectionClass::first; i <= (int)SelectionClass::last; i++) {
		for (int j = (int)AncestralAllele::first; j <= (int)AncestralAllele::last; j++) {
			Joint.v(i, j) = likelihood(Nallele, (SelectionClass)i, (AncestralAllele)j, ByClass);
		}
	}
}

// Beta functions are used to debug Beta value inference. Generalyl applied only to neutral polymorphic sites.
void	modelInsight::likelihoodBetaMatrix(ulong Nallele, bool ByClass, ezmatrix &Joint) const {
	double ps0, ps1, w, pw, norm = 0;
	Joint.resize((int)SelectionClass::size, (int)AncestralAllele::size);
	for (int i = (int)SelectionClass::first; i <= (int)SelectionClass::last; i++) {
		for (int j = (int)AncestralAllele::first; j <= (int)AncestralAllele::last; j++) {
			if ( i== (int)SelectionClass::Mono ) {
				Joint.v(i, j) = 0;
			} else { 
				likelihoodRaw(Nallele, (SelectionClass)i, (AncestralAllele)j, ps0, ps1, w);
				pw = ps0 * w; norm += pw;
				Joint.v(i, j) = (ByClass ? pw : ps0 );
			}
		}
	}
	if (norm != 0 ) Joint.mult( 1.0 / norm );
}

double	modelInsight::likelihoodBeta(ulong Nallele, SelectionClass Y, AncestralAllele Z, bool ByClass, parameters *Deravitive) const {
	// Conditioned on S=0 and Y!=Mono. To accomadate normalization, likelihoodBeta is implimented in terms of likelihoodBetaMatrix.
	//	Effective, but slow and primarily of use for debugging. Use accelerator for actual likelihood claculations...
	double ret = 0.0;
	ezmatrix joint;
	likelihoodBetaMatrix( Nallele, ByClass, joint );
	ret = joint.v( (int)  Y, (int) Z );
	if (Deravitive != NULL) {
		// SLOW^2, and numerical, if we use it much, make this an analytic deravitive.
		modelInsight mt = *this; double *pr = NULL, delta = 1e-6; parameters &d = (*Deravitive);
		// take a numerical deravitive......
		auto numd = [&](modelInsight &m, double *par, bool ubound) {
			double s = *par, x1 = s - delta, x2 = s + delta, y1, y2; if (x1<0) x1 = 0; if (ubound && (x2>1))x2 = 1;
			*par = x1; y1 = m.likelihoodBeta(Nallele, Y, Z, ByClass); *par = x2; y2 = m.likelihoodBeta(Nallele, Y, Z, ByClass); *par = s;
			return((y2 - y1) / (x2 - x1)); };
		pr = &mt.model.beta.b1;		d.beta.b1 = numd(mt, pr, true);
		pr = &mt.model.beta.b2;		d.beta.b2 = numd(mt, pr, true);
		pr = &mt.model.beta.b3;		d.beta.b3 = numd(mt, pr, true);
	}
	return ret;
}

double	modelInsight::likelihoodBetaDataLog(const BlockSet_t &bs, double PriorWeight, const modelInsight *Prior) const {
	double zxmaj, zxmin, l, ll = 0.0;
	ulong nallele = bs.numAlleles;
	ezmatrix like;								// data log likelihood for current block
	parameters tmp_param;						// parameters for current block
	inscomp::AncPriIdx_t idxMa = 0, idxMi = 0;	// index into probility table
	inscomp::PatternCounts_t nobs;				// count of observations of this pattern

	ll = 0.0; tmp_param = model;
	for (ulong iblock = 0; iblock < bs.numBlocks(); iblock++) {
		const inscomp::Block_t &bl = bs.blocks[iblock];
		bl.getHeader(NULL, NULL, NULL, &tmp_param.block.lambdaT, &tmp_param.block.theta);
		// This is log( expected [likelihood])
		likelihoodBetaMatrix(nallele, false, like); //like.logElem();
		// Ignore the mono case, likelihood is 0.
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
	};
	if (PriorWeight >0 && (Prior != NULL)) {
		ezmatrix t_lik, t_prob;
		Prior->likelihoodBetaMatrix(nallele, true, t_prob);				// calculate prior model sufficient property distribution
		likelihoodBetaMatrix(nallele, false, t_lik); t_lik.logElem();	// get current model log-likelihoods for each sufficient properties.
		for (int i = 0; i<t_prob.nrow(); i++) {
			for (int j = 0; i < t_prob.ncol(); j++) {
				if (i != (int)SelectionClass::Mono){
					double pseudocounts = (t_prob.v(i, j) * PriorWeight);
					if (pseudocounts>0)ll += pseudocounts * t_lik.v(i, j);
				}
			};
		};
	};
	return ll;
}

