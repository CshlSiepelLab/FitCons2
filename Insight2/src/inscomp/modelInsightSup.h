#pragma once

#include "modelInsight.h"
#include "modelInsightAccelerator.h"

// Utility functions for implimenting legacy calculations (from INSIGHT1), and supplimentary statistics, like
//	Dp, Pw, Alpha and Tau as well as uncertainties via Information Matrix / hessian and counts.

struct modelInsightSup : modelInsight {
	typedef inscomp::modelInsightAccelerator	accel_t;
	struct derived_t { // values derived from rho, eta,gamma over all positions in all blocks
		double alpha, tau, dp, pw;
		void clear() { alpha = tau = dp = pw = -1; };
		string fmt(double v, const char *f = "%12.8lf") const { static char b[32]; v=(v<0?0:v); sprintf(b, ((v == 0) || v >= 1e-4 ? f : "%.8le"), v); return(b); };
		string toStr() const {
			string s = "Dp:\t" + fmt(dp) + "\nPw:\t" + fmt(pw) + "\nAlpha:\t" + fmt(alpha) + "\nTau:\t" + fmt(tau) + "\n"; return s;
		};
	};	

	// Supplimentary tools for handeling the Insight model. Utilities and legacy code not central to the model.
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
		double probSel()	{ double c = count(); return((c <= 0 ? 0 : (cSelNodiv + cSelDiv + cSelLow) / c)); };	// S==1
		double probWeak()	{ double c = count(); return((c <= 0 ? 0 : cSelLow / c)); };							// Y==L, S==1
		double probAdap()	{ double c = count(); return((c <= 0 ? 0 : cSelDiv / c)); };							// Y==M, Z!=A,S==1
		double count()		{ return(cSelNodiv + cSelDiv + cSelLow + cNeutNodiv + cNeutDiv + cNeutLow + cNeutHigh + cHigh); };
		void clear()		{ cSelNodiv = cSelDiv = cSelLow = cNeutNodiv = cNeutDiv = cNeutLow = cNeutHigh = cHigh = cCnt = 0.0; return; }
		void renorm()		{ double c = count(); if (c<=0) return; cSelNodiv /= c; cSelDiv /= c; cSelLow /= c; cNeutNodiv /= c; cNeutDiv /= c; cNeutLow /= c; cNeutHigh /= c; cHigh /=c; cCnt = c; };
		posteriordist()		{ clear(); }
		inline void mult( double d ) { cSelNodiv*=d; cSelDiv*=d; cSelLow*=d; cNeutNodiv*=d; cNeutDiv*=d; cNeutLow*=d; cNeutHigh*=d; cHigh*=d; cCnt*=d; };
		//inline void accum(const posteriordist &p) { cSelNodiv += p.cSelNodiv; cSelDiv += p.cSelDiv; cSelLow += p.cSelLow; cNeutNodiv += p.cNeutNodiv; cNeutDiv += p.cNeutDiv ; cNeutLow += p.cNeutLow; cNeutHigh += p.cNeutHigh; cHigh += p.cHigh; cCnt += p.cCnt; };
		std::string header(int Detail = 40) { std::string t;
			if (Detail>=40) {			t="SelNoD\t\tSelDiv\t\t\tSelLow\t\t\tNeutNoD\t\t\tNeutDiv\t\t\tNeutLow\t\t\tNeutHi\t\t\tHi\t\t\tCount"; 
			} else if (Detail >= 30) {	t = "SelNoD\t\tSelDiv\t\t\tSelLow\t\t\tNeutNoD\t\t\tNeutDiv\t\t\tNeutLow\t\t\tNeutHi\t\t\tHi";
			} else if (Detail >= 20) {	t = " Sel \t\t Pw \t\t Pa ";
			} else {					t = " Sel ";	};
			return( t );
		};
		std::string toStr(const char *f = "%8.6lf", int Detail=40) { 
			string t;
			if (Detail >= 40) {
				t = std::string("") + fmt(cSelNodiv, f) + "\t" + fmt(cSelDiv, f) + "\t" + fmt(cSelLow, f) + "\t" + fmt(cNeutNodiv, f) + "\t" + fmt(cNeutDiv, f) + "\t" + fmt(cNeutLow, f) + "\t" + fmt(cNeutHigh, f) + "\t" + fmt(cHigh, f) + "\t" + fmt(cCnt, f);
			} else if (Detail >=30) {
				t = std::string("") + fmt(cSelNodiv, f) + "\t" + fmt(cSelDiv, f) + "\t" + fmt(cSelLow, f) + "\t" + fmt(cNeutNodiv, f) + "\t" + fmt(cNeutDiv, f) + "\t" + fmt(cNeutLow, f) + "\t" + fmt(cNeutHigh, f) + "\t" + fmt(cHigh, f);
			} else if (Detail >= 20) {
				t = std::string("") + fmt(probSel(), f) + "\t" + fmt(probWeak(), f) + "\t" + fmt(probAdap(), f);
			} else {
				t = std::string("") + fmt(probSel(), f);
			}
			return(t); };
		inline string fmt(double a, const char *f = "%8.6lf") { 
			static char buf[32];  sprintf(buf, (a < 1e-4 ? "%.2le" : f), (double) a); return((string) buf); };
		inline void accum( const posteriordist &p, double w=1.0) { 
			cSelNodiv+=p.cSelNodiv*w; cSelDiv+=p.cSelDiv*w; cSelLow+=p.cSelLow*w; cNeutNodiv+=p.cNeutNodiv*w; cNeutDiv+=p.cNeutDiv*w; cNeutLow+=p.cNeutLow*w; cNeutHigh+=p.cNeutHigh*w; cHigh+=p.cHigh*w; cCnt+=p.cCnt*w; }
	};

	bool initialized = false;
public:
	modelInsightSup( const modelInsight &M) : modelInsight( M ) { initialized = true; };
	// void init( const modelInsight &Model ) { model = Model.getparameters();  };
	void clear() { initialized = false; modelInsight::clear(); };

	// Set byBlock = true when applying to neutral polymorphism blocks, and false for weignted averages over data positions...
	static bool estimateLambdaTheta(const BlockSet_t &bs, double priorWeight, const modelInsight *Prior, double &lambda, double &theta, bool byBlock);

	bool hessian(const BlockSet_t &data, ezmatrix &Hess, const accel_t *accel = NULL, double PriorWeight = 0, const modelInsight *Prior = NULL, double delta = 1e-6) const;
	bool estimateErr(const BlockSet_t &data, parameters &Errs, const accel_t *accel = NULL, double PriorWeight = 0.0 , const modelInsight *Prior = NULL, ezmatrix *InvHess=NULL, double delta = 1e-6) const;
	void insight1Posterior(ulong Nallele, SelectionClass Y, double ZeqXmaj, double ZeqXmin, posteriordist &Post);
	void insight1PosteriorSum(const BlockSet_t &data, posteriordist &Post);
	// void insight1PosteriorDump(const BlockSet_t &data, const string &Fout);
	bool suppStats(const BlockSet_t &data, parameters &Errs, derived_t &deriv, derived_t &derivErrs, const accel_t *accel = NULL, double PriorWeight = 0, const modelInsight *Prior = NULL, ezmatrix *InvHess = NULL, double delta = 1e-6) const;
	bool suppStatsSimple(ulong NumAlleles, derived_t &deriv) const;
};

// calculate the second deravitive matrix using finite differences of analytic first deravitive...
bool modelInsightSup::hessian(const BlockSet_t &data, ezmatrix &Hess, const accel_t *accel, double PriorWeight, const modelInsight *Prior , double delta) const {
	Hess.clear();	if (!initialized) return false;

	const accel_t *acc = accel;	// allocate a likelihood/AP accelerator, if none is provided..
	if (acc == NULL || !(acc->initialized)) {
		accel_t *a = new accel_t( *this ); a->initialize( data, PriorWeight, Prior ); acc = a;
	};

	parameters t_model = model;	// make a copy of parameters, so this function is "const"

	Hess.resize(3,3);				// square matrix, rho, eta,gamma
	auto deriv = [&]( double *p, ezmatrix &m, int m_irow, bool max1 ) {									// numerical deravitive... finite differences of analytic first deriv...
		double s=*p,x1=*p-delta,x2=*p+delta; parameters y1,y2; if (x1<0) x1=0; if (max1 && x2>1) x2=1;	// get acceptable finite differences
		*p=x1; acc->calcNLL( &y1, &t_model ); *p = x2; acc->calcNLL(&y2, &t_model); *p=s;				// calculate first deravitives, NOTE: p points to an element of t_model....
		m.v(m_irow, 0) = (y2.rho - y1.rho ) / (x2 - x1); m.v(m_irow, 1) = (y2.eta - y1.eta ) / (x2 - x1); m.v(m_irow, 2) = (y2.gam - y1.gam ) / (x2 - x1);
	};
	deriv(&t_model.rho, Hess, 0, true );
	deriv(&t_model.eta, Hess, 1, false);
	deriv(&t_model.gam, Hess, 2, false);
	if (accel == NULL) { delete acc; acc = NULL; };
	return true;
}

// Get paramater error estimates using the hessian... optional parameters provide access to accelerated MAP (accel) or regular MAP (Prior/Weight).
//	InvHess Optionally returns the inverse of the hessian matris, which is useful for estimating a variety of confidence intervals.
bool modelInsightSup::estimateErr(const BlockSet_t &data, parameters &Errs, const accel_t *accel, double PriorWeight, const modelInsight *Prior, ezmatrix *InvHess, double delta ) const {
	ezmatrix a, ainv;
	bool ok = hessian(data, a, accel , PriorWeight, Prior , delta );
	if (!ok) return(ok);
	a.makeSym(); a.inv(ainv); // stabalize numerical error, should be arithmetically symmertric, disable for debugging...
	int i = 0; Errs.rho = sqrt(ainv.v(i, i)); i = 1; Errs.eta = sqrt(ainv.v(i, i)); i = 2; Errs.gam = sqrt(ainv.v(i, i));
	if (false) { // debugging
		ezmatrix m = a;
		std::cerr << "Inverse Hessian is  " << std::endl;
		std::cerr << ainv.toStr() << std::endl;
		a.mult(ainv, m);
		std::cerr << "Is Normal?: " << std::endl;
		std::cerr << m.toStr() << std::endl;
	};
	if (InvHess != NULL ) *InvHess = ainv;
	return ok;
}

// Get mean / weightes lambda and theta across blocks...
bool modelInsightSup::estimateLambdaTheta(const BlockSet_t &bs, double priorWeight, const modelInsight *Prior, double &lambda, double &theta, bool byBlock)  {
	typedef inscomp::Block_t::Coord Coord;

	Coord st, en;
	double nsite = 0.0, ncount = 0.0, cum_lam = 0.0, cum_thet = 0.0, tmp_lam = 0.0, tmp_thet = 0.0;
	for (ulong iblock = 0; iblock < bs.numBlocks(); iblock++) {
		auto &b = bs.blocks[iblock];
		b.getHeader(NULL, &st, &en, &tmp_lam, &tmp_thet);
		nsite = (double)(byBlock ? en - st : b.NumSites());	// weight by blocksize (neutrals) or number of informative sites
		ncount += nsite; cum_lam += nsite * tmp_lam; cum_thet += nsite * tmp_thet;
	};
	ncount *= bs.meanWeight; cum_lam *= bs.meanWeight; cum_thet *= bs.meanWeight;	// downweigh actual data according to mixture coefficent.

	// add priors, as pseudocounts, if requested.
	if ( (priorWeight > 0) && (Prior != NULL)) {
		const modelInsight::parameters p = Prior->getparameters();
		cum_lam		+= priorWeight * (p.block.lambdaT > 0 ? p.block.lambdaT : (ncount > 0 ? cum_lam  / ncount : 0));
		cum_thet	+= priorWeight * (p.block.theta   > 0 ? p.block.theta	: (ncount > 0 ? cum_thet / ncount : 0));
		ncount		+= priorWeight;
	}

	if (ncount <= 0) return false;
	lambda = cum_lam / ncount; theta = cum_thet / ncount;
	return true;
}
// Just calculate supp stats, but not errors, good for converting expected param values to estimates of expected derived values.
bool modelInsightSup::suppStatsSimple(ulong NumAlleles, derived_t &deriv) const {
	// Convienent aliases
	double rho = model.rho, eta = model.eta, gamma = model.gam, theta = model.block.theta*wattersonsA(NumAlleles), lambda = model.block.lambdaT;
	double lambdatheta = lambda * theta;
	if (!initialized) return(false);

	// Dp (sites under positive selection)
	double Dp = rho * eta * lambda;
	// Pw (sites udner weak selection)
	double Pw = rho * gamma * (theta - eta * lambdatheta);
	// Alpha
	double alpha = 0.0;
	if (rho*eta != 0) alpha = 1.0 / (1.0 + (1.0 - rho) / (rho*eta));
	// Tau
	double tau = 0.0;
	if (rho*eta != 0) tau = 1.0 / (1.0 + (1.0 - rho) / (rho*gamma));
	// return values
	deriv.dp = Dp * 1000;	deriv.pw = Pw * 1000;		deriv.alpha = alpha;		deriv.tau = tau;
	return true;
}
	// Calculate supplimentary stats
bool modelInsightSup::suppStats(const BlockSet_t &data, parameters &Errs, derived_t &deriv, derived_t &derivErrs, const accel_t *accel, double PriorWeight, const modelInsight *Prior, ezmatrix *InvHess, double delta) const {
	// Convienent aliases
	using butils::ezmatrix;
	//const BlockSet_t &bs = data;		
	double rho = model.rho, eta=model.eta, gamma = model.gam, theta=model.block.theta*wattersonsA(data.numAlleles), lambda = model.block.lambdaT;
	double lambdatheta = lambda * theta;

	if (!initialized) return( false );

	ezmatrix invhess;
	estimateErr( data, Errs, accel, PriorWeight, Prior, &invhess, delta );
	if (false)  std::cerr << "\nInverse Covariance Matrix:\n" << invhess.toStr() << std::endl;
	invhess.makeSym();
	if (InvHess != NULL ) *InvHess = invhess;

	// Dp (sites under positive selection)
	double eDp, Dp = rho * eta * lambda;
	if (rho == 0 || eta == 0) { Dp = 0; eDp = -1.0; } else {
		eDp = invhess.v(0, 0) / (rho*rho) + invhess.v(1, 1) / (eta*eta) + 2 * invhess.v(0, 1) / (eta*rho);
		eDp = rho*eta*lambda * sqrt(eDp);
	}
	// Pw (sites udner weak selection)
	double a, b, c;	// temp values
	double ePw, Pw = rho * gamma * (theta - eta * lambdatheta);
	if (rho == 0 || gamma == 0 || (theta - eta * lambdatheta) == 0) { Pw = 0; ePw = -1.0; } else {
		a = 1.0 / rho;  c = 1.0 / gamma;  b = -lambdatheta / (theta - eta * lambdatheta);
		ePw = a*a*invhess.v(0, 0) + b*b*invhess.v(1, 1) + c * c * invhess.v(2, 2);
		ePw += 2.0 * a * b *invhess.v(0, 1) + 2.0*a*invhess.v(0, 2) + 2.0 * b * c * invhess.v(1, 2);
		ePw = Pw * sqrt(ePw);
	}
	// Alpha
	double ealpha, alpha;
	if (rho*eta == 0) { alpha = 0.0; ealpha = -1.0; } else {
		alpha = 1.0 / (1.0 + (1.0 - rho) / (rho*eta));
		a = 1.0 / (rho * rho * eta); b = (1.0 - rho) / (rho * eta * eta);
		ealpha = a*a*invhess.v(0, 0) + b*b*invhess.v(1, 1) + 2.0 * a*b*invhess.v(0, 1);
		ealpha = alpha * alpha * sqrt(ealpha);
	}
	// Tau
	double etau, tau;
	if (rho*gamma == 0) { tau = 0.0; etau = -1.0; } else {
		tau = 1.0 / (1.0 + (1.0 - rho) / (rho*gamma));
		a = 1.0 / (rho * rho * gamma); b = (1.0 - rho) / (rho * gamma * gamma);
		etau = a * a * invhess.v(0, 0) + b * b * invhess.v(2, 2) + 2.0 * a * b * invhess.v(0, 2);
		etau = tau * tau * sqrt(etau);
	}
	// return values
	deriv.dp		= Dp * 1000;	deriv.pw = Pw * 1000;		deriv.alpha = alpha;		deriv.tau = tau;
	derivErrs.dp	= eDp * 1000;	derivErrs.pw= ePw * 1000;	derivErrs.alpha= ealpha;	derivErrs.tau = etau;
	return true;
}

void modelInsightSup::insight1Posterior(ulong Nallele, SelectionClass Y, double ZeqXmaj, double ZeqXmin, posteriordist &Post) {
	// The provide probabilties of UNOBSERVED (latent) variables (considered S & A). Yes, I know that the allele distribution over A is not explicitely 
	//	treated, and this is a little confusing (legacy). For efficency, we calculate the distribution over allowable S,A combinations for input
	//  observations. Input observation consists of the sufficient characterization of Selection Class Y {M,L,H} and the prob that Z=Xmaj & Z=Xmin.
	// informative alias for parameters
	double rho = model.rho, eta = model.eta, gamma = model.gam;
	double lambda = model.block.lambdaT, theta = model.block.theta * wattersonsA(Nallele);
	double beta1 = model.beta.b1, beta3 = model.beta.b3;
	double pzxmaj = ZeqXmaj, pzxmin = ZeqXmin;
	Post.clear();
	switch (Y) {
	case SelectionClass::Mono:
		Post.cNeutNodiv = (1 - rho)*(1 - theta)*pzxmaj*(1 - lambda);		// 3) S == 0, A == Z | Y=M,  (inf sites -> Z==Xma==A) - See Supp table S1.
		Post.cNeutDiv = (1 - rho)*(1 - theta)*(1 - pzxmaj)*lambda / 3;		// 4) S == 0, A != Z | Y=M,  (inf sites -> Z!=Xma==A)
		Post.cSelNodiv = rho*(1 - gamma*theta)*pzxmaj*(1 - eta*lambda);		// 0) S == 1, A == Z | Y=M,  (inf sites -> Z==Xma==A)
		Post.cSelDiv = rho*(1 - pzxmaj)*eta*lambda / 3;						// 1) S == 1, A != Z | Y=M,  (inf sites -> Z!=Xma==A)
		break;
	case SelectionClass::PolyL:
		Post.cSelLow = rho * gamma*theta / 3 * pzxmaj * (1 - eta*lambda);									// 2) S == 1, A == Xma | Y=L, Z==Xma==A
		Post.cNeutLow = (1 - rho) * beta1 * theta / 3 * ((1 - pzxmaj)*lambda / 3 + pzxmaj*(1 - lambda));	// 5) S == 0, A == Xma | Y=L, (Z==Xma==A)+(Z!=Xma==A) Low  Derived allele Freq 
		Post.cNeutHigh = (1 - rho) * beta3 * theta / 3 * ((1 - pzxmin)*lambda / 3 + pzxmin*(1 - lambda));	// 6) S == 0, A == Xmi | Y=L, (Z==Xmi==A)+(Z!=Xmi==A) High Derived allele Freq 				ptot      = pSelLow + pNeutLow + pNeutHigh;
		break;
	case SelectionClass::PolyH:
		Post.cHigh = 1.0;	break;											// 7) S=0, A=Xma or A=Xmi | Y = H (all HF pollys are Neut, by def)
    default: break;
	};
	Post.renorm(); Post.cCnt = 1.0;
}

void modelInsightSup::insight1PosteriorSum(const BlockSet_t &data, posteriordist &Post) {
	Post.clear();
	if (!initialized) return;

	parameters::block_t s_block = model.block;	// save these values
	posteriordist t_post, s_post;				// temporary and cumulative posterior counts...
	ulong nallele = data.numAlleles;

	inscomp::PatternCounts_t	nobs;
	inscomp::AncPriIdx_t		izxmin = 0, izxmaj = 0;
	double zxmaj = 0, zxmin = 0;
	SelectionClass sc;
	ulong nblocks = (ulong)data.numBlocks();
	for (ulong iblock = 0; iblock <= nblocks; iblock++) {
		auto &bl = data.blocks[iblock]; bl.getHeader(NULL, NULL, NULL, &(model.block.lambdaT), &(model.block.theta));
		sc = SelectionClass::Mono;
		for (ulong ipos = 0; ipos <= bl.monoNumPat(); ipos++) {
			bl.monoVal(ipos, izxmaj, nobs); zxmaj = *(bl.ancpri[izxmaj]); insight1Posterior(nallele, sc, zxmaj, 0, t_post); s_post.accum(t_post);
		};
		sc = SelectionClass::PolyL;
		for (ulong ipos = 0; ipos <= bl.polyLNumPat(); ipos++) {
			bl.polyLVal(ipos, izxmaj, izxmin, nobs); zxmaj = *(bl.ancpri[izxmaj]); zxmin = *(bl.ancpri[izxmin]); insight1Posterior(nallele, sc, zxmaj, zxmin, t_post); s_post.accum(t_post);
		};
		sc = SelectionClass::PolyH;
		for (ulong ipos = 0; ipos <= bl.polyHNumPat(); ipos++) {
			bl.polyHVal(ipos, izxmaj, izxmin, nobs); zxmaj = *(bl.ancpri[izxmaj]); zxmin = *(bl.ancpri[izxmin]); insight1Posterior(nallele, sc, zxmaj, zxmin, t_post); s_post.accum(t_post);
		};
	};
	model.block = s_block;	// recover saved parameters...
	return;
}
