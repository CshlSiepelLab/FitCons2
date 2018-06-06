#pragma once
#pragma once

// This class estiamtes the expected conditional information change when a Parent set of binary measurements (genomic positions)
//	is partitioned into two Gaussian child classes. Rho is the fraction of positions that are positive, the rrest are negavibe.
//	Information is about psotiivity. We assume that the total number of positives before and after partitioning is constant.
// We are give a centrality estimate for each child distribution (rho), as well an error (stdev) for the estimate, and the number of 
//	sites in the distriution (N). The parrent (unparititoned) distribution is taken as the union of the child distributions, and the 
//	measurements of the child distributiosn are taken as independent. We sampel a grid of NSamp points, over all admissible values
//	in each child, estiamte the parrental information and the difference, and then sum over all samples to get an expectation.
// Each sampel covers the SAME amnount of probability space and the admissible interval is limited to the intersection of 0<RHO<1, and 
//	RHO+/-(SampRange*StDevRho). Generally SampRange is set to 3 to cover > 99% of curve (+/- 3 stdev). 100=NSamp is probby too few
//	(10,000 samples) while 1000 is probably overkill. (1,000,000 samples). 300 is probably good. Rectuangular intervals are sampled at center.


#include <cstdint>		// uint32_t
#include <cassert>		// assert

#include <vector>		// - 
#include <algorithm>    // std::max

class expectedInfo {
private:
	typedef std::string string;
public:
	struct options {
		bool orderedEst;	// Only count as informative probbility masses that maintain the central variable ordering of p1>p2 
		void clear () { orderedEst=false; }
		options() { clear(); }
	};
	struct paramSet {
		double rho;
		double sigma;
		double n;
		double ns;
		paramSet() { clear(); }
		void clear() { rho = 0; sigma = -1; n = 0; ns = 0; }
		string fmt(double a, int w = 15, int p = 0) { char b[128]; sprintf(b, "%*.*lf", (int)w, (int)p, (double)a); return string(b); }
		string toStr() { return fmt(rho, 8, 5) + "  " + fmt(sigma, 8, 5) + "  " + fmt(n) + "  " + fmt(ns); }
		void set(double N, double Rho, double Sigma = -1.0) {
			n = N; rho = Rho; sigma = (Sigma <= 0 ? 3.0 : Sigma); ns = n*rho;
		}
	};
	struct resultSet {
		double expInf;
		double expRho;
		double pointParent;
		double expInfPart;
		double expRhoPart;
		double probPart;
		double stressN;		// should be 0, are more relevant under parental decompositio nstrategy, irrelevant udner child composition strategy.
		double stressS;		// SUS( ML_C1 ) + SUS( ML_C2 ) - SUS( ML_P ).... ideally this is small (near 0) but if ML values diverge from expectation, so can this term.
		double stressSR;	// (SUS( ML_C1 ) + SUS( ML_C2 ) - SUS( ML_P ))/ SUS(ML_P).... ideally this is small (near 0) but if ML values diverge from expectation, so can this term.
		double probTest;
		bool isSet() const { return !( expInf==0 && expRho == 0 && pointParent == 0 && probTest == 0 && stressS == 0); }
		void clear() { expInf = expRho = expInfPart = expRhoPart = pointParent = stressN = stressS = stressSR = probTest = probPart = 0.0; }
		resultSet() { clear(); }
		~resultSet() {};
		string fmt( double a, int w = 15, int p = 0) { char b[128]; sprintf(b, "%*.*lf", (int)w, (int)p, (double)a); return string(b); }
		string fmtE(double a, int w = 15, int p = 0) { char b[128]; sprintf(b, "%*.*le", (int)w, (int)p, (double)a); return string(b); }
		string toStr() { return fmt(expInf) + "  " + fmt(expRho, 8, 5) + "  " + fmt(expInfPart) + "  " + fmt(expRhoPart, 8, 5) + "  " + fmt(pointParent) + "  " + fmt(stressS) + "  " + fmt(stressSR, 8, 5) + "  " + fmt(stressN) + "  " + fmtE(probPart, 10, 2) + fmtE(probTest-1.0,10,2); }
	};
private:
	// paramSet p, c1, c2;
	struct probPr { double prob; double val; };
	struct probTr { double prob1; double prob2; double val; };
	typedef std::vector< double > vecD;
	typedef std::vector< probPr > vecP;
	typedef std::vector< probTr > vecT;
	class marginal : public vecP {
	public:
		static inline double bell(double x, double sig) { double d = (sig == 0.0 ? (x == 0.0 ? 1.0 : 0.0) : std::exp(-(x*x) / (sig*sig))); return d; }
		bool setMarg(const paramSet &Param, const uint32_t NSamp, const double SampRange) {
			if (NSamp < 1) return false;
			if (SampRange <= 0.0 ) return false;
			double rho_min, rho_max, rho_del, rho, norm;
			rho_min = std::max(0.0, Param.rho - SampRange*Param.sigma);
			rho_max = std::min(1.0, Param.rho + SampRange*Param.sigma);
			rho_del = (rho_max - rho_min) / NSamp;
			rho = rho_min + 0.5 * rho_del;	norm = 0.0;
			resize(NSamp);
			for (uint32_t i = 0; i < NSamp; i++) {
				(*this)[i].val = rho; (*this)[i].prob = bell(rho - Param.rho, Param.sigma);
				norm += (*this)[i].prob; rho += rho_del;
			}
			if (norm <=0) return false;
			norm = 1.0 / norm;
			for (auto & s : (*this)) s.prob *= norm;
			return true;
		}
	};
	static inline double Entrop(double P) { if (P <= 0.0 || P >= 1.0) return 0.0; return -P*log2(P); }
	static inline double Info(double N, double P) { if (N <= 0) return 0.0; return N * Entrop(P); }
	static inline double InfoBin(double N, double P) { if (N <= 0) return 0.0; return N * (Entrop(P) + Entrop(1.0 - P)); }
	static inline double deltaInfoP(double Pn, double Pp, double C1n, double C1p)  {
		double c2n = Pn - C1n; if (c2n <= 0) return 0.0;
		double c2sus = Pn*Pp - C1n*C1p; if (c2sus <  0) return 0.0;
		return InfoBin(Pn, Pp) - (InfoBin(C1n, C1p) + InfoBin(c2n, c2sus / c2n));
	}
	static inline double deltaInfoC(double C1n, double C1p, double C2n, double C2p, const double *Pp=NULL)  {
		double pn = C1n + C2n; if (pn <= 0) return 0.0;
		double prho = (C1n*C1p + C2n*C2p) / pn;					// Local consistency
		if (Pp!=NULL && *Pp > 0.0 && *Pp < 1.0 ) prho = *Pp;	// Override with global consistency, when available.
		return InfoBin(pn, prho) - (InfoBin(C1n, C1p) + InfoBin(C2n, C2p));
	}

public:
	expectedInfo() {  }
	~expectedInfo() { }


	static bool estimate(const paramSet &P, const paramSet &IC1, const paramSet &IC2, uint32_t NSamp, double SampRange, resultSet &res, const options *o=NULL )  {
		options opts; if (o!=NULL) opts = *o;	// use default options, if none provided....
		paramSet C1=IC1, C2=IC2;
		if (C1.rho > C2.rho ) std::swap( C1, C2 );	// Now C1 is the child with the LOWER value of rho.
		res.clear();
		res.stressN = (C1.n + C2.n) - P.n; res.stressS = (C1.ns + C2.ns) - P.ns;
		res.stressSR = res.stressS / P.ns; 	// However, it might be a good idea to monitor these stressors.... StressN shoudl be 0.
		marginal marg_c1, marg_c2;

		// sample proabilities of normally distributed values fro RHO centeres at C.p with stdev C.sigma and C.n observations.
		marg_c1.setMarg(C1, NSamp, SampRange);	
		marg_c2.setMarg(C2, NSamp, SampRange);

		double cum_inf = 0.0, cum_rho = 0.0, cum_prob = 0.0, cum_inf_part=0.0, cum_rho_part=0.0, cum_prob_part=0.0;
		for (uint32_t i = 0; i < NSamp; i++) {
			const auto &mc1 = marg_c1[i];
			for (uint32_t j = 0; j < NSamp; j++) {
				const auto &mc2=marg_c2[j];
				double prob = mc1.prob*mc2.prob;
				double einf = deltaInfoC(C1.n, mc1.val, C2.n, mc2.val, &(P.rho));
				double erho = ((C1.n * mc1.val + C2.n * mc2.val) / (C2.n + C1.n));
				cum_inf   += prob * einf;
				cum_rho   += prob * erho;
				cum_prob  += prob;
				if (mc1.val < mc2.val) {	// Only include samples whose central value is properly oriented.
					cum_inf_part  += prob * einf;
					cum_rho_part  += prob * erho;
					cum_prob_part += prob;
				}
			}
		}

		res.expInf		= cum_inf/cum_prob;					// essentially dividign by ~1.0
		res.expRho		= cum_rho/cum_prob;					// dividing by ~1.0...
		res.probTest	= cum_prob;							// should be very close to 1.0.
		// renormalize for accumulation over samples that meet admissability cosntraints.
		res.expInfPart	= cum_inf_part;						// dont renormalize, inadmissible samples are uninformative...
		res.expRhoPart  = cum_rho_part / cum_prob_part;		// do renorm to get estiamted parrental rho over admissable samples...
		res.probPart    = cum_prob_part;

		// Generate a point estiamte from childeren 
		res.pointParent = deltaInfoC(C1.n, C1.rho, C2.n, C2.rho, &(P.rho));
		return true;
	}
};
