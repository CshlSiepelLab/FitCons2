#pragma once

#include <vector>
#include <functional>	// greater template
#include <algorithm>	// sort

class mathplus {
public:
	static double normLL(std::vector<double> &lls) {	// reacale logprobs so the sum of the underlying probs is 1.
		double sumLL = sumOfLL( lls );
		//double ret = sumLL;
		//double sump = 0; for (long i=0; i<lls.size(); i++ ) sump += exp( lls[i] );
		//std::cerr << "NormLL prenorm - sumLog: " << sumLL << "  sumP: " << sump << "  log(sumP): " << mylog(sump) << std::endl;
		for (unsigned long i = 0; i<lls.size(); i++) lls[i] -= sumLL;
		//sumLL = sumOfLL(lls);
		//sump = 0; for (long i = 0; i<lls.size(); i++) sump += exp(lls[i]);
		//std::cerr << "NormLL postnorm - sumLog: " << sumLL << "  sumP: " << sump << "  log(sumP): " << mylog(sump) << std::endl;
		return( sumLL );
	}

	static double sumOfLL(const std::vector<double> &lls) {
		if (lls.size() < 1) return 0;
		std::vector<double> tmp = lls;
		std::sort(tmp.begin(), tmp.end(), std::greater<double>());			// sorted in descending order, improved stability
		double sum = 0.0, vmax = tmp[0];									// divide all by largest prob (largest LL)
		for (auto i = tmp.size(); i>0 ; --i) sum += exp(tmp[i-1] - vmax);	// cumulate from smallest to largest, this is the most stable...
		return(mylog(sum) + vmax);
	};
	static double sumOfNLL(const std::vector<double> &nlls) {
		if (nlls.size() < 1) return 0;
		std::vector<double> tmp = nlls;
		std::sort(tmp.begin(), tmp.end(), std::greater<double>());			// sorted in descending order, improved stability
		double sum = 0.0, vmin = tmp[tmp.size() - 1];						// divide all by largest prob (largest LL, i.e. smallest NLL)
		for (unsigned long i = 0; i < tmp.size(); ++i) sum += exp(-(tmp[i] - vmin));
		return(-log(sum) + vmin);
	};
	static double rat2frac(double r) { return (r <=0?0.0 : r / (1.0 + r)); };														// convert ratio a/b to fraction a/(a+b)	 [0-inf)
	static double frac2rat(double f)  { if (f < 0) f = 0; return (f >= 1 ? std::numeric_limits<double>::max() : f / (1.0 - f)); };	// convert fraction a/(a+b) to ratio a/b, if possible [0-1]
	inline static double mylog(double v) { return (v > .5 ? log1p(v - 1.0) : log(v)); };
	static constexpr double natToBit = 1.44269504088896;	// multiplicative factor that converts natural log to base 2 log (nats to bits)
	static constexpr double bitToNat = 1.0 / natToBit;		// multiplicative factor that converts base 2 log to natural log (bits to nats)
	static double entropy( double rho ) { if ((rho<=0.0) ||  (rho>=1.0)) return 0.0; return -rho*mylog(rho); }		// entropy for single variate w prob rho, nats.
	static double entropyB(double rho ) { return entropy(rho)+entropy(1.0-rho); }									// entropy for binary variable with Prob rho, nats.
};