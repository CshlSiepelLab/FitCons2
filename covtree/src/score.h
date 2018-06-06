#pragma once
// Classes used to encapsulate processing of covedefs , data and bitfields
//
// scoreSingelton - socres for a scoresets that subset a single covariate. This is used in two ways
//		1) To conditionally order nonmonotonic covaraites (Singleton is a subset consisting of a single value)
//		2) 1/2 of a partition of a covSet by a single covaraite. (Singleton is subset of covaraite values ans is paired with a second Singleton with the remaining values)
//
// scorePair - Two singletons, a nonoverlapping set and a coset that together span all the  parrentialcovariate values fro one covariate. All other
//		covariates inherit the parent's range of values.

// System includes
#include <cstdint>		// uint8_t ... etc
#include <map>
#include <vector>
#include <algorithm>    // std::min

// butils includes
#include "fastfile.h"
#include "mathplus.h"

//Local project includes
#include "covtreetypes.h"
#include "FmtSimple.h"
#include "cov.h"
#include "expectedInfo.h"
#include "covtreeArgs.h"

// A single covariateSet Subset, not a partition. Used both to get preliminary rho for ordering of individual nonmonotonic covaraite values
//	ans ALSO as as one of 2 parts in a binary partition of a covSet (with multiple cov values per cov)
class scoreSingelton : private FmtSimple {
public:
	typedef covDB_t::bitfieldDef_t BitEncoder_t;
	typedef std::vector<scoreSingelton> scoreSingelton_v;
private:
	const covsetDefSubset	*covSetMaster;				// Parrent Covariate set, contaisn pointer to Universal definition
	const BitEncoder_t		*covEncoderAnn, *covEncoderCTS;	// Encoder, converts covarate & value t bit field masks
	string covName, covValTag, fBaseName;				// String names for covariate, value and base for fileHandeling
	covsetDefSubset		covset;							// New Covset with a singleton value in the target Covarate
	CovBits_t			bitmaskCTS, bitmaskAnn;			// The daabase bitmask
														// Fields used to genereate Bed File
	string fnScan;
	fastfile::FILE *fscan;
	string scan_chr;
	LinAddr_t scan_st, scan_en;
	uint32_t scan_count;
	//
	uint64_t scan_count_sum, scan_count_sum_w;
	//
	insightShell::insightResults	scan_res;
	// Results String... / Set
	bool scanFlush() {
		// TODO check fore errors...
		if (!fscan) return false;
		if (scan_count != 0) {
			fastfile::putString(fscan, scan_chr.c_str()); fastfile::putString(fscan, "\t");
			fastfile::putInt(fscan, scan_st);     fastfile::putString(fscan, "\t");
			fastfile::putInt(fscan, scan_en);    fastfile::putString(fscan, "\t");
			fastfile::putInt(fscan, scan_count); fastfile::putString(fscan, "\n");
		}
		scan_chr = ""; scan_st = scan_en = scan_count = 0;
		return true;
	};
public:
	void clear() { covSetMaster =NULL; covEncoderAnn=NULL; covEncoderCTS=NULL; fscan = NULL; covName = covValTag=fBaseName=fnScan=scan_chr=""; covset.clear(); bitmaskAnn.clear(); bitmaskCTS.clear();scan_st = scan_en = 0; scan_count=0; scan_res.clear(); scan_count_sum = scan_count_sum_w=0; return;	}
	inline const string &CovName() const { return(covName); };
	inline const string &CovValTag() const { return(covValTag); };
	inline const covsetDefSubset &covsetDef() const { return covset; }
	inline uint64_t scanCount()  const { return scan_count_sum;}
	inline uint64_t scanCountW() const { return scan_count_sum_w; }
	inline bool  isCTS() const { return(covSetMaster->getCovByName(covName)->isCTS()); };
	inline bool  CanBePartitioned() const { return covset.CanBePartitioned(); };
	bool Init(const covsetDefSubset &MasterCovSet, const BitEncoder_t &EncoderAnn, const BitEncoder_t &EncoderCTS, const string &CovName, const string &CovValTag, const string &FileBaseName) {
		clear();
		// Stash the parameters for later use
		covSetMaster = &MasterCovSet;
		covEncoderAnn = &EncoderAnn; covEncoderCTS = &EncoderCTS;
		covName = CovName; covValTag = CovValTag; fBaseName = FileBaseName;

		// Getch the covaraiteSet from the master, and make a copy but replace one of the covariates with a signelton (jsut one of the valid values).
		const covDefSubset	*covOldP = MasterCovSet.getCovByName(CovName);
		if (!covOldP) return false;
		covDefSubset covNew; if (!covNew.Init(*covOldP, CovValTag)) return false;
		if (!covset.InitAndReplaceCov(MasterCovSet, covNew)) return false;

		// Now get the bitfield masks for the new covaraite set, with the singelton covariate value...
		for (int b = 0; b<2; b++) {
			bool iscts = (b == 1);
			CovBits_t &bitmask = (iscts ? bitmaskCTS : bitmaskAnn); bitmask.clear();
			const BitEncoder_t &Encoder = (iscts ? EncoderCTS : EncoderAnn);
			for (uint8_t i = 0; i < covset.getNumCovIDS(iscts); i++) {
				const covDefSubset *covdefsub = covset.getCovByID(iscts, i);
				for (uint8_t j = 0; j < covdefsub->numVals(); j++) {
					if (!Encoder.accumulateMask(i, covdefsub->vals[j]->value(), bitmask)) return false;
				};
			};
		};
		// TODO: Anything special with file name?
		// TODO: Anything special with resutls set?
		return true;
	};

	void const getMasks(CovBits_t &bmAnn, CovBits_t &bmCTS) const { bmAnn = bitmaskAnn; bmCTS = bitmaskCTS; return; };
	const string &getFname() const { return fnScan; };

	// Open file, clear temp variables
	bool scanInit( bool debug = false ) {
		if (fnScan != NULL) scanTerminate();
		{	// We should never have commas in the taglist, as it is just one tag, but sometimes a derived class will provide a list, so remove them, JIC.
			string covtagliststr = covValTag;	covtagliststr.erase(std::remove(covtagliststr.begin(), covtagliststr.end(), ','), covtagliststr.end());
			fnScan = fBaseName + "/" + covName + "-" + covtagliststr + ".bed";
		}
		scan_chr = ""; scan_st = scan_en = scan_count = 0; scan_count_sum = scan_count_sum_w = 0;
		if (!debug) {
			fscan = fastfile::ffopen(fnScan.c_str(), false);
			if (fscan == NULL) return false;
		};
		return true;
	}
	// Flush loci, clear temps, close file...
	bool scanTerminate() {
		if (fscan == NULL) return true;
		scanFlush(); fastfile::ffclose(&fscan);
		return true;
	};
	// bedfile format, 0 based, half open [start,end) 
	bool scanAddLocus(const string &chrom, uint32_t St, uint32_t En, uint32_t Count) {
		if (fscan == NULL) return false;
		// Just extend the current locus
		if (Count>0) {
			scan_count_sum += (En - St); scan_count_sum_w += (En - St)*Count; }
		if ((chrom == scan_chr) && (scan_en == St) && (Count == scan_count)) { scan_en = En; return true; }
		scanFlush();
		scan_chr = chrom; scan_st = St; scan_en = En; scan_count = Count;
		return true;
	}

	// Fetch Insight Results
	// insightString(); // generate Insight string ....
	bool resultsRead() { return scan_res.read(fnScan); }
	bool resultVal(const string &Group, const string &Field, double &Value) const	{ return(scan_res.getField(Group, Field, Value)); };
	bool resultVal(const string &Group, const string &Field, string &Value) const	{ return(scan_res.getField(Group, Field, Value)); };
	bool resultSet(const string &Group, const string &Field, const double &Value)	{ return(scan_res.setField(Group, Field, Value)); };	// manually add a result.
	bool resultSet(const string &Group, const string &Field, const string &Value)	{ return(scan_res.setField(Group, Field, Value)); };	// manually add a result.
	// Shortcut for frequently used value....
	double resultValRho()  const { double a; if (!resultVal("PARAMETERS",  "Rho",   a)) return -1.0; return a; };
	double resultValRhoU() const { string s; double a;	// Max observed stderr in test run is 1.57 (0<rho<1). Sometimes large stderr resunls in reported NAN values. Check for them.
		if (!resultVal("UNCERTAINTY", "Rho", s)) return -1.0; 
		if ( s.find_first_of("nan")!= string::npos) return -1.0;
		if (!resultVal("UNCERTAINTY", "Rho",   a)) return -1.0;
		return a; };
	double resultValEta()  const { double a; if (!resultVal("PARAMETERS",  "Eta",   a)) return -1.0; return a; };
	double resultValGam()  const { double a; if (!resultVal("PARAMETERS",  "Gamma", a)) return -1.0; return a; };
	double resultValNLL( bool WithPrior=false )  const { double a; if (!resultVal("STATS", (WithPrior?"NLPrior":"DataNLL"), a)) return -1.0; return a; };

	string toStrParams() const {
		string s = ""; double a;
		resultVal("PARAMETERS", "Rho", a);	s += fmt(a, 8, 6) + "\t";
		resultVal("UNCERTAINTY", "Rho", a);	s += fmt(a, 8, 6) + "\t";
		resultVal("PARAMETERS", "Dp", a); s += fmt(a, 8, 6) + "\t";
		resultVal("UNCERTAINTY", "Dp", a); s += fmt(a, 8, 6) + "\t";
		resultVal("PARAMETERS", "Pw", a); s += fmt(a, 8, 6) + "\t";
		resultVal("UNCERTAINTY", "Pw", a); s += fmt(a, 8, 6) + "\t";
		return s;
	}

	double resultValN(bool Inf = false )    const { 
		double b = 0.0, a = resultValNW( Inf ); 
		if ((a<0) || !resultVal("ARGS", "maxSampW", b) ) return -1.0; 
		if (b==0.0) return -1.0;
		return a / b;
	};
	double resultValNW(bool Inf=false)   const { double a=0.0; 
		if (!Inf) {
			if (!resultVal("STATS", "InPosW", a)) return -1.0; 
			return a; 
		}
		double t=0;
		if (!resultVal("STATS", "InfPosWH", t)) return -1.0; a += t;
		if (!resultVal("STATS", "InfPosWL", t)) return -1.0; a += t;
		if (!resultVal("STATS", "InfPosWM", t)) return -1.0; a += t;
		return a;
	}

	string toStrNsites(bool Inf = false) const {
		double a = resultValN(Inf);
		string s = fmt2(a, 13, true) + "\t";
		return s;
	}

	string toStrSummary(bool Inf = false) const {
		string s = ""; double a = 0;
		// resultVal("ARGS", "maxSampW", nCells);
		a = covset.getCovByName(CovName())->numVals();		s += fmt2(a, 2) + "\t";
		s += toStrNsites(true) + "\t";
		s += toStrNsites(false) + "\t";
		s += toStrParams() + "\t";
		s += fmt2( resultValNLL() * mathplus::natToBit, 13, true ) + "\t";
		s += fmt(  resultValNLL() * mathplus::natToBit/resultValN(true), 8, 6 ) + "\t";
		s += "|\t" + CovValTag();
		return s;
	}


	// Generate a list of covariate singletons. We score each and then order by scoring for binary partitioning.
	static bool GenerateSingeltons(const covsetDefSubset &MasterCovSet, const BitEncoder_t &EncAnn, const BitEncoder_t &EncCTS, const string &FileBase, scoreSingelton_v &SingCovSets) {
		string safeFileBase=FileBase + "/sing";
		fastfile::ffmkdir( safeFileBase.c_str());	// Try to make dir, returns 0 on success, but will fail if dir aready exists...
		SingCovSets.clear();
		for (int b = 0; b < 2; b++) {
			bool iscts = (b == 1);
			for (uint8_t icov = 0; icov<MasterCovSet.getNumCovIDS(iscts); icov++) {
				const covDefSubset * cov = MasterCovSet.getCovByID(iscts, icov); if (!cov) return false;
				// Ignore monotonic covariates, they to not need to be reordered at each step...
				if (!cov->univ->isMon()) {
					// If there is only 1 covaraite, we can't partition it, skip the pre-assersment
					if (cov->numVals() > 1) {
						for (uint8_t icoveV = 0; icoveV < cov->numVals(); icoveV++) {
							// Allocate & initialize the new singelton-replaces set
							SingCovSets.resize(SingCovSets.size() + 1);
							bool ok = SingCovSets[SingCovSets.size() - 1].Init(MasterCovSet, EncAnn, EncCTS, cov->Name(), cov->vals[icoveV]->tag(), safeFileBase);
							if (!ok) return false;
						}
					}
				}
			}
		}
		return true;
	};
};

// Represents a partition on a covsetSubset by partitionign a s ingle covariate from the parent. Implimented as 
//	a pair of scoreSingelton's 
class scorePair : private FmtSimple {
public:
	typedef scoreSingelton::BitEncoder_t					BitEncoder_t;
	typedef std::vector<double>								vecD;
	typedef std::vector<scorePair>							scorePair_v;
	//typedef std::pair<string,double>						covValWeight_t;
	typedef double											covValWeight_t;
	typedef std::vector<covValWeight_t> 					covValWeight_v;
	typedef std::map<string, covValWeight_v>				covsetCovOrder_t;
private:
	// OK this is slopy and wasteful, but I don't want to refactor scoreSingleton and this is not the main bottleneck, so do it simple.
	scoreSingelton				set, cos;	// will represent a partitioning of CovName into nonemoty, ordered subsets....
	vecD						covOrder;	// Ordering for covariate values of CovName
	expectedInfo::resultSet		infExp;		// expected information results
	double						cMinNpos;	// minimum child size
	double						cMaxPerr1;	// Max probability of ordering error
	double						cMaxPerr2;	// Max probability of distance error
	double						cMinDirInf;	// if >0 use expected directional DelInf rather than point DelInf, also report 0 unless > this ampount.
	double						cPerr1, cPerr2;	// Probability of an inversion (lower) or 2 fold distance error (higher probability)
	double						cParNLL;	// NLL of data from parent node, for calculation any nonnegative constant will do, but for reporting actual bits is useful
public:
	bool						qcOk;		// External agent sets this, true by default
	string						notes;		// processing notes for this split, for debugging...
	double pErr1()		{ return cPerr1; }
	double pErr2()		{ return cPerr2; }
	double maxpErr1()	{ return cMaxPerr1; }
	double maxpErr2()	{ return cMaxPerr2; }
	
	void clear() {	set.clear(); cos.clear(); covOrder.clear(); infExp.clear(); cMinNpos= cMinDirInf=0.0; cPerr1=cPerr2=cMaxPerr1 = cMaxPerr2 = 1.0; qcOk=true; cParNLL=-1.0; notes=""; return; }
	static void ParentalRho( double *RhoGet=NULL, double *RhoSet=NULL ) { 
		static double prho=0.0;
		if (RhoGet!=NULL) *RhoGet=prho; 
		if (RhoSet!=NULL) prho = *RhoSet; 
		if (prho<0.0) prho=0.0;
		if (prho>1.0) prho=1.0; 
		return; }
	inline const scoreSingelton &getPart(bool isSet) const { return (isSet ? set : cos); };
	inline bool  isCTS() const { return(set.isCTS()); };
	inline const string &CovName() const { return(set.CovName()); };
	inline const string &CovValTag(bool isSet) const { return((isSet ? set : cos).CovValTag()); };
	void const getMasks(CovBits_t &bmAnnSet, CovBits_t &bmCTSSet, CovBits_t &bmAnnCos, CovBits_t &bmCTSCos) const { set.getMasks(bmAnnSet, bmCTSSet); cos.getMasks(bmAnnCos, bmCTSCos); return; };
	const string &getFname(bool isSet) const { return (isSet ? set.getFname() : cos.getFname()); };

	bool scanInit(bool debug=false ) { return set.scanInit(debug) && cos.scanInit(debug); };
	bool scanTerminate() { return set.scanTerminate() && cos.scanTerminate(); };
	bool scanAddLocus(const string &chrom, uint32_t St, uint32_t En, uint32_t CountSet, uint32_t CountCos) {
		return set.scanAddLocus(chrom, St, En, CountSet) && cos.scanAddLocus(chrom, St, En, CountCos);
	}
	inline void scanGetSums( uint64_t &setsum, uint64_t &setsumW, uint64_t &cossum, uint64_t &cossumW)  const { 
		setsum=set.scanCount(); setsumW=set.scanCountW(); cossum=cos.scanCount(); cossumW=cos.scanCountW(); return; }
	inline bool scanValid() const { // did the scan return a partition with nonvoid set and nonvoid coset.
			return ( (set.scanCountW()!=0) && (cos.scanCountW()!=0)); }
	bool resultsRead() {
		bool a = set.resultsRead();  bool b = cos.resultsRead(); if (!a || !b) return false;
		// p1 - Probability that rho1 > rho2 | rho1_hat > rho2_hat , delta_rho1_hat, delta_rho2_hat (always > 50% for normal error). C++ has erf/erfc, now! Woo Hoo!
		// Gaussian convolution implies a normal distribution with variance as sum of convolved distributions and centrality at 0 (coincident centers)
		double a_rho = set.resultValRho(), a_rhou = set.resultValRhoU(); 
		double b_rho = cos.resultValRho(), b_rhou = cos.resultValRhoU();
		// make sure that a_rho<b_rho
		//if (a_rho > b_rho) { std::swap(a_rho,b_rho); std::swap(a_rhou, b_rhou);	} // actually this is unimportant...
		double sigma = std::sqrt(a_rhou*a_rhou + b_rhou * b_rhou);
		// The 0.5 multiplier on the outside is becasue erfc(1) returns 1.0 at x=0 and as high as 2.0 at -inf. In order to get prob mass to right of x, we need to divide by 2.
		cPerr1 = std::erfc(std::fabs(a_rho - b_rho) / sigma) * 0.5;				// OK this is the probability that rho2 > rho1 |  rho1_hat > rho2_hat  (probability of an ordering error)
		cPerr2 = std::erfc((std::fabs(a_rho - b_rho)*0.5) / sigma) * 0.5;		// OK this is the probability that ( rho1 - rho2 ) < (rho_hat1 - rho_hat2)/2 |  rho1_hat > rho2_hat  (probability of a "near" ordering error, always > p1)
		if (!setExpectedInf()) { qcOk = false; notes += " | QC: setExpInf failed"; return false; }					// Couldn't infer expected information

		if ( ( cMinNpos   > 0.0 ) && ( ( set.resultValN() < cMinNpos ) || ( cos.resultValN() < cMinNpos) )) { qcOk = false; notes += " | QC: TooFew Sites " + std::to_string(set.resultValN()) + " " + std::to_string(cos.resultValN()); }	// Too few sites in one child...
		if ( ( cMinDirInf > 0.0 ) && ( delInfBits() < cMinDirInf ) ) { qcOk = false; notes += " | QC: delInfBits() too low " + std::to_string(delInfBits()); }										// Too little information in seperation 

		return true;
	};

	bool resultVal(bool isSet, const string &Group, const string &Field, double &Value) const { return((isSet ? set : cos).resultVal(Group, Field, Value)); };
	bool resultVal(bool isSet, const string &Group, const string &Field, string &Value) const { return((isSet ? set : cos).resultVal(Group, Field, Value)); };
	bool resultSet(bool isSet, const string &Group, const string &Field, const double &Value) { return((isSet ? set : cos).resultSet(Group, Field, Value)); };	// manually add a result.
	bool resultSet(bool isSet, const string &Group, const string &Field, const string &Value) { return((isSet ? set : cos).resultSet(Group, Field, Value)); };	// manually add a result.

	double resultValRho(bool isSet)  const { return((isSet ? set : cos).resultValRho()); };
	double resultValInsNLL() const { double s=set.resultValNLL(),c=cos.resultValNLL(); if (s==-1.0 || c== -1.0) return -1.0; return s+c; };
	double resultValRhoU(bool isSet) const { return((isSet ? set : cos).resultValRhoU()); };

	// In gneral, use a conservative estimate of the number of bits saved by the split, take the minimum of the poitn estimate (using ML valeus of rho) and the 
	//	expectation (integrating saved information over the independent posterior distributions of rho in the childeren).
	// Note, if user provided NLL from INSIGHT2 in -ipriors argument, use Insight2NLL that (nats, not bits) instead of inframtion about S for (delInf).  To disable and use S, withold parrental NLL from -ipriors list
	double delInfBits() const { if (!qcOk) return 0.0; double d=(cParNLL>=0.0?cParNLL-(set.resultValNLL()+cos.resultValNLL()):std::min(delInfBitsP(),delInfBitsE())); return std::max(0.0,d);  };
	double delInfBitsP() const { double d; return(delInfBits(d) ? d : 0); };

	// Get Expected Value - Fixing MAX_STDEV, SAMP_RANGE and NUM_SAMNP here is a hack, and we are getting p calues fron C1/2 rather than from the actual parent. Another HACK.
	// But all in all, it shoud be good enough. 
	bool   setExpectedInf() {
		if (infExp.isSet()) return true;
		expectedInfo::paramSet	pr,c1,c2;
		const double   MAX_STDEV	= 3.0;	 // if stderr is nan, use a noninfromative value. Range is 0-1, so a stdev of 3 is huge. Largest observed in prelim run is, 1.57, smallest is 0.000514, the rest are nan.
		const double   SAMP_RANGE	= 3.0;   // spread samples over (up to) +/- SAMP_RANGE (3?) standard deviations, or support (0-1) whichever is less. (3 Stdev is 99.7% of gaussian, or 89% of Chebychev mass)
		const uint32_t NUM_SAMP     = 1000;  // Number of sampels per parmeter (1/gridsize), 100 is likely too few, and 1000 is likely generous (consider delta rho =.01 vs .001). Comp time is N^2 in this.
		{ double n, p, u; n = set.resultValN(); p = set.resultValRho(); u = set.resultValRhoU(); if (n<0 || p<0) return false; if (u <= 0) u = MAX_STDEV; c1.set(n, p, u); }
		{ double n, p, u; n = cos.resultValN(); p = cos.resultValRho(); u = cos.resultValRhoU(); if (n<0 || p<0) return false; if (u <= 0) u = MAX_STDEV; c2.set(n, p, u); }
		// Get the true parrential information from priors... 
		{	
			double n, p, u; n=c1.n + c2.n; p = (c1.ns +c2.ns) / n; 
			double prho=0.0; scorePair::ParentalRho( &prho ); if (prho>0.0 && prho < 1.0) p = prho;
			// Note, we dont have parrental information about error. It is not used, so fudge it....
			u = sqrt( (c1.sigma*c1.sigma*c1.n + c2.sigma*c2.sigma*c2.n) / n ); pr.set(n, p, u); 
			if (!expectedInfo::estimate(pr,c1,c2,NUM_SAMP,SAMP_RANGE,infExp)) return false;
		}
		return true;
	}
	bool delInfBitsE(expectedInfo::resultSet &Res) const { if (!infExp.isSet()) return false;  Res = infExp; return true; }
	double delInfBitsE() const  { if (!infExp.isSet()) return 0.0; return (cMinDirInf>0.0?infExp.expInfPart: infExp.expInf); }

	// Local consistency, splits are always nonnegative
	double delInfBitsLocal() const {
		double parn=set.resultValN() + cos.resultValN(); double parrho=(set.resultValN()*set.resultValRho() + cos.resultValN()*cos.resultValRho())/parn;
		double ucond = parn * mathplus::entropyB(parrho) * mathplus::natToBit;
		double cond  = (set.resultValN()*mathplus::entropyB(set.resultValRho()) + cos.resultValN()*mathplus::entropyB(cos.resultValRho()))* mathplus::natToBit;
		return ucond - cond;
	}
	// can be negative.... depending on consistency model and paretnal value.
	bool delInfBits(double &DeltaBits) const { double ucond = 0.0;	return(delInfBits(DeltaBits, ucond)); }
	bool delInfBits(double &DeltaBits, double &UcondBits) const { double cond = 0;	return(delInfBits(DeltaBits, UcondBits, cond)); }
	bool delInfBits(double &DeltaBits, double &UcondBits, double &CondBits) const {
		double rhoSet, rhoCos, nSet, nCos, rhoTot, nTot, maxW;
		// Can we get N (Numberof Samples) from Results file? We shoudl be able to.
		if (!resultVal(true, "PARAMETERS", "Rho", rhoSet))	return false;
		if (!resultVal(false, "PARAMETERS", "Rho", rhoCos))	return false;
		if (!resultVal(true, "STATS", "InPosW", nSet))		return false;
		if (!resultVal(false, "STATS", "InPosW", nCos))		return false;
		nTot = nSet + nCos;  if (nTot <= 0) return false;
		if (!resultVal(true, "ARGS", "maxSampW", maxW))		return false;
		if (maxW<1) return false;
		// Here default to local consistency,  DeltaBits is always >=0.
		rhoTot = ((nSet * rhoSet) + (nCos * rhoCos)) / (nSet + nCos);
		{
			// if requested, optain parrental estiamte of rho to get GLOBAL consistency. DeltaBits might negative.
			// While there is no "parrental rho" argument, we DO use the parrental RHO as a prior, so it is available.
			double rhoPar=0.0; scorePair::ParentalRho( &rhoPar );	// fetch the global parrental rho. If 0, use LOCAL consistance rather than global.
			if (rhoPar>0.0 && rhoPar<1.0) rhoTot = rhoPar;
		}
		CondBits = nSet * mathplus::entropyB(rhoSet) + nCos * mathplus::entropyB(rhoCos); CondBits *= mathplus::natToBit; CondBits /= maxW;
		UcondBits = nTot * mathplus::entropyB(rhoTot);  UcondBits *= mathplus::natToBit; UcondBits /= maxW;
		DeltaBits = UcondBits - CondBits;
		return true;
	}

	string toStrSum() const {
		// Covaraite Info
		string s = "", covname = set.CovName();
		//double numcells; resultVal(false, "ARGS", "maxSampW", numcells);
		double numvals = set.covsetDef().getCovByName(covname)->numVals(true);
		const scoreSingelton *a = &(getPart(true)), *b = &(getPart(false)); if (a->resultValRho() < b->resultValRho()) { std::swap(a, b); }
		// Info from a coset
		double a_numvals, a_rho, a_rhou, a_numpos;
		a_numvals = a->covsetDef().getCovByName(covname)->numVals();
		a_rho = a->resultValRho(); a_numpos = a->resultValN(); a_rhou = a->resultValRhoU();
		// Info from b coset
		double b_numvals, b_rho, b_rhou, b_numpos;
		b_numvals = b->covsetDef().getCovByName(covname)->numVals();
		b_rho = b->resultValRho();	b_numpos = b->resultValN(); b_rhou = b->resultValRhoU();
		// Info from joint analysis
		double condbits, ucondbits, expbits=-1.0, exprho=-1.0, delbits, prho;
		delInfBits(delbits, ucondbits, condbits);
		if (infExp.isSet()){	// this is set whe nresutls are read so baring error, it should always be present if anything is.
			expbits=infExp.expInf; exprho=infExp.expRho;
			if (cMinDirInf > 0.0) {		// use DIRECTED information
				expbits = infExp.expInfPart;  exprho=infExp.expRhoPart;
			}
		}
		// Format and return information
		double lprho = ((a_numpos * a_rho + b_numpos * b_rho) / (a_numpos + b_numpos));
		prho = lprho;
		{
			double pRho = 0.0; scorePair::ParentalRho(&pRho);
			if (pRho>0.0 && pRho < 1.0 ) prho = pRho;
		}
		double ins2_delnats=( cParNLL < 0.0 ? 0.0 : cParNLL-(set.resultValNLL()+cos.resultValNLL())), ins2_delbits = ins2_delnats * mathplus::natToBit;
		// p1 - Probability that rho1 > rho2 | rho1_hat > rho2_hat , delta_rho1_hat, delta_rho2_hat (always > 50% for normal error). C++ has erf/erfc, now! Woo Hoo!
		// Gaussian convolution implies a normal distribution with variance as sum of convolved distributions and centrality at 0 (coincident centers)
		//double sigma=std::sqrt( a_rhou*a_rhou + b_rhou * b_rhou );
		//double p1 = std::erfc( std::fabs(a_rho - b_rho ) / sigma ) * 0.5;			// OK this is the probability that rho2 > rho1 |  rho1_hat > rho2_hat  (probability of an ordering error)
		//double p2 = std::erfc( std::fabs((a_rho - b_rho)*0.5) / sigma) * 0.5;		// OK this is the probability that ( rho1 - rho2 ) < (rho_hat1 - rho_hat2)/2 |  rho1_hat > rho2_hat  (probability of a "near" ordering error, always > p1)
		s = fmt(covname, 8) + "\t" + fmt2(numvals, 2) + "\t" + fmt2(a_numvals + b_numvals, 2) + "\t" + fmt(prho, 8, 6) + "\t" + fmt(exprho, 8, 6) + "\t" + fmt(lprho, 8, 6) + "\t" + fmt2(cParNLL*mathplus::natToBit,9,true) + "\t";
		s += fmt2(a_numvals, 2) + "\t" + fmt(a_rho, 8, 6) + "\t" + fmt(a_rhou, 8, 6) + "\t" + fmt2(a_numpos, 13, true) + "\t" + fmt2(a->resultValN(true), 13, true) + "\t" + fmt2(a->resultValNLL()*mathplus::natToBit, 10, true) + "\t";
		s += fmt2(b_numvals, 2) + "\t" + fmt(b_rho, 8, 6) + "\t" + fmt(b_rhou, 8, 6) + "\t" + fmt2(b_numpos, 13, true) + "\t" + fmt2(b->resultValN(true), 13, true) + "\t" + fmt2(b->resultValNLL()*mathplus::natToBit, 10, true) + "\t";;
		s += fmt2(a_numpos + b_numpos, 13, true) + "\t" + fmt2(a->resultValN(true)+b->resultValN(true), 13, true) + "\t" + fmt2(delbits, 10, true) + "\t" + fmt2(expbits, 10, true) + "\t" + fmt2(delInfBitsLocal(),10,true) + "\t" + fmt2(ins2_delbits, 10, true) + "\t" + fmt2(ins2_delnats, 10, true) + "\t" + fmt(delbits / (a_numpos + b_numpos), 8, 6) + "\t";
		s += fmt(cPerr1, 10, 1,true) + "\t" + fmt(cPerr2, 10, 1, true) + "\t" + (qcOk?"1":"0") ;	// these values are NOT (bonferroni) corrected!
		return s;
	};

	bool Init(const covsetDefSubset &MasterCovSet, const BitEncoder_t &EncoderAnn, const BitEncoder_t &EncoderCTS, const string &CovName, const vecD &Order, const string &CovValTagsSet, const string &CovValTagsCos, const string &FileBaseName, const covtreeArgs *Args = NULL ) {
		// Stash the parameters for later use
		clear();
		covOrder = Order;
		if (!set.Init(MasterCovSet, EncoderAnn, EncoderCTS, CovName, CovValTagsSet, FileBaseName)) return false;
		if (!cos.Init(MasterCovSet, EncoderAnn, EncoderCTS, CovName, CovValTagsCos, FileBaseName)) return false;
		if (Args != NULL) {
			// apply heuristic limits in allowable partitions...
			cMaxPerr1=  Args->hMaxPerr1; cMaxPerr2 = Args->hMaxPerr2; cMinNpos = Args->hMinChildSz; cMinDirInf = Args->hMinDirInf; 
			cParNLL = (Args->splitOnRho?-1.0:std::stod(Args->vpriors[4]));	// NOTE insight works in nats..... dont mix NATS (from insight) with BITS (used internally)
		}
		return true;
	};


	static bool GenerateCosets(const covsetDefSubset &MasterCovSet, const covsetCovOrder_t &CovOrders, const BitEncoder_t &EncAnn, const BitEncoder_t &EncCTS, const string &FileBase, scorePair_v &PairCovSets, int v=0, const covtreeArgs *args=NULL ) {
		string safeFileBase = FileBase + "/part";
		fastfile::ffmkdir( safeFileBase.c_str());	// Returns 0 on success, but fails if dir exists, which makes check tough.
		PairCovSets.clear();
		for (int b = 0; b < 2; b++) {	// Loop over Covariet Classes: Annotations (0) and CTS covariates (1)
			bool iscts = (b == 1);
			if (v > 2 ) { std::cerr << "\t\t GenerateCosets: iscts - " << b << " NumCovIDS   " << ((int) MasterCovSet.getNumCovIDS(iscts)) << " ." << std::endl; }
			for (uint8_t icov = 0; icov<MasterCovSet.getNumCovIDS(iscts); icov++) {		// For each Class, loop over all covariates in that class
				const covDefSubset * cov = MasterCovSet.getCovByID(iscts, icov); if (!cov) return false;
				if (v > 0) 	std::cerr << "\t\t\t GenerateCosets: COV - " << cov->Name() << "  with " << cov->numVals() << " values." <<  std::endl;
				// If there is only 1 covariate value, we can't partition it, so skip that covariate...
				if (cov->numVals() > 1) {
					scorePair::covValWeight_v order;
					const scorePair::covValWeight_v *porder = NULL;
					if (!cov->univ->isMon()) {
						stringer a = cov->Name();
						order = CovOrders.at(a); porder = & order;
					}
					std::vector<covDefSubset> ValsSet, ValsCos;
					cov->AllPartitionSet(porder, ValsSet, ValsCos,(v>0?v-1:0));	// For the current ordered covariate, genreate a list of all covariate value sets and thier cosets...
					if (ValsSet.size() != ValsCos.size()) return false;
					if (ValsSet.size() > 0 ) {
						// should never happen, but if there is only one subset of values, then there's nothign to compare... jsut skip the covariate
						for (uint32_t ivalset = 0; ivalset <ValsSet.size(); ivalset++) {		// for the list of covariate value sets/ cosets
							PairCovSets.resize(PairCovSets.size() + 1);						// Allcoate a new scorepair object
																							// Then initialize the scorepair object....
							bool ok = PairCovSets[PairCovSets.size() - 1].Init(MasterCovSet, EncAnn, EncCTS, cov->Name(), order, ValsSet[ivalset].tagList(), ValsCos[ivalset].tagList(), safeFileBase,args);
							if (!ok) return false;
						}
					}
				}
			}
		}
		return true;
	};
};



