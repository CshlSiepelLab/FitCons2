#pragma once

#include <string>
#include "butils/butils.h"

// prior values determined genome-wide based on hg19
// Current numbers from ~/projects/fitcons2/baseline/hg19/wrk/2-hg19-posterior/hg19-20-parallel 10-Mar-2016 --Brad
// based on # Informative sites: 1,533,751,915 Mono, 6,033,316 PolyL, 1,514,816 PolyH 
// Neutrals based on # Sites: 4,110,726 PolyL and 1,062,043 PolyH
struct HG19 {
	// these values are posterior expectations (prior derived from "pcounts" observations drawn from distribution of sites udner ML values)
	static constexpr double rho    = 0.07342932116859754443;	//  0.072462  - 0.074395  95% symmetric credal interval
	static constexpr double eta    = 0.00021737936739721811;	//  5.686e-06 - 0.000943
	static constexpr double gam    = 0.64059032649056757425;	//  0.629478  - 0.651717
	// these values are genome wide data expecations (means), used only in calculating priors.
	static constexpr double lamda  = 0.00387775695968175290;
	static constexpr double theta  = 0.00094993578249541130;
	static constexpr double lamdaN = 0.00390664095888410348;	// based on sites with NEUTRAL polys only.
	static constexpr double thetaN = 0.00096222208636776438;
	// these values are genome wide maximum likelihood values
	static constexpr double rho_ml = 0.07342727348456748460;	// +/- .00052
	static constexpr double eta_ml = 6.899999750833038e-10;		// +/- .00823, actual value is 0. Non-0 value used for nuerical stability << 1 site in genome.
	static constexpr double gam_ml = 0.64058786208536111495;	// +/- .00431
	static constexpr double beta1  = 0.75666287190656289496;
	static constexpr double beta2  = 0.20531421372189487262;
	static constexpr double beta3  = 0.03802291437154223242;
	// from :~/projects/fitcons2/baseline/databaseRevision/prior-hg19/wrk/3-calcPrior.sh.log Mean of 4 cell type specific values.
	// gm 433.9512, hu 260.2172, h1 294.5597 and 398.1135 i6
	static constexpr double pcounts = 347.0;					// Actually 346.7104; pseudocounts 
	static std::string toStr() {
		std::ostringstream s; s.precision(15);
		s << "Expctn:\tRho:    " << rho    << "\tEta:   " << eta    << "\tGamma:   " << gam    << std::endl;
		s << "ML    :\tRho:    " << rho_ml << "\tEta:   " << eta_ml << "\tGamma:   " << gam_ml << std::endl;
		s << "Betas :\tBeta1:  " << beta1  << "\tBeta2: " << beta2  << "\tBeta3:   " << beta3  << std::endl;
		s << "Expctn:\tLambda: " << lamda  << "\tTheta: " << theta  << "\tLambdaN: " << lamdaN << "\tThetaN: " << thetaN << "\tPseudocounts: " << pcounts << std::endl;
		return( s.str() );
	};
};

struct Args {
	typedef std::string string;
	typedef inscomp::modelInsight model_t;
	struct Betas {
		struct bVal {
			double b1, b2, b3;
			bVal() { clear(); }
			void clear() { b1 = HG19::beta1 ; b2 = HG19::beta2; b3 = HG19::beta3; renorm(); };	// based on hg19 whole-genome values
			void renorm() { double b = b1 + b2 + b3; if (b > 0) { b2 /= b; b3 /= b; b1 = 1.0 - b2 - b3; }; };
			bool parse(const std::string &S) {
				butils::stringer s(S);  butils::stringer::vecS toks; s.tokenize(toks, ",");	auto n = toks.size();  clear();
				if (n < 3) return false;
				if (toks[0].length() > 0) b1 = std::stod(toks[0]);
				if (toks[1].length() > 0) b2 = std::stod(toks[1]);
				if (toks[2].length() > 0) b3 = std::stod(toks[2]);
				renorm(); return true;
			};
			std::string toStr() { std::ostringstream s; s << "Beta1: " << b1 << "\tBeta2: " << b2 << "\tBeta3: " << b3; return(s.str()); };
		};
		bVal init; bool refine; bVal prior;	// initial position for serarch, option to refine init values using ML/MaP, prior dist to use for MaP, prior weight determined elsewhere.
		Betas() { clear(); };
		// by default set prior weight to one observation of the rarest class.... This is fairly weak sauce.
		void clear() { init.clear(); refine = false; prior.clear(); };
		bool parse(const std::string &S) {
			butils::stringer s(S);  butils::stringer::vecS toks; s.tokenize(toks, ",");	auto n = toks.size(); for (auto &s : toks) s = s.trim(); clear();
			if (n < 3) return false;
			init.parse(toks[0] + "," + toks[1] + "," + toks[2]); prior = init;	// get initial value for betas
			if (n < 4) return true;
			if (toks[3].length() > 0) refine = (std::stoi(toks[3]) != 0);			// refine the provided betas using actual data?
			if (n < 7) return true;
			prior.parse(toks[4] + "," + toks[5] + "," + toks[6]);				// use different priors for MaP refinement....
			return true;
		}
		std::string toStr() { std::ostringstream s; s << "INIT:\t" << init.toStr() << "\tREFINE:\t" << refine << "\tPRIOR:\t" << prior.toStr(); return(s.str()); };
	};
	struct Param {
		double vinit, prior;
		bool refine;
		Param() { clear(); };
		void clear() { vinit = prior = 0.5; refine = false; }
		bool parse(const std::string &S) {
			butils::stringer s(S);  butils::stringer::vecS toks; s.tokenize(toks, ",");	auto n = toks.size();  for (auto &s : toks) s = s.trim(); clear();
			if (n < 1) return false;
			if (n >= 1 && toks[0].length()>0) prior = vinit = std::stod(toks[0]);
			if (n >= 2 && toks[1].length()>0) refine = (std::stoi(toks[1]) != 0);
			if (n >= 3 && toks[2].length()>0) prior = std::stod(toks[2]);
			return true;
		}
		std::string toStr() { std::ostringstream s; s << "Init:\t" << vinit << "\trefine:\t" << refine << "\tprior:\t" << prior; return(s.str()); };
	};
	Betas		betas;				// betas, if provided by user, one set of values for the entire input file...
	Param		lambdaT, theta;		// lambdaT and Theta are defined per block, these values are input-file (genome?) wide expectations, weighted across blocks...
	Param		lambdaTN, thetaN;	//   as above, but over neutral polymorphisms in blocks containing data, used for Beta refinement only.
	Param		rho, eta, gam;		// initial values, priors and requests to refine using MaP / ML
	int			v;					// verbosity level, 0 is standard, higher means more... 1->message per calc phase.. 2->message per block 3-> message per position
	string		dDB;				// directory name of database files
	string		dDBdef;				// default directory name of database files, before beign overridden by -ddb
	int			nAlleles;			// number of alleles, (usually twice the number of cell samples)
	int			nThreads;			// number of threads used in multithreading ops, if>1, should be >=4, only used in posterior -*parameter*- estimation (SLOW).
	bool		dataLikelihood;		// calculate the data likelihood.. (NLL, base e);
	bool		postSummary;		// generate sumary of posterior counts
	bool		postDetail;			// calculate posteriors (S & A) for each position, A is repressnted as sufficient identifying relationship between Z and X.
	bool		postParExpec;		// parameter Expectation from priors and likelihood, experimental, summary
	bool		postParExpecQ;		// Basket of features, refines model using expectation, produces supp stats based on expectation rather than MAP.
	bool		postParDist;		// parameter posteriors from priors and likelihood, experimental, detail
	int			postParGres;		// grid resolution for sampling parameter posteriors, good values are 7 (debugging) to 20 (lots of detail).
	bool		postAllele;			// Allele probability based on posterior model
	bool		printHelp;			// Print the help menu....
	int			postCntFac;			// Counterfactual allele output detail 10- S only, 20-S W A, 30 - Full Posterior
	int			postCntFacType;		// Counterfactual allele mode 1-simple, 2 more complex, 3 - not yet implimented....
	bool		dataSummary;		// Generate summary statistics from input data
	double		minEta;				// Allow Eta to go below 0, if desired by user... EXPERIMENTAL, may generate nonsense if set below -1.0...
	string		inlocifile;			// if "" then read database, then loop on cin looking for input files and arg lists..... otherwise just proces one file & exit.
	string		insfile;			// set to != "" to generate INSIGHT1 input files (X.beta.ins and X.ins)
	int			inscompat;			// make file compatible with Insight1, pick only neutral blocks containing POLYMORPHIC informative input sites
	double		priorWeight;		// weight of priors when calculating MaP and expectations, in # of observations (i.e. pseudo counts) use 0 to disable prior.
	double		maxSampleWeight;	// When Per-Positions weights are used in input .bed file, each position has an integral weight [0-Max]. This is the Max value,
									//		so dividing the input weight byt hnsi value yeilds a number inthe range [0-1.0]. Use for mixture models where weigth is
									//		number of cell types demonstrating this position and Max is the total number of cell-types.

	Args() { clear(); };

	void clear() {	// default behavior is do nothing.... quietly..
		betas.clear(); lambdaT.clear(); theta.clear(); rho.clear(); eta.clear(); gam.clear(); postCntFac = postCntFacType = inscompat = v = 0; dDB = dDBdef = inlocifile = insfile = ""; nAlleles = 108;
		postAllele= dataLikelihood = postSummary = postDetail = postParDist = postParExpec = dataSummary = postParExpecQ = printHelp = false; priorWeight = HG19::pcounts; minEta = 0.0;
		nThreads = 1; postParGres = 10;	// minimum grid resolution for useful calculation. 5 is OK for debugging, 20 is reccomended for std precision, 30 for high precision
		// revised versions from new  regularized lambda / theta..
		rho.vinit=rho.prior = HG19::rho ; eta.vinit = eta.prior = HG19::eta; gam.vinit = gam.prior = HG19::gam;
		lambdaT.vinit = lambdaT.prior = HG19::lamda; lambdaTN.vinit = lambdaTN.prior = HG19::lamdaN;
		theta.vinit   = theta.prior   = HG19::theta; thetaN.vinit   = thetaN.prior   = HG19::thetaN;
		maxSampleWeight = 1.0;
	}

	bool processArgs(int C, char **V);
	bool processArgs(const std::string &S);
	std::string toStr() {
		std::ostringstream s;
		s << "InsVer: " << INSIGHT_VER_STR << "\n";
		s << "DBdir:  " << dDB << "\n";
		s << "DBdirdf:" << dDBdef << "\n";
		s << "Infile: " << inlocifile << "\n";
		s << "Betas:  " << betas.toStr() << "\n";
		s << "Lambda: " << lambdaT.toStr() << "\n";
		s << "Theta:  " << theta.toStr() << "\n";
		s << "Rho:    " << rho.toStr() << "\n";
		s << "Eta:    " << eta.toStr() << "\n";
		s << "Gamma:  " << gam.toStr() << "\n";
		s << "Alls:   " << nAlleles << "\n";
		s << "nThreads:    " << nThreads << "\n";
		s << "DataLikelhd: " << dataLikelihood << "\n";
		s << "PostSum:     " << postSummary << "\n";
		s << "PostDetail:  " << postDetail << "\n";
		s << "PostParmExp: " << postParExpec << "\n";
		s << "PostParmDist:" << postParDist << "\n";
		s << "PostParmGres:" << postParGres << "\n";
		s << "PostAllele:  " << postAllele << "\n";
		s << "PostCntFacD: " << postCntFac << "\n";
		s << "PostCntFacT: " << postCntFacType << "\n";
		s << "PriorWt:     " << priorWeight << "\n";
		s << "DataSummary: " << dataSummary << "\n";
		s << "MinEta:      " << minEta << "\n";
		s << "Insight1File:" << insfile << "\n";
		s << "Indight1Cmpt:" << inscompat << "\n";
		s << "maxSampW    :" << maxSampleWeight << "\n";
		return(s.str());
	}

	void loadModels(model_t &init, model_t &pri, model_t &initN, model_t &priN);	// load paremeters from arg structure into models....

	static void displayArgs(int C, char *V[]) {
		for (int c = 1; c<C; c++) { std::cerr << " Arg: " << c << "\tval: " << V[c] << std::endl; };
	}

	void displayHelp(int C=0, char *V[]=NULL, bool DisplayArgs = false);

};

bool Args::processArgs(int C, char **V) {
	std::string sarg = "";
	for (int i = 1; i < C; i++) sarg += (i>1 ? " " : "") + std::string(V[i]);
	return(processArgs(sarg));
}

// Proces input user arguments
bool Args::processArgs(const std::string  &S) {
	using std::cerr;
	using std::endl;
	butils::stringer::vecS V;
	butils::stringer(S).tokenizeQuoted(V);
	for (auto &s : V) s = s.trim();
	V.push_back("-"); V.push_back("-");	// null arguments... this makes parsing easier
	auto C = V.size();
	if (C<4) return false;
	//cerr << "\n\n Found " << C << " args " << endl;
 	//for (ulong i=0;i<C;i++) cerr << i << "\t" << V[i] << endl;
	// process required arguments...
	dDB = V[0];		// Database directory. Must contain block.bedg, poly.bedg, polyn.bedg, monoDB.db, monoDB.cohroms, monoDB.tags
	unsigned int count = 1;
	while (count < C) {
		std::string p = V[count];
		// verbose
		if (p == "-" || p =="" ) {	// do nothing...
		} else if (p == "-v") {					// set verbosity
			v = 1;  p = V[count + 1];		if (p.front() != '-') { v = std::stoi(p); count++; };
		} else if (p == "-a") {			// change default allele count			
			if (++count >= C) return (false);
			p=V[count]; nAlleles = std::stoi(p);
		} else if (p == "-ddb") {			// override default DB directory with higher perfgormance working directory, if running multiple Insight isntances on one cluster node.
			if (++count >= C) return false;
			p = V[count];  if (p.front() == '-') return false;	// directory must not start 
			dDBdef = dDB; dDB=p;								// save old default dir, and override with working directory.
		} else if (p == "-fin") {
			p = V[count + 1];  if (p.front() != '-') { inlocifile = p; count++; 
				if ((inlocifile.length()>3) && (inlocifile.substr(inlocifile.length()-4)==".bed"))  inlocifile = inlocifile.substr(0,inlocifile.length()-4);
			};
		} else if (p == "-mineta") {	// require a number as the next field...
			minEta = 0.0; p = V[count + 1];  minEta = std::stod(p); count++;
		} else if (p == "-ins1comp") {
			p = V[count + 1];  inscompat = 1; if (p.front() != '-') { inscompat=std::stoi(p); count++; };
		} else if (p == "-mkins") {
			p = V[count + 1];  if (p.front() != '-') { insfile = p; count++; };
		} else if (p == "-datasum") {		// summarize input file, and intersection with daabase...
			dataSummary = true;  
		} else if (p == "-postsum") {		// summarize A/S posterior counts
			postSummary = true;  p = V[count + 1];  if (p.front() != '-') { postSummary = (std::stoi(p) != 0); count++; };
		} else if (p == "-postall") {		// summarize A/S posterior counts
			postAllele = true;  
		} else if (p == "-postcntfac") {		// summarize A/S posterior counts
			postCntFac = 21;  p = V[count + 1];  if (p.front() != '-')	  { postCntFac = std::stoi(p); postCntFacType = postCntFac%10; count++; };
		} else if (p == "-postdet") {		// geenrate A/S posteriors for each position
			postDetail = true;  p = V[count + 1];	if (p.front() != '-') { postDetail = (std::stoi(p) != 0); count++; };
		} else if (p == "-expval") {		// estimate expected values for parameters
			postParExpec = true;
		} else if (p == "-expdist") {		// Detailed dump of posterior probabilities for parameters
			postParDist = true; postParExpec = true;
		} else if (p == "-gres") {			// Sampleing grid resolution, careful runtime grows cubically with res. 5 is fast 30 is impossibly slow. 10-20 usual.
			p = V[count + 1];  if (p.front() != '-') { postParGres = std::stoi(p); count++; }
		} else if (p == "-nthread") {		// number if threads used in sampling posterior grid, inefficent, so if >1 should be >=4. <= # cores..., usually 4-15 are good numbers
			p = V[count + 1];  nThreads = 1; if (p.front() != '-') { nThreads = std::stoi(p); count++; }
		} else if (p == "-h") {
			printHelp = true;
		} else if (p == "-nll") {
			dataLikelihood = true;
		} else if (p == "-rho") {			// initial value, range, search override and prior (probability) for parameter
			if (++count >= C) return (false);
			if (!rho.parse(V[count])) return(false);
		} else if (p == "-eta") {
			if (++count >= C) return (false);
			if (!eta.parse(V[count])) return(false);
		} else if (p == "-gamma") {
			if (++count >= C) return (false);
			if (!gam.parse(V[count])) return(false);
		} else if (p == "-lambdaN") {
			if (++count >= C) return (false);
			if (!lambdaTN.parse(V[count])) return(false);
		} else if (p == "-thetaN") {
			if (++count >= C) return (false);
			if (!thetaN.parse(V[count])) return(false);
		} else if (p == "-lambda") {
			if (++count >= C) return (false);
			if (!lambdaT.parse(V[count])) return(false);
		} else if (p == "-theta") {
			if (++count >= C) return (false);
			if (!theta.parse(V[count])) return(false);
		} else if (p == "-betas") {		// override values for beta, rather than find using EM...
			if (++count >= C) return (false);
			if (!betas.parse(V[count])) return(false);
		} else if (p == "-maxsampw") {		// Max weight for any observation (even if not encountered in data set)
			if (++count >= C) return (false);
			p = V[count];   maxSampleWeight = std::stof(p);
		} else if (p == "-priorwt") {		// Weight of prior in observations, or pseudocounts... 
			p = V[count + 1];  if (p.front() != '-') { priorWeight = std::stod(p); count++; };
		} else if (p == "-qmap" || p == "-qexp" ) {			// QUICK MAP/ML extimate of parameters and uncertanties... ML by default, ad posterior couts to get MAP w/ all of hg19 as prior.
			dataLikelihood = rho.refine = eta.refine = gam.refine = betas.refine = theta.refine = thetaN.refine = lambdaT.refine = lambdaTN.refine =  true; priorWeight = HG19::pcounts;
			if (p == "-qexp") { postParExpec = postParExpecQ = true; }
			p = V[count + 1];  if (p.front() != '-') { priorWeight = std::stod(p); count++; };
		} else {									// unknown argument...
			cerr << "Unable to parse argument " << count << " :" << p << ":, exiting." << endl;
			return false;
		}
		count++;
	};
	return(true);
}


void Args::loadModels(model_t &initM, model_t &priM, model_t &initNM, model_t &priNM) {
	auto copyBeta = [](const Betas::bVal &bsrc, model_t::parameters::betas_t &bdst) { bdst.b1 = bsrc.b1; bdst.b2 = bsrc.b2; bdst.b3 = bsrc.b3;  };
	model_t::parameters init=initM.getparameters(), pri=priM.getparameters(), initN = initNM.getparameters(), priN=priNM.getparameters();
	// data set parameters
	copyBeta(betas.init, init.beta); copyBeta(betas.prior, pri.beta );
	// model parameters
	init.rho = rho.vinit; pri.rho = rho.prior; init.eta = eta.vinit; pri.eta = eta.prior; init.gam = gam.vinit; pri.gam = gam.prior;
	// block param
	init.block.lambdaT = lambdaT.vinit;	pri.block.lambdaT = lambdaT.prior;
	init.block.theta = theta.vinit;		pri.block.theta = theta.prior;
	// copy general parameters to neutral model, then override with neutral specific values
	initN = init; priN = pri;
	initN.block.lambdaT = lambdaTN.vinit;	priN.block.lambdaT	= lambdaTN.prior;
	initN.block.theta	= thetaN.vinit;		priN.block.theta	= thetaN.prior;
	initN.rho			= 0;				priN.rho			= 0;	// neutral models ALWAYS have S=0, thus rho =0;
	// init.sup -- these are derived only, and never provided directly by user.
	initM.setparameters(init); priM.setparameters(pri); initNM.setparameters(initN); priNM.setparameters(priN);
	return;
}

// block.bedg, poly.bedg, polyn.bedg, monoDB.db, monoDB.cohroms, monoDB.tags
void Args::displayHelp(int C, char *V[], bool DisplayArgs) {
	using std::cerr;
	using std::endl;
	cerr << endl << endl << "Insight2 Version: " << INSIGHT_VER_STR << endl;
	cerr << endl << "Insight2 DirDB [-fin Fname] [args, become defaults in server mode]" << endl;
	cerr << endl << "May be invoked in 2 modes - single use and server." << endl;
	cerr << "\tSingle Use - requires specification of -fin argument. Reads db, processes command line and exits." << endl;
	cerr << "\tServer     - without -fin, reads db, then awaits series of lines from stdin." << endl;
	cerr << "\t\t\tthe word \"done\" on a line by itself terminates server mode and exits." << endl;
	cerr << "\t\t\tLine format:" << endl;
	cerr << "\t\t\t\tInFileName # [list of command line args]" << endl << endl;
	cerr << "DirDB dierctory must contain 6 files" << endl;
	cerr << "\tblock.bedg- string_chrom int_start int_end float_theta float_lambda" << endl;
	cerr << "\tpoly.bedg - string_chrom int_start int_end char_freqLH float_priMaj float_priMin" << endl;
	cerr << "\tpolyN.bedg- Poly file, but filtered to only contain positions in neutral loci." << endl;
	cerr << "\tMonomorphism database, consistign of .tags .chroms .db " << endl;
	cerr << "\t\tmonoDB.tags      int_ind str_tag  - a list of up to 255 tags, ind is in range 1-255, tag is any string (usually a string formatted float)" << endl;
	cerr << "\t\tmonoDB.chroms    sorted lsit of chromosome extents, bed format. string_chromid 0 int_chromLen" << endl;
	cerr << "\t\tmonoDB.db        binary file, one byte per position in .chroms, each byte is the index of an entry in .tags. 0 for missing data." << endl;
	cerr << endl << "QuickCommands: " << endl;
	cerr << "\t-h            - Print this help menu." << endl;
	cerr << "\t-v [level]    - Verbosity flag. 0 is default without flag. 1 is default with flag increasingly verbose up to ~5." << endl;
	cerr << "\t-ddb dname    - override DirDB, but in the argument list. Has no effect when provided in server mode as DB is read before server loop is entered." << endl;
	cerr << "\t-fin fname    - source file for positions of interest. Bed format. Do not include .bed suffix." << endl;
	cerr << "\t-qmap [N]     - ML/MaP estimate for beta,rho,eta,gamma parameters. N must be an integer>=0. N=0 performs ML. N>0 (default:" << HG19::pcounts << ") provides N pseudocoutns of hg19 distrib as prior." << endl;
	cerr << "\t-qexp [N]     - Expectation for rho,eta,gamma parameters. N, if specified Must be >0, (default:" << HG19::pcounts << ") hg19 pseudocount prior. Betas estimated via ML." << endl;
	cerr << "\t-mkins fname  - Generate Insight1 input files from database and estimate betas." << endl;
	cerr << endl << "DetailedCommands: " << endl;
	cerr << "\t-a int        - number of alleles, echoed in first line of output file. Default is 108." << endl;
	cerr << "\t-mineta D     - Minimum value for eta in ML asessment. Generally 0. May be set lower to investigate ML stability. Less than -5 not reccomended." << endl;
	cerr << "\t-datasum      - Primt summary of nucleotide positions and intersection with database." << endl;
	cerr << "\t-postsum      - Summary of Insight1 Posterior information, over all sites innput file from refined model." << endl;
	cerr << "\t-postdet      - Summary of Insight1 Posterior for each site with database information." << endl;
	cerr << "\t-postall      - Insight2 Posterior allele distributions for each site with database information." << endl;
	cerr << "\t-postcntfac I - Insight2 Counterfactual analysis of impact of novel allele at each site w/ database information." << endl;
	cerr << "\t                I is a 2 digit number XY, with Y=[1,2] representing type of analysis and Y=[0,1,2,3] verbosity of output." << endl;
	cerr << "\t-expval       - Generate expected values for parameters rho, eta, gamma." << endl;
	cerr << "\t-expdist      - Generate prior, likelihood and posterior probabilities at point sampled grid and marginals for each parameter." << endl;
	cerr << "\t-gres         - Grid resolution for posterior parameter asessment. Runtime grows cubically. 5 is fast, 30 is impossibly slow. 10-20 typ. default 10." << endl;
	cerr << "\t-nthread      - Number of threads used in posterior parameter asessment. Inefficent, so if >1, use >=4. 4-12 typ." << endl;
	cerr << "\t-nll          - Print total data likelihood under model for positions in database, in nats." << endl << endl;
	cerr << "\t-maxsampw D   - Maximum input position weight (number of samples). Defaults to 1 (unweighted case)." << endl << endl;
	cerr << "\t-priorwt D    - Strength of the prior measured as D pasudocounts of data from whole autosome hg19 distribution." << endl << endl;
	cerr << "\t-ins1comp [N] - Insight1 compatibilty mode def N=1: 0=Insight2, 1=use only blocks w/ informative poly sites in fin (Insight1)." << endl;
	cerr << endl << "Parameter initialization:" << endl;
	cerr << "\tFor each parameter listed below, provide default values, as well as a flag indicating if that parameter should be refined." << endl;
	cerr << "\tGenerally Initial values are the hg19, whole autosome expectations. Priors default to Init values when not provided." << endl;
	cerr << "\tFormat is [Init[,Refine[,Prior]]" << endl;
	cerr << "\t\tInit   - real value for initial search value in ML gradient descent." << endl;
	cerr << "\t\tRefine - 0 - supress refinement of this variable, 1 - allow." << endl;
	cerr << "\t\tPrior  - Defailts to Init. Overrides init value as prior value for this parameter." << endl << endl;
	cerr << "\t\tExample: -rho .65,1,.077" << endl;
	cerr << "\tAvailable parameters are:" << endl;
	cerr << "\t\t-rho, -eta, -gamma, -lambda, -lambdaN, -theta, -thetaN" << endl;
	cerr << "\t\tlambda and theta are only used in prior model, and are generally be inferred from actual data. N referes to Neutral Sites." << endl << endl;
	cerr << "\t-betas b1,b2,b3 - override default values for beta when not inferred from data. Values renormalized to 1. Defaults are based on whole autosome hg19." << endl << endl;
	cerr << endl;
	if (DisplayArgs) {
		cerr << "Args:" << endl;
		Args::displayArgs(C, V);
		cerr << endl << endl;
	}
	return;
}
