// Insight2.cpp : Defines the entry point for the console application.
//
#define INSIGHT_VER_STR "0.16e"
// 
// 12.b Updated for new prior values based on whole genome ML values with regularized lambda / theta dabase.
// .c Fixed several issues in posterior
// .d Added cheap multithreading for posteriors - use nthreads=1 to use single threading (old code) tested up to 16 threads, 
//		for all of hg 19 requries 1GB of mem / thread plus about 3.5 GB for the database. Less for smaller data sets.
// 14 Updated defaults eith expected values and prior counts (see Insight2args::HG19 for constants & intervals...).
//	all basd on new posterior labda/ theta database using pseudocoutns and at least 1 count per window....
// 14a - fixed some bugs in server mode interation and termination. Allow -fin file to have a .bed suffix (removed upon entry)
//		fixed early termination of input file reads past database blocks that truncated count of # of physical positions.
// 15 - Support user argument for setting max weights. This allows for proper interpertation as a mixture model. 7 Observations
//		of position A, in a sample size of 8, is 7/8 of an observations, not 7 observations. ML/MAP/Posterior & unweighted calcus
//		are uneffected, but this does effect asolute NLL, MAP/ML confidence and posterior calculations...
// 15a- Bugfix, 1-Insight1 ins file generation bug (both regular and compabilility). Also enhancement, setting -maxsampw 1 (default)
//		disables reading and interpertation of 4'th field in bed file. Uses default weight is 1. So if tag is noninteger it is
//		ignored rather than generating an error.
// 15b - Bugfix, argument processing (-a argument)
// 15c - Bugfix, reduced search range for Eta,Gamma from 0-20,0-20 to 0-2,0-5 to handle a stability issue in CalcNLL (p<0 -> p=MinDbl, but Deriv -> inf)
// 15d - Bugfix, replaced command line tokenizer with tokenizeQuoted so suplicated delimiters (multiple spaes) are treated as 1
// 16  - Compat, added FastFile::fflush() call to flush .model IO to disk before .done semaphore is written, preventing premature read by extenral program.
// 16a - Add -ddb argument to override dirDB in argument to utilize local working directory for database reads.
// 16b - update optimizar failure behavior to use the best found value even if technical covneregence fails (rather than the prior). Fixed 2 new compiler warnings.
// 16c - Loosten hard limit on eta from 2.0 to 6 and gamma from 5.0 to 12.
// 16d - Fix a test in the Locus.h SkipPrevious() function that inaproproaitely missed positions in the last block.
// 16e - Tighten initial tolerance for beta calculation to improve consistency. Apply adaptive tolerance backoff in case of convergence failure.

#include <string>		// slow string operations
#include <iostream>		// slow IO for debugging and status messages
#include <fstream>		// slow IO for short output files
#include <sstream>		// slow IO for debugging and status messages
#include <ctime>		// slow, for debugging
#include <cassert>

#include "butils/butils.h"		// namespace: general utilities
#include "insdb/insdb.h"		// namespace: fast read and extract of database files
#include "inscomp/inscomp.h"	// namespace: compressed storage and calculation on database subsets.

#include "Insight2args.h"		// Interpret user args, print help and debugging info if needed

using inscomp::BlockSet_t; using inscomp::optInsight; using inscomp::optBeta;
using inscomp::modelInsight; using inscomp::modelInsightPosterior; using inscomp::modelInsightSup;


struct suppStats {
	typedef butils::ulong	ulong;
	typedef std::string		string;
	ulong		inLoci,   inPos;				// pertain to input bed file from extractInputLoci()
	uint64_t	inLociW,  inPosW;
	ulong		infPosM,  infPosL,  infPosH;	// pertain to number of informative positions (intersection of input with database)
	uint64_t	infPosWM, infPosWL, infPosWH;	// pertain to weighted number of informative positions (intersection of input with database)
	double	nll, nlpp;							// nll of Data and nll of Prior, from calcLikelihood(), if requested
	double  theta, thetaN, lambda, lambdaN;		// calculated if refine lambda/theta is requested... is requested...
	double	maxWeight;							// from user input, is important for itnerperting weighted valeus like  inPosW...
	void clear() { theta = thetaN = lambda = lambdaN = inLoci = inPos = infPosH = infPosL = infPosM = 0; inLociW = inPosW = infPosWH = infPosWL = infPosWM=0; nll = nlpp = 0.0; maxWeight = 1.0; };
	string fmt(double v, const char *f = "%.0lf") { static char b[32]; sprintf(b, ( (v==0) || (fabs(v) >= 1e-4) ? f : "%.6le"), v); return(b); };
	string toStr() {
		string s = "";
		s += "InLoci:  \t" + fmt(inLoci)      + "\nInPos:   \t"  + fmt(inPos) + "\nInLociW: \t" + fmt((double)inLociW) + "\nInPosW:  \t" + fmt((double)inPosW) + "\n";
		s += "InfPosM: \t" + fmt(infPosM)     + "\nInfPosWM: \t" + fmt((double)infPosWM) + "\n";
		s += "InfPosL: \t" + fmt(infPosL)     + "\nInfPosWL: \t" + fmt((double)infPosWL) + "\n";
		s += "InfPosH: \t" + fmt(infPosH)     + "\nInfPosWH: \t" + fmt((double)infPosWH) + "\n";
		s += "DataNLL: \t" + fmt(nll,"%.6lf") + "\nNLPrior: \t"  + fmt(nlpp,"%.6lf") + "\n";
		s += "MaxWeight:\t" + fmt(maxWeight, "%.2f") + "\n";
		return(s);
	};
	suppStats() { clear(); };
};


typedef inscomp::optInsight::vecD vecD;				// std::vector of doubles...
typedef vecD::size_type vecDi;

bool extractInputLoci(insdb::DB &db, const Args &a, BlockSet_t &bs, BlockSet_t &bsn, suppStats &stats, suppStats &statsN);
bool addpriorCounts(const inscomp::optInsight::parameters &priorparam, BlockSet_t &bs );
bool estimateLambdaTheta(const BlockSet_t &bs, bool byBlock, double &LambdaT, double &Theta, double priorWeight, const inscomp::modelInsight *priorModel);
bool refineBeta(const BlockSet_t &bsn, modelInsight &model, double PriorWeight = 0, const modelInsight *mPrior = NULL);
bool refineREG(const BlockSet_t &bs, modelInsight &model, modelInsight &modelE, const Args &a, modelInsightSup::derived_t &deriv, modelInsightSup::derived_t &derivE, double PriorWeight = 0.0, modelInsight *PriorModel=NULL);
bool posteriorDetail(insdb::DB &db, optInsight::parameters &model, const Args &a, std::string &fnin, std::string &fnout);
bool calcPostAllele(insdb::DB &db, const modelInsight &model, const Args &a, std::string &fnin, std::string &fnout);
bool calcPostCounterFac(insdb::DB &db, const modelInsight &model, const Args &a, int Detail, std::string &fnin, std::string &fnout);
bool calcLikelihood(const modelInsight &model, const BlockSet_t &bs, double &nll, double &nlpp, double PriorWeight = 0, const modelInsight *PriorModel = NULL );
bool calcParamPost(const BlockSet_t &bs, const modelInsight &model, double PriorWeight, const modelInsight &PriorModel, const Args &a, modelInsightPosterior::postSet &ppost, int GridRes=10);
bool dataStatistics(insdb::DB &db, const modelInsight &model, const Args &a, std::string &fnin, std::string &fnout);
bool makeIns(insdb::DB &db, const Args &a, modelInsight::parameters &model, const std::string &fnout, const std::string &fnoutb);
std::string strtime() { std::time_t ctt = std::time(0);	char *savetime = std::asctime(std::localtime(&ctt)); savetime[strlen(savetime)-1]='\0'; return( std::string(savetime) ); }

// Compile a database block to a high performance pattern block.
inline bool compileFastBlock(insdb::DB::fastblock &cumulant, BlockSet_t &bs );

using inscomp::modelInsight; using inscomp::modelInsightSup;

int
main(int argc, char *argv[])
{
	using std::cerr;
	using std::cout;
	using std::endl;
	using std::string;

	Args			adef;	// hold default user args...
	insdb::DB		db;		// read and access insight database files.
	BlockSet_t		bs,bsn;	// a dense, computation friendly, storage format for position-weighted database subsets. full (bs), and neutral polys (bsn).
	double			atime;	// double used to monitor wall-runtime, denominates in CLOCK_TICKS (ala clock() )

 	if (!adef.processArgs(argc, argv)) {
		cerr << "\nINSIGHT2: Error Processing Input args:\n\n";
		adef.displayHelp(argc, argv); cerr << "\n\nExiting." << endl;
		return(-1); }
	if (adef.printHelp || adef.v > 1 ) { 
		adef.displayHelp(argc, argv); 
		if (adef.printHelp) { 
			if (adef.v>0) { cerr << "\nPrior Parameter Values:\n" << HG19::toStr() << std::endl; };
			cerr << "\n\nExiting.." << endl; return( -2 ); 
		}
	};
	bool server_mode = adef.inlocifile.size() < 1;	// server mode or one-shot mode?

	if (adef.v > 0) { cout << strtime() << "\tProcessing database files from " << adef.dDB << "..."; cout.flush(); }

	if (!server_mode) { cout << strtime() << "\tReading binary Insight databse from " << adef.dDB << "/" << endl; };
	atime = (double)clock();
	int r = db.ReadDB(adef.dDB, adef.v); bool ok = (r == 0);
	if (!ok) { cerr << "Error reading database - exiting. Err Code:" << r << endl; cerr.flush();  return(-1); }
	if (adef.v > 0) { cout << "Database read compelte." << endl; };

	// Now loop over input files
	{
		using butils::ulong;
		using insdb::DB;
		Args a;						// argumens for this batch of loci
		string	bfname;				// source file defining this batch of loci
		std::string savetime = "";	// timign string, for debuggign purposes...
		//bool process_ok = false;	// did we generate this input file without error?
		bool ok = true, server_done = false;				// status of current operation
		inscomp::optInsight::parameters model;			// estimated parameters, as we calculate them
		inscomp::optInsight::parameters modelE;			// Uncertainty in estimated parameters,
		inscomp::optInsight::parameters prior_param;	// Priors an expected priors for parameters of population...

		suppStats	stats, statsN;						// spplimentry quantities of itnerest about the input data...

		// a model contains all relevant parameters, and expectations for block-specific parameters. 
		//		N means neutral (polys in blocks covered by data). Init is initial, Prior is population fitted, E means error (parameter uncertainty) 
		modelInsight modelPrior, modelInit, modelRefined, modelRefinedE, modelInitN, modelPriorN, modelRefinedN;
		modelInsightSup::derived_t deriv, derivE;
		ulong server_request = 0;						// count the nubmer of server requests so far.
		std::string donefile="";						// sync file for server-communication...

		if (adef.v>0) { cerr << "Entering event loop with server mode = " << server_mode << endl; }
		while ( !feof(stdin) && (  (server_mode && !server_done) || (!server_mode && (server_request<1)) ) ) {
			if (!server_mode && (server_request > 0 )) break;							// if interractive (non-server) mode, terminate after one loop
			if (!server_mode) server_request = 1;

			// record time for benchmarking purposes
			a = adef;																	// accept  defaults for arguments

			if (server_mode) {
				if (donefile.length()>0) { std::ofstream fo(donefile); fo.close(); donefile = ""; };	// in the event of an error, notify client of completion..
				butils::stringer linein; getline(std::cin, linein); linein = linein.trim();		// we are in server mode, read an input line from stdin
				if (a.v>0) { cerr << "\n" << server_request << "Server Mode Read Input Line: " << linein << endl; };
				// Skip empty lines
				if (linein.length() < 1) { cerr << "Empty Line in Input, skipping." << endl; continue; };
				// exit the server loop ?
				{	butils::stringer uc = linein; uc.toupper();
					if (uc == "DONE") { cerr << "Read DONE from input. Terminating Server mode." << endl; server_done = true; continue; };
				}
				std::vector<butils::stringer> toks;
				linein.tokenize(toks, "#"); for (auto &s : toks) s = s.trim();					// split input line, trim arguments
				a.inlocifile = toks[0].trim();
				if (a.inlocifile.length() < 1) { cerr << "Null input file in server mode, ignoring." << endl; continue; };
				if ((a.inlocifile.length()>3) && (a.inlocifile.substr(a.inlocifile.length() - 4) == ".bed"))  a.inlocifile = a.inlocifile.substr(0, a.inlocifile.length() - 4);
				donefile= a.inlocifile + ".done";
				if (toks.size() > 1) ok = a.processArgs("-\t-\t-\t-\t" + toks[1]);	// process arguments, if any. Requires 4 arguments, "-" does nothing, but increments arg count
				if (!ok) { cerr << "Error processing input line: " << linein << "\n\tSkipping.\n" << endl; continue; }
				if (a.printHelp) { a.displayHelp(); cerr << endl << "Arguments:" << endl << a.toStr() << endl << endl; continue; };
			}
			// { std::time_t ctt = std::time(0);  std::cout << std::endl << savetime << std::endl << std::asctime(std::localtime(&ctt)) << std::endl; }
			{	std::time_t ctt = std::time(0);	savetime = std::asctime(std::localtime(&ctt)); };
			if (server_mode) atime = (double) clock();

			// 1- process arguments, clear data, load defaults 
			stats.clear(); statsN.clear();
			modelInit.clear(); modelPrior.clear(); modelRefined.clear(); modelRefinedE.clear(); modelInitN.clear(); modelPriorN.clear(); modelRefinedN.clear();
			deriv.clear(); derivE.clear();

			// collect the priors from the arguments into cohesive models
			a.loadModels( modelInit, modelPrior, modelInitN, modelPriorN);
			modelRefined = modelInit; modelRefinedN = modelInitN;	// start refined models off as initial... then refine if requested...

			// 2 - EXTRACT DATA SUBSET from database and compile into computationally efficent representation
			bool do_refine = a.betas.refine || a.lambdaT.refine || a.theta.refine || a.rho.refine || a.eta.refine || a.gam.refine;
			bool do_posterior = a.postDetail || a.postParDist || a.postParExpec || a.postSummary;
			if (do_refine || do_posterior || a.dataLikelihood ) {
				if (!server_mode) { cout << strtime() << "\tIntersecting database with positive loci from " << a.inlocifile << ".bed" << endl; };
				if (a.inlocifile == "") {
					cerr << "Posterior values, likelihood or parameter refinement requested, but no source loci provided. Invalid combination. Skipping request." << endl; continue;
				};
				ok = extractInputLoci(db, a, bs, bsn, stats, statsN);	// take subset of database (db) identified by input file (a.inlocfile) and put remaidner in BlocksetBS.
				if (!ok) {
					cerr << "Insight2: extractInputLoci() failed to process input file : " << a.inlocifile << " :  Skipping request." << endl; continue;
				};
			};

			// {std::time_t ctt=std::time(0);  savetime=std::asctime(std::localtime(&ctt)); }
			// 3A Get Expected Lambda / Theta from data set... used in some prior calcualtions, where block-averages replace individual block values.
			model.block.lambdaT = a.lambdaT.vinit; model.block.theta = a.theta.vinit; 
			if (a.lambdaT.refine || a.theta.refine) {
				double lambda = 0, lambdaN = 0, theta = 0, thetaN = 0;
				double priorW = a.priorWeight;
				ok = estimateLambdaTheta(bs,  false, lambda,  theta,  priorW, &modelPrior );	// Weight by number of informative positions. (positives)
				stats.lambda = lambda; stats.theta = theta;
				ok = estimateLambdaTheta(bsn, true,  lambdaN, thetaN, priorW, &modelPriorN);	// Weight by number of total positions in block (neutrals)
				stats.lambdaN = lambdaN; stats.thetaN = thetaN;
				if (!ok) { cerr << "Insight2: estimateLambdaTheta failed. Skipping request." << endl; continue; };
				{
					modelInsight::parameters p = modelRefined.getparameters(), pn=modelRefinedN.getparameters();
					if (a.lambdaT.refine)	{ p.block.lambdaT = lambda; pn.block.lambdaT = lambdaN; };
					if (a.theta.refine)		{ p.block.theta = theta;	pn.block.theta = thetaN; };
					modelRefined.setparameters(p); modelRefinedN.setparameters(pn);
					if (a.v>0) { cerr << "\nrefined LT :\t" << lambda << "\t" << theta << "\nrefined LTN:\t" << lambdaN << "\t" << thetaN << "\n" << endl; };
				}
			}

			// 4A - Infer / refine BETAs if requested - Not weighted by input position weight.... perhaps should be .... (porportinate to Weight of input positions in each block...)
			if (a.betas.refine) {
				//double beta1 = 0, beta2 = 0, beta3 = 0, lambdan = 0, thetan=0;
				double priorW = a.priorWeight;
				if (!server_mode) { cout << strtime() << "\tEstimating Beta values " << endl; };
				ok = refineBeta( bsn, modelRefinedN, priorW, &modelPriorN );
				if (!ok) {
					cerr << "Insight2: Beta optimization failed, using initial values for model." << endl;
				}
				// copy betas from neutral model to active model...
				auto p = modelRefined.getparameters(); p.beta = modelRefinedN.getparameters().beta; modelRefined.setparameters(p);
				if (a.v>0) { cerr << "refined Betas:\t" << p.beta.toStr() << endl; };
			};

			// 4B - Infer/refine INSIGHT parameters rho / eta / gamma, derived values, and their uncertainty.
			model.rho = a.rho.vinit; model.eta = a.eta.vinit; model.gam = a.gam.vinit; modelE.clear();
			modelInsightSup::derived_t deriv, derivE;
			if (a.rho.refine || a.eta.refine || a.gam.refine) {
				// double priorW = (a.postParExpec ? 0 : a.priorWeight);	// if we are calculating posteriors, only use ML value to find grid center.
				double priorW = a.priorWeight;	
				if (a.v > 0) cerr << "INSIGHT2: Refining REG (RhoEtaGamma) parameters." << endl;
				if (!server_mode) { cout << strtime() << "\tEstimating MAP/ML parameters" << endl; };
				ok = refineREG( bs, modelRefined, modelRefinedE, a, deriv, derivE, priorW, &modelPrior );
				if (!ok) {
					cerr << "INSIGHT2: ERROR Refining REG (RhoEtaGamma) parameters. Skipping task." << endl; continue;
				};
				if (a.postParExpec) {
					string fout = a.inlocifile + ".model.map"; std::ofstream fo(fout);
					fo << "PARAMETERS:\n"  << modelRefined.toStr()  << endl << endl;
					fo << "UNCERTAINTY:\n" << modelRefinedE.toStr() << endl << endl;
					fo.close();
				}
			};
			
			// 5 - Generate Full Joint Data & Parameter Parameter distributon, replace parameter estimates with expectations, if requested.
			if (a.postParDist || a.postParExpec) {
				modelInsightPosterior::postSet ppost; 
				if (a.v > 0) { cout << "Insight2: calculating full joing parameter distribution and expectations." << endl; };
				if (!server_mode) { cout << strtime() << "\tEstimating posterior parameters distributions." << endl; };
				ok = calcParamPost(bs, modelRefined, a.priorWeight, modelPrior, a, ppost, (ulong) a.postParGres );
				if (!ok) { cerr << "Insight2: calcParamPost() returned false, unable to complete data sampling, skipping." << endl;  continue;
				} else {
					//cout << "EXPECTATIONS:\t" + ppost.toStrExp() << endl;
					if (a.postParExpec) {
						string fout = a.inlocifile + ".ppost.marg"; std::ofstream fo(fout);
						fo << ppost.toStr() << endl;
						fo.close();
					}
					if (a.postParDist) {
						string fout = a.inlocifile + ".ppost.full"; std::ofstream fo(fout);
						fo << ppost.elements[0].header() << "\n";
						for (ulong i=0; i<ppost.elements.size(); i++) fo << ppost.elements[i].toStr() << "\n";
						fo.close();
					}
				}
				if (a.postParExpecQ) { // use expectations as point estimates for further calculations
					auto p = modelRefined.getparameters();
					p.rho = ppost.exp.rho; p.eta = ppost.exp.eta; p.gam = ppost.exp.gam; do_refine = true;
					modelRefined.setparameters(p);
					modelInsightSup sup(modelRefined);
					sup.suppStatsSimple(a.nAlleles,deriv);	// recalculate dervied quantities based on expectations.... warning derivE still referes to ML values! 
					// use credibal interval data from ppost if you are itnerested in error quantities...
				}
			}

			// 6 - Generate INSIGHT1 input files and beta2, if that's all that is wanted. User must manually provide betas, perhaps
			//	using INSIGHT1 to estimate Betas, then appending these values to the .ins file...
			if (a.insfile != "") {
				modelInsight::parameters p = modelRefined.getparameters();
				if (a.v > 0) cerr << "INSIGHT2: Generating INSIGHT1 input files with basename " << a.insfile << " ." << endl;
				if (!server_mode) { cout << strtime() << "\tGenerating iInsight1 EM input file(s) " << a.insfile + ".ins" << " , " << a.insfile + ".beta.ins" << endl; };
				ok = makeIns(db, a, p, a.insfile + ".ins", a.insfile + ".beta.ins");
				if (!ok) { cerr << "INSIGHT2: makeIns() failed to generate files based at\n\t" + a.insfile + "\nSkipping request." << endl << endl; continue; };
			}

			// 7 - We can get likelihood from the refinement proceedure, but if we just want data likelihood, we need another pass...
			if (a.dataLikelihood) {
				if (a.v > 0) { cout << "Insight2: calculating data likelihood." << endl; };
				if (!server_mode) { cout << strtime() << "\tRecalculating data likelihood." << endl; };
				ok = calcLikelihood(modelRefined, bs, stats.nll, stats.nlpp, a.priorWeight, &modelPrior);
				stats.nll = stats.nll; stats.nlpp = stats.nlpp;
				if (!ok) {
					cerr << "Insight2: ERROR: problem claculating data likelihood.... Skipping." << endl;
				};
				if (a.v > 0) {
					cout << "Insight2: nll  = " << stats.nll << " (" << stats.nll / butils::mathplus::mylog(2) << " bits)" << endl;
					cout << "Insight2: nlpp = " << stats.nlpp << " (" << stats.nlpp / butils::mathplus::mylog(2) << " bits)" << endl;
				};
			};

			// 8 - Data posterior and supplimentary statistics, much like the -posterior function from insight1, but calculated in terms of sufficient properties, not A&Z
			if (a.dataSummary || a.postSummary || a.postDetail) {
				if (a.v > 0) { cout << "Insight2: calculating data summary stats." << endl; };
				if (!server_mode) { cout << strtime() << "\tCalculating supplimentary statistics - may be slow if per-site latent posteriors are requested." << endl; };
				string fnin = a.inlocifile + ".bed"; string fnout = a.inlocifile + ".lpos";	// latent class posteriros...
				ok = dataStatistics(db, modelRefined, a, fnin, fnout);
				if (!ok) {
					cerr << "Insight2: ERROR: problem claculating summary statts.... Skipping." << endl;
				};
			}

			// 9 - Dump the refined Model, if refinement occured.
			if (do_refine) {
				string fout = a.inlocifile + ".model";
#ifdef OLD_CODE
				std::ofstream fo(fout);
				fo << "Insight2 Version:\t" << INSIGHT_VER_STR << endl;
				fo << "PARAMETERS:\n"		<< modelRefined.toStr();
				if (deriv.alpha>=0) fo << deriv.toStr() << endl;
				fo << endl;
				fo << "UNCERTAINTY:\n"		<< modelRefinedE.toStr();
				if (deriv.alpha >= 0) fo << derivE.toStr() << endl;
				fo << endl;
				fo << "STATS:\n"			<< stats.toStr()			<< endl << endl;
				fo << "STATSN:\n"			<< statsN.toStr()			<< endl << endl;
				fo << "NEUTRAL_MODEL:\n"	<< modelRefinedN.toStr()	<< endl << endl;
				fo << "ARGS:\n"				<< a.toStr()				<< endl << endl;
				fo << "PRIOR_MODEL:\n"		<< modelPrior.toStr()		<< endl << endl;
				fo << "PRIORN_MODEL:\n"		<< modelPriorN.toStr()		<< endl << endl;
				fo.close();
#else
				FILE *fo = fopen(fout.c_str(),"wb"); string tmp;
				tmp = "Insight2 Version:\t" + string(INSIGHT_VER_STR) + "\n";	fputs( tmp.c_str(),fo);
				tmp = "PARAMETERS:\n" + modelRefined.toStr();					fputs(tmp.c_str(), fo);
				if (deriv.alpha >= 0) { tmp = deriv.toStr() + "\n"; 			fputs(tmp.c_str(), fo); }
				tmp = "\n";														fputs(tmp.c_str(), fo);
				tmp = "UNCERTAINTY:\n" + modelRefinedE.toStr();					fputs(tmp.c_str(), fo);
				if (deriv.alpha >= 0) { tmp = derivE.toStr() + "\n";			fputs(tmp.c_str(), fo); }
				tmp = "\n";														fputs(tmp.c_str(), fo);
				tmp = "STATS:\n" + stats.toStr() + "\n\n";						fputs(tmp.c_str(), fo);
				tmp = "STATSN:\n" + statsN.toStr() + "\n\n";					fputs(tmp.c_str(), fo);
				tmp = "NEUTRAL_MODEL:\n" + modelRefinedN.toStr() + "\n\n";		fputs(tmp.c_str(), fo);
				tmp = "ARGS:\n" + a.toStr() + "\n\n";							fputs(tmp.c_str(), fo);
				tmp = "PRIOR_MODEL:\n" + modelPrior.toStr() + "\n\n";			fputs(tmp.c_str(), fo);
				tmp = "PRIORN_MODEL:\n" + modelPriorN.toStr() + "\n\n";			fputs(tmp.c_str(), fo);
				butils::fastfile::fflush( fo ); fclose( fo ); fo = NULL; butils::fastfile::fflush();
#endif
			};

			// 10A - Posterior Allele probabilty analysis
			if (a.postAllele) {
				string fnin= a.inlocifile+".bed", fnout=a.inlocifile+".dpos.alp";
				if (!server_mode) { cout << strtime() << "\tCalculating Posterior Allele probabilities to " << fnout << "." << endl; };
				calcPostAllele( db, modelRefined, a, fnin, fnout );
			}

			// 10B - Posterior Allele counterfactual analysis
			if (a.postCntFac > 0) {
				string fnin= a.inlocifile+".bed", fnout= a.inlocifile+".dpos.cfa" ;
				if (!server_mode) { cout << strtime() << "\tCalculating Posterior Allele counterfactuals " << fnout << "." << endl; };
				calcPostCounterFac(db, modelRefined, a, a.postCntFac, fnin, fnout );
			}

			// 11 - Dump Simple output to .insres file.... 
			{
				// define twice to deal with "No default arg in lambda expression" warning....
				auto f  = [&] (double v) { const char *fm = "%.10lf";	static char b[64]; v=(v<0?0:v); sprintf(b, ((v == 0) || (fabs(v) >= 1e-4) ? fm : "%.10le"), v); return(string(b)); };
				auto f2 = [&](double v, const char *fm ) { 	static char b[64]; v =(v<0?0:v); sprintf(b, ((v == 0) || (fabs(v) >= 1e-4) ? fm : "%.10le"), v); return(string(b)); };
				inscomp::optInsight::parameters p = modelRefined.getparameters(), pe = modelRefinedE.getparameters();
				string so("");
				so += "Mode:";
				if (a.postParExpecQ) { so += "\tExpectation"; }
				else if (do_refine) { so += ( a.priorWeight <=0 ? "\tMaxLikelihood" : "\tMaxAPosteriori" ); }
				else so += "\tOther1";
				so += "\n";

				if ( (do_refine || do_posterior) && (deriv.alpha<0)) {	//deriv.alpha<0 checks to see sup stats have been calcualted yet (>=0), if not, calcualte them.
					modelInsightSup sup(modelRefined); sup.suppStatsSimple(a.nAlleles, deriv); };
				so += "ModelP:\tRho\t\tEta\t\tGamma\t\tDp\t\tPw\t\tAlpha\t\tTau\t\tBeta1\t\tBeta2\t\tBeta3\n";
				so += "Model:\t" + f(p.rho) + "\t" + f(p.eta) + "\t" + f(p.gam) + "\t" + f(deriv.dp) + "\t" + f(deriv.pw) + "\t" + f(deriv.alpha) + "\t" + f(deriv.tau) + "\t";
				so += f(p.beta.b1) + "\t" + f(p.beta.b2) + "\t" + f(p.beta.b3) + "\n";
				// if error values were calculated for MAP / ML, use them. But not if we are reporting expectations. User must find credible intervals in ppost.sum file
				if ((pe.rho>0) && (!a.postParExpecQ)) {
					so += "Errs:\t" + f(pe.rho) + "\t" + f(pe.eta) + "\t" + f(pe.gam) ;
					if (derivE.tau>=0) so += "\t" + f(derivE.dp) + "\t" + f(derivE.pw) + "\t" + f(derivE.alpha) + "\t" + f(derivE.tau);
					so += "\n";
				}
				if (a.dataLikelihood) so += "NLL:\t" + f2(stats.nll,"%.4lf") + "\t" + f2(stats.nlpp, "%.4lf") + "\n";
				double etime = (((double) clock() )- atime)/CLOCKS_PER_SEC;
				so += "Runtime:\t" + f2(etime,"%.4lf") + "\tsecs\n";
				string fon = a.inlocifile + ".insres"; std::ofstream fo(fon);
				fo << so << endl; fo.close();
				if (a.v>0) cerr << so << endl;
				// dump one liner for each input line processed.
				if (a.v<1) cout << strtime() << "\tSourcefile: " << a.inlocifile << "\t Runtime: " << f2(etime, "%.4lf") + "\tsecs" << endl;
			}

			if (server_mode) {
				// if we are in server mode, create a new file indicating that all processing is done... This allows externinal monitors to feed us a new task.
				if (donefile.length()>0) { std::ofstream fo(donefile); fo.close(); donefile = ""; };	// in the event of an error, notify client of completion..
			}
			if (!server_mode) { cout << strtime() << "\tDone." << endl; };

			server_request++;	// increment the server request count... 
		};
		// force a flush of all outstanding IO to disk.... then write ".done" semaphore
		if (donefile.length()>0) { butils::fastfile::fflush(); std::ofstream fo(donefile); fo.close(); donefile = ""; };	// in the event of an error, notify client of completion..
	};
	return 0;
}



bool posteriorDetail(insdb::DB &db, inscomp::optInsight::parameters &model, const Args &a, std::string &fnin, std::string &fnout) {
	using std::cout; using std::cerr; using std::endl;
	using butils::fastfile::SEOL; using butils::fastfile::SERR;
	using insdb::DB;
	using inscomp::WLocus; using inscomp::optInsight;

	FILE *fin = NULL, *fout = NULL;
	bool ok = true;
	if (a.v > 0) { cout << "\tPosteriorDetail(): Opening input file  " << fnin << endl; }; fin = fopen(fnin.c_str(), "rb"); if (!fin)  return false;
	if (a.v > 0) { cout << "\tPosteriorDetail(): Opening output file " << fnout << endl; }; fin = fopen(fnout.c_str(), "w");  if (!fout) return false;

	WLocus			lin;					// locus read from BED file
	DB::bigblock	cumulant;				// Saves intersection of block and locus.
	DB::cumestatus	cstat = DB::NewBlock;	// 3 values NewBlock, NewLoci, BlocksDone...
											//optInsight		opt;
	modelInsight::SelectionClass cls = modelInsight::SelectionClass::Mono;
	double zeqx_maj = 0.0, zeqx_min = 0.0;
	const char *chrom = NULL;
	uint32_t	pos = 0, poscount = 0;
	modelInsightSup::posteriordist sitepost;
	fprintf(fout, "%s\n", sitepost.header().c_str());
	modelInsight	tmp_mod;
	modelInsightSup tmp_mod_sup(tmp_mod);
	modelInsight::parameters tmp_model = model;	// save current model, replace "block" expectations with actual block values...
	db.locStreamInit(DB::AllData);
	bool ignoreweight = (a.maxSampleWeight<=1);
	while (ok && (lin.readln(fin, ignoreweight) == SEOL)) {
		if (a.v > 2) { cout << "Got Locus: "; lin.writeln(stdout); }
		// cumulate DB sites overlapping the input locus. This call alters lin as it consumes positions, untill they are all processed
		do {
			cstat = db.locStreamBig(lin, cumulant);	// a cumulant is a set of positions within a block...
			// if we have exhausted block, process any positions in the block and clear the cumulant
			if ( (cstat == DB::NewBlock || cstat == DB::BlocksDone) && (cumulant.head != NULL)) {
				tmp_model.block.theta = cumulant.head->theta; tmp_model.block.lambdaT = cumulant.head->lambda; cumulant.db = &db;
				tmp_mod.setparameters( tmp_model );
				chrom = &(cumulant.head->chrom[0]);
				for (vecDi ipos = 0; ipos < cumulant.sites.size(); ipos++) {
					if (cumulant.sites[ipos].ismono) {
						cls = modelInsight::SelectionClass::Mono;
						pos = (uint32_t)cumulant.sites[ipos].monopos;
						zeqx_maj = db.monoMajDouble(cumulant.sites[ipos].monoval); zeqx_min = 0.0;
					} else {
						const  insdb::Poly *s = cumulant.sites[ipos].polyptr;
						cls = (s->freq == 'L' ? modelInsight::SelectionClass::PolyL : modelInsight::SelectionClass::PolyH);
						zeqx_maj = s->majDouble(); zeqx_min = s->minDouble(); pos = s->st;
					}
					tmp_mod_sup.setparameters( tmp_mod.getparameters() );
					tmp_mod_sup.insight1Posterior(  a.nAlleles, cls, zeqx_maj, zeqx_min, sitepost);	// TODO maywe shoudl get nAlelles from somewhere else...
					fprintf(fout, "%s\t%lu\t%s\n", chrom, (unsigned long)pos, sitepost.toStr().c_str()); poscount++;
				}
			};
		} while (cstat == DB::NewBlock);	// if one locus spans multiple blocks, keep going till we've exhausted the locus, or all blocks...
		ok = (cstat == DB::NewLocus);		// fetch the next locus!
	};
	ok = (feof(fin) || (cstat == DB::BlocksDone));
	fclose(fin); fin = NULL; fclose(fout); fout = NULL;
	if (a.v > 0) { std::cout << "\tPosteriorDetail(): Processed " << poscount << " positions. Returning with status: " << ok << " ." << endl; };
	return ok;
}

bool calcPostAllele(insdb::DB &db, const modelInsight &model, const Args &a, std::string &fnin, std::string &fnout) {
	using std::cout; using std::cerr; using std::endl;
	using butils::fastfile::SEOL; using butils::fastfile::SERR;
	using insdb::DB;
	using inscomp::WLocus; using inscomp::optInsight;

	FILE *fin = NULL, *fout = NULL;
	bool ok = true;
	modelInsight tmp_model = model;
	modelInsight::parameters tmp_param = model.getparameters();

	if (a.v > 0) { cout << "\tPosteriorAllele(): Opening input file  " << fnin << endl; }; fin = fopen(fnin.c_str(), "rb"); if (!fin)  return false;
	if (a.v > 0) { cout << "\tPosteriorAllele(): Opening output file " << fnout << endl; }; fout = fopen(fnout.c_str(), "w");  if (!fout) return false;

	WLocus			lin;					// locus read from BED file
	DB::bigblock	cumulant;				// Saves intersection of block and locus.
	DB::cumestatus	cstat = DB::NewBlock;	// 3 values NewBlock, NewLoci, BlocksDone...
											//optInsight		opt;
	modelInsight::SelectionClass cls = modelInsight::SelectionClass::Mono;
	//double zeqx_maj = 0.0, zeqx_min = 0.0;
	const char *chrom = NULL;
	uint32_t	pos = 0, poscount = 0;
	modelInsightSup::posteriordist sitepost;
	butils::ezmatrix dist;
	char c_cls='N';

	fprintf(fout, "#Chrom\tPos\tSelCls\t   [S]   \tP(X=Xa)\tP(X=Xi)\tP(X=Xo)\tP(S0&Xa)\tP(S0&Xi)\tP(S0&Xo)\tP(S1&Xa)\tP(S1&Xi)\tP(S1&Xo)\tWeight\n");
	db.locStreamInit(DB::AllData);
	bool ignoreweight = (a.maxSampleWeight<=1);
	while (ok && (lin.readln(fin,ignoreweight) == SEOL)) {
		if (a.v > 2) { cout << "Got Locus: "; lin.writeln(stdout); }
		// cumulate DB sites overlapping the input locus. This call alters lin as it consumes positions, untill they are all processed
		do {
			cstat = db.locStreamBig(lin, cumulant);	// a cumulant is a set of positions within a block...
			if ((cstat == DB::NewBlock || cstat == DB::BlocksDone) && (cumulant.head != NULL) ) {
				tmp_param.block.theta = cumulant.head->theta; tmp_param.block.lambdaT = cumulant.head->lambda; cumulant.db = &db;
				tmp_model.setparameters(tmp_param);
				chrom = &(cumulant.head->chrom[0]);
				for (vecDi ipos = 0; ipos < cumulant.sites.size(); ipos++) {
					if (cumulant.sites[ipos].ismono) {
						cls = modelInsight::SelectionClass::Mono; c_cls = 'M';
						pos = (uint32_t)cumulant.sites[ipos].monopos;
						//zeqx_maj = db.monoMajDouble(cumulant.sites[ipos].monoval); zeqx_min = 0.0;
					} else {
						const  insdb::Poly *s = cumulant.sites[ipos].polyptr; c_cls = s->freq;
						cls = (s->freq == 'L' ? modelInsight::SelectionClass::PolyL : modelInsight::SelectionClass::PolyH);
						pos = s->st;
						//zeqx_maj = s->majDouble(); zeqx_min = s->minDouble(); pos = s->st;
					}
					tmp_model.posteriorJointAlleleAndS( a.nAlleles, cls, dist ); dist.norm1();
					// now print: chrom pos class [S] P(Z=Xmaj) P(Z=Xmin) P(Z=Xoth) full joint distribution of S,Z=X... (6 values)
					{	
						auto fmt = [&] (double d) { char b[32]; sprintf(b,(d!=0 && fabs(d)<1e-4?"%.2le":"%.6lf"),d); return( (std::string) b); };
						double s1=0;
						for (int j=0;j<dist.ncol();++j) s1 += dist.v(1,j);	// calculate posterior expectation value of S
						fprintf(fout, "%s\t%lu\t%c\t%s", chrom, (unsigned long)pos,c_cls,fmt(s1).c_str());					// header and [S]
						for (int j = 0; j<dist.ncol(); ++j) fprintf(fout, "\t%s", fmt(dist.v(0,j)+ dist.v(1, j)).c_str());	// Expected allele dsitrbution
						for (int i = 0; i<dist.nrow(); ++i) 
							for (int j = 0; j<dist.ncol(); ++j) fprintf(fout, "\t%s", fmt(dist.v(i, j) ).c_str());			// full joint of S & allele dist | Class
						fprintf(fout, "\t%.0lf\n",(double)lin.counts);
					};
				};
				cumulant.clear();
			};
		} while (cstat == DB::NewBlock);	// if one locus spans multiple blocks, keep going till we've exhausted the locus, or all blocks...
		ok = (cstat == DB::NewLocus);		// fetch the next locus!
	};
	ok = (feof(fin) || (cstat == DB::BlocksDone));
	fclose(fin); fin = NULL; fclose(fout); fout = NULL;
	if (a.v > 0) { std::cout << "\tPosteriorAlleleDist (): Processed " << poscount << " positions. Returning with status: " << ok << " ." << endl; };
	return ok;
}

// Processing a counterfactual may involve remapping the sufficient properties of a virtual observation (counterfactual) to some other
//	set of sufficient properties. (IE M,X!=Xmaj -> L,X=Xmin) to estimate the posterior contignent on the counterfactual. More advances
//	assessment may mix probabilities of known sets of sufficient properties to estimate counterfactual posterior.
bool calcCounterFacCore( modelInsightSup &M, int Mode, butils::ulong Nallele, modelInsight::SelectionClass modCls, const double &ZeqXma, const double &ZeqXmi, modelInsightSup::posteriordist &postZeqXma, modelInsightSup::posteriordist &postZeqXmi, modelInsightSup::posteriordist &postZeqXot) {
	typedef modelInsight::SelectionClass  SelectionClass;
	modelInsightSup::posteriordist post;			// needed when mixing posteriors...
	double ZeqXot = (1.0 - ZeqXma - ZeqXmi );

	//double nobs = Nallele / (Nallele + 1.0);	// change in probability because we have an additional observation...
	//auto obdn = [&] (double frac) { return ((frac * Nallele)    / (Nallele + 1)); };
	//auto obup = [&] (double frac) { return ((frac * Nallele +1) / (Nallele + 1)); };

	postZeqXma.clear(); postZeqXmi.clear(); postZeqXot.clear();

	// mode 1, consider primarly M->L... this is a very simple analysis, simply remap classes to accommodate counterfactual, do not apply expected class changes
	if (Mode==1) {
		if (modCls == SelectionClass::Mono) {
			M.insight1Posterior(Nallele+1, modCls, ZeqXma, ZeqXmi, postZeqXma);							// if class is mono and obs is Xma we need no remapping, ZeqXmi should == 0
			M.insight1Posterior(Nallele+1, SelectionClass::PolyL, ZeqXma, ZeqXot/3.0, postZeqXot);		// if class is mono and obs is X!=Xma we need to treat this as an observation of L, with Xmi<-X, 
			// The observed "Other" allele becomes the minor. In a more fulsome analysis we'd have to propigate the whole Z distrib to get the right counterfactual value for ZeqXmi
		} else {
			M.insight1Posterior(Nallele + 1, modCls, ZeqXma, ZeqXmi, postZeqXma);				// one more major allele, ok within model
			M.insight1Posterior(Nallele + 1, modCls, ZeqXma, ZeqXmi, postZeqXmi);				// one more minor allele, ok within model
			if (modCls == SelectionClass::PolyL) {												// an Xoth observation in a LF poly is special... forces S->0
				auto a = M.getparameters(), save=a;
				a.rho=0; M.setparameters(a);
				M.insight1Posterior(Nallele + 1, modCls, ZeqXma, ZeqXmi, postZeqXot);
				M.setparameters(save);
			} else  {
				M.insight1Posterior(Nallele + 1, modCls, ZeqXma, ZeqXmi, postZeqXot);
			};
		}
		return true;
	}

	// TODO mode 2, consider possibility that new observation will induce a change in class type...
	//	That is, the fractional possibility that adding an XMaj to an HF class will turn it into a LF class
	//	and the possibility adding an Xmin to a LF class turns it to HF. In expectation, just multiply by estimated probability of being at
	//	the appropriate threshold when the observation comes in.
	//	value.... 
	if (Mode==2) {
		if (modCls == SelectionClass::Mono) {
			M.insight1Posterior(Nallele + 1, modCls, ZeqXma, ZeqXmi, postZeqXma);									// if class is mono and obs is Xma we need no remapping, ZeqXmi should == 0
			M.insight1Posterior(Nallele + 1, SelectionClass::PolyL, ZeqXma, ZeqXot / 3.0, postZeqXot);				// if class is mono and obs is X!=Xma we need to treat this as an observation of L, with Xmi<-X, 																											// The observed "Other" allele becomes the minor. In a more fulsome analysis we'd have to propigate the whole Z distrib to get the right counterfactual value for ZeqXmi
		} else if ( modCls == SelectionClass::PolyL ) {
			M.insight1Posterior(Nallele + 1, modCls, ZeqXma, ZeqXmi, postZeqXma);				// Observe major allele, ok within model no class change
			{																					// Observe Minor allele....
				double pclasschange = 1.0 / floor( (.15 * Nallele - 1.0) );						// Chance we were at the upper bound of the LF class, before obs
				M.insight1Posterior(Nallele + 1, SelectionClass::PolyL, ZeqXma, ZeqXmi, postZeqXmi);		// Didnt change class
				M.insight1Posterior(Nallele + 1, SelectionClass::PolyH, ZeqXma, ZeqXmi, post);	// Now a HF class
				postZeqXmi.mult( 1.0 - pclasschange); postZeqXmi.accum( post, pclasschange);
			}
			{																					// Observe Xoth observation in a LF poly is special... forces S->0
				auto a = M.getparameters(), save = a;
				a.rho = 0; M.setparameters(a);
				M.insight1Posterior(Nallele + 1, modCls, ZeqXma, ZeqXmi, postZeqXot);
				M.setparameters(save);
			} 
		} else if (modCls == SelectionClass::PolyH ){
			{																						// Observe Major allele....
				double pclasschange = 1.0 / floor(((0.5 - .15) * Nallele));							// Chance we were at the lower bound of the HF class, before obs
				M.insight1Posterior(Nallele + 1, SelectionClass::PolyH, ZeqXma, ZeqXmi, postZeqXma);	// Didnt change class
				M.insight1Posterior(Nallele + 1, SelectionClass::PolyL, ZeqXma, ZeqXmi, post);			// Now an LF class
				postZeqXma.mult(1.0 - pclasschange); postZeqXma.accum(post, pclasschange);
			}
			M.insight1Posterior(Nallele + 1, modCls, ZeqXma, ZeqXmi, postZeqXmi);		// Observe minor allele, never changes HF class, may swap Xmi and Xma  but no influence est of S/A/W
			M.insight1Posterior(Nallele + 1, modCls, ZeqXma, ZeqXmi, postZeqXot);		// Observe OTHER allele, never changes HF class, may swap Xmi and Xoth but no influence est of S/A/W
		};
		return true;
	}

	// TODO - much more sophisticated analysis, considers polymorphism Class (U) as a latent variable informed by observed poly class Y, this was foreshadowed in Adam's
	//	Original INSIGHT writeup, but never made it into the code. Treat actual poly as distribution over U and calculate probabilities that a new observation
	//	alters distribution over U.... then take expested posterior given new observation.... not implimented... allows observation of XMaj to *increase* rho for LF poly.
	//	A quick and dirty approximation might be made by taking the expectation of a -*replacelemt*- of 1 actual observation with the counterfactual, then recomputing
	//	expected difference in posterior....
	if (Mode == 3) {
	};

	return false;
}

bool calcPostCounterFac(insdb::DB &db, const modelInsight &model, const Args &a, int Detail, std::string &fnin, std::string &fnout) {
	using std::cout; using std::cerr; using std::endl;
	using butils::fastfile::SEOL; using butils::fastfile::SERR;
	using insdb::DB;
	using inscomp::WLocus; using inscomp::optInsight;

	FILE *fin = NULL, *fout = NULL;
	bool ok = true;
	modelInsightSup tmp_model = model;
	modelInsight::parameters tmp_param = model.getparameters();

	if (a.v > 0) { cout << "\tPosteriorAllele(): Opening input file  " << fnin << endl; }; fin = fopen(fnin.c_str(), "rb"); if (!fin)  return false;
	if (a.v > 0) { cout << "\tPosteriorAllele(): Opening output file " << fnout << endl; }; fout = fopen(fnout.c_str(), "w");  if (!fout) return false;

	WLocus			lin;					// locus read from BED file
	DB::bigblock	cumulant;				// Saves intersection of block and locus.
	DB::cumestatus	cstat = DB::NewBlock;	// 3 values NewBlock, NewLoci, BlocksDone...
											//optInsight		opt;
	modelInsight::SelectionClass cls = modelInsight::SelectionClass::Mono;
	double zeqx_maj = 0.0, zeqx_min = 0.0;
	const char *chrom = NULL;
	uint32_t	pos = 0, poscount = 0;
	modelInsightSup::posteriordist sitepost, postZeqXma, postZeqXmi, postZeqXot;
	butils::ezmatrix dist;
	char c_cls = 'N';

	fprintf(fout, "#Chrom\tPos\tSelCls\t");
	fprintf(fout, "Ac:%s\t", sitepost.header(Detail).c_str());
	fprintf(fout, "Xa:%s\t", sitepost.header(Detail).c_str());
	fprintf(fout, "Xi:%s\t", sitepost.header(Detail).c_str());
	fprintf(fout, "Xo:%s\t", sitepost.header(Detail).c_str());
	fprintf(fout, "Weight\n");
	db.locStreamInit(DB::AllData);
	bool ignoreweight = (a.maxSampleWeight<=1);
	while (ok && (lin.readln(fin,ignoreweight) == SEOL)) {
		if (a.v > 2) { cout << "Got Locus: "; lin.writeln(stdout); }
		// cumulate DB sites overlapping the input locus. This call alters lin as it consumes positions, untill they are all processed
		do {
			cstat = db.locStreamBig(lin, cumulant);	// a cumulant is a set of positions within a block...
			if ((cstat == DB::NewBlock || cstat == DB::BlocksDone) && (cumulant.head != NULL)) {
				tmp_param.block.theta = cumulant.head->theta; tmp_param.block.lambdaT = cumulant.head->lambda; cumulant.db = &db;
				tmp_model.setparameters(tmp_param);
				chrom = &(cumulant.head->chrom[0]);
				for (vecDi ipos = 0; ipos < cumulant.sites.size(); ipos++) {
					if (cumulant.sites[ipos].ismono) {
						cls = modelInsight::SelectionClass::Mono; c_cls = 'M';
						pos = (uint32_t)cumulant.sites[ipos].monopos;
						zeqx_maj = db.monoMajDouble(cumulant.sites[ipos].monoval); zeqx_min = 0.0;
					} else {
						const  insdb::Poly *s = cumulant.sites[ipos].polyptr; c_cls = s->freq;
						cls = (s->freq == 'L' ? modelInsight::SelectionClass::PolyL : modelInsight::SelectionClass::PolyH);
						zeqx_maj = s->majDouble(); zeqx_min = s->minDouble(); pos = s->st;
					}
					// now print: chrom pos class [S] P(Z=Xmaj) P(Z=Xmin) P(Z=Xoth) full joint distribution of S,Z=X... (6 values)
					{	// we do this by positing an additional counterfactual observation of X, as one of Xmaj, Xmin or Xoth...
						// then observing the impct of the observation on expected distributioins, under the model. THIS IS ONLY AN ESTIMATE
						//	as certain coutnerfactual observations invalidate model invariants of class designations
						//	and require a more nuanced analysis to fully infer. P(H->L), P(L->H), P(M->L) with new observation?
						//	As an alternative counterfactual as -*replacing*- an observation, P(L->M)? Very experimental.
						std::string t;

						fprintf(fout, "%s\t%lu\t%c", chrom, (unsigned long)pos, c_cls);				// header 
		
						tmp_model.insight1Posterior(a.nAlleles, cls, zeqx_maj, zeqx_min, sitepost); // Actual observed data
						t = sitepost.toStr("%.6lf", Detail);
						fprintf(fout, "\t%s", t.c_str());											// header and posterior

						// counterfactual, observe Xmaj allele 
						calcCounterFacCore(tmp_model, a.postCntFacType, a.nAlleles, cls, zeqx_maj, zeqx_min, postZeqXma, postZeqXmi, postZeqXot);
						t = postZeqXma.toStr("%.6lf", Detail); fprintf(fout, "\t%s", t.c_str());
						// counterfactual, observe Xmin allele 
						t = postZeqXmi.toStr("%.6lf", Detail); fprintf(fout, "\t%s", t.c_str());
						// counterfactual, observe Xoth allele 
						t = postZeqXot.toStr("%.6lf", Detail); fprintf(fout, "\t%s", t.c_str());
						fprintf(fout, "\t%.0lf\n",(double)lin.counts);	// weight for this poistion
					};
				}
				cumulant.clear();
			};
		} while (cstat == DB::NewBlock);	// if one locus spans multiple blocks, keep going till we've exhausted the locus, or all blocks...
		ok = (cstat == DB::NewLocus);		// fetch the next locus!
	};
	ok = (feof(fin) || (cstat == DB::BlocksDone));
	fclose(fin); fin = NULL; fclose(fout); fout = NULL;
	if (a.v > 0) { std::cout << "\tPosteriorAlleleDist (): Processed " << poscount << " positions. Returning with status: " << ok << " ." << endl; };
	return ok;
}

inline bool compileFastBlock(insdb::DB::fastblock &cumulant, BlockSet_t &bs) {
	if (cumulant.items < 1) return true;

	auto &bm = bs.blockMaker;
	bm.clear();

	// get block header 
	auto ch = cumulant.head;
	bm.set(ch->chrom.c(), ch->st, ch->en, ch->lambda, ch->theta);

	// Get monomorphic positions 
	auto &cm = cumulant.mono;
	for (vecDi i = 0; i < cm.size(); i++)
		if (cm[i] > 0) bm.addMono((inscomp::AncPriIdx_t)i, cm[i]);
	// FYI, inscomp::AncPriIdx_t is small, typically a 1 ubyte type, or at most 2. CM has at most 255/64K entries
	// So its OK to downcast the larger vecDi, (an 8 byte ulong) to the smaller AncPriIdx_t (a 1 ubyte).

	// get polymorphic positions
	auto &cp = cumulant.poly;
	for (vecDi i = 0; i < cp.size(); i++) {
		insdb::Poly *tp = cp[i].poly;
		bm.addPolyX(tp->freq, tp->aPriMaj, tp->aPriMin, cp[i].count );
	}

	bm.compact();
	if (bm.hasSiteData()) bs.addBlock(bm);	// Base type of BlockMaker is Block_t, we we add to BlockMaler, then compact to the base.
	return true;
}

// {std::time_t ctt = std::time(0);  savetime = std::asctime(std::localtime(&ctt)); }

// The code to extract the blocks and crate an INS file is nearly the same, so rather
//	than write two routines, add a littl logic to dump INS files, if needed.
bool extractInputLoci(insdb::DB &db, const Args &a, BlockSet_t &bs, BlockSet_t &bsn, suppStats &stats, suppStats &statsn) {
	using std::cin;	using std::cout; using std::cerr; using std::endl; using std::string;
	using insdb::DB;
	using inscomp::WLocus;
	using butils::fastfile::SEOL;
	using butils::fastfile::SERR;
	// Used for generating .ins file

	// used for all purposes
	FILE *fin = NULL;
	std::string fnin = a.inlocifile + ".bed";

	if (a.v > 0) { std::cout << "\tOpening input file " << fnin << endl; }; fin = fopen(fnin.c_str(), "rb"); if (!fin) return false;
	if (a.v > 1) { std::cout << "\tFile opened for binary read " << endl; }; 

	db.locStreamInit(DB::AllData);
	DB::fastblock	cumulant;
	DB::cumestatus	cstat = DB::NewBlock;	// 3 values NewBlock, NewLoci, BlocksDone...
	bool ok = true; bs.clear(); db.copyAncPri(NULL, &bs.ancPri);
	WLocus lin;
	bool ignore_weights = (a.maxSampleWeight <= 1.0 );
	while (ok  && (lin.readln(fin,ignore_weights) == SEOL)) {
		if (a.v > 2 ) { cout << "Got Locus: "; lin.writeln(stdout); }
		// NOTE: Conversion from # of observations to fraction of bservation occurs in the accelerator.. not here
		//	Across 8 cell types seeing the same pattern 1 in exactly 1 cell type at each of 8 positions in a block,
		//	is the SAME as seeing that pattern in all 8 cell typbs, but only at one position. SumCounts is a sufficient statistic.
		stats.inLoci++; stats.inLociW += lin.counts;
		stats.inPos += (lin.en - lin.st); stats.inPosW += (lin.en - lin.st)*lin.counts;
		// cumulate DB sites overlapping the input locus. This call alters lin as it consumes positions, untill they are all processed
		while (ok && ((cstat = db.locStreamFast(lin, cumulant, lin.counts)) == DB::NewBlock)) {
			ok = compileFastBlock(cumulant, bs);	// compile the sites we have been accumulating and store the cumulant, signals start of new block.
			cumulant.clear();
		};
		//if (cstat == DB::BlocksDone) break; // this short circuit truncates valeus of stats.inpos... so remove it. --Brad
	};
	if (ok) { ok = compileFastBlock(cumulant, bs); cumulant.clear(); };	// compile block we were workign on when input data was exhausted. Compiling empty cumulant is not an error and produces no output.
	fclose(fin); fin = NULL;
	bs.blockSetValues(a.nAlleles, a.maxSampleWeight );
	stats.infPosM  = db.infSiteCount().M;  stats.infPosL  = db.infSiteCount().L;  stats.infPosH  = db.infSiteCount().H;
	stats.infPosWM = db.infSiteCount().wM; stats.infPosWL = db.infSiteCount().wL; stats.infPosWH = db.infSiteCount().wH;
	stats.maxWeight = a.maxSampleWeight;
	if (a.v > 0 && ok ) { std::cout << "\tSucessfully read input file. " << endl; };

	WLocus linn;
	bsn.clear();
	// fetch all Neutral Polys in DB blocks covered by BlockSet.
	db.locStreamInit(DB::NeutralPolysOnly); db.copyAncPri(NULL, &bsn.ancPri); cumulant.clear();
	for (uint32_t iblock = 0; (ok && (iblock < bs.blocks.size())); iblock++) {
		if (a.inscompat > 0 && !bs.blocks[iblock].hasPolys()) continue;
		// get the block extents
		bs.blocks[iblock].getHeader(&(linn.chrom), &(linn.st), &(linn.en));
		// cumulate neutral polys from the block. - NOTE return value should ALWAYS be DB::NewLocus bs.blocks is a subset of db.
		db.locStreamFast(linn, cumulant); ok = compileFastBlock(cumulant, bsn);	cumulant.clear();
	};
	if (ok) { ok = compileFastBlock(cumulant, bsn); cumulant.clear(); };	// compiling an empty cumulant is not an error and produces no output
	bsn.blockSetValues( bs.numAlleles);
	statsn.infPosM  = db.infSiteCount().M;  statsn.infPosL  = db.infSiteCount().L;  statsn.infPosH  = db.infSiteCount().H;
	statsn.infPosWM = db.infSiteCount().wM; statsn.infPosWL = db.infSiteCount().wL; statsn.infPosWH = db.infSiteCount().wH;
	statsn.maxWeight = a.maxSampleWeight;

	return ok;
}

bool estimateLambdaTheta(const BlockSet_t &bs, bool byBlock, double &LambdaT, double &Theta, double priorWeight, const inscomp::modelInsight *priorModel ) {
	bool ok = inscomp::modelInsightSup::estimateLambdaTheta(bs, priorWeight, priorModel, LambdaT, Theta, byBlock );
	return ok;
}

bool refineBeta(const BlockSet_t &bsn, modelInsight &model, double PriorWeight , const modelInsight *mPrior ) {
	// initial search values should be loaded into "model". "model" will be updated on return, if successful.
	inscomp::optBeta opter( model );
	bool ok = opter.optimize( bsn, model, PriorWeight, mPrior );
	return ok;
}

bool refineREG(const BlockSet_t &bs, modelInsight &model, modelInsight &modelE, const Args &a, modelInsightSup::derived_t &deriv, modelInsightSup::derived_t &derivE, double PriorWeight, modelInsight *PriorModel) {
	using std::cerr; using std::cout; using std::endl;
	using butils::ulong;
	//typedef inscomp::optInsight optInsight;

	optInsight opt( model );								// optimizer object
	optInsight::options opts;						// options for optimizer object
	optInsight::insightProvonance	prov;

	// Identify which parameters are to be refined - Error estimates for unrefined, or boundary-intersecting parameter values may be unreliable!!
	opts.fixRho = !a.rho.refine; opts.fixEta = !a.eta.refine; opts.fixGam = !a.gam.refine; opts.minEta = a.minEta; opts.v = (a.v<1?0:a.v-1);

	// get ML/MaP model parameters...
	bool ok = opt.optimize( bs, model, opts, PriorWeight, PriorModel, &prov);
	if (!ok) { // sometimes numerical issues near the optimum in a data set that closely matches the prior prevents converegence. Just use the best value we did find.
		cerr << "Insight2: ERROR: Rho/Eta/Gamma parameter refinement failed. Using last estimated values." << endl; /* return false; */ };
	if (a.v > 1) 	cerr << "Insight optimizer returns status " << ok << " and value " << model.getparameters().rho << endl;
	if (a.v > 0) 	cerr << "\tRho, eta, gamma estimated as " << model.getparameters().rho << "\t" << model.getparameters().eta << "\t" << model.getparameters().gam << endl;

	// Get uncertainties in model parameters
	modelInsight::parameters	p = model.getparameters();
	modelInsightSup c_sup(  model );
	ok = c_sup.estimateErr( bs, p, &opt.mAccel );	// provide ML accelerator...  more efficent then providing the priors and weight, again....
	if (!ok) {
		cerr << "Insight2: ERROR: Rho/Eta/Gamma parameter refinement error coaculation failed. " << endl; return false;	};
	modelE.setparameters(p);
	if (a.v > 0) 	cerr << "\tRho, eta, gamma err as " << p.rho << "\t" << p.eta << "\t" << p.gam << endl;

	// Get derived parameters adn their uncertainties....
	p = modelE.getparameters( );
	ok = c_sup.suppStats( bs, p, deriv, derivE, &opt.mAccel);
	if (!ok) {
		cerr << "Insight2: ERROR: Rho/Eta/Gamma derived values calculation error. " << endl; return false; 	};
	// TODO return derived parameters (pass in suterctures for this)

	return true;
}

bool dataStatistics(insdb::DB &db, const modelInsight &model, const Args &a, std::string &fnin, std::string &fnout ) {
	using std::cout; using std::cerr; using std::endl; using std::string;
	using butils::fastfile::SEOL; using butils::fastfile::SERR; using butils::ezmatrix; using butils::ulong;
	using insdb::DB;
	using inscomp::WLocus; //  using inscomp::optInsight;
	typedef modelInsight::SelectionClass  SelectionClass;
	typedef modelInsight::AncestralAllele AncestralAllele;

	// input file characteristics
	double pos_w=0, pos_u=0, loc_w=0, loc_u=0;

	// count the number of data elements overlapping loci
	ezmatrix counts_w(3), counts_u(3), counts_nw(3), counts_nu(3);		// 3x3 matrix for counts, weighted, unweighted
	double   countsn_w=0, countsn_u=0, countsn_nw=0, countsn_nu=0;		// Total counts, for normalization

	// Files for input and results..
	FILE *fin = NULL, *fouts = NULL, *foutd = NULL;
	if (a.v > 0) { cout << "\tPosteriorDetail(): Opening input file  " << fnin << endl; };  fin = fopen(fnin.c_str(), "rb");  if (!fin)  return false;
	if (a.dataSummary) {	string fn = fnout + ".sum";
		if (a.v > 0) { cout << "\tPosteriorDetail(): Opening output file " << fn << endl; }; fouts = fopen(fn.c_str(), "w");  if (!fouts) return false; }
	if (a.postDetail) {		string fn = fnout + ".det";
		if (a.v > 0) { cout << "\tPosteriorDetail(): Opening output file " << fn << endl; }; foutd = fopen(fn.c_str(), "w");  if (!foutd) return false; }

	WLocus			lin;					// locus read from BED file
	insdb::vBlock	block_list;				// list of blocks with valid TEST data in them, used to identify set of NEUTRAL blocks
	DB::bigblock	cumulant(&db);			// Saves intersection of block and locus. BigBlock maintains position information as well..
	DB::cumestatus	cstat = DB::NewBlock;	// 3 values NewBlock, NewLoci, BlocksDone...

	// assorted temporary variables
	bool ok = true;
	SelectionClass cls = SelectionClass::Mono;
	double zeqx_maj = 0.0, zeqx_min = 0.0, zeqx_oth = 0.0;
	const char *chrom = NULL;
	uint32_t	pos = 0, poscount = 0;

	//NOTE:    summaries of posteriors, like counts_w are adjusted for fratctional positions (ie 0.0+0.9+1.0+0.1+0.7 = 2.7 HF polys over 5 positions)
	//HOWEVER: summaries of INPUT data should not be adjusted (stats.inposw), 0 + 9 + 10 + 1 + 7 = weight of 27 with maxWeight=10.
	//			Sumaries of informative positions (countsn_w) bridge this. Raw values are more useful for debugging input, but weighted values are more
	//			useful for debuggign calculations.... calcs win.. so these summaries are normalized against a.maxSampleWeight.
	//			To Convert use WeightedSum/(a.MaxSampleWeight) to get output calibrated in # of positions (foreach position, W/maxW is the fraction
	//			of each position and is in [0,1]. lin.counts is unnormalzied, while  lincounts is the normalized version

	// ACTIVE Positions: Collect Posteriors and summary data about positions in active blocks!
	modelInsightSup mod_sup( model ); 
	modelInsight::parameters tmp_model = mod_sup.getparameters();
	modelInsightSup::posteriordist posdist, posdist_cum, posdist_cumw; posdist_cum.clear(); posdist_cumw.clear();
	db.locStreamInit(DB::AllData); block_list.clear();
	if (foutd!=NULL) fprintf(foutd,"Chrom\tPos\t%s\n",posdist.header().c_str() );
	std::vector<bool>	block_has_polys;	// used in insight1 compatibility mode, we only draw neutral blocks for blockdss with informative poly sites.
	bool ignore_weights = (a.maxSampleWeight<=1.0);
	while (ok && (lin.readln(fin,ignore_weights) == SEOL)) {
		if (a.v > 2) { cout << "Got Locus: "; lin.writeln(stdout); }
		if (ignore_weights) lin.counts = 1;
		// gather input file statistics at each new locus read....
		double lincounts= lin.counts / a.maxSampleWeight;
		pos_u += lin.en-lin.st; pos_w += (lin.en - lin.st)*(lin.counts); loc_u++; loc_w += (lin.counts);
		// cumulate DB sites overlapping the input locus. This call alters lin as it consumes positions, untill they are all processed
		do {
			cstat = db.locStreamBig(lin, cumulant);	// a cumulant is a set of positions within a block...
			if (cstat == DB::NewBlock ) { // process sites stoded in cumulant, than clear the cumulant...
				block_list.push_back( *(cumulant.head) ); block_has_polys.push_back( false );
				tmp_model.block.theta = cumulant.head->theta; tmp_model.block.lambdaT = cumulant.head->lambda; chrom = &(cumulant.head->chrom[0]);
				mod_sup.setparameters(tmp_model);			
				for (vecDi ipos = 0; ipos < cumulant.sites.size(); ipos++) {
					if (cumulant.sites[ipos].ismono) {
						cls = SelectionClass::Mono;
						pos = (uint32_t)cumulant.sites[ipos].monopos;
						zeqx_maj = db.monoMajDouble(cumulant.sites[ipos].monoval); zeqx_min = 0.0; zeqx_oth = 1.0 - zeqx_maj - zeqx_min;
						counts_u.v((int)cls, (int)AncestralAllele::Xmaj) += zeqx_maj;				counts_u.v((int)cls, (int)AncestralAllele::Xoth) += (zeqx_oth);					countsn_u++;
						counts_w.v((int)cls, (int)AncestralAllele::Xmaj) += (zeqx_maj*lincounts);	counts_w.v((int)cls, (int)AncestralAllele::Xoth) += (zeqx_oth)*(lincounts);		countsn_w+= lincounts;
					}
					else {
						const insdb::Poly *s = cumulant.sites[ipos].polyptr;
						cls = (s->freq == 'L' ? SelectionClass::PolyL : SelectionClass::PolyH);
						zeqx_maj = s->majDouble(); zeqx_min = s->minDouble(); zeqx_oth = 1.0 - zeqx_maj - zeqx_min; pos = s->st;
						counts_u.v((int)cls, (int)AncestralAllele::Xmaj) += zeqx_maj;				counts_u.v((int)cls, (int)AncestralAllele::Xmin) += zeqx_min;					counts_u.v((int)cls, (int)AncestralAllele::Xoth) += (zeqx_oth);					countsn_u++;
						counts_w.v((int)cls, (int)AncestralAllele::Xmaj) += (zeqx_maj*lincounts);	counts_w.v((int)cls, (int)AncestralAllele::Xmin) += (zeqx_min*lincounts);		counts_w.v((int)cls, (int)AncestralAllele::Xoth) += (zeqx_oth)*(lincounts);		countsn_w += lincounts;
						block_has_polys[block_has_polys.size()-1] = true;
					}
					mod_sup.insight1Posterior(a.nAlleles, cls, zeqx_maj, zeqx_min, posdist);  posdist.cCnt = lincounts;
					if (foutd != NULL ) {
						string t = posdist.toStr();
						fprintf(foutd, "%s\t%lu\t%s\n", chrom, (unsigned long)pos, t.c_str()); };
						posdist_cumw.accum(posdist, lincounts); posdist.cCnt = 1.0; posdist_cum.accum(posdist);
				};	// done processing all the sites in the cumulant
				cumulant.clear();
			};
		} while (cstat == DB::NewBlock);	// if one locus spans multiple blocks, keep going till we've exhausted the locus, or all blocks...
		ok = (cstat == DB::NewLocus);		// fetch the next locus!
	};
	ok = (feof(fin) || (cstat == DB::BlocksDone));
	fclose(fin); fin = NULL; if (foutd != NULL) { fclose( foutd ); foutd = NULL; };

	if (fouts != NULL) {
		// NEUTRAL Positions : Collect Block / Neutral Poly summary Data!
		//  loop over the blocks associated with active positions, but look at all neutral polys in each block..
		//unsigned long c_m, c_l, c_h;
		db.locStreamInit(DB::NeutralPolysOnly ); cstat = DB::NewBlock; WLocus tmp_loc; cumulant.clear();
		for (ulong iblock=0; iblock< block_list.size(); ++iblock ) {
			insdb::Block &bl = block_list[iblock];
			if (a.inscompat > 0 && block_has_polys[iblock]) continue;	// Ins1 COmpatibility mode, only draw neutrals from blocks with informative poly sites...
			tmp_loc.chrom = bl.chrom; tmp_loc.st = bl.st; tmp_loc.en = bl.en; tmp_loc.counts = 1;
			cumulant.clear(); db.locStreamBig(tmp_loc, cumulant); //cumulant.numPos( c_m, c_l, c_h );
			double szi = 1.0 / cumulant.sites.size();	// weight by block vs position...
			for (ulong ipos =0; ipos < cumulant.sites.size(); ++ipos ) {
				assert(!cumulant.sites[ipos].ismono );
				const insdb::Poly *s = cumulant.sites[ipos].polyptr; 
				cls = (s->freq == 'L' ? SelectionClass::PolyL : SelectionClass::PolyH);
				zeqx_maj = s->majDouble(); zeqx_min = s->minDouble(); zeqx_oth = 1.0 - zeqx_maj - zeqx_min; pos = s->st;
				counts_nu.v((int)cls, (int)AncestralAllele::Xmaj) += zeqx_maj;				counts_nu.v((int)cls, (int)AncestralAllele::Xmin) += zeqx_min;						counts_nu.v((int)cls, (int)AncestralAllele::Xoth) += (zeqx_oth);		countsn_nu++;
				counts_nw.v((int)cls, (int)AncestralAllele::Xmaj) += (zeqx_maj*szi);		counts_nw.v((int)cls, (int)AncestralAllele::Xmin) += (zeqx_min*szi);				counts_nw.v((int)cls, (int)AncestralAllele::Xoth) += (zeqx_oth)*(szi);	countsn_nw += szi;
			};
		};

		// DUMP the summary Data!
		// 1-Requested Position Data
		fprintf(fouts, "REQUESTED POSITIONS\n");
		fprintf(fouts, "\t        \t     #Loci     \t#Positions\n");
		fprintf(fouts, "\tRaw     \t%15.2lf\t%15.2lf\n",(double)loc_w, (double)pos_w);
		fprintf(fouts, "\tWeighted\t%15.2lf\t%15.2lf\n",(double)loc_u, (double)pos_u);

		// 2-Active  & Neutral Position Summary, for each of weighted adn unweighted data... 
		//	(weighted neutral data is normalized by BLOCK, unweighted just sums positions actross all blocks
		//	(Active data is weighted by locus weight, unweighted active data simply counts each position once.
		//	the most representative data for Neutrals is clearly Weighted. Personally, I'd argue this for Active data too...
		std::vector<string> sel_nam = { "Mono", "PolyL", "PolyH" }, anc_nam = { "ZeqXmaj", "ZeqXmin", "ZeqXoth" };
		fprintf(fouts, "\nPOSITIONS WITH DATA\n");
		for (int b_neut = 0; b_neut < 2; b_neut++) {
			ezmatrix &wgt=(b_neut ? counts_nw  : counts_w ), &unw = (b_neut ? counts_nu  : counts_u);
			double   wgtn=(b_neut ? countsn_nw : countsn_w), unwn = (b_neut ? countsn_nu : countsn_u);
			string t_set = (b_neut ? "NEUTRAL" : "TESTED");
			for (int b_weighted = 0; b_weighted < 2; b_weighted++) {
				ezmatrix &val = (b_weighted ? wgt   : unw ); double   valn = (b_weighted ? wgtn : unwn );
				fprintf(fouts, "\n%s - %s\n", t_set.c_str(), (b_weighted ? "WEIGHTED" : "UNWEIGHTED"));
				for ( ulong isel = (ulong) SelectionClass::first; isel <= (ulong) SelectionClass::last; isel++ ) {
					if (isel== (ulong)SelectionClass::first) 
						fprintf(fouts, "\t\t%s\t\t%s\t\t%s", anc_nam[0].c_str(), anc_nam[1].c_str(), anc_nam[2].c_str() );
					fprintf( fouts , "\n%s",sel_nam[isel].c_str());
					for (ulong ianc = (ulong)AncestralAllele::first; ianc <= (ulong)AncestralAllele::last; ianc++) {
						fprintf(fouts,"\t%15.2lf",(double) val.v(isel,ianc)); }
				};
				fprintf(fouts, "\nTotal: %18.4lf\n", (double)valn);
			};
		};

		if (a.postSummary) {
			// 3-Posterior Summary, simple sum of posterior countes, unweighted...
			fprintf(fouts, "\nPOSTERIOR SUMMARY - unweighted \n");
			fprintf(fouts, "\t\t%s\n",posdist_cum.header().c_str());
			fprintf(fouts, "\t%s\n", posdist_cum.toStr("%18.4lf").c_str());
			fprintf(fouts, "\nPOSTERIOR SUMMARY - weighted \n");
			fprintf(fouts, "\t\t%s\n", posdist_cumw.header().c_str());
			fprintf(fouts, "\t%s\n", posdist_cumw.toStr("%18.4lf").c_str());
		};

		fclose(fouts); fouts = NULL;
	};

	if (a.v > 0) { std::cout << "\tPosteriorDetail(): Processed " << poscount << " positions. Returning with status: " << ok << " ." << endl; };
	return ok;
}

bool calcLikelihood(const modelInsight &model, const BlockSet_t &bs, double &nll, double &nlpp, double PriorWeight, const modelInsight *PriorModel ) {
	inscomp::modelInsightAccelerator ma( model );
	ma.initialize( bs, PriorWeight, PriorModel ); nlpp = ma.calcNLL(); nll = nlpp;
	if ( PriorWeight != 0 && PriorModel != NULL ) {
		ma.initialize(bs); nll = ma.calcNLL();
	};
	return true;
}

bool calcParamPost(const BlockSet_t &bs, const modelInsight &model, double PriorWeight, const modelInsight &PriorModel, const Args &a, modelInsightPosterior::postSet &ppost, int GridRes) {
	modelInsightPosterior pos( model );
	ppost.clear();
 	pos.calcParamPost( bs, PriorWeight, PriorModel, ppost, GridRes, .05, a.v, a.nThreads );	// TODO set default grid, 20 migth be better... last arg is Credal Thresh .01? .001?
	// TODO dump this to file spcified in Args - A...
	return true;
}

bool makeIns(insdb::DB &db, const Args &a, modelInsight::parameters &model, const std::string &fnout, const std::string &fnoutb ) {
	using std::cin;	using std::cout; using std::cerr; using std::endl; using std::string;
	using insdb::DB;
	using inscomp::WLocus;
	using butils::fastfile::SEOL; using butils::fastfile::SERR; using butils::ulong;

	std::vector<DB::bigblock> bs;	// blockset, really we only keep headers, so this is not to big...

	FILE *fout = NULL, *foutb = NULL;
	FILE *fin = NULL;
	std::string fnin = a.inlocifile + ".bed";

	if (a.v > 0) { std::cout << "\tOpening input file " << fnin << endl; }; fin = fopen(fnin.c_str(), "rb"); if (!fin) return false;
	if (a.v > 1) { std::cout << "\tFile opened for binary read " << endl; };

	fout = fopen(fnout.c_str(), "wb"); foutb = fopen(fnoutb.c_str(), "wb");
	if ((fout == NULL) || (foutb == NULL)) {
		if (fout) fclose(fout); if (foutb == NULL)  fclose(foutb);
		cerr << "INSIGHT2: extractLoci, unable to open .ins/.beta.ins files based at : \n\t" << a.insfile << "\nSkipping .ins file generation\n" << endl;
		fout = foutb = NULL; return false;
	}
	if (a.v > 1) { std::cout << "\tInsight1: output files opened for binary write.\n\t" << fnout << "\n\t" << fnout << "\n" << endl; };

	db.locStreamInit(DB::AllData); bs.clear();
	DB::bigblock cumulant; cumulant.db = &db;
	DB::cumestatus	cstat = DB::NewBlock;	// 3 values NewBlock, NewLoci, BlocksDone...
	bool ok = true, ignoreweight=true; 
	WLocus lin;
	fprintf(fout, "samples %d\n",(int)a.nAlleles);
	while (ok && (lin.readln(fin,ignoreweight) == SEOL)) {
		if (a.v > 2) { cout << "Got Locus: "; lin.writeln(stdout); }
		// cumulate DB sites overlapping the input locus. This call alters lin as it consumes positions, untill they are all processed
		while (ok && ((cstat = db.locStreamBig(lin, cumulant)) == DB::NewBlock)) {
			cumulant.dump(fout);																// dump prints nothing if there are no sites
			// if (cumulant.sites.size() >= 0) { cumulant.clear(true); bs.push_back(cumulant); };	// just in case we need poly counts for compatibility, save full block
			if (cumulant.sites.size() > 0) { bs.push_back(cumulant); };	// just in case we need poly counts for compatibility, save full block...
			cumulant.clear(); cumulant.db = &db;
		};
		if (cstat == DB::BlocksDone) break;
	};
	if (ok) {
		cumulant.dump(fout);																// dump prints nothing if there are no sites
		// if (cumulant.sites.size() >= 0) { cumulant.clear(true); bs.push_back(cumulant); };	// if >0 sites, save the header. But the header only.
		if (cumulant.sites.size() >  0) {					    bs.push_back(cumulant); };	// if >0 sites, save the header. But the header only.
		cumulant.clear(); cumulant.db = &db;
	};
	// if we calculated the betas, dump them
	if (a.betas.refine) fprintf(fout,"beta\t%.10lf\t%.10lf\t%.10lf\n",(double)model.beta.b1, (double) model.beta.b2, (double)model.beta.b3);
	fclose(fin); fin = NULL; fclose(fout); fout = NULL;
	if (a.v > 1) { std::cout << "\tWrote Insight1 .ins file, attempting .beta.ins file.\n"  << endl; };

	ulong countH = 0, countL = 0, m = 0, l = 0, h = 0;
	db.locStreamInit(DB::NeutralPolysOnly); 
	fprintf(foutb, "samples %d\n", (int)a.nAlleles);
	cumulant.clear(); cumulant.db = &db;  cstat = DB::NewBlock;	ok = true; lin.clear(); // 3 values NewBlock, NewLoci, BlocksDone...
	for ( ulong iblock = 0; iblock < bs.size(); iblock++ ) {
		DB::bigblock &b = bs[iblock];	// alias
		if (a.inscompat > 0 ) {
			b.numPos(m, l, h); if ((h+l)<1) continue;	// compabilility mode: use only blocks with polymorphic informative sites to calculate beta...
		};
		lin.chrom = b.head->chrom; lin.st = b.head->st; lin.en = b.head->en; lin.counts = 1;
		db.locStreamBig(lin, cumulant);	// we are only using known blocks, so don't need to worry about status, or block boundaries.
		cumulant.numPos(m, l, h); countH += h; countL += l;
		if (l + h > 0) { cumulant.dump(foutb, true); };
		cumulant.clear(); cumulant.db = &db;
	};
	printf("Beta2:\t%.12lf\t%.0lf\t%.0lf\n", ((double)countH) / (countH + countL), (double) countH, (double) countL);
	fclose(foutb); fout = NULL;
	if (a.v > 1) { std::cout << "\tWrote Insight1 .ins file, attempting .beta.ins file.\n" << endl; };

	return ok;
}


#ifdef OLD_CODE
void suppStats::gather(const BlockSet_t & bs) {
	// gather simple stats on input database ovrlap with input loci....
	using butils::ulong;



	return;
}

// create a grid with good overall coverage, and more detail near a focus....
struct gridder {
	struct focus { double foc, disp; };
	inline static double ratToFrac(double rat) { return rat / (rat + 1.0); };		// convert between ratio (a/b) [unbounded] and fraction (a/(a+b)) [0,1]
	inline static double fracToRat(double frac) { return (frac == 1.0 ? std::numeric_limits<double>::max() : frac / (1.0 - frac)); };
	inline static double elemSize(butils::ulong I, const vecD &grid) {
		double r = 0; auto imax = grid.size(); if (imax <= 0 || I >= imax) return r; --imax;
		if (I > 0)		r += 0.5*(grid[I] - grid[I - 1]);
		if (I < imax)	r += 0.5*(grid[I + 1] - grid[I]);
		return r;
	}
	static bool makeGrid(double gmin, double gmax, const std::vector<focus> &gfoc, unsigned long resolution, vecD &grid) {
		double delta = (gmax - gmin) / resolution; grid.clear();
		for (unsigned long i = 0; i <= resolution; i++) grid.push_back(gmin + i * delta);									// add coarse grid
		for (auto ilev = 0; ilev < gfoc.size(); ilev++) {															// loop over each focus
			double foc = gfoc[ilev].foc, disp = gfoc[ilev].disp;
			double v = 0, base = foc - 0.5*disp;  delta = disp / resolution;
			for (unsigned long i = 0; i <= resolution; i++) v = base + i * delta;  if (v >= gmin && v <= gmax) grid.push_back(v);	// add fine grids
		}
		std::sort(grid.begin(), grid.end()); grid.erase(std::unique(grid.begin(), grid.end()), grid.end());			// sort and remove duplicates
		return true;
	};
};


bool posteriorDetail(insdb::DB &db, inscomp::optInsight::parameters &model, const Args &a, std::string &fnin, std::string &fnout) {
	using std::cout; using std::cerr; using std::endl;
	using butils::fastfile::SEOL; using butils::fastfile::SERR;
	using insdb::DB;
	using inscomp::WLocus; using inscomp::optInsight;

	FILE *fin = NULL, *fout = NULL;
	bool ok = true;
	if (a.v > 0) { cout << "\tPosteriorDetail(): Opening input file  " << fnin << endl; }; fin = fopen(fnin.c_str(), "rb"); if (!fin)  return false;
	if (a.v > 0) { cout << "\tPosteriorDetail(): Opening output file " << fnout << endl; }; fin = fopen(fnout.c_str(), "w");  if (!fout) return false;

	WLocus			lin;					// locus read from BED file
	DB::bigblock	cumulant;				// Saves intersection of block and locus.
	DB::cumestatus	cstat = DB::NewBlock;	// 3 values NewBlock, NewLoci, BlocksDone...
											//optInsight		opt;
	optInsight::SelectionClass cls = optInsight::SelectionClass::Mono;
	double zeqx_maj = 0.0, zeqx_min = 0.0;
	const char *chrom = NULL;
	uint32_t	pos = 0, poscount = 0;
	optInsight::posteriordist sitepost;
	fprintf(fout, "%s\n", sitepost.header().c_str());

	optInsight::parameters tmp_model = model;	// save current model, replace "block" expectations with actual block values...
	db.locStreamInit(DB::AllData);
	while (ok && (lin.readln(fin) == SEOL)) {
		if (a.v > 2) { cout << "Got Locus: "; lin.writeln(stdout); }
		// cumulate DB sites overlapping the input locus. This call alters lin as it consumes positions, untill they are all processed
		do {
			cumulant.clear(); cstat = db.locStreamBig(lin, cumulant);	// a cumulant is a set of positions within a block...
																		// TODO need to multiply by wattersons A?
			tmp_model.block.theta = cumulant.head->theta; tmp_model.block.lambdaT = cumulant.head->lambda; cumulant.db = &db;
			chrom = &(cumulant.head->chrom[0]);
			for (auto ipos = 0; ipos < cumulant.sites.size(); ipos++) {
				if (cumulant.sites[ipos].ismono) {
					cls = optInsight::SelectionClass::Mono;
					pos = (uint32_t)cumulant.sites[ipos].monopos;
					zeqx_maj = db.monoMajDouble(cumulant.sites[ipos].monoval); zeqx_min = 0.0;
				} else {
					const  insdb::Poly *s = cumulant.sites[ipos].polyptr;
					cls = (s->freq == 'L' ? optInsight::SelectionClass::PolyL : optInsight::SelectionClass::PolyH);
					zeqx_maj = s->majDouble(); zeqx_min = s->minDouble(); pos = s->st;
				}
				optInsight::posteriorIndividual(tmp_model, a.nAlleles, cls, zeqx_maj, zeqx_min, sitepost);	// TODO maywe shoudl get nAlelles from somewhere else...
				fprintf(fout, "%s\t%lu\t%s\n", chrom, (unsigned long)pos, sitepost.toStr().c_str()); poscount++;
			}
		} while (cstat == DB::NewBlock);	// if one locus spans multiple blocks, keep going till we've exhausted the locus, or all blocks...
		ok = (cstat == DB::NewLocus);		// fetch the next locus!
	};
	ok = (feof(fin) || (cstat == DB::BlocksDone));
	fclose(fin); fin = NULL; fclose(fout); fout = NULL;
	if (a.v > 0) { std::cout << "\tPosteriorDetail(): Processed " << poscount << " positions. Returning with status: " << ok << " ." << endl; };
	return ok;
};


// P(M) -> P(priMod|priD) -> P(priM|priE)^priN. Under a uniform hyper-prior
//	P(priM|priE) \porp P(priE|priM). So we calculate P(M) from the values of 
//	P(priE|priM)^priN, normalized to 1 across a computationally feasible sample 
//	grid that covers the space of plausable parameter values.....
struct paramPosteriors {
	paramPosteriors() { clear(); };
	typedef butils::ulong		ulong;
	typedef std::vector<double>	vecD;
	typedef std::string			string;
	struct mSample {
		double val, size, prior, nll, joint, post, densPr, densPost;
		mSample() { clear(); };
		void clear() { val = size = prior = nll = joint = post = densPr = densPost = 0.0; return; };
		const char *fmt(const double &val, const char *f = "%.6lf") { static char buf[32]; sprintf(buf, (val < 1e-3 ? "%7.4e" : f), val); return(buf); }
		static string strHead() { string s = "Rho\tPrior\tDataNLL\tJointNLL\tPostr\tCellSz\tPriorDens\tPostDens"; return(s); };
		string toStr(const char *f = "%.6lf") { return(string("") + fmt(val, f) + "\t" + fmt(prior, f) + "\t" + fmt(nll, f) + "\t" + fmt(joint, f) + "\t" + fmt(post, f) + "\t" + fmt(size, f) + "\t" + fmt(densPr, f) + "\t" + fmt(densPost, f)); };
	};
	struct paramPosterior {
		double rho, eta, gam;
		double prior1, nll, cellSz, priorU, priorN, jointN, postN, postD;		// 1 obs of expected prior, data log likelihood, normalized posterior probability (and its density).
		paramPosterior() { clear(); };
		const char *fmt(const double &val, const char *f) { static char buf[32]; sprintf(buf, (val < 1e-3 ? "%7.4e" : f), val); return(buf); }
		static string strHead() { return("Rho\tEta\tGamma\tPriorN\tDataNLL\tJointN\tPostN\tCellSize\tPrior1\tPriorU\tPostD"); };
		string toStr(const char *f = "%.6lf") { return(string("") + fmt(rho, f) + "\t" + fmt(eta, f) + "\t" + fmt(gam, f) + "\t" + fmt(priorN, f) + "\t" + fmt(nll, f) + "\t" + fmt(jointN, f) + "\t" + fmt(postN, f) + "\t" + fmt(cellSz, f) + "\t" + fmt(prior1, f) + "\t" + fmt(priorU, f) + "\t" + fmt(postD, f)); };
		void clear() { rho = eta = gam = prior1 = nll = cellSz = priorU = priorN = jointN = postN = postD = 0.0; return; };
	};
	ulong	gridDensity;						// number of probability samples
	vecD	gridRho, gridEta, gridGam;			// grid of sample points
	std::vector<paramPosterior> posts;			// set of posterior records, one for each sample point
	std::vector<mSample>		mRho;			// marginalized values fro Rho...
	double expRho, expEta, expGam;				// posterior expected values....
	double priorWeight;							// PRIOR: priN number of pseudocounts (observations of proExpec) used in calculating prior,
	butils::ezmatrix	priExpec;				// PRIOR: priD distribution of observations used to genrate prior P(M) \porp P(priD|PriM)^priN
	inscomp::optInsight::parameters priModel;	// PRIOR: priM model parameters used to generate prior. Generally drawn from whole population ML/expectation values.
	const char *fmt(const double &val, const char *f = "%.6lf") { static char buf[32]; sprintf(buf, f, val); return(buf); }
	void clear() {
		gridDensity = 0; gridRho.clear(); gridEta.clear(); gridGam.clear(); posts.clear(); mRho.clear(); expRho = expEta = expGam = 0.0; priorWeight = 0; priExpec.clear(); priModel.clear(); return;
	};
	string toStrExp(const char *f = "%.6lf") { return(string("") + "Expectation-Rho,Eta,Gamma:\t" + fmt(expRho, f) + "\t" + fmt(expEta, f) + "\t" + fmt(expGam, f)); };
	string toStrRhoDist(const char *f = "%.6lf") {
		string s = "ExpectationOverEtaGamma-:\t" + mSample::strHead() + "\n";
		for (auto i = 0; i < mRho.size(); ++i) s += mRho[i].toStr() + "\n";
		return(s);
	}
	string toStrFull(const char *f = "%.6lf") {
		string s = "";
		s += "GridDensity:\t" + (string)fmt(gridDensity, "%.0lf") + "\n";
		s += "GridRho:"; for (auto i = 0; i < gridRho.size(); i++) s += "\t" + (string)fmt(gridRho[i]); s += "\n";
		s += "GridEta:"; for (auto i = 0; i < gridEta.size(); i++) s += "\t" + (string)fmt(gridEta[i]); s += "\n";
		s += "GridGam:"; for (auto i = 0; i < gridGam.size(); i++) s += "\t" + (string)fmt(gridGam[i]); s += "\n\n";
		s += "PriorWeight:\t"; s += fmt(priorWeight); s += "\n\n";
		s += "PriorExpec:\n"; s += priExpec.toStr(); s += "\n\n";
		s += "PriorModel:\n"; s += priModel.toStr(); s += "\n\n";
		// this is the long table... put it last...
		s += "SampleGrid:\t"; s += paramPosterior::strHead() + "\n"; for (auto i = 0; i < posts.size(); i++) s += posts[i].toStr() + "\n"; s += "\n";
		return(s);
	};
};

bool calcLikelihood(const modelInsight &model, const BlockSet_t &bs, double &nll, int V) {
	nll = model.likelihoodDataLog(bs);
	using butils::ulong;
	using inscomp::optInsight;
	modelInsight::parameters params = model.getparameters();
	modelInsight mod_tmp; mod_tmp.setparameters(model.getparameters());
	double prob, zeqxmaj, zeqxmin;
	inscomp::AncPriIdx_t izeqxmaj, izeqxmin;
	inscomp::PatternCounts_t nobs;
	nll = 0;
	for (auto iblock = 0; iblock < bs.blocks.size(); iblock++) {
		auto &b = bs.blocks[iblock];
		// TODO: multipley by watterson A here?
		b.getHeader(NULL, NULL, NULL, &params.block.lambdaT, &params.block.lambdaT);
		zeqxmin = 0; izeqxmin = 0;
		for (ulong i = 0; i < b.monoNumPat(); i++) {
			prob = zeqxmaj = 0; nobs = 0;
			b.monoVal(i, izeqxmaj, nobs); zeqxmaj = 0;  if (izeqxmaj != 0) zeqxmaj = *(bs.ancPri[izeqxmaj]);
			prob = model.likelihood(bs.numAlleles, modelInsight::SelectionClass::Mono, zeqxmaj, zeqxmin, false);
			if (prob>0) nll += nobs * log(prob);
		};
		for (ulong i = 0; i < b.polyLNumPat(); i++) {
			prob = zeqxmaj = zeqxmin = 0; nobs = 0;
			b.polyLVal(i, izeqxmaj, izeqxmin, nobs); zeqxmaj = zeqxmin = 0;
			if (izeqxmaj != 0) zeqxmaj = *(bs.ancPri[izeqxmaj]); if (izeqxmin != 0) zeqxmin = *(bs.ancPri[izeqxmin]);
			prob = optInsight::likelihoodExpected(params, bs.numAlleles, optInsight::SelectionClass::PolyL, zeqxmaj, zeqxmin, false);
			if (prob>0) nll += nobs * log(prob);
		};
		for (ulong i = 0; i < b.polyHNumPat(); i++) {
			prob = zeqxmaj = zeqxmin = 0; nobs = 0;
			b.polyHVal(i, izeqxmaj, izeqxmin, nobs); zeqxmaj = zeqxmin = 0;
			if (izeqxmaj != 0) zeqxmaj = *(bs.ancPri[izeqxmaj]); if (izeqxmin != 0) zeqxmin = *(bs.ancPri[izeqxmin]);
			prob = optInsight::likelihoodExpected(params, bs.numAlleles, optInsight::SelectionClass::PolyH, zeqxmaj, zeqxmin, false);
			if (prob>0) nll += nobs * log(prob);
		};
	};
	return true;
};
#endif
