// covtree.cpp : Defines the entry point for the console application.

#define COVTREE_VER_STR "00.00.12" " (" __DATE__ " " __TIME__ ")"

// 00.00.01 - 2016-08-23 0355 ET. Compiles, runs, and iterates. All the basics seem here. Time to committ to alpha
// 00.00.02 - 2016-08-24 XXXX ET. Added linux compatible changes C++1y (2013 std). Added -restart option, fixed prior
//			parsing. Variety of small fixes.
// 00.00.03 - 2016-08-25 1855 Fix some IO bugs, add -sarg serialized argument support to command line.  
// 00.00.04 - 2016-08-27 XXXX Add more checking to asynch read of insight .model files. Add -dwrk to assign local working directory for cluster usage. Fix problem when there are no nonmonotonic splits to evalaute.
// 00.00.04a- 2016-08-29 
// 00.00.04b- fix bug in non-persistable arguments like dWrk, that must be passsed to deseralized arrg,
// 00.00.05 - add multithreading for Insight2.
// 00.00.06 - minor bugfixes ad support for expected information gain.
// 00.00.07 - Allow imposition of global Rho consistency by using parrental RHO value from priors, rather than locally consistent rho value. 
//			- Implimented via the -implprho arg. Otherwise get Parental Insight2 Rho estiamtes from Priors if provided >0 and <1. Use implied if no valid prior provided.
// 00.00.07a- New -idb flag allows increased cluster performance by loading Insight2 DB from local storage.
// 00.00.08 - 2016-09-15 Big update. Improved expectation handeling. Calculate prove rho1 > rho2 from rho_hat & rho_delta's (w/ bonferroni correction). 
//				Compute conditional expected inf gain | cond on rho1>rho2. Lots of cleanup to output, added err1 & err2 terms (prob rho1>rho2, prob rho1 > (rho2 + ((rhohat1-rhohat2)/2)).).
//				Errors fixed, extra outPart columns..
// 00.00.09 - 2016-09-17 Small but profound update. Changed default INFo measurement to InsightNLL, rather than [S]. many changes to support this. Added fields to output format.
// 00.00.10 - 2016-09-28 Several very minor updates. Documentation.
// 00.00.11 - 2016-12-20 Handeling of improper partitions (those where one coset has no genomic positions). Passing an empty bed file (or one identical to parent) to Insight2 can result in problems.
//				Does NOT effect standard covaraites (V1-3), but does effect nested covarites such as celltype integration by covaraite-rho (V4).
// 00.00.12 - 2018-06-15 Minor update. When no prioers are provided at the root node we get a crash because we can't idntify deifferneces in NLL withotu a parent NLL
//				Workaround is to provide HG19 NLL as prior. Fix is to run INSIGHT2 on the parent (union of all positions) before testing splits. Draw prior from parent. Perhaps in 00.13...

// system includes - C first, then C++
#include <cassert>
#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <fstream>

// butils modules
#include "stringer.h"
#include "fastfile.h"
#include "mathplus.h"

// global type definitions for this project (template specifications)
#include "covtreetypes.h"

// Project specific includes
#include "covtreeArgs.h"
#include "insightShell.h"
#include "cov.h"			// Covaraite & CovaraiteSet definition, and database directories
#include "bitfield.h"		// Covaraite DB to bitfield mapping...
#include "posMapper.h"		// map a contiguous linear range (ie 0-2881033285) to a set of loci (ie hg19 autosome)
#include "score.h"			// scoreSingleton and scorePair objexcts, used to encapsulate scoring of sovsetSubsets.


// These are the key objects employed by main()
struct runtimeElements {
	covtreeArgs			args;				// argument object, uset to parse, access, serialize and deserlize arguments to the program.
	string::vecS		argumentStack;		// vector of serrialized covtreeArgs, each representing a pending binary subset decomposition... this is vehicle for local recursion.
	covsetDef			covDef;				// The universial set of all covaraites and all thier values. This is static once loaded.	
	covsetDefSubset		covSubsetMaster;	// This is a complete (improper) subset of covDef, essentially all of CovDef as a "Subset" object
	posMapperHG19_t		mapHG19;			// Map to/from sorted chromosome space to linear space.
	insightShellSet		insight;			// Object representation of external / async Insight2 process...
	covDB_t				covdb;				// The covariate database. For each genomic psoition in each cell type, identifies a set of anntoations and a (CTS) covariate value for each CTS covariate.
	bool init( int argc, char **argv, string &status );	// Initialize the components that this program uses from the arguments provide by the user. True if AOK, otherwise check Status.
	bool shutdown( string &status );
};

bool runtimeElements::init( int argc, char **argv, string &status ){
	status = "runtimeInit: ";

	// Process the arguments provided by the user
	bool ok = args.processArgs(argc, argv);
	if (args.v>1) { std::cerr << FmtSimple::strNow() << " Init: Processed input arguments.\n"; };
	if (args.sargExit) return true;

	if (!ok) args.displayHelp();
	if (args.v>0) std::cerr << args.toStr() << std::endl;
	if (!ok) { status += "Unable to parse user arguments"; return(false); }
	if (args.restart == "") {
		// no restart is requested, just process the user provided args
		argumentStack.push_back(args.seralize());
	} else {
		// restart requested, read  the requested args from the restart file
		string line_hi, line_lo, test = args.restart.aslower();
		if (args.v>2) { std::cerr << FmtSimple::strNow() << "\t Init: processing restart from " << args.dirWrk() << "/covtree.out  test = " << test << "\n"; };
		std::ifstream t(args.dirWrk() + "/covtree.out");
		std::getline(t,line_hi); std::getline(t, line_lo);	// read searlized arg arguments from restart file.
		// remove comments, if needed
		if (line_hi[0] == '#') line_hi = line_hi.substr(1); line_hi.trim();
		if (line_lo[0] == '#') line_lo = line_lo.substr(1); line_lo.trim();
		if (args.v>3) { std::cerr << FmtSimple::strNow() << "\t\t Init: line_hi = " << line_hi << "\n"; }
		if (args.v>3) { std::cerr << FmtSimple::strNow() << "\t\t Init: line_lo = " << line_lo << "\n"; }
		if ((test == "hi") || (test == "pr")) { argumentStack.push_back(line_hi); fastfile::ffmkdir(string(args.dirWrk() + "/hi").c_str()); }
		if ((test == "lo") || (test == "pr")) { argumentStack.push_back(line_lo); fastfile::ffmkdir(string(args.dirWrk() + "/lo").c_str()); }
		if ((test == "hi") || (test == "pr")) {
			status += "\nattempting to deseralize line_hi";
			if (!args.deseralize( line_hi )) return false;
		} else if (test == "lo") {
			status += "\nattempting to deseralize line_lo";
			if (!args.deseralize(line_lo)) return false;
		} else { status += "\nUnable to parse test case not hi|lo|pr."; return false; }
	}

	// Get Covariate definition, that is the universe of all covaraites and values that we are working with
	if (args.v>2) { std::cerr << FmtSimple::strNow() << "\tInit: Reading Covaraite definitions from: " << args.fnCovDef << ".\n"; };
	if (! covDef.read(args.fnCovDef) ) {
		status += "Unable to parse covaraite definitions from: " + args.fnCovDef; return false ; };

	// Allocate a linear mapping, defaults to HG19 autosome... but you can override from a bedfile... coordinates must be synced to tagsets, so not reccomended.
	if (args.v>2) { std::cerr << FmtSimple::strNow() << "\tInit: Constructing linear mapping to refrence genome.\n"; };
	if (!mapHG19.read()) {
		status += "Unable to initialize linear address space mapper from hg19."; return false; };

	// Start up the INSIGHT2 process
	{
		string flog = args.dirWrk() + "/insight2.log";
		if (args.v>2) { 
			std::cerr << FmtSimple::strNow() << "\tInit: starting INSIGHT2 process. \n";
			std::cerr << FmtSimple::strNow() << "\t\texe     :" << args.insightExe << "\n";
			std::cerr << FmtSimple::strNow() << "\t\tdb      :" << args.insightDB << "\n";
			std::cerr << FmtSimple::strNow() << "\t\tThreads :" << args.insightThreads << "\n";
			std::cerr << FmtSimple::strNow() << "\t\tLog     :" << flog << "\n";
			std::cerr << FmtSimple::strNow() << "\t\tArg     :" << args.getInsightArgs() << "\n\n";
		}
		// Identify the maximum number of samples, that is the number of cell types...
		{	
			string insightdb=(args.insightDBwrk == "" ? args.insightDB : args.insightDBwrk );
  			if (!insight.serverInitialize(args.insightExe, insightdb, flog, args.getInsightArgs(), args.insightThreads )) {
				status += "Unable to initialize INSIGHT2."; return false; };
		}
	};

	// Initialize the bitfield structure of the database, and load database .
	{
		if (args.v>2) { std::cerr << FmtSimple::strNow() << "\tInit: Mapping covarate field definitions to bitfield representation.\n"; };
		// Create a mapping between the defined fields and the bitfields... this is requried before loading the DB...
		covDB_t::bitfieldDef_t map_ann, map_cts;
		for (const auto & s : covDef.covsAnn()) {
			map_ann.addField(s.tagW(), (uint8_t)s.val().size() - 1);
		};

		for (const auto & s : covDef.covsCTS()) {
			map_cts.addField(s.tagW(), (uint8_t)s.val().size() - 1);
		};

		// Identify the name and file location of each tagset.. One for Annotations, and one per cell type for CTS data.
		covDB_t::dbElemSrc		srcTmp, srcAnn; srcAnn.name = "E999"; srcAnn.fname = covDef.dirAnn() + "/" + srcAnn.name + "/covs.bedg";
		covDB_t::vdbElemSrc_t	srcCts;
		for (const auto & s : args.vcts) {
			srcTmp.name = s; srcTmp.fname = covDef.dirCTSb() + "/" + srcTmp.name + "/covs.bedg"; srcCts.push_back(srcTmp); srcTmp.clear(); };

		// Now we can load the DB! - keeps coppies of maps internally...
		string dirbase="";
		if (args.v>2) { std::cerr << FmtSimple::strNow() << "\tInit: Loading covaraite database from: " << dirbase << " | " << srcAnn.fname << "\n"; };
		if (!covdb.loadDB(dirbase, srcAnn, map_ann, srcCts, map_cts, mapHG19.lastPos() + 1,(args.v<=0?0: args.v-1))) {
			status += "Unable to load database into bitfield."; return false;
		};
	}

	// Generate Universal subset, that is the subset that contains all elements in the covDef.
	if (args.v>2) { std::cerr << FmtSimple::strNow() << "\tInit: Generating master cvariate subset \n"; };
	if (!covSubsetMaster.Init(covDef)) {
		status += "Unable to generate complete master subset from covaraite definitions."; return false;
	};

	if (args.v>1) { std::cerr << FmtSimple::strNow() << " Init: Initialization complete.\n"; };
	status +="AOK"; return true;
}

bool runtimeElements::shutdown(string &status) {
	status = "";
	covdb.clear();
	covSubsetMaster.clear();
	covDef.clear();
	mapHG19.clear();
	insight.serverTerminate(); insight.msleep(3000); insight.serverTerminate(true);
	argumentStack.clear();
	args.clear();
	return true;
}

// 0 - Init runtime, process args, load daabase, start INSIGHT2,

// 1 - Get Ordering / Insight Scores for nonmonotonic covariates, conditioned on covDefST
bool generateOrderings(runtimeElements &rt, const covtreeArgs &args, const covsetDefSubset &covDefST, scorePair::covsetCovOrder_t &covset_scores, std::vector<scoreSingelton> &sings);

// 2 - Generate adn getinsight score for each partition
bool evaluatePartitions(runtimeElements &rt, const covtreeArgs &args, const covsetDefSubset &covDefST, const scorePair::covsetCovOrder_t &covset_scores, scorePair::scorePair_v &partitions);

// 3 - sort partitionigns first by covaraite, then covariates by best score.
typedef std::vector<const scorePair *> scorePair_vp;
//bool sortPartitions(const scorePair::scorePair_v &partitions, std::vector<scorePair_vp> &covParts, int v = 1);
bool sortPartitions(scorePair::scorePair_v &partitions, std::vector<scorePair_vp> &covParts, int v = 1);

// 4 - dump results in viewable sumamry format
bool saveResults(const runtimeElements &rt, const covtreeArgs &args, const covsetDefSubset &covDefST, const std::vector<scoreSingelton> &sings, const std::vector<scorePair_vp> &covParts);

// 5 - test best partition and split again, if as needed.
bool processRecursion(const runtimeElements &rt, const covtreeArgs &args, const std::vector<scorePair_vp> &covParts, string::vecS &argStack);

// this is used to change values of strings duing debugging...
void setstring( std::string &s, const char *c) { 
	s = (string)c; }

int main(int argc, char **argv ) {
	using std::cerr;
	using std::endl;
	using std::cout;
	using std::flush;

	std::cerr << FmtSimple::strNow() << " Main: Entering initalization.\n";
	// Create and initialize runtime elements
	runtimeElements rt;
	{
		string status;
		if (!rt.init(argc, argv, status)) {
			cerr << "CovTree: Unable to initialize runtime elements.\n\tInit() returns: " << status << "\n\n";
			return -1;
		};
		//  just serialize the input argument list and exit....
		if (rt.args.sargExit) {
			std::cout << rt.args.seralize() << "\n" << std::endl;
			return 0;
		}
	}

	std::cerr << FmtSimple::strNow() <<  " Main: Initalization complete, entering recursive subdivision.\n";
	// Process data
	while (rt.argumentStack.size()>0) {

		// Get arguments for this iteration
		covtreeArgs args;
		{ 
			// string sarg= rt.argumentStack.back();			// this would give us depth first processeing
			string sarg = rt.argumentStack[0];					// breadth first... I think this is what we want
			bool ok;
			args = rt.args;										// Copy trasient elements (not persisted through serialziation) from user input.
			ok = args.deseralize(sarg);							// Overwrite persistable arguments with those from string. TODO check for error here!
			// rt.argumentStack.pop_back();						// depth first
			rt.argumentStack.erase(rt.argumentStack.begin());	// breadth first, erase the first element
			if (!ok) {
				std::cerr << "\n\nCovTree-main: FATAL-ERROR Unable to parse arguments, skipping.\n\t" << sarg << "\n\n";
				continue;
			}
			if (args.v >0) { std::cerr << "\n\n" << FmtSimple::strNow() << " Main: Got Arguments " << sarg << "\n"; }
		}
		try {
			{ 
			// TODO: if we have no parental values (sometimes provided as priors) we should run INSIGHT2 on the parrental (universal) set, and then use ther results
			// as priors... perhaos for 00.13. If we impliment this, remove default priors from covtreeArgs! --Brad

			double d_zero=0.0; scorePair::ParentalRho( NULL, &d_zero ); };
			// Genereate working subset
			if (args.v > 0) { std::cerr << FmtSimple::strNow() << " Main: Initializing master covaraiteSet structure.\n"; }
			covsetDefSubset covDefST;								// Master Subset, that is the subset defined by the user input
			covDefST.Init(rt.covSubsetMaster, args.covStart );		// TODO check for error here... 

			// Develop Conditional (based on covdefST) ordering for covaraits that are nonmonotonic
			if (args.v > 0) { std::cerr << FmtSimple::strNow() << " Main: Generating conditional ordering for non-monotonic covaraites.\n"; }
			scorePair::covsetCovOrder_t covset_scores;
			std::vector<scoreSingelton> sings;				// we'll want this temp value later when we dump sumamry stats
			generateOrderings(rt, args, covDefST, covset_scores, sings );

			// The Meat of the algorithm, generate and evaluate all ordered covariates
			if (args.v > 0) { std::cerr << FmtSimple::strNow() << " Main: Processing all possible single, binary, covarate splits.\n"; }
			{	// ACK. We need the parental estimate of RHO to get global consistency. This is ONLY passed in as the prior on rho.
				double d_pRho=0.0;
				if (!args.implicitParRho ) {
					// for global consistency use parrental RHO estimate from INSIGHT in parrent data set. This allows total SUS to change
					//		as the sum of the child SUS need not = parrent SUS. Global in the sense that ths sum of all deltaInf = root Inf.
					// for local consistency assume parrental rho for two chidleren is (p1*n1 + p2*n2) / (n1+n2), this keeps the SUS locally consistent
					//		at each split, and favors more conerent classes over those with fewer (more extreme) combined SUS.
					if (args.vpriors.size()>0 && args.vpriors[0].trim() != "") d_pRho = std::stod( args.vpriors[0].c_str());
				}
				scorePair::ParentalRho( NULL, &d_pRho ); 
			}
			scorePair::scorePair_v partitions;
			evaluatePartitions(rt, args, covDefST, covset_scores, partitions);

			// Sort the partitions within each covaraits, then by best score for each covariate
			// The result is a vector of sd
			if (args.v > 0) { std::cerr << FmtSimple::strNow() << " Main: Sorting partitions by informative content.\n"; }
			std::vector<scorePair_vp> covParts;
			sortPartitions(partitions, covParts, args.v );	// covParts contains a vector (one for each covartes) of vectors (one for each split) of poitners to elements in partitions.
			// Celar covParts BEFORE cearing partitions, it concains pointers to elements in partitions...

			// now covParts[0][0] contains a pointer to the scorePair with the highest score!!
			if (args.v > 1 ) {
				for (const auto & scov : covParts) {
				std::cerr << scov[0]->CovName() << "\n";
				for (const auto & scovp : scov) {
					std::cerr << "\t" << scovp->delInfBits() << "\t" << scovp->CovValTag(true) << "\t|\t" << scovp->CovValTag(false) << "\n";
				}
				std::cerr << "\n";
			}
			}

			// Format and print results
			if (args.v > 0) { std::cerr << FmtSimple::strNow() << " Main: Saving Partition Results.\n"; }
			saveResults( rt, args, covDefST, sings, covParts);				// TODO process error 

			// Writes recursion request to shell, for forked processing, or adds jobs to argument stack, as requested by user.
			if (args.v > 0) { std::cerr << FmtSimple::strNow() << " Main: Checkign and processing needed reciursion.\n"; }
			processRecursion(rt, args, covParts, rt.argumentStack);			// TODO process error 
		} catch (const std::exception &e) {
			std::cerr << "\n\n" << FmtSimple::strNow() << " CovTree-MAIN: FATAL - Caught exception : \n\t" << e.what() << " \n\tThis is bad. Writing .done semaphore and attempting to continue.\n\n";
		} catch (...) {
			std::cerr << "\n\n" << FmtSimple::strNow() << " CovTree-MAIN: FATAL - Caught an unknown exception. This is very bad. Writing .done semaphore and attempting to continue.\n\n";
		}
		// write .done semaphore - This must be the LAST thing in the loop!
		{
			string fout = args.dirWrk() + "/" + args.fnOut + ".done";
			if (args.v>1) { std::cerr << FmtSimple::strNow() << "\n Main: Writing completion semaphore : " << fout << "\n"; }
			std::fstream fo(fout.c_str(), std::ios_base::out); fo.close();
		};
		
	}; // WHILE rt.argumentStack.size()>0
	
	std::cerr << FmtSimple::strNow() << " Main: Recursion complete, cleaning up.\n"; 
	// Shutdown runtime elements
	{
		string status;
		if (!rt.shutdown(status)) {
			cerr << "CovTree: Unable to shutdown runtime elements.\n\tshutdown() returns: " << status << "\n\n"; return -2; };
	};

	std::cerr << FmtSimple::strNow() << " Main: EXITING.\n"; 
	return(0);
}

// the first three are const inputs, the last is the output
// rt.covdb.scanpat requires altering state, which is a non-const operation, so we have to allow non-cosnt ref..
bool generateOrderings( runtimeElements &rt, const covtreeArgs &args, const covsetDefSubset &covDefST, scorePair::covsetCovOrder_t &covset_scores, std::vector<scoreSingelton> &sings ) {

	//scoreSingelton sing;

	if (args.v>2) { std::cerr << FmtSimple::strNow() << "\t EvalOrdering - Genreating conditional single-covaraite-value subsets.\n"; };
	scoreSingelton::GenerateSingeltons(covDefST, rt.covdb.mapAnn(), rt.covdb.mapCTS(), args.dirWrk() , sings);
	covDB_t::bitfieldDB::vscanTarget_t targets;
	for (uint32_t i = 0; i < sings.size(); i++) {
		targets.resize(i + 1);
		targets[i].clear();
		if (args.v>3) { std::cerr << sings[i].CovName() << "-" << sings[i].CovValTag() << "    "; };
		sings[i].getMasks(targets[i].patA1, targets[i].patC1); // only hte first pattern is used (1, not 2)
		sings[i].getMasks(targets[i].patA2, targets[i].patC2);
	};
	if (sings.size() < 1) {
		if (args.v>2 ) { std::cerr << FmtSimple::strNow() << "\t EvalOrdering - No nonmonotonic splits found, terminating ordering inference.\n"; };
		return true;
	}
	if (args.v>3) { std::cerr << "\n"; };

	// Initialize scanning and position variables.
	if (args.v>2) { std::cerr << FmtSimple::strNow()<< "\t EvalOrdering - Initializing database scan loop.\n"; };
	rt.covdb.scanpatInit();
	LinAddr_t lin_loc_start = 0, lin_loc_len = 0, lin_chr_first = 0, lin_chr_last = 0, lin_LAST = rt.mapHG19.lastPos();
	uint8_t chrid;
	uint32_t chr_first, chr_start, chr_len; stringer chr_name;

	// Skip to the desired chromosome, if desired for debugging purposes.
	if (args.debugOneChrom.length()>0) {
		chrid=0;
		if (args.v>2) { std::cerr << FmtSimple::strNow() << "\t EvalOrdering - Skipping to chromosome " << args.debugOneChrom << "\n"; };
		if (!rt.mapHG19.chrToID(args.debugOneChrom,chrid)) return false;
		lin_loc_start= rt.mapHG19.firstPos( chrid ); lin_LAST = rt.mapHG19.lastPos(chrid);
	}

	// Map the starting position to chromosomeID, offset etc...
	if (args.v>2) { std::cerr << FmtSimple::strNow() << "\t EvalOrdering - Initializing database scan.\n"; };
	rt.mapHG19.posToID(lin_loc_start, chrid);
	rt.mapHG19.idInfo(chrid, chr_first, chr_name, chr_start, chr_len); lin_chr_first = chr_first; lin_chr_last = lin_chr_first + chr_len - 1;

	// If debugSkipInsight activated, we assume INSIGHT has alredy been run on this data set,
	//	dont regenerate the loci & rerun Insight, just read the existing scores.
	bool do_debug = ((args.debugSkipInsight & 1) > 0 );
	for (uint32_t i = 0; i < sings.size(); i++) { 
		// this ALSO initializes file names, important for retrieving resutls even if we are not recalculating them
		sings[i].scanInit(do_debug); 
		if (args.v>3) { std::cerr << FmtSimple::strNow() << "\t\t EvalOrdering - " << sings[i].getFname() << "\n"; };
	}

	if (!do_debug) {
		// process locus scan and generate weighted bed files for Insight2 input.
		uint32_t dbg_cnt = 0;
		while (lin_loc_start <= lin_LAST) {
			rt.covdb.scanpatGetNext(targets, lin_loc_start, lin_loc_len);
			if (lin_loc_start > lin_chr_last) {
				rt.mapHG19.posToID(lin_loc_start, chrid);
				rt.mapHG19.idInfo(chrid, chr_first, chr_name, chr_start, chr_len);
				lin_chr_first = chr_first; lin_chr_last = lin_chr_first + chr_len - 1;
			}
			// Make sure that the locus returned from GetNext does not straddle a chromosome boundary (should never happen)
			assert(lin_chr_first <= lin_loc_start);
			assert(lin_chr_last >= (lin_loc_start + lin_loc_len - 1));
			// Save the info to each scan range... onlu the first pattern is used 
			uint32_t t_start = lin_loc_start - lin_chr_first, t_end = t_start + lin_loc_len;// bed format, half open...
			for (uint32_t i = 0; i < targets.size(); i++) {
				sings[i].scanAddLocus(chr_name, t_start, t_end, targets[i].hits1);
			}
			lin_loc_start += lin_loc_len;
			if (args.v>3) {
				if ((lin_loc_start / 100000000) > dbg_cnt) {
					std::cerr << FmtSimple::strNow() << "\t\t EvalOrdering - " << FmtSimple::fmt(lin_loc_start, 10) << "\t" << chr_name << "\t" << FmtSimple::fmt(chr_start, 10) << "\t" << FmtSimple::fmt(lin_loc_len, 6) << "\n";
					dbg_cnt = lin_loc_start / 100000000;
				};
			};
		};
		for (uint32_t i = 0; i < sings.size(); i++) { 
			sings[i].scanTerminate(); } // close and flush files

		// Queue calculation of INSIGHT2 scores for each element in the range
		if (args.v>2) { std::cerr << FmtSimple::strNow() << "\t EvalOrdering - Requesting INSIGHT2 evaluation of single covaraite-value data sets.\n"; };
		for (uint32_t i = 0; i < sings.size(); i++) {
			string insight_job = sings[i].getFname() + " # " + args.getInsightArgs();
			if (args.v>3) { std::cerr << FmtSimple::strNow() << "\t\t EvalOrdering - InsInit " << insight_job << "\n"; };
			rt.insight.jobInit(insight_job);
		};
	};

	if (args.v>2) { std::cerr << FmtSimple::strNow() << "\t EvalOrdering - Awaiting completion of INSIGHT2 calculation.\n"; };
	// Await completion semapoure and read INSIGHT2 results
	for (uint32_t i = 0; i < sings.size(); i++) {
		if (args.v>3) { std::cerr << FmtSimple::strNow() << "\t\t EvalOrdering - InsWAIT " << sings[i].getFname() << "\n"; };
		rt.insight.jobWait(sings[i].getFname());
		sings[i].resultsRead();
	};
	if (args.v>2) { std::cerr << FmtSimple::strNow() << "\t EvalOrdering - Reordering non-monotonic covaraits based on INSIGHT2 results.\n"; };

	// process sorting of covaraites according to Insight scores...
	string cov_name = "";
	scorePair::covValWeight_v	cov_val_scores; cov_val_scores.clear();
	// Create Ordering for each non-mono covariate 
	for (uint32_t i = 0; i < sings.size(); i++) {
		if (i==0) {cov_name = sings[i].CovName(); };	// careful, sings.size() may be 0.
		if (cov_name != sings[i].CovName()) {
			covset_scores[cov_name] = cov_val_scores;					// Store ordering
			cov_val_scores.clear();	cov_name = sings[i].CovName();		// Clear orderling list...
		}
		double v_res = sings[i].resultValRho();
		// null score for a null value
		if (( cov_val_scores.size() == 0 ) && (covDefST.getCovByName(sings[i].CovName())->hasNULL()) ) v_res = 0.0;
		cov_val_scores.push_back(v_res);
		assert(sings[i].CovValTag() == covDefST.getCovByName(sings[i].CovName())->vals[cov_val_scores.size()-1]->tag());		// make sure we are pulling the right value...
	};
	if (sings.size()>0) {
		covset_scores[cov_name] = cov_val_scores;					// Store ordering
		cov_val_scores.clear();	cov_name = "";						// Clear ordering list...
	};

	if (args.v>2) { std::cerr << FmtSimple::strNow() << "\t EvalOrdering - Ordering completed sucessfully.\n"; };

	return true;
}

// Evaluate partitions - get raw data. Not concerned with QC.
bool evaluatePartitions(runtimeElements &rt, const covtreeArgs &args, const covsetDefSubset &covDefST, const scorePair::covsetCovOrder_t &covset_scores, scorePair::scorePair_v &partitions) {

	// Generate all ordered binary partitions....
	if (args.v>2) {	std::cerr << FmtSimple::strNow().trim() << "\t EvalPart - Generating list of covarate-wise Partitions of Master Covariate set.\n"; };
	scorePair::GenerateCosets(covDefST, covset_scores, rt.covdb.mapAnn(), rt.covdb.mapCTS(), args.dirWrk(), partitions, (args.v>0?args.v-1:0), &args);
	if (partitions.size()<1) {
		if (args.v>0) { std::cerr << FmtSimple::strNow().trim() << "\t EvalPart - No partitions found to test, exiting.\n"; };
		return true;
	}

	// Generate database encoded bitmasks associated with each set/coset in the partition
	if (args.v>2) {	std::cerr << FmtSimple::strNow().trim() << "\t EvalPart - Generating database bitmask for each partition coset. \n"; };
	covDB_t::bitfieldDB::vscanTarget_t targets;
	targets.clear();
	for (uint32_t i = 0; i < partitions.size(); i++) {
		targets.resize(i + 1); targets[i].clear();
		partitions[i].getMasks(targets[i].patA1, targets[i].patC1, targets[i].patA2, targets[i].patC2);
	}

	// Initialize scanning and position variables.
	if (args.v>2) {	std::cerr << FmtSimple::strNow().trim() << "\t EvalPart - Initializing database Scan\n"; };
	rt.covdb.scanpatInit();
	LinAddr_t lin_loc_start = 0, lin_loc_len = 0, lin_chr_first = 0, lin_chr_last = 0, lin_LAST = rt.mapHG19.lastPos();
	uint8_t chrid;
	uint32_t chr_first, chr_start, chr_len; stringer chr_name;

	// Skip to the desired chromosome, if desired for debugging purposes.
	if (args.debugOneChrom.length()>0) {
		chrid = 0;
		if (args.v>2) { std::cerr << FmtSimple::strNow() << "\t EvalPart - Skipping to " + args.debugOneChrom + "\n"; };
		if (!rt.mapHG19.chrToID(args.debugOneChrom, chrid)) return false;
		lin_loc_start = rt.mapHG19.firstPos(chrid); lin_LAST = rt.mapHG19.lastPos(chrid);
	}
	
	// Map the starting position to chromosomeID, offset etc...
	rt.mapHG19.posToID(lin_loc_start, chrid);
	rt.mapHG19.idInfo(chrid, chr_first, chr_name, chr_start, chr_len); lin_chr_first = chr_first; lin_chr_last = lin_chr_first + chr_len - 1;

	bool do_debug = ((args.debugSkipInsight & 2) > 0);
	for (uint32_t i = 0; i < partitions.size(); i++) {
		// this ALSO initializes file names, important for retrieving resutls even if we are not recalculating them
		partitions[i].scanInit( do_debug ); 
		if (args.v>4) { std::cerr << FmtSimple::strNow() << "\t\t\t EvalPart - Initializing S: " + partitions[i].getFname(true) + "\n"; };
		if (args.v>4) { std::cerr << FmtSimple::strNow() << "\t\t\t EvalPart - Initializing C: " + partitions[i].getFname(false) + "\n"; };
	}

	if (!do_debug ){
		if (args.v>2) {	std::cerr << FmtSimple::strNow().trim() << "\t EvalPart - Entering scan loop of covarate database \n"; };
		uint64_t dbg_cnt = 0;
		// process locus scan and generate weighted bed files for Insight2 input.
		while (lin_loc_start <= lin_LAST) {
			rt.covdb.scanpatGetNext(targets, lin_loc_start, lin_loc_len);
			if (lin_loc_start > lin_chr_last) {
				rt.mapHG19.posToID(lin_loc_start, chrid);
				rt.mapHG19.idInfo(chrid, chr_first, chr_name, chr_start, chr_len);
				lin_chr_first = chr_first; lin_chr_last = lin_chr_first + chr_len - 1;
			}
			// Make sure that the locus returned from GetNext does not straddle a chromosome boundary (should never happen)
			assert(lin_chr_first <= lin_loc_start);
			assert(lin_chr_last >= (lin_loc_start + lin_loc_len - 1));
			// Save the info to each scan range... onlu the first pattern is used 
			uint32_t t_start = lin_loc_start - lin_chr_first, t_end = t_start + lin_loc_len;// bed format, half open...
			for (uint32_t i = 0; i < partitions.size(); i++) {
				partitions[i].scanAddLocus(chr_name, t_start, t_end, targets[i].hits1, targets[i].hits2);
			}
			lin_loc_start += lin_loc_len;
			if (args.v > 3) {
				if ((lin_loc_start / 100000000) > dbg_cnt) {
					std::cerr << FmtSimple::strNow().trim() << "\t\tEvalPart \t" << FmtSimple::fmt(lin_loc_start,10) << "\t" << chr_name << "\t" << FmtSimple::fmt(chr_start,10) << "\t" << FmtSimple::fmt(lin_loc_len, 6 ) << "\n";
					dbg_cnt = lin_loc_start / 100000000;
				}
			}
		};

		if (args.v>2) {	std::cerr << FmtSimple::strNow() << "\tEvalPart - Scan loop completed - closing " << partitions.size() << " Insight Input files.\n"; };
		for (uint32_t i = 0; i < partitions.size(); i++) {
			partitions[i].scanTerminate(); } // close and flush files

		// make sure both the set and the coset had positions, and both are different from the parent. if not, delete the improper partition.
		if (args.v>2) {	std::cerr << FmtSimple::strNow() << "\tEvalPart - Checking for proper partition cosets (each size>0).\n"; };
		for (uint32_t ix = (uint32_t) partitions.size(); ix > 0; ix--) {
			uint32_t i = ix - 1;
			if (!partitions[i].scanValid()) {
				uint64_t s = 0, sw = 0, c = 0, cw = 0;
				partitions[i].scanGetSums(s, sw, c, cw);
				if (args.v>2) { std::cerr << FmtSimple::strNow() << "\t\tEvalPart - Invalid Partition " << s << " " << sw << " " << c << " " << cw << " removing files " << partitions[i].getFname(true) << " " << partitions[i].getFname(false) << "\n"; }
				fastfile::ffrmf(partitions[i].getFname(true)); 				fastfile::ffrmf(partitions[i].getFname(false));
				// to do this, we need a copy constructor, too much hassle.. jsut mark it as invalid and move on...
				// partitions.erase(std::remove(partitions.begin(), partitions.end(), i), partitions.end());	// slow, but this is not a bottleneck (rarely done).
				partitions[i].qcOk = false;	// yes I know we arent "itnerested" in qc, but an improper partition is not really a partition at all. And can cause Insight2 crash.
				partitions[i].notes += "| QC failed improper partition - 1 : " + std::to_string( sw ) + "  2 : " + std::to_string(cw);
			}
		}

		// Initiate clculation of INSIGHT2 scores for each valid element in the range
		if (args.v>2) { std::cerr << FmtSimple::strNow() << "\tEvalPart - Scan loop completed - requesting Insight scores.\n"; };
		for (uint32_t i = 0; i< partitions.size(); i++ ) {
			if (partitions[i].scanValid()) {
				string insight_job;
				insight_job = partitions[i].getFname(true) + " # " + args.getInsightArgs();
				if (args.v>3) { std::cerr << FmtSimple::strNow() << "\t\tEvalPart - InitINS2 - s - " << insight_job << "\n";}
				rt.insight.jobInit(insight_job, do_debug);
				insight_job = partitions[i].getFname(false) + " # " + args.getInsightArgs();
				if (args.v>3) { std::cerr << FmtSimple::strNow() << "\t\tEvalPart - InitINS2 - c - " << insight_job << "\n"; }
				rt.insight.jobInit(insight_job, do_debug);
			}
		}
	};

	// Get INSIGHT2 results
	if (args.v>2) { std::cerr << FmtSimple::strNow().trim() << "\tEvalPart - Awaiting Insight scores and reading results.\n"; };

	for (uint32_t i = 0; i < partitions.size(); i++) {
		if (partitions[i].scanValid()) {
			if (args.v>3) { std::cerr << FmtSimple::strNow().trim() << "\t\t DONE-" <<  partitions[i].getFname(true) << "\n"; };
			rt.insight.jobWait(partitions[i].getFname(true),  0);		// wait forever, if need be
			if (args.v>3) { std::cerr << FmtSimple::strNow().trim() << "\t\t DONE-" <<  partitions[i].getFname(false) << "\n"; };
			rt.insight.jobWait(partitions[i].getFname(false), 0);		// wait forever, if need be...
			partitions[i].resultsRead();
		};
	};

	if (args.v>2) { std::cerr << FmtSimple::strNow().trim() << "\tEvalPart - Exiting Sucessfully.\n"; };
	return true;
}

// partitions should really be const scorePair::scorePair_v &partitions, but we need to perform QC.....
bool sortPartitions( scorePair::scorePair_v &partitions, std::vector<scorePair_vp> &covParts, int v ) {

	// We don't really operate on partitions.. some are invalid.. instead create a vector of pointers and operate on that
	std::vector<scorePair *> partitions_pv(partitions.size());
	for (size_t  i=0; i<partitions.size(); i++) partitions_pv[i] = &( partitions[i]);

	// How many are proper / valid? Filter...
	if (v>2) std::cerr << FmtSimple::strNow() << "\t SortPartitions() : found " << partitions_pv.size() << " total partitions.\n";

	for (size_t ix = partitions_pv.size(); ix>0; ix--) {
		size_t i = ix -1;
		if (!partitions_pv[i]->qcOk && (partitions_pv[i]->notes.find("improper") != string::npos)) 
			partitions_pv.erase(partitions_pv.begin() + i); // Linear! but rarely used...
			//partitions_pv.erase(std::remove(partitions_pv.begin(), partitions_pv.end(), partitions_pv.begin()+i), partitions_pv.end());
	}
	if (v>2) std::cerr << FmtSimple::strNow() << "\t SortPartitions() : found " << partitions_pv.size() << " proper partitions.\n";

	// Test each partition for the advanced QC requirements / heuristics requirted by user.
	uint32_t numok=0;
	numok = 0; for ( const auto &s : partitions_pv ) { if (s->qcOk) numok++; }
	if (v>2) std::cerr << FmtSimple::strNow() << "\t SortPartitions() : partitions passing qc1 (min size & min inf) = " << numok << " .\n";
	for ( scorePair &s : partitions ) {
		if (s.qcOk) {
			if (s.maxpErr1()<1.0 && ((s.pErr1()*numok)>s.maxpErr1())) { s.qcOk = false;	s.notes += " | QC: Corrected P(rho Inversion) too hi = " + std::to_string(s.pErr1()*numok);	} //
			if (s.maxpErr2()<1.0 && ((s.pErr2()*numok)>s.maxpErr2())) { s.qcOk = false;	s.notes += " | QC: Corrected P(rho <minDist)  too hi = " + std::to_string(s.pErr2()*numok);	} //
		}
	}
	numok = 0; for (const auto &s : partitions_pv) { if (s->qcOk) numok++; }
	if (v>2) std::cerr << FmtSimple::strNow() << "\t SortPartitions() : partitions passing qc2 (mac prob rho1 > rho2) = " << numok << " .\n";
	// QC failures all report delInfBits of 0.

	// Sort partitions, the result is a vector (1 elem for each cov type), of vectors of partition pointers.
	//	Within a cell type all partitions are sorted by decreasing delInfBits (delta Info, more is better), and then covs are sorted
	//	in descending order by bestpartition. To the very best is in part[0][0], the best for covariate 3 woudl be part[2][0], the
	//		third best for covaraite 4 would be part[3][2]...
	scorePair_vp covValParts;
	string covname = "";
	for (uint32_t i = 0; i < partitions_pv.size(); i++) {
		if ((i != 0) && (covname != partitions_pv[i]->CovName())) {
			if (covValParts.size()>0){
				std::sort(covValParts.begin(), covValParts.end(), [](const auto & a, const auto & b) { return(a->delInfBits() > b->delInfBits()); });	// sort in DESCENDING order
				covParts.push_back(covValParts);
			}
			covValParts.clear(); covname = "";
		} 
		//covValParts.push_back(&(partitions[i])); covname = partitions[i].CovName();
		covValParts.push_back(partitions_pv[i]); covname = partitions_pv[i]->CovName();
	}
	if (covValParts.size()>0) {
		std::sort(covValParts.begin(), covValParts.end(), [](const auto & a, const auto & b) { return(a->delInfBits() > b->delInfBits()); });			// sort in DESCENDING order
		covParts.push_back(covValParts); covValParts.clear();
	}
	// Now sort covaraites 
	std::sort(covParts.begin(), covParts.end(), [](const auto & a, const auto & b) { return(a[0]->delInfBits() > b[0]->delInfBits()); });				// sort in DESCENDING order
	return true;
}

bool saveResults(const runtimeElements &rt, const covtreeArgs &args, const covsetDefSubset &covDefST, const std::vector<scoreSingelton> &sings, const std::vector<scorePair_vp> &covParts) {

	if (args.v>2) { 
		std::cerr << FmtSimple::strNow() << "\t Save Results - beginning result dump.\n"; };

	// Convienent object for formatting result summaries.
	struct CovSingSum : public FmtSimple {
		string name;
		uint32_t numVals, numValsMax;
		double   numPos, infCond2, infUcond2, infDel2, infDelEnt2;
		void clear() { name = ""; numVals = numValsMax = 0; numPos = infCond2 = infUcond2 = infDel2 = infDelEnt2 = 0; return; };
		uint32_t covSummary(const covsetDefSubset &covs, uint32_t NumCells, const std::vector<scoreSingelton> &sings, const string &CovName) {
			clear();
			name = CovName; numValsMax = (uint32_t) covs.getCovByName(name)->univ->val().size();
			//numVals = covs.getCovByName(name)->numVals();	// Can include NULL
			double pos = 0, rho = 0, sumrho = 0;
			for (uint32_t i = 0; i<sings.size(); i++) {
				const auto &asing = sings[i];
				if (asing.CovName() == name) {
					numVals++;	// Jsut number of values for which we retrieved results...
					asing.resultVal("STATS", "InPosW", pos);	numPos += (pos / NumCells);		// Unweight the number of positions
					rho = asing.resultValRho(); // .resultVal("PARAMETERS", "Rho", rho);
					sumrho += rho * numPos;
					infCond2 += mathplus::entropy(rho) * numPos * mathplus::natToBit;
				}
			}
			sumrho /= numPos;
			{
				double pRho=0.0;
				scorePair::ParentalRho(&pRho);
				if (pRho > 0.0 &&  pRho < 1.0) sumrho = pRho;
			}

			infUcond2 = mathplus::entropyB(sumrho) * numPos * mathplus::natToBit;
			infDel2 = infUcond2 - infCond2; infDelEnt2 = infDel2 / numPos;
			return numVals;
		}
		string toStr() {
			string s = fmt(name, 8) + "\t" + fmt(numVals, 2) + "\t" + fmt(numValsMax, 2) + "\t" + fmt2(numPos, 13, true) + "\t" + fmt2(infUcond2, 13, true) + "\t" + fmt2(infCond2, 13, true) + "\t" + fmt2(infDel2, 10, true) + "\t" + fmt(infDelEnt2, 8, 6);
			return(s);
		}
	};

	// Dump Singleton Analysis Summary...
	{
		string fnout = args.dirWrk() + "/" + args.fnOut + "sing.txt";
		if (args.v>2) { std::cerr << FmtSimple::strNow() << "\t Save Results - Saving singelton data to :" << fnout << ":\n"; };
		std::ofstream fo; fo.open(fnout, std::ios_base::out);
		fo << "#Nonmonotonic Covariate Value Asessment\n";
		string old_cov_name = "";
		for (uint32_t i = 0; i<sings.size(); i++) {
			auto &asing = sings[i];
			if (old_cov_name != asing.CovName()) {
				CovSingSum sum; sum.covSummary(covDefST, rt.covdb.NumCells(), sings, asing.CovName());
				fo << "\nCOV:\n"; fo << "COV1\t" + sum.toStr() + "\n";
				old_cov_name = asing.CovName();
				while ((i < sings.size()) && (old_cov_name == asing.CovName())) {
					auto &bsing = sings[i]; string s;
					s = FmtSimple::fmt(bsing.CovValTag(), 2) + "\t" + FmtSimple::fmt(bsing.CovName(), 8) + "\t";
					s += bsing.toStrNsites(true); s += bsing.toStrParams();
					fo << s << "\n"; i++;
				}
			}
		}
		fo.close();
	}

	// Dump Partitioning Analysis Summary...
	string fnout = args.dirWrk() + "/" + args.fnOut + "part.txt"; 
	if (args.v>2) { std::cerr << FmtSimple::strNow() << "\t Save Results - Saving partition data to :" << fnout << ":\n"; };
	std::ofstream fo; fo.open(fnout, std::ios_base::out);
	fo << "#Covariate Partitioning Summary\n\n";
	fo << "BestSplit:\n#Most informative partition overall\n";
	if(( covParts.size() < 1) || (covParts[0].size()<1)) {
		if (args.v>0) { std::cerr << FmtSimple::strNow() << "\t Save Results - no partitions found, or no best partition found, skipping.\n"; };
		fo << "#No Splits found \n";
		fo.close(); return true; 
	}


	const scorePair &sbest = *(covParts[0][0]); const scoreSingelton *a=NULL, *b=NULL;
	fo << "SplitCov:\t";
	fo << sbest.toStrSum() << "\n";
	a = &(sbest.getPart(true)); b=&(sbest.getPart(false)); if (a->resultValRho() < b->resultValRho()) {std::swap(a,b); }
	fo << "A:\t" << a->toStrSummary() << "\n";
	fo << "B:\t" << b->toStrSummary() << "\n";
	fo << "\n";
	if (args.v>3) {
		std::cerr << FmtSimple::strNow() << "\t\t COVARIATE: " << sbest.CovName() << "\n";
		a = &(sbest.getPart(true)); b = &(sbest.getPart(false)); if (a->resultValRho() < b->resultValRho()) { std::swap(a, b); }
		std::cerr << FmtSimple::strNow() << "\t\t BestPart A: " << a->toStrSummary().trim()  << "\n";
		std::cerr << FmtSimple::strNow() << "\t\t BestPart B: " << b->toStrSummary().trim() << "\n\n";
	};

	fo << "CovList:\n#Most informative partition for every covariate\n";
	for (const auto & covpart : covParts) {
		//if (first) { first = false; continue; }
		if (covpart.size()<1) continue;
		const scorePair &sbest = *(covpart[0]);
		fo << "SplitCov:\t";
		fo << sbest.toStrSum() << "\n";
		a = &(sbest.getPart(true)); b = &(sbest.getPart(false)); if (a->resultValRho() < b->resultValRho()) { std::swap(a, b); }
		fo << "A:\t" << a->toStrSummary() << "\n";
		fo << "B:\t" << b->toStrSummary() << "\n";
		fo << "\n";
		if (args.v>4) {
			std::cerr << FmtSimple::strNow() << "\t\t\t " << FmtSimple::fmt(sbest.CovName(), 10) << " CovBest A: " << a->toStrSummary() << "\n";
			std::cerr << FmtSimple::strNow() << "\t\t\t " << FmtSimple::fmt(sbest.CovName(), 10) << " CovBest B: " << b->toStrSummary() << "\n\n";
		};
	};
	fo << "\n";

	fo << "AllSplits:\n#Every partition in order of informativeness\n";
	//Get a vector of pointers to partitions
	// std::vector< decltype(covParts[0][0]) > all_parts;
	std::vector< const scorePair  * > all_parts;
	for (const auto & acov : covParts) {
		for (auto apart : acov) all_parts.push_back(apart);
	};

	// Remover improper partitions... these aren't really partitions at all..
	for (size_t ix = all_parts.size(); ix>0; ix--) {
		size_t i = ix - 1;
		if (!all_parts[i]->qcOk && (all_parts[i]->notes.find("improper") != string::npos))  
			all_parts.erase(all_parts.begin() + i);// Linear! but rarely used
			// all_parts.erase(std::remove(all_parts.begin(), all_parts.end(), all_parts.begin()+i), all_parts.end());
	}

	// Sort largest to smallest
	std::sort(all_parts.begin(), all_parts.end(), [](const auto & a, const auto & b) { return(a->delInfBits() > b->delInfBits()); } );
	// Dump in one line each...
	fo << "AllSplits:\t" << all_parts.size() << "\n";
	for (const auto & covpart : all_parts) {
		a = &(covpart->getPart(true)); b = &(covpart->getPart(false)); if (a->resultValRho() < b->resultValRho()) { std::swap(a, b); }
		fo << covpart->toStrSum() << "\t|\t" << a->CovValTag() << "\t|\t" << b->CovValTag() << "\t" << covpart->notes << "\n";
		// Notes is the QC info, if the split failed and is not being considered, this epxlains why.
	}
	fo << "\n"; fo.close();

	if (args.v>2) { std::cerr << FmtSimple::strNow() << "\t Save Results - Done. \n"; };
	return true;
}

// Write recursion output, if needed so shell cah fork childeren, or add new childeren to argStack in current thread for processing, as
//		requsted by user.
bool processRecursion(const runtimeElements &rt, const covtreeArgs &args, const std::vector<scorePair_vp> &covParts, string::vecS &argStack ) {

	if (args.v>2) { std::cerr << FmtSimple::strNow() << "\t Recursion - entering recursion test\n"; };

	// if none of the partitions pass QC than do not recurse..
	int numqcparts = 0;
	for ( const auto &s : covParts) { for ( const auto &r : s) { if (r->qcOk) numqcparts++; } };

	// Check stopping criteria and initiate appropriate recursive form, internal (1 thread), mixed or external.
	bool has_parts = (covParts.size()>0) && (covParts[0].size()>0) && (numqcparts>0);
	bool do_recurse = ((args.maxDepth>1) && has_parts && (covParts[0][0]->delInfBits() > args.minBits));	// if args.minbits is unset it is 0, we want to require SOME valid split, must have delInf>0
	if (args.v>2) {
		if (do_recurse) {
			std::cerr << FmtSimple::strNow() << "\t Recursion - qualifies for recursion. MaxDepth=" << args.maxDepth << " best DelInf = " << covParts[0][0]->delInfBits() << " (min " << args.minBits << "). \n";
		} else {
			std::cerr << FmtSimple::strNow() << "\t Recursion - split fails recursion test. MaxDepth=" << args.maxDepth;
			if (has_parts) {
				std::cerr << " best DelInf = " << covParts[0][0]->delInfBits() << " (min " << args.minBits << "). \n";
			} else {
	 			std::cerr << " no partitions found, or no positive best partition.\n";
			}
		}
	}

	string fout = args.dirWrk() + "/" + args.fnOut;
	std::fstream fo(fout.c_str(), std::ios_base::out);
	if (!has_parts) {
		// if there are no more splits to test.... just gracefully indicate that we are done...
		fo << "#$\t No More Partitions to test HI\n";
		fo << "#$\t No More Partitions to test LO\n";
		if (args.v>3) {
			std::cerr << FmtSimple::strNow() << "\t\t HI: No Partitions to test\n";
			std::cerr << FmtSimple::strNow() << "\t\t LO: No partitions to test\n\n";
		}
	} else {
		const scoreSingelton *a = &(covParts[0][0]->getPart(true)), *b = &(covParts[0][0]->getPart(false));
		covtreeArgs part_a(args), part_b(args);
		// place the set/coset with higher rho first... simplifies things later...
		if (a->resultValRho() < b->resultValRho() ) std::swap(a,b);
		// start by coppying parrental args, then modify needed fields
		part_a.maxDepth--; part_b.maxDepth--;
		part_a.baseDir	+= "/" + string("hi");	part_b.baseDir += "/" + string("lo");
		if (part_a.dWrk != "") { part_a.dWrk += "/" + string("hi"); }
		if (part_b.dWrk != "") { part_b.dWrk += "/" + string("lo"); }
		fastfile::ffmkdir(part_a.dirWrk().c_str());			// TODO check for success...
		fastfile::ffmkdir(part_b.dirWrk().c_str());			// TODO check for success...
		// set up new priors
		part_a.vpriors[0] = std::to_string(a->resultValRho()); part_b.vpriors[0] = std::to_string(b->resultValRho());
		// From insight2 posterior estimates/expectations, dont let MAP values of 0 drive priors to 0. Default to genome wide priors when MAP valeus are too small...
		if (a->resultValEta()> 0.0002) part_a.vpriors[1] = std::to_string(a->resultValEta());
		if (b->resultValEta()> 0.0002) part_b.vpriors[1] = std::to_string(b->resultValEta());
		part_a.vpriors[2] = std::to_string(a->resultValGam()); part_b.vpriors[2] = std::to_string(b->resultValGam());
		part_a.vpriors[4] = std::to_string(a->resultValNLL(false)); part_b.vpriors[4] = std::to_string(b->resultValNLL(false));	// Insight2 NLL - ML
		part_a.vpriors[5] = std::to_string(a->resultValNLL(true));  part_b.vpriors[5] = std::to_string(b->resultValNLL(true));	// Insight2 NLL - MaP

		// set up new initial covariate subset
		part_a.covStart = a->covsetDef().toTemplate();	part_b.covStart = b->covsetDef().toTemplate();

    	// Can the partition element be further subdivided?
	    bool can_part_a = a->CanBePartitioned(), can_part_b = b->CanBePartitioned();

		// perform recursion by pushing args back into the procesign stack (singel threaded) or write them to an external file
		//	(for externalized parallel processing).
		{
			// Lines starting with # are comments and ignored, so push back both arg strings, but comment out the ones we are processing
			//	in this thread.. (and not in the shell) - genrally process HI coset before LO one (mroe interesting). Now that we have
			//	moved to breadth first processing, it matters a lot less....
			if (!do_recurse) {
				// even if we are NOT recursing, save the parameters we WOULD use to recurse, so we can pick it up
				//  again later, if we want to...
				fo << "#$\t" + part_a.seralize() + "\n";				fo << "#$\t" + part_b.seralize() + "\n";
			} else {
				switch (args.extRecurse) {
					case 0:
						if (can_part_a) argStack.push_back(part_a.seralize());
						if (can_part_b) argStack.push_back(part_b.seralize());
						fo << (can_part_a ? "#\t" : "#$\t") << part_a.seralize() + "\n";
						fo << (can_part_b ? "#\t" : "#$\t") << part_b.seralize() + "\n";			break;
					case 1: 
						if (can_part_a) argStack.push_back(part_a.seralize());
						fo << (can_part_a ? "#\t" : "#$\t") << part_a.seralize() + "\n";
						fo << (can_part_b ? "" : "#$\t")    << part_b.seralize() + "\n";			break;
					case 2: 
						fo << (can_part_a ? "" : "#$\t") << part_a.seralize() + "\n";
						fo << (can_part_b ? "" : "#$\t") << part_b.seralize() + "\n";				break;
						// default:  // TODO... handle this.....
				}
			}
			if (args.v>3) {
				std::cerr << FmtSimple::strNow() << "\t\t HI: " << (can_part_a ? "" : "#$\t") << part_a.seralize() + "\n";
				std::cerr << FmtSimple::strNow() << "\t\t LO: " << (can_part_b ? "" : "#$\t") << part_b.seralize() + "\n\n";
			}
		}
		fo.close();

		// Dump status message
	};

	// Done.
	if (args.v>2) { std::cerr << FmtSimple::strNow() << "\t Recursion - exiting recursion processing\n"; };
	return true;
}

#ifdef OLD_CODE
int test_Insight(int argc, char **argv) {
	using std::cerr;
	using std::endl;
	using std::cout;
	using std::string;
	using std::flush;
	covtreeArgs args;
	{
		bool ok = args.processArgs(argc, argv);
		if (!ok) args.displayHelp();
		if (args.v>0) {
			cerr << args.toStr() << endl;
		}
		if (!ok) return(-1);
	};
	// Get Args

	string astring;
	insightShell insight;
	cout << "\nExists? : " << args.baseDir << " : " << fastfile::ffexists(args.baseDir.c_str()) << "\n";
	string flog = args.baseDir + "/insight2.log";
	cout << "\nExists? : " << flog << " : " << fastfile::ffexists(flog.c_str()) << "\n";
	cout << "\n\nInitializing server \n" << std::flush;
	bool ok = insight.serverInitialize(args.insightExe, args.insightDB, flog, args.insightArgs);
	cout << "\n\nServer Initialized \n" << std::flush;
	stringer jin1 = "C:/Temp/t/covtree/tfbs.bed";
	cout << "\n\nStarting job: " << jin1 << flush;
	ok = insight.jobInit(jin1);
	stringer jin2 = "C:/Temp/t/covtree/hg19.bed";
	cout << "\n\nStarting job: " << jin2 << flush;
	ok = insight.jobInit(jin2);

	bool done = false; int i = 0;
	while (!done) {
		cout << "\nWaiting for job  " << jin1 << " " << i;
		done = insight.jobWait(jin1, 2000); i++;
	}
	cout << "\nJob Done " << jin1 << " " << i << "\n";
	insightShell::insightResults res1;
	cout << "\nAttempting to read results for jin1 : " << insight.jobResultsGet(jin1, res1) << "\n";
	for (auto & cl : res1.iter()) {
		for (const auto & cli : cl.second) {
			cout << cl.first << "\t" << cli.first << "\t" << cli.second.s << "\n";
		};
	};

	cout << "\n\n";
	done = false; i = 0;
	while (!done) {
		cout << "\nWaiting for job  " << jin2 << " " << i;
		done = insight.jobWait(jin2, 2000); i++;
	}
	cout << "\nJob Done " << jin2 << " " << i;
	insightShell::insightResults res2;
	cout << "\nAttempting to read results for jin2 : " << insight.jobResultsGet(jin2, res2) << "\n";
	for (auto & cl : res2.iter()) {
		for (const auto & cli : cl.second) {
			cout << cl.first << "\t" << cli.first << "\t" << cli.second.s << "\n";
		};
	};

	double v;
	cout << "\n\n";
	cout << "Rho1 = " << (res1.getField("PARAMETERS", "Rho", v) ? v : 0) << "  delta = " << (res1.getField("UNCERTAINTY", "Rho", v) ? v : 0) << "\n";
	cout << "Rho2 = " << (res2.getField("PARAMETERS", "Rho", v) ? v : 0) << "  delta = " << (res2.getField("UNCERTAINTY", "Rho", v) ? v : 0) << "\n";

	stringer dres = "C:/Temp/t/covtree/res/1";
	cout << "Making directory : " << dres << "  :  " << (fastfile::ffmkdir(dres) == 0) << "\n";
	cout << "Moving results   : " << jin1 << "  :  " << insight.jobResultsMove(jin1, dres) << "\n";
	cout << "Moving results   : " << jin2 << "  :  " << insight.jobResultsMove(jin2, dres) << "\n";

	cout << "\nRequesting server termination " << flush;
	done = false; insight.serverTerminate();
	for (i = 0; i<10 && (!done); i++) {
		done = !insight.serverReady();
		if (!done) {
			cout << "\nWaiting for term  " << i;
			insightShell::msleep(1000);
			done = !insight.serverReady(); i++;
		};
	};

	if (insight.serverReady()) {
		cout << "\nForcing server termination " << flush;
		insight.serverTerminate(true);
		cout << "\nServer Ready Status  : " << insight.serverReady() << "\n" << flush;
	}
	cout << "\n\nServer terminated\n" << flush;

	cout << "\n\nReady to exit..\n\tPress return to continue:" << flush;
	getline(std::cin, astring);
	// Read Cov Defs

	// Read Split state (for init or restart)

	// find best split
	//	for each covatiate.... 
	//		if monotonic, make list of splits
	//		if not monotonic, get insigth score, sort by score make lsit of splits
	//		Find best split
	//	Find covaraite who;sd best split is best...
	// execute split
	//	Create childeren alogn split lines
	//	Housekeeping
	// recurse...

	// Save & dump results.

	return 0;
}
#endif
