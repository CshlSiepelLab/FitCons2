#pragma once
// #define COVTREE_VER_STR "00.00.00" // Pull from executale

#include <iostream>
#include <sstream>
#include "stringer.h"
#include "FmtSimple.h"

struct covtreeArgs {
	typedef stringer string;
	string fnCovDef;				// 1			File containing definition of covariate and covariate sets.
	string insightExe;				// 2			Full path to Insight Executable....
	string insightDB;				// 3			Directory containing Insight DB.
	string baseDir;					// 4			Directory, source of tree decomposition recursion (Working Directory)
	int		v;						// -v			Verbosity. 1 is normal, 0 is quiet.. 2,3,4 provode increasing amounts of info...
	string insightArgs;				// -iargs		any special Insight2 arguments?
	string insightDBwrk;			// -idb         local, nonpersistable override for insight DB doe use on cluster /high eprformance systems with lcoa lmemory.
	int	   insightThreads;			// -iprocs		Number of instances of insight2 to run. 1 is min (default), 3 or 5 are generally most useful
	string fRefGen;					// -refgen		reference genome, usually HG19, all covariates must fully partition this space.
	string fnBedPosMask;			// -bedmask		bed file containing position mask, only these positions will be considered (presently ignored)
	string covStart;				// -covstart	serialzied covsetDefSubset object, start with these values or "" for all
	string sCTS;					// -scts		Tags for each cell type (comma seperated, no whitespace)
	stringer::vecS vcts;			//				vector of cell types from sCTS...
	string spriors;					// -ipriors		string encoded, "," seperated priors, or "", rho eta gamma N (generally N is empty and if eta is 0, default eta is used)
	stringer::vecS vpriors;			//				vectors of string encoded priors, or "", rho eta gamma N (generally N is empty and if eta is 0, default eta is used)
	string fnOut;					// -fout		file name for output recursion, "fOut.done" is hte completion semaphour, written to baseDir.
	int		extRecurse;				// -erecurse	0, 1 or 2 - write 0, 1, or 2 serialzied arg strigns to fnOut and terminate. Other combos are handeled using internal recursion...
	double	minBits;				// -minbits		Termination condition. Dont recurse if best split yeilds less than this many bits.
	int		maxDepth;				// -maxdepth	Termination condition. Don't recurse beyond this many levels. If <1, then program invocation does nothing...
	bool	implicitParRho;			// -implprho	Use implicit parent rho (true = sue weighted sum of child rhos, false = use iprior rho for parent, if provided.
	bool	splitOnRho;				// -rhosplit	Split on entropy of RHO rather than NLL of Insight data. Usign Insight removes the problem of inconsistent estimators and is included only for historical consistency
	string	debugS; 				// -debug		Opts: Comma sep list of options Ins,Chrom... etc..
	int		debugSkipInsight;		// Internal		0 skip nothing, 1 skip Singeltons, 2 Slip  Partitions (3 - skip both)
	string	debugOneChrom;			// Internal		Only process a singler chromosome, name
	string	restart;				// -restart     "Hi", "Lo", "Pr" rear resutls of previous partition and continue iteration wiht Hi (first) Lo (second) or both lines in the output (covtree.out) file- NOT SERIALIZED!!									
	string  dWrk;					// -dwrk		Working directory, if not "", replaces basdeDir for file storags but NOT for logical recursion (ie /hi /lo directories)m external shell muyst move files... only useful when erecurse=2
	bool	sargExit;				// -sarg        If only one more parameter, tht is taken as a serialized args object, if more arguments, Generates a serlization of args, and exits...
	double  hMinChildSz;		    // -hcpos		Heuristic constriant - dont consider splits if one child has fewer than hMinChildSz positions. defaults to 0.
	double  hMaxPerr1;			    // -hcerr1		Heuristic constriant - dont consider splits if the probability of a child ordering error is > then this. Defaults to 1 to disable. Values > .5 (or <=0) are meaningless.
	double  hMaxPerr2;			    // -hcerr2		Heuristic constriant - dont consider splits if the probability of that rhoc1-rhoc2 < (rhoc1_hat-rhoc2_hat)/2 is greater than this value, more stringent than hMaxPerr2...
	double  hMinDirInf;				// -hcmndinf    Minimum Directed Information, for rhoc1_hat>rhoc2_hat what is the expected inf gain for all rhoc1 > rhoc2, defaults to 0.
	// These are best handled by the caller in a script...
	// X - Script - icleanupStyle
	// X - Script - ... scleanupMethod (mv or rsync|args|src|dest)

	// The working base is Temp directory (perhaps in SSD) that is used for processing. 
	//	The shell msut rebase any needed files (like reatart source), and must cleanup and move the results 
	//	files from this temp dir to the baseDir. Output will reference BaseDir and its chidleren, but filed will
	//	be written to dirWrk() isntead, if it is non-null. dWrk is transient, it is not persisted via serilization.
	string dirBase() const { return baseDir; };
	string dirWrk() const { if (dWrk!="") return dWrk; return dirBase(); };

	string getInsightArgs() const {
		string iarg = insightArgs + " -maxsampw " + std::to_string(vcts.size());
		// Priors, from the parent distribution. Weight is fixed at around 300 pseudocounts. Args are MLStart,Refine?,Prior
		if ((vpriors.size()>0) && (vpriors[0] != "") && (FmtSimple::toDouble(vpriors[0]) > 0.0))  iarg += " -rho "   + vpriors[0] + ",1," + vpriors[0];
		if ((vpriors.size()>1) && (vpriors[1] != "") && (FmtSimple::toDouble(vpriors[1]) > 1e-6)) iarg += " -eta "   + vpriors[1] + ",1," + vpriors[1];
		if ((vpriors.size()>2) && (vpriors[2] != "") && (FmtSimple::toDouble(vpriors[2]) > 0.0))  iarg += " -gamma " + vpriors[2] + ",1," + vpriors[2];
		return iarg;
	};
		
	void setDebugOpts(const string &argstr )  {
		stringer::vecS toks; 
		debugSkipInsight=0; debugOneChrom="";
		if (argstr.length()<1) return;
		argstr.tokenize( toks, ",");
		if (toks.size()>0)	debugSkipInsight= std::atoi( toks[0].c_str());
		if (toks.size()>1)	debugOneChrom = toks[1].trim();
		return;
	}

	void clear() {
		v = 0;								// Default to no debugging
		covStart="";						// special value, ignore nothing.
		sCTS=""; vcts.clear();				// default is NO cts data.
		// reinstate these blank priors weh ner fix parrental values, in the mean time default to hg19 when no parrental priors are provided
		//spriors=""; 
		//vpriors.clear(); vpriors.resize(6);	// 4 blank priors?, use Insight2 default....
		spriors="0.07342666985087427378,0.00021737936739721811,0.64058097515779444109,347,125851001.911685,125851025.131036";	// HG19 valeusfrom INSIGHT1 0.16e, etas is earlier value for expecation
		vpriors.clear(); spriors.tokenize(vpriors, ","); vpriors.resize(6); for (auto &s : vpriors) s.trim();					// draw from arg processing code.....
		fnOut="covtree.out";
		extRecurse=0;						// use internal recursion only
		minBits = 100000.0;					// Fairly conservative...
		maxDepth=4;							// very shallow, for testing only.
		implicitParRho=false;				// If false & user provides rho<>0 prior, use it for parent. Else force Prho = (C1n*C1rho + C2n*C2rho)/(C1n +C2n).
		debugS="";							// Comma sep list of options Ins,Chrom... etc..
		setDebugOpts( debugS );
		restart="";
		sargExit=false;
		dWrk="";							// Not persistable
		insightThreads=1;					// Not persistable
		// debugSkipInsight=0;					// 0 skip nothing, 1 skip Singeltons, 2 Slip  Partitions (3 - skip both)
		// debugOneChrom="";					// Only process a singler chromosome, namel
		insightDBwrk="";					// no default workign direcotr, use insightDB
		hMinChildSz=0.0;
		hMaxPerr1=1.0;						// -hcerr1		Heuristic constriant - dont consider splits if the probability of a child ordering error is > then this. Defaults to 1 to disable. Values > .5 (or <=0) are meaningless.
		hMaxPerr2=1.0;						// -hcerr2		Heuristic constriant - dont consider splits if the probability of that rhoc1-rhoc2 < (rhoc1_hat-rhoc2_hat)/2 is greater than this value, more stringent than hMaxPerr2...
		hMinDirInf=0.0;						// -hcmndinf    Minimum Directed Information, for rhoc1_hat>rhoc2_hat what is the expected inf gain for all rhoc1 > rhoc2, defaults to 0.
		splitOnRho=false;					// -rhosplit	Override split on Insight2 NLL with Split on Rho.
#if defined(_WIN32) || defined(_WIN64)
		fnCovDef	= "C:/Temp/t/covtree/covedef.txt";
		insightExe  = "C:/Users/bg279.CORNELL/Source/Git/FitCons/INSIGHT-2/Insight2/x64/Release/Insight2.exe";
		insightDB   = "C:/Temp/t/db4";
		insightArgs = "-qmap";		// use posterior MAP estimate
		baseDir     = "C:/Temp/t/covtree/base";
		fRefGen	    = "C:/Temp/t/hg19.bed";
		fnBedPosMask="";
#else
		// unix values not yet defined
		insightExe	= "";
		insightDB	= "";
		insightArgs = "-qmap";		// use posterior MAP estimate
		baseDir		= "";
		fnCovDef	= "";
		fRefGen		= "";
		fnBedPosMask = "";
#endif
	};

	// machine readible list
	// Convert set of arguments into a single line string (no newlines)..
	string seralize() const {
		string t_priors="";
		// when serializing, take values from vector, NOT original string...
		for (const auto & s : vpriors ) t_priors += (t_priors.size()>0?",":"") + s;
		string q="\"", c="," , cq = c + q, qc = q + c, qcq=q+c+q, s="";
		s += q + fnCovDef + qcq +insightExe + qcq +insightDB+qcq + baseDir+ qc + "-v" + c + std::to_string(v) + (insightArgs.length()>0?c + "-iargs" + cq + insightArgs + q:"");
		s += c + (splitOnRho ? "-rhosplit" + c : "" ) + "-refgen" + cq + fRefGen + qc + "-bedmask"  + cq + fnBedPosMask + qc + "-covstart" + cq + covStart + qc + "-scts" + cq + sCTS + q;
		s += (t_priors.length()>0?c + "-ipriors" + cq + t_priors + q : "" ) + c + "-implprho" + c + (implicitParRho?"1":"0") + c +"-fout" + cq + fnOut + qc + "-erecurse" + c + std::to_string(extRecurse) + c + "-minbits" + c + std::to_string(minBits);
		s += c + "-maxdepth" + c + std::to_string( maxDepth) + c + "-debug" + cq + debugS + q ;
		s += (hMinChildSz	<= 0.0	?   "" : c + "-hcpos"	 + c + std::to_string(hMinChildSz));
		s += (hMaxPerr1		>= 1.0	?	"" : c + "-hcerr1"	 + c + std::to_string(hMaxPerr1));
		s += (hMaxPerr2		>= 1.0	?	"" : c + "-hcerr2"	 + c + std::to_string(hMaxPerr2));
		s += (hMinDirInf	<= 0.0	?	"" : c + "-hcmndinf" + c + std::to_string(hMinDirInf));
		return s;
	};
	// Convert a single line string into values...
	bool deseralize( const string &sblob) { 
		stringer::vecS toks; sblob.tokenizeQuoted(toks,", \t");
		return( processArgs( toks ) );
	};

	// human readable list
	string toStr() const { 
		std::ostringstream s; 
		s << "covtreeVersion:\tV" << COVTREE_VER_STR << "\n\n";
		s << "covtreeArgs:\n" ;
		s << "\tCovDefFile    \t" << fnCovDef << "\n";
		s << "\tInsight2      \t" << insightExe << "\n";
		s << "\tInsight2DB    \t" << insightDB << "\n";
		s << "\tBaseDir       \t" << baseDir << "\n\n";
		s << "\t-scts         \t" << sCTS << "\n\n";
		s << "covtreeOpts:\n";
		s << "\t-v            \t" << v << "\n";
		s << "\t-iargs        \t" << insightArgs << "\n";
		s << "\t-idb          \t" << insightDBwrk << "\n";
		s << "\t-refgen       \t" << fRefGen << "\n";
		s << "\t-bedmask*     \t" << fnBedPosMask << "\n";
		s << "\t-covstart     \t" << covStart << "\n";
		s << "\t-ipriors      \t" << spriors << "\n";
		s << "\t-implprho     \t" << (implicitParRho?"1":"0") << "\n";
		s << "\t-rhosplit     \t" << (splitOnRho?"1":"0") << "\n";
		s << "\t-iprocs       \t" << insightThreads << "\n";
		s << "\t-fout         \t" << fnOut << "\n";
		s << "\t-erecurse     \t" << extRecurse << "\n";
		s << "\t-minbits      \t" << minBits<< "\n";
		s << "\t-maxdepth     \t" << maxDepth << "\n";
		s << "\t-hcpos        \t" << hMinChildSz << "\n";
		s << "\t-hcerr1		  \t" << hMaxPerr1 << "\n";
		s << "\t-hcerr2		  \t" << hMaxPerr2 << "\n";
		s << "\t-hcmndinf     \t" << hMinDirInf << "\n";
		s << "\t-debug        \t" << debugS << "\n";
		s << "\t-restart      \t" << restart << "\n";
		s << "\t-dwrk         \t" << dWrk << "\n";
		s << "\t-sarg         \t" << sargExit << "\n";
		s << "\n\n";
		return(s.str());
	};

	covtreeArgs() { clear(); };

	bool processArgs(int C, char **V) {
		stringer::vecS vs;
		for (int i = 1; i < C; i++) vs.push_back( V[i] );
		return(processArgs(vs));
	}

	bool processArgs(stringer::vecS &S);

	static void displayArgs(int C, char *V[]) {
		for (int c = 1; c<C; c++) { std::cerr << " Arg: " << c << "\tval: " << V[c] << std::endl; };
	}

	void displayHelp(int C = 0, char *V[] = NULL, bool DisplayArgs = false);

};

// Proces input user arguments
bool covtreeArgs::processArgs(stringer::vecS  &V) {
	using std::cerr;
	using std::endl;
	for (auto &s : V) s = s.trim();
	V.push_back("-"); V.push_back("-");	// null arguments... this makes parsing easier
	auto C = V.size();
	if (C<5) return false;
	fnCovDef	= V[0];	// Covaraite definition file...
	insightExe	= V[1];	// Insight Executable
	insightDB	= V[2];	// Insight Database directory
	baseDir		= V[3];	// Base operating  directory
	//cerr << "\n\n Found " << C << " args " << endl;
	//for (ulong i=0;i<C;i++) cerr << i << "\t" << V[i] << endl;
	// process required arguments...
	unsigned int count = 4;
	while (count < C) {
		stringer p = V[count];
		// verbose
		if (p == "-" || p == "") {			// do nothing...
		} else if (p == "-v") {				// set verbosity
			v = 1;  if ((count+1) >= C) return (true);
			p = V[count + 1];		if (p.front() != '-') { v = std::stoi(p); count++; };
		} else if (p == "-h") {				// set verbosity
			return false;					// just print args and exit...
		} else if (p == "-iargs") {				// Any special Insight2 args? Simple -qmap is the usual default.
			if (++count >= C) return (false);
			insightArgs = V[count];
		} else if (p == "-rhosplit") {			// By default use data NLL under Insight. This overrides to use RHO entropy.
			splitOnRho = true;
		} else if (p == "-idb") {				// Override default insight db with local one, for multipele cluster threads. NOT PERSISTABLE
			if (++count >= C) return (false);
			p=V[count]; if ( p.front() == '-' ) return false;
			insightDBwrk = p;
		} else if (p == "-iprocs") {			// Any special Insight2 args? Simple -qmap is the usual default.
			if (++count >= C) return (false);
			insightThreads = std::stoi(V[count]);
		} else if (p == "-debug") {				// Debugging options are used to effect optimized code structure to identify problems...
			if (++count >= C) return (false);
			debugS = V[count]; setDebugOpts( debugS );
		} else if (p == "-refgen") {			// reference genome if different from hg19... afain file name of sorted bed file
			if (++count >= C) return (false);
			fRefGen = V[count];
		} else if (p == "-bedmask") {			// Limit analysis to these positions... string is name of a sorted, nonoverlapping,  bed file....
			if (++count >= C) return (false);
			fnBedPosMask = V[count];
		} else if (p == "-restart") {			// Restart recursion from previous output file "Hi" "Lo" "Pr" (pr = both) or "" for inactive..
			if (++count >= C) return (false);
			restart = V[count];
		} else if (p == "-covstart") {			// Initial covariates, Serialization of covsetDefSubsetObject, with "" means all. ANN:{CovA1;0,1,2},{CovA2:3,4,5}:CTS:{CovS1;0,1},{CovS2;1,2,3}
			if (++count >= C) return (false);
			covStart = V[count];
		} else if (p == "-dwrk") {				// Working directory, only suefel when erecurse = 2, physical files are written here, but paths in output & recursion are still referenced from baseDir
			if (++count >= C) return (false);
			dWrk = V[count];
		} else if (p == "-sarg") {				// serialize / deserialize argument list from command line. If next arg is - serialize args and exit, else accept next arg as serlaized arglist.
			if (++count >= C) return (false);
			if (V[count]=="-") {
				sargExit = true;
			} else {
				if (!deseralize(V[count]) ) return false;
			};
		} else if (p == "-scts") {				// Cell types for cell type sensative covariates
			if (++count >= C) return (false);
			sCTS= V[count]; vcts.clear(); if (sCTS.trim()!="") sCTS.tokenize(vcts,",");
		} else if (p == "-ipriors") {			// 1-Rho,2-Eta,3-Gamma,4-Pcounts,5-DataNLL,6-PosteriorNLL
			if (++count >= C) return (false);
			// NOTE any prior can be left blank by provideign a null CSV (ie 1,2,,4). Generally 4 is left blank. If 5/6 are missing, then DeltaInfof S is used to tedermine best split (bits), otherwise DeltaNLL from Insight is used...
			spriors = V[count]; spriors.tokenize(vpriors, ","); vpriors.resize(6); for (auto &s : vpriors) s.trim();
		} else if (p == "-implprho") {			// Cell types for cell type sensative covariates
			implicitParRho = true;				// optional second argument
			if ( (C>(count+1)) && ((V[count+1]=="1") || (V[count+1]=="0"))) { count++; implicitParRho = (V[count] == "1"); }
		} else if (p == "-fout") {				// Output file. If external recursion, serialzed args list goes here. (fnOut.done) is semaphour for completion.
			if (++count >= C) return (false);
			fnOut = V[count];
		} else if (p == "-erecurse") {			// Extternal Recursion Type write 0, 1, or 2 files to shell for external recursion, others handeled in thsi thread.
			if (++count >= C) return (false);
			extRecurse = std::atoi(V[count].c_str());
		} else if (p == "-hcpos") {				// Min number of sites in a each child needed to consider split
			if (++count >= C) return (false); p=V[count];
			if (p.startswith("-")) return false;
			hMinChildSz = std::stod( p );
		} else if (p == "-hcerr1") {			// MaxProb of ordering error in child rhos for allowable split - bonf corrected
			if (++count >= C) return (false); p = V[count];
			if (p.startswith("-")) return false;
			hMaxPerr1 = std::stod(p);
		} else if (p == "-hcerr2") {			// Max prob of distance error in chold rhos (is distance < 1/2 of this) - bonf corrected
			if (++count >= C) return (false); p = V[count];
			if (p.startswith("-")) return false;
			hMaxPerr2 = std::stod(p);
		} else if (p == "-hcmndinf") {			// Minimum expected directed info gain needed to to consider this split...
			if (++count >= C) return (false); p = V[count];
			if (p.startswith("-")) return false;
			hMinDirInf = std::stod(p);
		} else if (p == "-minbits") {			// Min number of its saved in order to recurse on childeren
			if (++count >= C) return (false);
			minBits = 0.0; minBits = std::stod(V[count]);
			if (minBits <= 0) {
				cerr << "-minbits X must have X > 0. Try 10.0 (for very small) to 1e6 (for very big).\n";
				return false; }
		} else if (p == "-maxdepth") {			// Maximum recursion depth
			if (++count >= C) return (false);
			maxDepth = 0; maxDepth = (int) std::stoul(V[count]);
			if (maxDepth == 0) {
				cerr << "-maxdepth must be > 0. Try 5 (for very small) to 20 (for very big).\n"; return false;
				return (false);
			}
		} else {									// unknown argument...
			cerr << "Unable to parse argument " << count << " :" << p << ":, exiting." << endl;
			return false;
		}
		count++;
	};
	return(true);
}

// block.bedg, poly.bedg, polyn.bedg, monoDB.db, monoDB.cohroms, monoDB.tags
void covtreeArgs::displayHelp(int C, char *V[], bool DisplayArgs) {
	using std::cerr;
	using std::endl;
	cerr << endl << endl << "covtree Version: " << COVTREE_VER_STR << endl;
	cerr << endl << "covtree fCovDef fInsight2 dInsight2DB dOutBase [args]" << endl << endl;
	cerr << "Required Arguments: " << endl;
	cerr << "\tfnCovDef      - Path to covariate definition file, see documentation for details." << endl;
	cerr << "\tinsightExe    - Path to Insight2 executable." << endl;
	cerr << "\tdInsight2DB   - Directory containing Insight2 database." << endl;
	cerr << "\tdOutBase      - Directory at top of hierarchy containing results." << endl;
	cerr << endl;
	cerr << "Options: " << endl;
	cerr << "\t-h                - Print this help menu and exit." << endl;
	cerr << "\t-v [level]        - Verbosity flag. 0 is quiet. 1 is default. Higher values are increasingly verbose up to ~5." << endl;
	cerr << "\t-iargs    [s]     - s is a quoted string containing any special Insight2 arguments. See Docs for Insight2. Defult is -qmap." << endl;
	cerr << "\t-idb      [d]     - Local Insight2Database Directory. Overrides Insight2DB argument, but is not persistable." << endl;
	cerr << "\t-iprocs   [u]     - Number of instances of Insight to run in parallel. Each uses 3-6 GB of memory. Defult is 1." << endl;
	cerr << "\t-refgen   [f]     - f is a psth to a sorted bedfile defining refernce genome limits. Hg19 is default. Covaraites must fully parition this reference space." << endl;
	cerr << "\t-bedmask  [f]     - f is pathname to a sorted bedfile defining subset of genomic positions. Only operate on these positions (not yet implimented)." << endl;
	cerr << "\t-covstart [s]     - Begin with covaraite subset. Serialized covsetDefSubset object. See Docs. Generally like ...Anno:{CvnA1;1,2};{CvnA2;3,4,5}:CTS:{CvnS1;1,2};{CvnS2;4,5,8}." << endl;
	cerr << "\t-scts     [s]     - Comma seeperated list of cell-type names." << endl;
	cerr << "\t-ipriors  [s]     - Comma seeperated list of floats rho,eta,gamma,N,DatNLL,PosteriorNLL any can be mising and will default to hg19. N is pseudocounts." << endl;
	cerr << "\t-implprho [d]     - Alone or followed by 1 sets implicit-parent-Rho mode. ParRho = weighted sum of child rho. If 0 (default) use -iprior rho as parent rho, if it is provided, >0, and <1." << endl;
	cerr << "\t-rhosplit         - if present, split on rho entropy, rather than default of Insight2NLL. " << endl;
	cerr << "\t-fout     [s]     - Nale of output file for external recursion. Written to dOutBase/[s]. [s].done ic completion semaphor file." << endl;
	cerr << "\t-erecurse [u]     - write [u]=0, 1 or 2 children to fout for external recursions. Others are calculated serially in this process." << endl;
	cerr << "\t-minbits  [d]     - Terminate recursion when best split yeilds fewer than this many bits. Def is 100,000. Must be positive. 100 is small, 1e6 is big." << endl;
	cerr << "\t-hcpos    [u]     - Minimum number of positions in each child class needed to consider a split. Defaults to 0. Generally 300-10,000 is good. 1000 is niave est." << endl;
	cerr << "\t-hcerr1   [d]     - Maximum Bonferroni Corrected probability of ordering error in childs. Defaults to 1.0 (disabled) 0.05 is good est." << endl;
	cerr << "\t-hcerr2   [d]     - Maximum Bonferroni Corrected probability of distance error in childs (rho1-rho2)<(rho1_hat-rho2_hat)/2. Defaults to 1.0 (disabled) 0.05 is good est." << endl;
	cerr << "\t-hcmndinf [d]     - If >0 use expected directedInformation to determine ordering: E[DeltaInf(rho1,rho2)&(rho1>rho2)] Instead of point estimate DeltaInf(rho1_hat,rho2_hat). Value is min inf." << endl;
	cerr << "\t-maxdepth [u]     - Terminate recursion no later than this depth. If 0, program does nothing and exits." << endl;
	cerr << "\t-restart  [s]     - Restart recursion from previous covtree.out file in OutBase. Creats Hi/Lo sub dires. String Vals = Hi, LO or PR (both) children." << endl;
	cerr << "\t-dwrk     [s]     - Write all physical files in this directory, however recursion logically proceeds on dirbase. Shell muse move files, so only useful when erecurse = 2." << endl;
	cerr << "\t-debug    [d,[s]] - Debugging options d=1-skip singelton,2-skip partiton,3-skip both,0-skip none; s=process only single chromosme with tag s." << endl;
	cerr << "\t-sarg     s       - If s=- print serialized args and exit, else deseralize next arg as arglist." << endl; 
	cerr << endl;
	if (DisplayArgs) {
		cerr << "Args:" << endl;
		covtreeArgs::displayArgs(C, V);
		cerr << endl << endl;
	}
	return;
}

