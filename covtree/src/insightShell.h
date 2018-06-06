#pragma once
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <stdio.h>
#include <fcntl.h>
#include <string>
#include <chrono>
#include <thread>
#include <map>
#include <vector>

#include "fastfile.h"
#include "stringer.h"

// Open an instance of Insight2, feed it commands. Wait for results etc.... 

class insightShell
{
public:
	typedef stringer string;
	class insightResults {
	public:
		struct dualdat { double d; stringer s; dualdat() {d=0;s="";}; };
		typedef std::map<string,dualdat> elems_t;
		typedef std::map<string,elems_t> clases_t;
	private:
		clases_t vals;
		double safestod( const string &s, bool *ok=NULL) { double d=0; try {d=std::stod(s);} catch (...) { if (ok!=NULL) *ok=false; }; return(d); };
	public:
		void clear() { vals.clear(); };
		insightResults() {clear(); };
		~insightResults() {clear(); };
		bool getField(const string &group, const string &field, string &st) const {
			st="";
			if (vals.count(group) < 1) return false;
			if (vals.at(group).count(field)<1) return false;
			st = vals.at(group).at(field).s;
			return true;
		}
		bool getField(const string &group, const string &field, double &value) const {
			value=0;
			if (vals.count(group) < 1) return false;
			if ( vals.at(group).count(field)<1) return false;
			value=vals.at(group).at(field).d;
			return true;
		}
		const clases_t &iter() { return( vals ); };
	private:
		bool read_try(const string &fnbed ) {
			clear(); string fin; string aline; elems_t elems;  string elemsname="";
			fin = fnbed.substr(0,fnbed.find_last_of("."))+".model";
			// std::cerr << "Reading results from " << fin << "\n";
			std::fstream sin; sin.open(fin, std::ios_base::in);
			if (!sin.good()) { sin.close(); return false; };
			while (sin.good()) {
				std::getline(sin,aline); if (aline.length()==0) continue;
				if (aline == aline.asupper()) {
					// std::cerr << "\tHeader" << aline << "\n";
					if (elems.size()>0) { vals[elemsname]=elems;}
					elems.clear(); elemsname=aline.substr(0,aline.length()-1);
				} else {
					stringer::vecS toks;
					// std::cerr << "\t\tRead   : " << aline << "\n";
					aline.tokenize(toks,":"); if (toks.size()<2) {sin.close(); return false; };
					{
						stringer key, val;
						key = toks[0].trim();
						val = toks[1].trim();
						// std::cerr << "\t\tParsed : " << toks[0] << " :: " << toks[1] << "\n";
						elems[key].d = safestod( val );
						elems[key].s = stringer( aline.substr(aline.find_first_of(":")+1)).trim();
						// std::cerr << "\t\tStored : " << elems[toks[0].trim()].d << " :: " << elems[toks[0].trim()].s << "\n";
					}
				}
			};
			if (elems.size()>0) { vals[elemsname] = elems; }
			sin.close(); return true;
		}

		bool read_ok() const {
			double a;
			//Last field in the file - this should insure a good read...
			if (!getField("PRIORN_MODEL", "Beta3", a)) return false; if (a<0) return false;
			// Parameters we use in this program
			if (!getField("PARAMETERS", "Rho", a)) return false;	if (a<0) return false;
			if (!getField("PARAMETERS", "Eta", a)) return false;	if (a<0) return false;
			if (!getField("PARAMETERS", "Gamma", a)) return false;	if (a<0) return false;
			// Other values used in theis program
			if (!getField("ARGS", "maxSampW", a)) return false;		if (a<0) return false;
			if (!getField("STATS", "InPosW", a)) return false;		if (a<0) return false;
			if (!getField("STATS", "InfPosWM", a)) return false;		if (a<0) return false;
			if (!getField("STATS", "InfPosWH", a)) return false;		if (a<0) return false;
			if (!getField("STATS", "InfPosWL", a)) return false;		if (a<0) return false;
			return true;
		}
	public:
		bool setField(const string &group, const string &field, const double &v ) {
			// the map [] operator inserts objects as necessary.
			vals[group][field].d = v; vals[group][field].s = std::to_string(v);
			return true;
		}
		bool setField(const string &group, const string &field, const string &s ) {
			// the map [] operator inserts objects as necessary.
			vals[group][field].d = 0; vals[group][field].s = s;
			return true;
		}
		bool read(const string &fnbed) {
			// Sadly, the network filesystem does not serialzie reads, with writes from other processes (say Insight2),
			//	even under sync(). So we have to resort to retries and heuristic timeouts... sigh... apologies.
			int tries = 0; const int MAX_TRYS = 60;
			if (read_try(fnbed) && read_ok()) return true;
			while (++tries<=MAX_TRYS) {
				fastfile::fflush();
				fastfile::sleepMillis(1000);
				if (read_try(fnbed) && read_ok()) return true;
			}
			return false;
		}

	};

private:
	enum fstatus { FOK, FFAIL, FERROR, FEXISTS };
	bool fileExists( const string &fname ) const { return( fastfile::ffexists(fname.c_str()) ); };
	fstatus makeDir(const string &fname) const { return( (fastfile::ffmkdir( fname.c_str() ) == 0)? FOK : FERROR); };
	fstatus moveFile( const string &src, const string &dst, bool force=false) const {
		if (!fileExists(src)) return( FERROR );
		if (force && fileExists(dst)) { fastfile::ffrmf(dst); };
		return(fastfile::ffmv(src.c_str(), dst.c_str())==0? FOK : FERROR ); 
	};

	// convert a job init strign to a semaphore fie name...
	string jobToSem(const string &job) const {
		stringer fsem = job.substr(0, job.find_first_of("#"));	// get rid of job arguments
		fsem = fsem.trim();										// get rid of excess whitespace
		fsem = fsem.substr(0, fsem.find_last_of(".")) + ".done";//change suffix from .bed to .done
		return fsem;
	}

	string logFile="";
	FILE *serverfin = NULL;
	// Hide OS specific details for popen / _popen
	// If we really need a hard KILL on this, open _pipes, use _dup, use _spawn, then swap pipes back. This gives the PID
	//	of the child AND provides a pipe for writing (and reading, if desired), while staying within the POSIX subsystem.
	FILE *lpopen(const string &fcmd, const string &flog) {
#ifdef _WIN32 
		string acmd=fcmd+" > " + '"' + flog + '"' + " 2>&1";	// sindows requires sderr redirection AFTER command.
		serverfin = _popen(acmd.c_str(),"w");
#else
		string acmd="( " + fcmd + " 2>&1 ) > " + '"' + flog + '"';
		serverfin  = popen(acmd.c_str(), "w");
#endif
		return  serverfin;
	}
	void lclose(FILE **handle) {
		if (handle==NULL) return;
		if (*handle == NULL) return;
#ifdef _WIN32 
		_pclose(*handle); *handle=NULL;
#else
		pclose(*handle); *handle = NULL;
#endif
		return;
	}
public:
	static void msleep(unsigned int msecs) { std::this_thread::sleep_for(std::chrono::milliseconds(msecs)); };

	insightShell() { serverfin=NULL; logFile=""; };
	~insightShell() { serverTerminate( true );	};

	// Start the server OS independent code
	bool serverInitialize( const string &ExePath, const string &DBdir, const string &LogFile, const string &DefArgs = "-qmap") {
		string fcmd=ExePath + " " + '"' + DBdir + '"' + " " + DefArgs;
		if (serverStarted()) return false;
		serverfin = lpopen( fcmd, LogFile ); 
		return serverStarted();
	};	

	bool serverStarted() const { 
		if (serverfin==NULL) return false;
		if (ferror(serverfin) || feof(serverfin)) return false; ;
		if (fputs("",serverfin)==EOF) return false;
		return true;
	};
	bool serverReady() const	{ return (serverStarted()); };		// Check the log file to see if we are ready to go...
	const string & serverGetLogName() const { return logFile; };	// Grab the processign log file name

	// request termination, or kill the server process (force = true)
	bool serverTerminate(bool force = false) {
		if (!serverStarted()) { lclose(&serverfin); return true; }
		fprintf(serverfin,"done\n"); fflush(serverfin); 
		if ( force ) {
			msleep(2000); lclose(&serverfin); msleep(2000);
		}
		return( ! serverStarted() );
	}; 

	// Initiate a request
	bool jobInit(const string &request, bool debug = false) {
		if (!serverStarted()) return false;
		// Identify the input filename
		stringer fnflag = jobToSem( request );
		// Clear the "done" semaphour, if it exists...
		if (!debug){ 
			int stat=0;
			// std::cerr << "COVTREE-INSIGHT DEBUG - Semaphore " << fnflag << "\n";
			fastfile::ffrmf(fnflag);
			// Write he input file to the job Queue pipe, extended format to chance default args is is "Fname # args"
			//	ie "hg119.bed # -qmap -v 4"
			// std::cerr << "COVTREE-INSIGHT DEBUG - Command   " << request << "\n";
			stat = fputs(request.c_str(), serverfin);	if (stat == EOF) return false;
			// std::cerr << "COVTREE-INSIGHT DEBUG - status    " << stat << " ";
			stat = fputs(request.c_str(),serverfin);	if (stat == EOF) return false;
			// std::cerr << stat << " ";
			stat = fputs("\n", serverfin);				if (stat == EOF) return false;
			// std::cerr << stat << " ";
			stat = fflush(serverfin);					if (stat != 0) return false;
			// std::cerr << stat << "\n";
			if (stat != 0) return false;
		};
		return true;
	};  

	// Asynch check to see if request is done
	bool jobPoll(const string &finName) const {
		string fdone=jobToSem(finName);
		return( fileExists(fdone));
	};	

	// Wait for request to complete, or timeout. True if done, false on not done (including timeout)
	bool jobWait(const string &finName, uint32_t timeoutMS=3000 ) const {
		static const uint32_t delta = 500; // miliseconds
		string fsem = jobToSem(finName);
		uint32_t waited=0; 
		while ((waited < timeoutMS) || (timeoutMS==0)) {
			if (jobPoll(fsem)) return true;
			msleep( delta ); waited+=delta; };
		return jobPoll(fsem);
	};	

	// Read results from .model file
	bool jobResultsGet(const string &finName, insightResults &Res ){
		string modfname=jobToSem(finName);
		modfname = modfname.substr(0,finName.find_last_of("."))+".model";	// replace .bed with .model to get input file name
		return( Res.read( modfname ) );
	};		

	// Move all results to a target directory
	bool jobResultsMove(const string &finName, const string &destDir, bool force = true ) {
		bool ok = true;
		makeDir( destDir );	// Ignore results, this fails quietly if dest exists...
		string pbase = jobToSem(finName);
		pbase = pbase.substr( 0, pbase.find_last_of("."));					// remove the trailing suffix, .bed from input file
		string fbase = pbase.substr( pbase.find_last_of("/")+1);			// get basename for data set
		ok &= ( moveFile( pbase + ".bed",    destDir + "/" + fbase + ".bed",    force ) == FOK );		// move source file
		ok &= ( moveFile( pbase + ".insres", destDir + "/" + fbase + ".insres", force ) == FOK );		// move results summaryfile
		ok &= ( moveFile( pbase + ".model",  destDir + "/" + fbase + ".model",  force ) == FOK );		// move full model
		ok &= ( moveFile( pbase + ".done",   destDir + "/" + fbase + ".done",   force ) == FOK );		// move semaphoure
		return ok;
	};		
};

class insightShellSet : public insightShell {
private:
	std::vector< insightShell * > inst;
	int last_inst=0;
public:
	insightShellSet() {}
	~insightShellSet() { serverTerminate(true); }

	bool serverInitialize(const string &ExePath, const string &DBdir, const string &LogFile, const string &DefArgs = "-qmap", const int NumThreads=1 ) { 
		if (inst.size() > 0 ) return false;
		if (NumThreads<1) return false;
		inst.resize( NumThreads ); for (int i=0;i< (int) inst.size();i++) inst[i]=NULL;
		inst[0]=this;
		for (int i=1; i<NumThreads;i++) { 
			inst[i] = new insightShell(); };
		for (int i=0; i<NumThreads; i++) {
			string lf =  LogFile + "." + std::to_string(i);
			if (!inst[i]->insightShell::serverInitialize(ExePath, DBdir, lf, DefArgs)) {
				serverTerminate(true); 	return false;	}
		}
		last_inst=(NumThreads-1);
		return true;
	}

	bool serverTerminate(bool force = false) {
		bool ok = true, b;
		for (int i = 0; i<numThreads(); i++) {
			if (inst[i] != NULL) { 
				b=inst[i]->insightShell::serverTerminate(force); 
				ok = ok && b; }
		}
		if (force) {
			for (int i = 0; i<numThreads(); i++) {
				if (i!=0) delete inst[i];
				inst[i] = NULL;
			}
			inst.clear();
		};
		return ok;
	}

	bool jobInit(const string &request, bool debug = false, const int Thread = -1) {
		int tid=Thread;	// round robin scheduling, by default
		// If user selects a specific thread, dispatch job to it, otherwise use round robin
		if (Thread<0 || Thread>=numThreads()) tid=-1;
		if (tid < 0) tid=(++last_inst % numThreads() );
		return( inst[tid]->jobInit( request, debug) );
	}

	int numThreads() const  { return((int)inst.size()); };

	const string serverGetLogName( int ThreadID = 0 ) const { if (ThreadID<0 || ThreadID>=numThreads()) return(string("")); return( inst[ThreadID]->serverGetLogName() ); };
};
