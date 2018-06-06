#pragma once

#include "modelInsight.h"
#include "modelInsightAccelerator.h"

// Class for calculating posterior distributions over parameters rho, eta, gamma
//	Based on data likelihood, and priors drawn from global population distributions.

struct modelInsightPosterior : modelInsight {
	typedef std::vector<double>	vecD;
	typedef butils::ulong ulong;
	typedef vecD::size_type vecDi;

	// precomputed cofactors for each term. These could likely be converted to floats
	struct postElem {
		double rho, eta, gam;					// parameter values
		double drho, deta, dgam, sz, szlog;		// dimensions of cell (differential)
		double priorRaw, priorDens, priorAbs;	// the raw, density (for graphing) and absolute (normalized) prior for this cell
		double lik;								// data negative log likelihood at rho, eta, gamma.
		double posterior, posteriorDens;		// normalized posterior, and posterior density (for graphing)
		void clear() { rho=eta=gam=drho=deta=dgam=sz=szlog=priorRaw=priorDens=priorAbs=lik=posterior=posteriorDens=0; }
		postElem() { clear(); };
		inline string f( double v, const char *fmt="%.6lf" ) { static char buf[32]; sprintf(buf,(v>=1E-4?fmt:"%.2le"),(double)v); return( string(buf) ); };
		string header() { return( (string) "rho\teta\tgamma\tprior\tlogLik\tposterior\tcellZz\tpriorRaw\tpriorDens\tpostDens\tdRho\tdEta\tdGamma"); };
		inline string toStr(const char *s = "%.6lf") {
			string a = f(rho,s) + "\t" + f(eta,s) + "\t" + f(gam,s) + "\t" +  f(priorAbs,s) + "\t" + f(lik, s) + "\t" + f(posterior, s) + "\t";
			a += f(sz,s) + "\t" + f(priorRaw, s) + "\t" + f(priorDens, s) + "\t" + f(posteriorDens, s) + "\t";
			a += f(drho, s) + "\t" + f(deta, s) + "\t" + f(dgam, s); return(a);
		};
	};

	struct postSet {
		typedef std::vector<vecD> vvD;
		std::vector<postElem> elements;
		vecD gridRho, gridEta, gridGam;				// sampling grids, finer near ML/MAP values - actual parameter values at center of cell, not log.
		vecD margRhoPri, margEtaPri, margGamPri;	// Marginal parameters priors, NLL
		vecD margRhoPos, margEtaPos, margGamPos;	// Marginal posterior probability, NLL
		vvD  margRhoLik, margEtaLik, margGamLik;	// Marginal posterior probability, NLL
		vecD margRhoD,	 margEtaD,	 margGamD;		// Marginal probability element size, for calculating densities...
		double ci = .05;							// cutoff for credal interval for expectations, Ok to set to .001, .01, .05....
		parameters	exp;							// expected values
		parameters	expCIm;							// lower credal interval bound
		parameters	expCI;							// median value
		parameters	expCIp;							// upper credal interval bound
		ezmatrix priorObs;							// population distribution of input states, so P(param) \porp [P(X|param)]^N, for N pseudo counts.
		void clear() { elements.clear(); exp.clear(); priorObs.clear(); margRhoD.clear();	// individual objects...
			gridRho = gridEta = gridGam = margEtaD = margGamD = margRhoPos = margEtaPos = margGamPos = margRhoPri = margEtaPri = margGamPri = margRhoD;	// all the vecD
			margRhoLik.clear(); margEtaLik.clear(); margGamLik.clear();			// the vecvecD
		};
		inline string f(double v, const char *s = "%.6lf")       const { static char buf[32]; const char *s2 = ((v!=0)&&fabs(v)<1e-4?"%.3le":s); sprintf(buf, s2, (double)v); return(string(buf)); };
		inline string fv(const vecD &v, const char *s = "%.6lf") const { string a; for (ulong i=0;i<v.size();i++) a += (i>0?"\t":"") + f(v[i],s); return a;  };
		inline static vecD vv2v( const vvD &vv) { vecD v(vv.size()); for (ulong i=0;i<vv.size();i++) v[i]=vv[i][0]; return(v); };
		inline string toStr(const char *s = "%.6lf") { string a=""; vecD tmp;
			a += "CredalInterval:\t" + f(1.0-ci) + "\n";
			a += "Expec   Rho:   \t" + f(exp.rho)    + "\tEta:\t" + f(exp.eta)    + "\tGamma:\t" + f(exp.gam)    + "\n";
			a += "CredMin Rho:   \t" + f(expCIm.rho) + "\tEta:\t" + f(expCIm.eta) + "\tGamma:\t" + f(expCIm.gam) + "\n";
			a += "Median  Rho:   \t" + f(expCI.rho ) + "\tEta:\t" + f(expCI.eta ) + "\tGamma:\t" + f(expCI.gam)  + "\n";
			a += "CredMax Rho:   \t" + f(expCIp.rho) + "\tEta:\t" + f(expCIp.eta) + "\tGamma:\t" + f(expCIp.gam) + "\n";
			a += "RHO\nGrid:\t"   + fv(gridRho) + "\n"; a += "Prior:\t" + fv(margRhoPri) + "\n";  a += "LigLik:\t" + fv(vv2v(margRhoLik)) + "\n";  a += "Post:\t" + fv(margRhoPos) + "\n"; a += "Dens:\t" + fv(margRhoD) + "\n\n";
			a += "ETA\nGrid:\t"   + fv(gridEta) + "\n"; a += "Prior:\t" + fv(margEtaPri) + "\n";  a += "LigLik:\t" + fv(vv2v(margEtaLik)) + "\n";  a += "Post:\t" + fv(margEtaPos) + "\n"; a += "Dens:\t" + fv(margEtaD) + "\n\n";
			a += "GAMMA\nGrid:\t" + fv(gridGam) + "\n"; a += "Prior:\t" + fv(margGamPri) + "\n";  a += "LigLik:\t" + fv(vv2v(margGamLik)) + "\n";  a += "Post:\t" + fv(margGamPos) + "\n"; a += "Dens:\t" + fv(margGamD) + "\n\n";
			return( a );
		}


	};
	// create a grid with good overall coverage, and more detail near a focus....
	struct gridder {
		typedef std::vector<double>	vecD;
		typedef butils::ulong		ulong;
		struct focus { 
			double fmin, fmax; ulong res; bool logspace;
			void loadFocus( double focus, double width, ulong Fres, bool DoLog=false)	{ fmin = focus - (width/2); fmax= focus + (width/2); res = Fres; logspace=DoLog;};
			void loadRange(double RMin, double RMax, ulong Rres, bool DoLog=false)		{ fmin = RMin; fmax = RMax; res = Rres; logspace=DoLog;};
			void clear() { fmin = fmax =0.0; res = 0; logspace=false;};
		};
		inline static double ratToFrac(double rat) { return rat / (rat + 1.0); };		// convert between ratio (a/b) [unbounded] and fraction (a/(a+b)) [0,1]
		inline static double fracToRat(double frac) { return (frac == 1.0 ? std::numeric_limits<double>::max() : frac / (1.0 - frac)); };
		inline static double elemSize(ulong I, const vecD &grid) {
			double r = 0; auto imax = grid.size(); if (imax <= 0 || I >= imax) return r; --imax;
			if (I > 0)		r += 0.5*(grid[I] - grid[I - 1]);
			if (I < imax)	r += 0.5*(grid[I + 1] - grid[I]);
			return r;
		}
		static bool makeGrid(const std::vector<focus> &gfoc, vecD &grid) {
			double gmin = gfoc[0].fmin, gmax = gfoc[0].fmax, v=0; grid.clear();
			for (vecDi  ilev = 0; ilev < gfoc.size(); ilev++) {																		// loop over each focus
				const focus &f = gfoc[ilev]; double delta = (f.fmax - f.fmin) / f.res;
				if (!f.logspace) {	// linear spacing
					for (ulong i = 0; i <= f.res; i++) { v = f.fmin + i * delta;  if (v >= gmin && v <= gmax) grid.push_back(v); };		// add finer grids
				} else {			// lograthmic (base 2)  sizing...
					unsigned long res=(f.res+1)/2; double cent=(f.fmax+f.fmin)/2.0, delt = fabs(f.fmax-cent);
					for (ulong i=0; i<res;i++) { 
						v = cent - delt; if (v >= gmin && v <= gmax) grid.push_back(v); 
						v = cent + delt; if (v >= gmin && v <= gmax) grid.push_back(v);
						delt*=0.5;	// halv the interval, could use .75 if .5 is too extreme...
					};
				}
			}
			std::sort(grid.begin(), grid.end()); grid.erase(std::unique(grid.begin(), grid.end()), grid.end());						// sort and remove duplicates
			return true;
		};
	};

	bool initialized; 

	modelInsightPosterior(const modelInsight &M) : modelInsight( M ) { initialized = true; };

	void clear() { initialized = false;  modelInsight::clear();  };

	void generateSamplingGrid( postSet &ppost, ulong NRes = 20 ) const ;
	bool calcParamPost(const BlockSet_t &bs, double PriorWeight, const modelInsight &Priors, postSet &ppost, ulong gres = 20, double CredalThresh = .05, ulong verbose = 0, int NThreads = 0);
	static std::string strtime() { std::time_t ctt = std::time(0);	char *savetime = std::asctime(std::localtime(&ctt)); savetime[strlen(savetime) - 1] = '\0'; return(std::string(savetime)); }
};

void modelInsightPosterior::generateSamplingGrid(postSet &ppost, ulong NRes ) const {
	std::vector<gridder::focus>	focs;
	gridder::focus			f;

	//  Generate Sampling Grid: 
	//   WARNING calculation time grows cubically with grid res, so be careful. <10 is not reccomended, nor is > 100. 10 or 20 are likely good values.
	// Rho....
	double eps = 1.0e-9, st = 0.0 + eps, en = 1.0 - eps;	// keep us just off the boundaries....
	f.loadRange(st,en,NRes);					focs.push_back(f);		// Span range 0-1          at resolution 0.050
	f.loadFocus(model.rho, 0.200, NRes);		focs.push_back(f);		// span range ML +/- 0.200 at resolution 0.010
	f.loadFocus(model.rho, 0.020, NRes);		focs.push_back(f);		// span range ML +/- 0.020 at resolution 0.001
	gridder::makeGrid( focs, ppost.gridRho);	focs.clear();

	en = 10.0 - eps;	// keep us just off the boundaries....
	{
		bool lg = (model.eta < .01); ulong nres2=NRes*(lg?2:1); 
		f.loadRange(st, en, NRes);		 				focs.push_back(f);		// Span range 0-10 at resolution .50
		f.loadFocus(model.eta, 2.000, NRes,lg); 		focs.push_back(f);		// span range ML +/- 1.000 at resolution .100
		f.loadFocus(model.eta, 0.200, nres2,lg);		focs.push_back(f);		// span range ML +/- 0.100 at resolution .010
		f.loadFocus(model.eta, 0.020, nres2,lg);		focs.push_back(f);		// span range ML +/- 0.010 at resolution .001
	}
	gridder::makeGrid(focs, ppost.gridEta);	focs.clear();

	f.loadRange(st, en, NRes);					focs.push_back(f);
	f.loadFocus(model.gam, 2.000, NRes);		focs.push_back(f);
	f.loadFocus(model.gam, 0.200, NRes);		focs.push_back(f);
	f.loadFocus(model.gam, 0.020, NRes);		focs.push_back(f);
	gridder::makeGrid(focs, ppost.gridGam);	focs.clear();

	return;
}

bool modelInsightPosterior::calcParamPost(const BlockSet_t &bs, double PriorWeight, const modelInsight &Priors, postSet &ppost, ulong gres, double CredalThresh, ulong verbose, int NThreads ) {
//bool calcParamPost(const inscomp::optInsight::parameters &model, const inscomp::optInsight::parameters &priorparam, const BlockSet_t &bs, const Args &a, butils::ulong gres, paramPosteriors &ppost) {
	//using inscomp::optInsight;

	// if grid has been set manually, set gres =0. Then we don't clear the grid, and use what is alread ythere...
	if (gres > 1) { 
		if ( verbose > 0 ) std::cout << strtime() << "\t\tCalcparampost - Generating sampling Grid" << std::endl;
		ppost.clear();  
		generateSamplingGrid(ppost,gres); }
	
	if (verbose > 0) std::cout << strtime() << "\t\tCalcparampost - Generating Unnormalized priors " << std::endl;

	// Get the population-wide expected input distribution, use as basis for prior...
	Priors.likelihoodMatrix(bs.numAlleles, true, ppost.priorObs);	// TODO: this may be redundant, but it is relatively low overhead...
	// int nthreads = a.threads;
	postElem e;
	if (NThreads <2) {
		// The fastest way to do this is just to use accelerator to generate posterior... but we don't get to seperate prior and posterior...
		//	So we do it twice and gerneateexplicit priors and posteriors...
		if (verbose > 0) std::cout << strtime() << "\t\tCalcparampost - Initializing likelihood accelerator" << std::endl;
		modelInsightAccelerator a_prior(*this), a_like(*this);		// accelerators for prior and likelihood... note init prior accel with THIS. Prior model just used to generate pseudocounts...
		BlockSet_t no_data;											// allocate a null data set, so we can get at priors
		no_data.ancPri = bs.ancPri; no_data.numAlleles = bs.numAlleles;	// TODO these should be requried for constructor!!
		a_prior.initialize(no_data, PriorWeight, &Priors);			// just get priors. PriorWeight Pseudocounts of posterior distrib based on Priors...
		a_like.initialize( bs );									// just get likelihood of data, manually combine Prior with Like to get Posterior

		// generate unnormalized priors, element sizes and likelihoods
		if (verbose > 0) std::cout << strtime() << "\t\tCalcparampost - Generating likelihoods." << std::endl;
		parameters params; params.clear();
		double drho, deta,dgam;
		for (vecDi irho = 0; irho < ppost.gridRho.size(); irho++) {
			params.rho = e.rho = ppost.gridRho[irho]; e.drho = drho = gridder::elemSize((ulong)irho, ppost.gridRho);
			if (verbose > 1) std::cout << "\r" << strtime() << "\t\t\tCalcparampost - Generating Likelihood. Irho= " << irho << " of " << ppost.gridRho.size() << std::flush;
			for (vecDi ieta = 0; ieta < ppost.gridEta.size(); ieta++) {
				if (verbose > 2) std::cout << "\r" << strtime() << "\t\t\t\tCalcparampost - Generating Likelihood. Irho= " << irho << " / " << ppost.gridRho.size() << "  Ieta= " << ieta << " / " << ppost.gridEta.size() << std::flush;;
				params.eta = e.eta = ppost.gridEta[ieta]; e.deta = deta = gridder::elemSize((ulong)ieta, ppost.gridEta);
				for (vecDi igam = 0; igam < ppost.gridGam.size(); igam++) {
					params.gam = e.gam =ppost.gridGam[igam]; e.dgam = dgam = gridder::elemSize((ulong)igam, ppost.gridGam);
					a_prior.setparameters( params ); a_like.setparameters( params );
					e.priorRaw = -a_prior.calcNLL( ); e.lik = -a_like.calcNLL(); e.sz = e.drho * e.deta * e.dgam; e.szlog = mylog( e.sz ); e.priorAbs = e.priorRaw + e.szlog;
					ppost.elements.push_back(e);
				}
			}
		}
	} else {
		// The fastest way to do this is just to use accelerator to generate posterior... but we don't get to seperate prior and posterior...
		//	So we do it twice and gerneateexplicit priors and posteriors...
		if (verbose > 0) std::cout << strtime() << "\t\tCalcparampost - Initializing MT likelihood accelerator" << std::endl;

		// Contains all working data for a single thread
		struct t_OpUnit {
			t_OpUnit(modelInsightPosterior &t, postSet &tppost, double PriorWeight, const modelInsight &Priors, const BlockSet_t &bs) : a_prior(t), a_like(t), ppost( tppost) {
				irho=ieta=igam=0; drho=deta=dgam=0.0; params.clear(); e.clear(); 
				params = t.getparameters();									// important for setting beta, lamda theta etc...
				BlockSet_t no_data;											// allocate a null data set, so we can get at priors
				no_data.ancPri = bs.ancPri; no_data.numAlleles = bs.numAlleles;	// TODO these should be requried for constructor!!
				a_prior.initialize(no_data, PriorWeight, &Priors);			// just get priors. PriorWeight Pseudocounts of posterior distrib based on Priors...
				a_like.initialize(bs);									// just get likelihood of data, manually combine Prior with Like to get Posterior
			};
			modelInsightAccelerator a_prior, a_like;
			vecDi irho, ieta, igam;
			parameters params;
			double drho, deta, dgam;
			postSet &ppost;
			postElem e;
			static void CalcOpS(t_OpUnit *u) { u->CalcOp(); };
			void CalcOp( void ) {
				params.rho = e.rho = ppost.gridRho[irho]; e.drho = drho = gridder::elemSize((ulong)irho, ppost.gridRho);
				params.eta = e.eta = ppost.gridEta[ieta]; e.deta = deta = gridder::elemSize((ulong)ieta, ppost.gridEta);
				params.gam = e.gam = ppost.gridGam[igam]; e.dgam = dgam = gridder::elemSize((ulong)igam, ppost.gridGam);
				a_prior.setparameters(params); a_like.setparameters(params);
				e.priorRaw = -a_prior.calcNLL(); e.lik = -a_like.calcNLL(); e.sz = e.drho * e.deta * e.dgam; e.szlog = mylog(e.sz); e.priorAbs = e.priorRaw + e.szlog;
			};

		};

		if (verbose > 0) std::cout << strtime() << "\t\tCalcparampost - Initializing MT thread structures" << std::endl;
		std::vector< t_OpUnit *> vUnits;
		std::vector< std::thread *> vThreads;
		for (int ithread = 0; ithread < NThreads; ithread++) { vUnits.push_back( new t_OpUnit( *this, ppost, PriorWeight, Priors, bs )); vThreads.push_back(nullptr); };
		int thread_cnt = 0;

		// generate unnormalized priors, element sizes and likelihoods
		if (verbose > 0) std::cout << strtime() << "\t\tCalcparampost - Generating likelihoods." << std::endl;
		for (vecDi irho = 0; irho < ppost.gridRho.size(); irho++) {
			if (verbose > 1) std::cout << "\r" << strtime() << "\t\t\tCalcparampost - Generating Likelihood. Irho= " << irho << " of " << ppost.gridRho.size() << std::flush;
			for (vecDi ieta = 0; ieta < ppost.gridEta.size(); ieta++) {
				if (verbose > 2) std::cout << "\r" << strtime() << "\t\t\t\tCalcparampost - Generating Likelihood. Irho= " << irho << " / " << ppost.gridRho.size() << "  Ieta= " << ieta << " / " << ppost.gridEta.size() << std::flush;;
				for (vecDi igam = 0; igam < ppost.gridGam.size(); igam++) {
					if (thread_cnt >= NThreads) {
						// quick and dirty, if we reached a thread limit, wait for threads to complete, in in order...
						for (int ithread=0; ithread <NThreads; ithread++)  { 
							// save results and delete thread
							vThreads[ithread]->join(); ppost.elements.push_back( vUnits[ithread]->e ); delete vThreads[ithread]; vThreads[ithread] = nullptr;
						};
						thread_cnt = 0;
					}
					// create up to thread_cnt threads, each with its own object to work on.... not as efficent as an async thread pool, but much easier to code.
					t_OpUnit &o=*(vUnits[thread_cnt]); o.irho=irho; o.ieta=ieta; o.igam = igam;
					vThreads[thread_cnt] = new std::thread( &(t_OpUnit::CalcOpS), &o );
					thread_cnt++;
				}
			}
		}
		// clean up remaining threads.....
		for (int ithread = 0; ithread <thread_cnt; ithread++) {
			vThreads[ithread]->join(); ppost.elements.push_back(vUnits[ithread]->e); delete vThreads[ithread]; vThreads[ithread] = nullptr;
		};
		thread_cnt = 0;
	};

	// normalize the prior P(r,e,g), and calculate prior density...  then calculate raw posterior
	if (verbose > 0) std::cout << strtime() << "\tRenormalizing Priors " << std::endl;
	vecD tmp_ll; tmp_ll.resize( ppost.elements.size()); for (ulong i = 0; i < ppost.elements.size(); i++) tmp_ll[i] = ppost.elements[i].priorAbs;
	//double norm = butils::mathplus::normLL(tmp_ll); 
	butils::mathplus::normLL(tmp_ll); 
	for (ulong i = 0; i < ppost.elements.size(); i++) {
		auto &p = ppost.elements[i];
		p.priorAbs	 = tmp_ll[i];				// save normalized prior (logspace) 
		p.priorDens	 = p.priorAbs - e.szlog;	// genreate prior density as prior / cellsz (logspace)
		p.posterior  = p.priorAbs + p.lik;		// posterior is porportionate to prior times likelihood (logspace)
	}
	tmp_ll.clear();

	// normalize posterior and calculate posterior density.
	if (verbose > 0) std::cout << strtime() << "\tRenormalizing Posteriors" << std::endl;
	tmp_ll.resize(ppost.elements.size()); for (ulong i = 0; i < ppost.elements.size(); i++) tmp_ll[i] = ppost.elements[i].posterior;
	// norm= butils::mathplus::normLL(tmp_ll); 
	butils::mathplus::normLL(tmp_ll); 
	for (ulong i = 0; i < ppost.elements.size(); i++) {
		auto &p = ppost.elements[i];
		p.posterior		 = tmp_ll[i];				// save normalized posterior, logspace
		p.posteriorDens  = p.posterior - e.szlog;	// generate psoterior density as posterior / cellsz... logspace..
	}
	tmp_ll.clear();

	// Calculate expectations, and marginalize to get single variable distributions....
	// Start by cl;earing accumulators......
	ppost.margRhoPri = ppost.margRhoPos =ppost.gridRho; ppost.margEtaPri = ppost.margEtaPos = ppost.gridEta; ppost.margGamPri = ppost.margGamPos = ppost.gridGam;
	{ for (ulong i = 0; i < ppost.margRhoPri.size(); i++) ppost.margRhoPri[i] = 0; }; ppost.margRhoD = ppost.margRhoPos = ppost.margRhoPri;
	{ for (ulong i = 0; i < ppost.margEtaPri.size(); i++) ppost.margEtaPri[i] = 0; }; ppost.margEtaD = ppost.margEtaPos = ppost.margEtaPri;
	{ for (ulong i = 0; i < ppost.margGamPri.size(); i++) ppost.margGamPri[i] = 0; }; ppost.margGamD = ppost.margGamPos = ppost.margGamPri;

	// TODO clculate log(expected likelihood), rather than expected(log likelihood)...
	// Now loop over cells and sum to generate martinals. 
	if (verbose > 0) std::cout << strtime() << "\tCalculating marginal distributions" << std::endl;
	auto &p=ppost;	ulong ind = 0; // gather marginals..
	p.margRhoLik.clear(); p.margRhoLik.resize(p.gridRho.size());
	p.margEtaLik.clear(); p.margEtaLik.resize(p.gridEta.size());
	p.margGamLik.clear(); p.margGamLik.resize(p.gridGam.size());

// Use LEL, that is Log( expected Likelihood) rather than ELL <expected( log likelihood)>... What we really want is expected likelihood and log(EL) is justa way to represent it.
#define POSTERIOR_LEL
#ifndef POSTERIOR_LEL
	// for ELL set these to size 1, adn jsut accumulate weighted (log likelihoods), otherwise we have to accumulate log(weighted likelihoods)
	for (vecDi i = 0; i< p.margRhoLik.size(); ++i) p.margRhoLik[i].push_back(0);
	for (vecDi i = 0; i< p.margEtaLik.size(); ++i) p.margEtaLik[i].push_back(0);
	for (vecDi i = 0; i< p.margGamLik.size(); ++i) p.margGamLik[i].push_back(0);
#endif
	for (vecDi irho = 0; irho < ppost.gridRho.size(); irho++) {
		for (vecDi ieta = 0; ieta < ppost.gridEta.size(); ieta++) {
			for (vecDi igam = 0; igam < ppost.gridGam.size(); igam++) {
				const postElem &pe = p.elements[ind];
				double post = exp(pe.posterior), prio = exp(pe.priorAbs), llik = pe.lik;
#ifdef POSTERIOR_LEL
				// log (expected likelihood)
				p.margRhoPos[irho] += post; p.margRhoPri[irho] += prio; p.margRhoLik[irho].push_back(llik + pe.szlog); p.margRhoD[irho] += pe.sz;
				p.margEtaPos[ieta] += post; p.margEtaPri[ieta] += prio; p.margEtaLik[ieta].push_back(llik + pe.szlog); p.margEtaD[ieta] += pe.sz;
				p.margGamPos[igam] += post; p.margGamPri[igam] += prio; p.margGamLik[igam].push_back(llik + pe.szlog); p.margGamD[igam] += pe.sz;
#else
				// expected (log likelihood)
				p.margRhoPos[irho] += post; p.margRhoPri[irho] += prio; p.margRhoLik[irho][0] += llik * pe.sz; p.margRhoD[irho] += pe.sz;
				p.margEtaPos[ieta] += post; p.margEtaPri[ieta] += prio; p.margEtaLik[ieta][0] += llik * pe.sz; p.margEtaD[ieta] += pe.sz;
				p.margGamPos[igam] += post; p.margGamPri[igam] += prio; p.margGamLik[igam][0] += llik * pe.sz; p.margGamD[igam] += pe.sz;
#endif
				ind++;
			}
		}
	}

	vecDi i;
#ifdef POSTERIOR_LEL
	// Convert collection of weighted log likelihoods into a single value representig the log of sum of the weighted likelihoods.
	double t;
	for (i = 0; i<p.gridRho.size(); ++i) { t = butils::mathplus::sumOfLL(p.margRhoLik[i]); p.margRhoLik[i].resize(1); p.margRhoLik[i][0] = t; };
	for (i = 0; i<p.gridEta.size(); ++i) { t = butils::mathplus::sumOfLL(p.margEtaLik[i]); p.margEtaLik[i].resize(1); p.margEtaLik[i][0] = t; };
	for (i = 0; i<p.gridGam.size(); ++i) { t = butils::mathplus::sumOfLL(p.margGamLik[i]); p.margGamLik[i].resize(1); p.margGamLik[i][0] = t; };

	// Priors and Posteriors should ALREADY be normalized, now add up the volume of parameter space associaetd with each bin in the marginal
	
	// divide the expected likelihood by the sum of the weights (bu subtracting the log of the weight_sum)
	for (i = 0; i<p.margRhoLik.size(); i++) p.margRhoLik[i][0] -= mylog(p.margRhoD[i]); 
	for (i = 0; i<p.margEtaLik.size(); i++) p.margEtaLik[i][0] -= mylog(p.margEtaD[i]);
	for (i = 0; i<p.margGamLik.size(); i++) p.margGamLik[i][0] -= mylog(p.margGamD[i]);
#else
	// divide by the sum of the weights applied to each value.
	for (i = 0; i<p.margRhoLik.size(); i++) p.margRhoLik[i][0] /= p.margRhoD[i];
	for (i = 0; i<p.margEtaLik.size(); i++) p.margEtaLik[i][0] /= p.margEtaD[i];
	for (i = 0; i<p.margGamLik.size(); i++) p.margGamLik[i][0] /= p.margGamD[i];
#endif
	if ( false ) {
		// test for normalization, note that marginal priors and posteriors are in probspace, not logprobspace
		double t1,t2;
		t1 = t2 = 0.0; for (ulong i = 0; i<p.gridRho.size(); i++) { t1 += p.margRhoPri[i]; t2 += p.margRhoPos[i]; }; fprintf(stdout, "\nRHO: %.16lf \t %.16lf\n", (double)t1, (double)t2);
		t1 = t2 = 0.0; for (ulong i = 0; i<p.gridEta.size(); i++) { t1 += p.margEtaPri[i]; t2 += p.margEtaPos[i]; }; fprintf(stdout, "\nETA: %.16lf \t %.16lf\n", (double)t1, (double)t2);
		t1 = t2 = 0.0; for (ulong i = 0; i<p.gridGam.size(); i++) { t1 += p.margGamPri[i]; t2 += p.margGamPos[i]; }; fprintf(stdout, "\nGAM: %.16lf \t %.16lf\n", (double)t1, (double)t2);
	}

	// calculate expectations...
	if (verbose > 0) std::cout << strtime() << "\tCalculating expectations" << std::endl;
	double sum_v, sum_w, w; p.exp.clear();
	sum_v = sum_w = 0; for (i = 0; i<p.margRhoLik.size(); i++) { w = p.margRhoPos[i]; sum_v += (p.gridRho[i] * w); sum_w += w; }; p.exp.rho = sum_v / sum_w;
	sum_v = sum_w = 0; for (i = 0; i<p.margEtaLik.size(); i++) { w = p.margEtaPos[i]; sum_v += (p.gridEta[i] * w); sum_w += w; }; p.exp.eta = sum_v / sum_w;
	sum_v = sum_w = 0; for (i = 0; i<p.margGamLik.size(); i++) { w = p.margGamPos[i]; sum_v += (p.gridGam[i] * w); sum_w += w; }; p.exp.gam = sum_v / sum_w;

	// calculate credal boundaries,
	// note first and last cells are 1/2 cells with grid point at exterior boundary. Interior cells are "centered" on the 
	//	grid sample point and have extents [ (grid[i-1]+grid[i])/2, (grid[i]+grid[i+1])/2], this is ASSYMETRICAL, thus the complexity
	if (verbose > 0) std::cout << strtime() << "\tCalculating credial intervals" << std::endl;
	auto gridpos = [&] ( const vecD &marg, const vecD &grid, double margPt, double &gridPt ) {
		if ( grid.size() != marg.size() || grid.size() < 2 ) return false;
		vecD pbound( grid.size() ); double prob_prev=1.0, prob_next; pbound[0]=0;
		for (vecDi i = 1; (i+1)<grid.size(); ++i) {	// calculate probability mass below each grid point (thus first value is 0, and last value is 1.0)
			prob_next = (0.5*(grid[i] - grid[i-1])) / (0.5*(grid[i+1] - grid[i]) + 0.5*(grid[i] - grid[i-1]));
			pbound[i] = pbound[i-1] + prob_prev * marg[i-1] + prob_next * marg[i];	// Increment by add trailing fraction of previous sample to leading fractuon of current one.
			prob_prev = 1.0 - prob_next; 
		}; pbound[pbound.size()-1] = 1.0;
		for (vecDi i=1;i<pbound.size();++i){
			if (margPt<pbound[i]) { // linear interpolation...
				gridPt = grid[i-1] + (grid[i]-grid[i-1])*(margPt-pbound[i-1])/(pbound[i]-pbound[i-1]);
				return true; }; 
		}
		return false;
	};
	p.ci = CredalThresh;
	gridpos(p.margRhoPos, p.gridRho, CredalThresh/2, p.expCIm.rho); gridpos(p.margRhoPos, p.gridRho, 0.5, p.expCI.rho); gridpos(p.margRhoPos, p.gridRho, 1.0 - CredalThresh/2, p.expCIp.rho);
	gridpos(p.margEtaPos, p.gridEta, CredalThresh/2, p.expCIm.eta); gridpos(p.margEtaPos, p.gridEta, 0.5, p.expCI.eta); gridpos(p.margEtaPos, p.gridEta, 1.0 - CredalThresh/2, p.expCIp.eta);
	gridpos(p.margGamPos, p.gridGam, CredalThresh/2, p.expCIm.gam); gridpos(p.margGamPos, p.gridGam, 0.5, p.expCI.gam); gridpos(p.margGamPos, p.gridGam, 1.0 - CredalThresh/2, p.expCIp.gam);
	return true;
}


