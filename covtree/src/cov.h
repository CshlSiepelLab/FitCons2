#pragma once
#include <string>
#include <vector>
#include <sstream>	// low performance convienence IO
#include <iostream>
#include <fstream>
#include <tuple>

#include "stringer.h"


// Covariate values are enumerated from 1, with 0 being a missing value 
// Each Covariate has value has an positive integer ID number (generally <32 or so)
//	Together covariate values must compactly and uniquely fill the range (no holes)
//	Each value has 
//		ID a 1 based index....
//		Tag (generally a 0 padded representation of the value in 1-2 chars), All tags for a single covarate
//			should have the same width.... ie "1"-"5" or "01"-"22", NOT 1-22
//		Nmemonic (<8 chars, no spaces) and a desription (Smallish <20-200 chars is nice, spaces OK)
//
class covVal
{
public:
	typedef stringer string;
private:
	uint8_t		ivalue;
	string		stag;
	string		snemonic;
	string		sdescription;
public:
	void clear() { ivalue = 0; stag = snemonic = sdescription = ""; return; }
	covVal() { clear(); }
	// covVal( const std::string &Init) { parse(Init);	}
	~covVal() { clear(); };
	inline uint8_t value() const { return ivalue; };
	inline const string &tag() const { return stag; };
	inline const string &nemonic() const { return snemonic; };
	inline const string &desc() const { return sdescription; };
	// OK, err
	bool parse(const stringer &Init, int wtag = 1, int wnem = 8) {
		stringer::vecS toks;
		Init.tokenizeQuoted( toks );
		if (toks.size()<2) return false;
		stag=toks[0];											// Tag as string
		ivalue = (uint8_t)std::stoi(stag);						// Tag as number
		snemonic = toks[1];										// Get Nmeonic
		sdescription = (toks.size()<3 ? snemonic : toks[2]);	// get long description (quoted)
		if (stag.length() != (uint32_t) wtag) return false;				// Check values....
		if (snemonic.length() > (uint32_t) wnem) return false;
		return true;
	};
	static void skipSP(std::istream &fin) {
		while ((fin.peek() == ' ') || (fin.peek() == '\t') || (fin.peek() == '\r')) fin.get();
	}
	static bool getTok(std::istream &fin, std::string &tok) {
		tok = "";
		if (fin.peek() == (int) '\n') return true;
		skipSP(fin);						// skip leading spaces or tabs, if any
		if (fin.peek() == (int) '\n') return true;
		if (fin.peek() == (int) '"') {		// is token quoted?
			std::getline(std::getline(fin, tok, '"'), tok, '"');
		} else {
			fin >> tok;
		}
		return true;
	}
	static bool getLine(std::istream &fin, std::string &line) {
		line = ""; std::getline(fin, line);
		while (line.substr(0, 1) == "#") std::getline(fin, line);
		return true;
	}
};

// A single covariate, with a set of values, potentially some subset of the universe of all possible values.
class covDef
{
public:
	typedef stringer string;
	typedef std::vector<covVal> vecVals_t;
private:
	string nemonic;
	string description;
	string dirbase;
	uint8_t tagw;				// covaraite value tag width
	bool bcts;					// Cell Type Sensative? alternative is "annotation" ( 1 set for all cell types)
	bool bmonotonic;			// Monotonic? Increasing class has uniform impact on conservation (ie RnaSeq). Alternative is ChromHMM or MeltingTemp.
	vecVals_t	vals;			// collection of values present in this subset, ordered, sparse.
public:
	void clear() { nemonic = description = dirbase = ""; bcts = bmonotonic = false; vals.clear(); tagw = 1; };
	inline bool			isCTS()		const { return( bcts ); };		// is Cell Type Sensitive - that is, is there a covariate value for each cell type, or do all cells share 1 value per position
	inline bool			isMon()		const { return(bmonotonic); };	// is momontonc, that is is selective pressure supposed to be monotonic (up or down) function of this  variable  or is it contect dependent?
	inline string		tag()		const { return nemonic; };
	inline uint8_t		tagW()		const { return tagw; };
	const  vecVals_t &	val()		const { return vals ;};
	inline string		Name()		const { return nemonic ;};
	inline uint8_t		NumVals()	const { return (uint8_t) vals.size(); };
	covDef() {};
	~covDef() {};
	const covVal *valByTag(const string &Tag)  const {
		for ( auto & v : vals) { // Slow linear search, but done infrequently.....
			if (v.tag()==Tag) return( & v); };	
		return( NULL );
	}
	bool read(std::istream &fin, int wnem = 8) {
		string stmp; stmp=""; string::vecS toks;
		clear();
		if (!fin.good()) return false;
		stmp=""; std::getline( fin, stmp ); if (stmp.length()<1) return false;
		stmp.tokenizeQuoted( toks ); if (toks.size()<6) return false;
		/* 1 */ nemonic=toks[0];	if (nemonic.length() > (uint32_t) wnem) return(false);
		/* 2 */ description=toks[1];
		/* 3 */ bcts = (toks[2].aslower()== "cts");
		/* 4 */ bmonotonic = (toks[3].aslower()== "mono");
		/* 5 */ tagw = std::stoi(toks[4]); if (tagw<1) return false;
		/* 6 */ dirbase=toks[5];
		{ char ctmp[1024]; ctmp[0] = 0; sprintf(ctmp, "%0*d NULL MissingValue", tagw, 0); stmp = ctmp; }
		covVal vtmp; 
		while (stmp.length() > 0) {
			if (!vtmp.parse(stmp, tagw, wnem)) return false; 
			vals.push_back(vtmp);
			covVal::getLine(fin, stmp);stmp.trim();
		}
		if (vals.size()<2) return false;					// need at least 2 values to define a covaraite
		for (uint8_t i = 0; i<(uint8_t)vals.size(); i++) {		// index of value must equal tag.
			if (vals[i].value() != i) return false;
		}
		return true;
	};
};


class covsetDef {
public:
	typedef stringer string;
	typedef std::vector< covDef > vecCovDef_t;
private:
	string dirann, dirctsbase, name;
	vecCovDef_t	covAnno;
	vecCovDef_t	covCTS;
public:
	void clear() { dirann = dirctsbase = name = ""; covAnno.clear(); covCTS.clear(); };
	covsetDef() { clear(); };
	const vecCovDef_t &covsAnn( ) const		{ return covAnno;};
	const vecCovDef_t &covsCTS( ) const		{ return covCTS; };
	string dirAnn() const					{ return dirann; };
	string dirCTSb() const					{ return dirctsbase; };
	string Name() const						{ return name; };
	bool read(const string &fname )			{
		clear();
		string stmp, tag; std::vector<string> toks;

		// get covset tag.
		std::ifstream fin( fname );
		covVal::getLine( fin, stmp ); while (!fin.bad() && stmp.length()<1) covVal::getLine(fin, stmp);
		tag= stmp.aslower(); if (tag!= "covset") return false;
		//std::cerr << "Found Covset\n";//debugging

		// Get Covset definition line
		covVal::getLine(fin, stmp); while (!fin.bad() && stmp.length()<1) covVal::getLine(fin, stmp);
		stmp.tokenizeQuoted(toks," \t");
		if (toks.size()<3) return false;
		name=toks[0].trim(); dirann=toks[1].trim(); dirctsbase = toks[2].trim();
		//std::cerr << "Found Covset Def Line:\n\t" << stmp << "\n";// debugging

		// get individual covaraites....
		covDef cov_tmp;
		// get a line
		do { covVal::getLine(fin, stmp); } while (stmp.length()<1);
		/// while there is more data, process line and get next one
		do {
			// process line
			// std::cerr << "Read Cov Line:\n\t" << stmp << "\n"; // debugging
			if ( stmp.trim().aslower() == "end")		break;	// this is the normal end of file
			if ( stmp.trim().aslower() != "covariate")	return false;	// anything else but hte start of a covariate block is an error
			if (!cov_tmp.read( fin )) return false;
			// save processed data
			(cov_tmp.isCTS() ? covCTS : covAnno ).push_back( cov_tmp );
			// read next line
			do { covVal::getLine(fin, stmp); } while (fin.good() && stmp.length()<1);
		} while (stmp.length()>0);

		return true;
	};
};

// Some subset of the universe of all possible values.
class covDefSubset {
public:
	typedef stringer string;

	std::vector<const covVal *>	vals;		// collection of values present in this subset, ordered, sparse.
	const covDef		*univ;				// Class that defines all possible values. NULL if this instance is the source.
	const covDefSubset	*parent;			// Superset of possible values, can be source. Null if this is the source.
	stringer Name()				const		{ return univ->Name(); };
	uint32_t numVals(bool Univ=false) const	{ return((uint32_t)(Univ?univ->NumVals() : vals.size())); };								// number of values in this subset..
	uint32_t numValPartitions() const		{ auto sz = numVals(); return( sz > 1 ? sz - 1 : 0); };			// number of nonempty, nonuniversal, linear partitions.

	void clear( bool ClearUniv = false )	{ if (ClearUniv) univ=NULL; parent=NULL; vals.clear(); return; };
	covDefSubset() { clear(true); };
	~covDefSubset() { clear(true); };

	inline bool isCTS() const {		return( univ && univ->isCTS() );	};
	inline bool hasNULL() const {	for (const auto &s : vals ) {if (s->value()==0) return true; }; return false; };	// Check to see contains a NULL (missing data) element

	// Useful interfaces to Init feature...
	bool Init(const covDef &U)				{	return Init( &U );				};	// make a copy from the Unviersal / compelte set. Base Case, self parrenting.
	bool Init(const covDefSubset &P, const string &Template = "")					// Take a subset of the Parent, most common form. Default is copy.
											{ return Init(NULL,&P,Template); };	

	// Most generic form....
	bool Init( const covDef *U = NULL, const covDefSubset *P= NULL, const string &Template ="" ) {	// General Form...
		if ((U==NULL) && ( P==NULL)) return false;
		if (U==NULL) U = P->univ; if (U==NULL) return false;
		if (P==NULL) P = this;
		clear( true );
		univ = U; parent = P;
		stringer::vecS toks;  stringer tmp=Template; tmp.trim();
		if (tmp.substr(0, 1) == "{") {
			if (tmp.length()<3) return false;
			if (tmp.substr(tmp.length()-1,1)!="}") return false;
			tmp = tmp.substr(1, tmp.size() - 2);			// delete brackets
		}
		if (tmp.find_first_of(";") != string::npos) {
			tmp.tokenize(toks, ";"); if (toks.size() > 2) return false;
			if (toks[0] != univ->tag() ) return false;		// check tag name, is this the right covariate?
			tmp = toks[1];									// Trim name from values
		}
		if ((tmp=="") || (tmp=="*")) { 
			if (parent == this ) {
				for ( auto & s : univ->val()) vals.push_back( &s ); 
			} else { 
				vals = parent->vals;
			}
		} else {
			tmp.tokenize(toks, ","); if (toks.size()<1) return false;
			std::sort(toks.begin(),toks.end());
			for (auto & s : toks) {
				auto t = univ->valByTag( s );
				if (t==NULL) return false;
				vals.push_back( t );
			}
		}
		return true;
	}

	// Convert list of covaite values into a comma seperated list of value tags (not Nemonics!)
	string tagList() const {
		string sout; for ( auto & s : vals ) { sout += (sout.size()>0?",":"") + s->tag(); }
		return sout;
	}

	string toTemplate() const {
		return( string("{") + Name() + ";" + tagList() + "}" );
	};

	// Return a vector of CovariateValue sets, each with a single value. Generally used to identify conditional ordering of covaraite values for non-mono covariates.
	void AllValueSet( std::vector<covDefSubset> &AllVals ) const { 
		uint32_t index = 0;
		for (auto v : vals) {
			// if (v->value != 0 ) {	// Ignore the NULL class... hardwired .. sigh...
				AllVals.resize(index+1);
				AllVals[index].Init( *this, v->tag() );	
				index++;
			// }
		}
		return;
	}

	// Generate all monotonic partitions of ordered vales. For N values this produces N-1 {Set, coset} pairs. Pairs with emptu Set or Coset are ignored.
	//	of PORDER is null, then assume valeus are already (gnerally natively) ordered....
	// This is inefficent, but run rarely.
	bool AllPartitionSet( const std::vector<double> *POrder, std::vector<covDefSubset> &ValsSet , std::vector<covDefSubset> &ValsCoSet, int v=0) const {
		if (vals.size()<2) return false;								// can't split a singelton
		if (POrder && (POrder->size() != vals.size())) return false;	// order must be same size as vals, if POrder is not provided asume set is already ordered.

		// Make a copy and remove NULL value if it is part of this covariate. Missing data positions (NULL) go into NO subclass...
		std::vector<const covVal *>	local_vals = vals;
		std::vector<double> local_order;
		if (POrder!=NULL) { 
			local_order = *POrder; 
		} else {
			local_order.resize( vals.size());
			for (uint32_t i=0;i<vals.size(); i++) local_order[i]=i;	// generate monotonic ordering...
		}
		if (local_vals[0]->value() == 0) {
			local_vals.erase(local_vals.begin() + 0); local_order.erase(local_order.begin() + 0); };
		if ( local_vals.size() <2 ) return false;						// again, singlets can't be split... why waste time...

		if (v > 0) 	std::cerr << "\nAllPartition Set:  " << this->Name() << std::endl;

		// Generae tuples for sorting by Order Value
		typedef std::tuple<double, const covVal *> ord_t;
		std::vector<ord_t> tmp_vals;
		if (v>1) std::cerr << "\t Presort:\n";
		for ( uint32_t i=0; i<local_vals.size(); i++) {
			if (v>1) std::cerr << "\t\t" << i << "\t" << local_order[i] << "\t" << local_vals[i]->tag() << std::endl;
			tmp_vals.push_back( ord_t(local_order[i],local_vals[i]) ); }
		std::sort( tmp_vals.begin(), tmp_vals.end(), [](const ord_t & a, const ord_t & b) { return( std::get<0>(a) < std::get<0>(b) ) ;});
		if (v>2) {
			std::cerr << "\t Postsort:\n";
			for (uint32_t i = 0; i<tmp_vals.size(); i++) {
				std::cerr << "\t\t" << i << "\t" << std::get<0>(tmp_vals[i]) << "\t" << std::get<1>(tmp_vals[i])->tag() << std::endl;
			}
		}
		// Allocate the Set and Coset elements of each patition
		uint32_t num_splits = (uint32_t) tmp_vals.size();
		if (num_splits<1) return false;
		num_splits--;
		ValsSet.resize(num_splits); ValsCoSet.resize(num_splits);

		if (v>1) std::cerr << "\nPostSort - subsets:  " << this->Name() << std::endl;
		// Generate Partitions on ordered covariate values
		for (uint32_t i_splitafter = 0; i_splitafter < num_splits; i_splitafter++ ) {
			string tags_set, tags_coset;
			if (v>1)  std::cerr << "\t\t" << i_splitafter << "\t";
			for (uint32_t j = 0;				j <= i_splitafter;		j++) tags_set   += (tags_set.length()  >0 ? "," : "") + std::get<1>(tmp_vals[j])->tag();
			for (uint32_t j = i_splitafter+1;	j <  tmp_vals.size();	j++) tags_coset += (tags_coset.length()>0 ? "," : "") + std::get<1>(tmp_vals[j])->tag();
			if (v>1) std::cerr << tags_set << "\t|\t" << tags_coset << std::endl;
			ValsSet[  i_splitafter ].Init( *this, tags_set );	// Init will sort the tags into lexi order...
			ValsCoSet[i_splitafter ].Init( *this, tags_coset);
		};
		if (v>0) std::cerr << "\n";
		return true;
	}

};

// .. a full cov string might look like this
// Annotation:{CovNam;v1,v2,v3,v4},{CovNam;v1,v2,v3,v4},...:CTS:{...
class covsetDefSubset {
	typedef stringer string;
	typedef std::vector<covDefSubset>			covDefSubset_v;
	covDefSubset_v			covAnnotations;		// collection of values present in this subset, ordered, sparse.
	covDefSubset_v			covCTS;				// collection of values present in this subset, ordered, sparse.
	const covsetDef			*univ;				// Class that defines all possible values. 
	const covsetDefSubset	*parent;			// Superset of possible values, can be source. *this if this is the root.
public:
	covsetDefSubset()	{ univ = NULL; parent = NULL; clear(true); };
	~covsetDefSubset()	{ clear(true); };
	void clear( bool ClearUniv = false ) { 
		covCTS.clear(); covAnnotations.clear();
		if (ClearUniv) univ = NULL; 
		parent = NULL;  return; 
	};

	bool CanBePartitioned() const {		// is is possible to subtartition this covsetSubset. It is possible of SOME cov in the set has more than one value.
		for (const auto & s : covAnnotations ) if (s.numVals()>1) return true;
		for (const auto & s : covCTS         ) if (s.numVals()>1) return true;
		return false;
	}

	// Use these
	bool Init(const covsetDef &U ) {									// Create Initial subset from universe. Contains all values
		return( Init( &U, NULL , "") ); };								// Subset another subset, really, template should never == ""
	bool Init(const covsetDefSubset &P, const string &Template = "") { 
		return( Init( NULL, &P, Template ));	};
	bool InitAndReplaceCov(const covsetDefSubset &P, const covDefSubset &Repl ) {	// Copy parent, but replace appropriate covaraite with Repl
		return( Init( NULL, &P, "", &Repl )); 
	};

	// Most general routine..... dont use this outside class (make private?)
	bool Init(const covsetDef *U = NULL, const covsetDefSubset *P = NULL, const string &Template = "", const covDefSubset *Repl=NULL ) {
		if ((U==NULL) && (P==NULL)) return false;
		if (U==NULL) U=P->univ;
		if (P==NULL) P = this;
		clear( true );
		univ = U; parent = P;
		if (Template == "") { // Copy all states from Parent or Univ
			if (parent == this) {	// Copy from univ... 
				for (const auto & s : univ->covsAnn()) { covDefSubset sub; sub.Init(s); covAnnotations.push_back(sub); }
				for (const auto & s : univ->covsCTS()) { covDefSubset sub; sub.Init(s); covCTS.push_back(sub); }
			} else {				// Just copy parent values, possibly replacign one covariate def....
				bool do_replace;
				for (const auto & s : parent->covAnnotations)	{ 
					do_replace = Repl && !Repl->univ->isCTS() && (Repl->Name() == s.Name());
					covAnnotations.push_back((do_replace?*Repl:s)); }
				for (const auto & s : parent->covCTS)			{ 
					do_replace = Repl && Repl->univ->isCTS() && (Repl->Name() == s.Name());
					covCTS.push_back((do_replace ? *Repl : s)); }
			}
		} else {
			// Copy selected set from parent... requires valid parent.....
			if (P==this) return false;
			string::vecS toks, toks1, toks3; 
			Template.tokenize(toks,":");
			if (toks.size()!=4) return false;

			// Handle annotations
			if (toks[0].asupper().substr(0,3)!="ANN") return false;
			toks[1].tokenizeQuoted( toks1,",","{}");	// extract individual covaraites
			if (toks1.size()!=parent->covAnnotations.size()) return false;
			uint32_t index=0;
			for (index = 0; index < parent->covAnnotations.size(); index++) {
				covDefSubset sub;
				if (!sub.Init(parent->covAnnotations[index], toks1[index])) return false;
				covAnnotations.push_back( sub );
			}

			// Handel Cell Type Specific Data..
			if (toks[2].asupper().substr(0, 3) != "CTS") return false;
			toks[3].tokenizeQuoted(toks3, ",", "{}");	// extract individual covaraites
			for (index = 0; index < parent->covCTS.size(); index++) {
				covDefSubset sub;
				if (!sub.Init(parent->covCTS[index], toks3[index])) return false;
				covCTS.push_back(sub);
			}
		}
		return true;
	};

	const uint8_t  getNumCovIDS(bool IsCTS) const {
		return (uint8_t) (IsCTS? covCTS.size() : covAnnotations.size() );
	}

	const covDefSubset *getCovByID(bool IsCTS, uint8_t ID ) const {
		const covDefSubset_v &covs = (IsCTS? covCTS : covAnnotations );
		if (ID>=covs.size()) return NULL;
		return &covs[ID];
	}

	const covDefSubset *getCovByName(const string &CovName) const {
		auto s = getCovByName( CovName, true );
		if (s==NULL) s = getCovByName(CovName, false);
		return s;
	}

	const covDefSubset *getCovByName(const string &CovName, bool isCTS ) const {
		for (const auto & s : (isCTS? covCTS : covAnnotations))
			if (s.Name() == CovName ) return( &s );
		return NULL;
	}

	const string toTemplate() const {
		string out="", tmp="";
		for ( const auto & s : covAnnotations ) { tmp += (tmp.size()>0?",":"") + s.toTemplate(); };
		out += "ANN:" + tmp; tmp = "";
		for (const auto & s : covCTS) { tmp += (tmp.size()>0 ? "," : "") + s.toTemplate(); };
		out += ":CTS:" + tmp; tmp = "";
		return( out );
	}
	// BinaryPartitions( Covariate, Ordering )	// Get all monotinic splits of chosen covaraite based ordering, alogn with current lask for all otehr covariates
};

