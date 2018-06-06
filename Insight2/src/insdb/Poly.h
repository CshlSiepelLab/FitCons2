#pragma once

#include "util.h"

#include "butils/butils.h"
#include "Locus.h"
#include "Mono.h"


class Poly : public Locus {
public:
	typedef butils::stringS stringS;					// 16 chars
	typedef butils::fastfile::FFStatus FFStatus;		// MSVC diesnt like using in class scope....
	typedef Mono::dbElem dbElem;
	typedef butils::RTable<stringS, dbElem>	RTableS;	// high performance map from index to value, one for string one for double...
	typedef butils::RTable<double, dbElem>	RTableD;
private:
	static RTableS mapString;
	static RTableD mapDouble;
public:
	dbElem	aPriMaj;
	dbElem	aPriMin;
	char	freq;		// Valid valeus are N (uninit), L (Low Freq), H (Hi Freq)
	inline double		majDouble() const { return(*(Poly::mapDouble[aPriMaj])); };
	inline double		minDouble() const { return(*(Poly::mapDouble[aPriMin])); };
	inline const char *	majString() const { return(& ((*(Poly::mapString[aPriMaj]))[0])); };
	inline const char *	minString() const { return(& ((*(Poly::mapString[aPriMin]))[0])); };

	Poly() { freq = 'N'; aPriMaj = aPriMin = (dbElem)0; }
	virtual ~Poly(void) { };

	static void setTables( void *tString, void *tDouble) {
		mapString = *((RTableS *) tString); mapDouble = *((RTableD *) tDouble);
	}

	FFStatus  read( FILE *f ) {
		using namespace butils;
		FFStatus s = Locus::read(f);
		ulong i=0;
		if (s==fastfile::SEOF) return( s );
		if (s!= fastfile::SEOR) return(fastfile::SERR );
		freq=(char)fastfile::getChar( f );
		if ((freq!='L') && (freq!='H')) return(fastfile::SERR );
		if (fastfile::termChar(fastfile::getChar( f ) ) != fastfile::SEOR ) return(fastfile::SERR );
		i = 0;  if (fastfile::getInt(f, &i) != fastfile::SEOR) return(fastfile::SERR);
		aPriMaj = (dbElem)i;
		i = 0;  s = fastfile::getInt(f, &i);
		aPriMin = (dbElem)i; 
		if (s!= fastfile::SEOL) s= fastfile::SERR;
		return( s );
	};
	FFStatus readln( FILE *f ){
		using namespace butils;
		FFStatus s=read(f);
		if (s==fastfile::SEOF) return( s );
		if (s==fastfile::SEOR) s=fastfile::clearLine(f);
		return( s );
	}
	// low performance conversion for debugging..
	std::string toString( void ) {
		std::string s = Locus::toString();
		return(  s + " " + majString( ) + " " + minString( ) );
	}

	bool writeInp( FILE *f ) const { 
		bool res = false;
		if (!f) return(false); 
		res = (fprintf(f,"site\t%s:%lu\t%c\t%s\t%s\n",chrom.c(),(ulong)st,freq,majString(),minString())>0);
		return(res ); }

	// bool set() {};
};

// definition of static variables...
Poly::RTableS Poly::mapString;
Poly::RTableD Poly::mapDouble;
