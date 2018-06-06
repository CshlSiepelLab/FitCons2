#pragma once

#include <cstring>	// strcnp
#include <cstdio>
#include <vector>

#include "butils/butils.h"

#include "Locus.h"


class Block : public Locus {
public:
	typedef butils::fastfile::FFStatus FFStatus;
	float		theta;
	float		lambda;
	void		clear( void ) { theta = lambda = 0.0; }
	// TODO: SLOW, but rarely used, OK....
	FFStatus	writeInp( FILE * f)		{ if ( fprintf(f,"block\t%s:%lu-%lu\ttheta\t%.9e\tlambda\t%.9e\n",chrom.c(),(ulong)st,(ulong)en,(float)theta,(float)lambda) <1) return( butils::fastfile::SERR ); return(butils::fastfile::SEOL) ; }
	FFStatus	readInpS( char *s) { 
		using namespace butils::fastfile;
		ulong t_st=0, t_en=0; double t_theta=0.0, t_lambda = 0.0;
		// replace first : and first - with tabs.... slow but works...
		char *p = s; bool done = false;
		for (done = false; !done && *p != '\0'; p++) { if (*p == ':') { *p = '\t'; done = true; }; };
		for (done = false; !done && *p != '\0'; p++) { if (*p == '-') { *p = '\t'; done = true; }; };
		// TODO: SLOW but Ok for now, rarely used and when used, is used only on blocks which are sparse...
		int n = sscanf(s, "block %s %lu %lu theta %le lambda %le", chrom.c(), &t_st, &t_en, &t_theta, &t_lambda);
		if (n != 5) return(SERR);
		st = (Coord) t_st; en = (Coord) t_en; theta = (float) t_theta; lambda = (float) t_lambda;
		return(SEOL);
	}
	Block() { clear(); }
	virtual ~Block(void) {};

	FFStatus read( FILE *f ) {
		using namespace butils;
		clear();
		FFStatus s = Locus::read(f);
		if (s== fastfile::SEOF) return( s );
		if (s != fastfile::SEOR) return(fastfile::SERR); theta = lambda = 0;
		if (fscanf(f," %f",&theta)  !=1) return( fastfile::SERR );
		if (fscanf(f," %f",&lambda) !=1) return( fastfile::SERR );
		return( fastfile::termChar( fastfile::getChar( f ) )); 
	}

	FFStatus  readln( FILE *f ) {
		using namespace butils;
		FFStatus s=read(f);
		if (s==fastfile::SEOF) return(fastfile::SEOF);
		if (s==fastfile::SEOR) s= fastfile::clearLine(f);
		return( s );
	}

	// low performance conversion for debugging..
	std::string toString( void ) {
		return(Locus::toString() + " " + std::to_string(theta) +" " + std::to_string(lambda));
	}

};

typedef  std::vector<Block> vBlock;
