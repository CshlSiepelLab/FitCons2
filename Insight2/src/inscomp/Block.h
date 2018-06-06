#pragma once

#include <vector>

#include "butils/butils.h"

// This is just a STORAGE class with convienet access functions. 
//	use the BlockCompiler to create the class


// to save memory convert Tcount from uint32, down to uint16, less than uint16 is likely a bad idea....
// when running with 1 cell type uint16 is safe.... with multiple it we might conceivablly hit 100K (unlikely),
//	just over the 64K limit. Test if this is a concern....


template< typename Tval, typename Tind, typename Tcount>
class Block 
{
public:
	typedef butils::ulong	ulong;
	typedef insdb::Coord	Coord;
	typedef butils::StringWrapper<char, 8>	stringT;
	static butils::RTable<Tval, Tind>	ancpri;	// we don't really need this here, we only have one table for now, make this static..
protected:
	// Major and Minor allele ancenstral priors, two indices into an Rtable of real probability values [0.0-1.0]
	struct TPolyA		{ Tind		maj; Tind min; 
		bool operator< ( const TPolyA &R) const  { return(maj != R.maj ? (maj < R.maj) : (min < R.min)); }; };	
	struct TMonoCount	{ Tind		ind; Tcount	cnt; };	// Allele Ancestral Prior for mono sites, and the number of tiems it was observed in block
	struct TPolyCount	{ TPolyA	ind; Tcount cnt; };	// Allele Ancestral Prior pair, and the number of tiems it was observed in block

private:
	stringT	chrom;		// Name of chromosome, 8 bytes
	Coord	st, en;		// 4 (or 8) bytes each
	double	lambda;		// blocks are rare, so strore full values here, 8 bytes each
	double	theta;
	// we get chrom, start and end positions from base class, locus...

protected:
	// if need be, these could be replaced with arrays.... dropping overhead by about 100 bytes / block.
	std::vector<TMonoCount> mono;		// mono and poly are very small geenrally 1 or 2 ancestral prioris, each is a 1 uchar/ushort index into a table of reals (Rtable)
	std::vector<TPolyCount> polyL;		// 
	std::vector<TPolyCount> polyH;		//

	inline void polyVal( const std::vector<TPolyCount> &p, ulong I, Tind &IAncPriMaj, Tind &IAncPriMin, Tcount &Nobs ) const { 
		IAncPriMaj=p[I].ind.maj; IAncPriMin=p[I].ind.min;Nobs=p[I].cnt;
		};

public:

	Block(void) { clear(); }
	~Block(void) { clear(); };								// keep this non-virtual to avoid vtable allocation

	inline void clear() { lambda = theta = -1.0; st = en = 0; chrom[0] = '\0'; mono.clear(); polyL.clear(); polyH.clear(); }

	inline void set( const char *Ch, Coord St, Coord En, double Lam, double Thet ) {
		en=En; st=St; lambda=Lam; theta=Thet; strncpy(chrom.c(),Ch,(size_t) chrom.siz()); }

	// pat repreent the number of unique patterns in the block, each pattern may be observed ar several sties
	inline ulong monoNumPat( void )		const { return( (ulong) mono.size() ); };		// 64 -> 32 bits
	inline ulong polyHNumPat(void)		const { return((ulong)polyH.size()); };
	inline ulong polyLNumPat(void)		const { return((ulong)polyL.size()); };
	inline ulong polyNumPat(void)		const { return(polyHNumPat() + polyLNumPat()); };
	// number of sites, that is the sum of the counts of each observed pattern..
	inline ulong monoNumSites(void)		const { ulong c = 0; for (ulong i = 0; i < monoNumPat(); i++) c += mono[i].cnt;  return c; };		// 64 -> 32 bits
	inline ulong polyHNumSites( void )	const { ulong c = 0; for (ulong i = 0; i < polyHNumPat(); i++)  c += polyH[i].cnt; return c; };
	inline ulong polyLNumSites( void )	const { ulong c = 0; for (ulong i = 0; i < polyLNumPat(); i++)  c += polyL[i].cnt; return c; };
	inline ulong polyNumSites(void)		const { return(polyHNumSites() + polyLNumSites()); };
	inline ulong NumSites(void)			const { return(polyNumSites() + monoNumSites()); };
	inline bool  hasPolys(void)			const { if (polyLNumPat()>0) return true; if (polyHNumPat()>0) return true; return false; };
	inline bool  hasSiteData(void)		const { return(monoNumPat() > 0 || hasPolys() ); };

	// get indices into Rtable of ancestral priors, as werll as 
	inline void monoVal( ulong I, Tind &IAncPri, Tcount &Nobs ) const { IAncPri=mono[I].ind; Nobs=mono[I].cnt;  };
	inline void polyLVal( ulong I, Tind &IAncPriMaj, Tind &IAncPriMin, Tcount &Nobs ) const { polyVal(polyL,  I, IAncPriMaj,IAncPriMin, Nobs); }
	inline void polyHVal( ulong I, Tind &IAncPriMaj, Tind &IAncPriMin, Tcount &Nobs ) const { polyVal(polyH,  I, IAncPriMaj,IAncPriMin, Nobs); }
	inline void getHeader(stringT *Chrom, Coord *St, Coord *En = NULL, double *Lambda = NULL, double *Theta = NULL) const {
		if (Chrom)	*Chrom = chrom;
		if (St)		*St= st;
		if (En)		*En = en;
		if (Lambda) *Lambda = lambda;
		if (Theta)	*Theta = theta;
	}

};

// sigh static template inits are not so terribly easy.
//   TEMPLATE                                              TYPE                        NAME
template< typename Tval, typename Tind, typename Tcount> butils::RTable<Tval, Tind>  Block<Tval,Tind,Tcount>::ancpri;

