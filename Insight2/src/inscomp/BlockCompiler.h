#pragma once

#include <vector>
#include <map>

#include "Block.h"
//#include <utility>

// this class allows a user to dynamically load a block, one element at a time, then
//	compress the important elements into a compact Block object that is memory efficent

// Tval   - working mathematicval type, typically double, could be float.
// Tind   - unsigned type typcially uint8 or uint16, number of unuique Tvals used
// Tcount - type to hold max count of observed Tvals, typically this is a uint16 or uint32

template< class Tval, class Tind, class Tcount>
class BlockCompiler : public Block<Tval,Tind,Tcount>
{
	typedef Block<Tval, Tind, Tcount>	Block_t;
	typedef  typename Block_t::TPolyA TPolyA;
private:
	//Rtable<Tval,Tind> &ancPri;	// reference to rtable, an indexed table that turns a short int index (1 or 2 bytes) into a long float value (4 or 8 bytes)
									// this doesn't change and coudl really be static...

	TPolyA	tmpPoly;		// temporary value, used when loading poly elements into a map...

	// we accumulate counts of each ancestral prior, for mono, and poly sites. Mono sites have 1 anc prior, Poly have 2 (Maj & Min).
	//	for bervity we only store index into table of prior probs, because priors are 8 bites each and indices are 1 or 2 bytes.
	//	we also compact using RLL, storing an idnex (or index pair) and the number of times we encouter it...

	std::vector<Tcount>		tMono;	// mono have a few (256/65536) possible vlaues, just allcoate an array and increment index...
	std::map<TPolyA,Tcount>	tPolyL;	// the space of pairs of values can be much bigger. Array is OK if Tind is uint8, but NOT if Tind is uint16.
	std::map<TPolyA,Tcount>	tPolyH;
	//
	inline void addPoly(std::map<TPolyA, Tcount> &p, TPolyA &APriPr, uint32_t counts = 1) { 	auto it = p.find(APriPr); if (it != p.end()) { it->second += counts; } else { p[APriPr] = counts; }; };
public:
	BlockCompiler( ) { tMono.resize( std::numeric_limits<Tind>::max() ); clear(); }
	~BlockCompiler() {};

	void clear( bool ClearBase=true ) { for (uint32_t i=0;i<tMono.size();i++) { tMono[i]=0; }; tPolyL.clear(); tPolyH.clear(); if (ClearBase) Block_t::clear(); }

	void addMono( Tind IAncPri, uint32_t counts = 1) { tMono[IAncPri]+=counts; };							// increase count of this observed prior
	//void addMono( Tval  AncPri ) { Tind i=ancPri.map(AncPri); tMono[i]++; };	// increase count of this observed prior

	inline void addPolyX(char F, Tind IAncMaj, Tind IAncMin, uint32_t counts = 1)	{ tmpPoly.maj = IAncMaj; tmpPoly.min = IAncMin; addPolyX(F, tmpPoly, counts); };	// from 2 indices
	inline void addPolyX(char F, TPolyA &IAnc, uint32_t counts = 1)					{ addPoly((F=='H'?tPolyH:tPolyL), IAnc, counts); };												// fron am index pair struct

	void compact( void ) {
		// coleases slow and verbose map formats into a list (vector/array) of value/observation pairs...
		//	compact the verbose wrapepr classes into the base class and return a reference to the compact base 
		//	the base can then be coppied to an element in a vector... (BlockList)

		uint32_t cnt=0, i=0, pos=0;
		// how many entries in mono table
		for (cnt=i=0; i<tMono.size();i++) { if (tMono[i]>0)cnt++; }
		Block_t::mono.resize(cnt);
		// load entries into mono table (val, cnt)
		for (pos=0, i=0; i<tMono.size();i++) { if (tMono[i]>0) {Block_t::mono[pos].ind=i;Block_t::mono[pos].cnt=tMono[i];pos++;} }

		// now load values into poly vectors... again compact format...
		Block_t::polyL.resize( tPolyL.size() ); cnt=0;
		for (auto it=tPolyL.begin(); it !=tPolyL.end(); it++) { Block_t::polyL[cnt].ind=it->first;Block_t::polyL[cnt].cnt=it->second;cnt++; }
		Block_t::polyH.resize( tPolyH.size() ); cnt=0;
		for (auto it=tPolyH.begin(); it !=tPolyH.end(); it++) { Block_t::polyH[cnt].ind=it->first;Block_t::polyH[cnt].cnt=it->second;cnt++; }
	}

	void expand( void ) {
		clear( false );
		// expands compacted form into map form for easy accumulation
		ulong i;
		for (i=0;i<Block_t::mono.size();  i++)	{ tMono[Block_t::mono[i].ind]=Block_t::mono[i].cnt; };
		for (i=0;i<Block_t::polyL.size(); i++)	{ tPolyL[Block_t::polyL[i].ind]=Block_t::polyL[i].cnt; };
		for (i=0;i<Block_t::polyH.size(); i++)	{ tPolyH[Block_t::polyH[i].ind]=Block_t::polyH[i].cnt; };
	}

	void expand(Block_t &aBlock) { (*this) = aBlock; expand(); }

};
