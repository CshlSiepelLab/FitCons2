#pragma once

#include <vector>

#include "butils/butils.h"

#include "Block.h"
#include "BlockCompiler.h"

template< class Tval, class Tind, class Tcount>
class BlockSet
{
public:
	// types
	typedef Block<Tval, Tind, Tcount>		Block_t;
	typedef Tind							AncPriInd_t;
	typedef Tval							AncPriVal_t;
	typedef Tcount							AncPriCount_t;
	// Locals
	butils::RTable<Tval, Tind>				ancPri;			// small set (<255 or < 65K) of unique values for of ancestral priors spanning (0.0,1.0).
	std::vector<Block_t>					blocks;			// memory efficent block structures
	BlockCompiler<Tval, Tind, Tcount>		blockMaker;		// temp object for reading values and compressing them into calc-friendly blocks

	// global values for blockset
	uint32_t	numAlleles;				// This IS a proeprty of the data (blockset)
	double		meanWeight;				// Mean positional weight of input positions but ONLY FOR informative sites
										//	Sasly, this is set by outside code Insight2::extractInputLoci [Sum_nsites(W_nsite/MaxWeught)]/ |Nsites|
	// double		beta1, beta2, beta3; // this is not a property of the data, retahr it is part of the model.

	BlockSet(void) { clear(); ancPri.clear(); };
	~BlockSet(void) {};
	
	void clear() { blocks.clear();  blockMaker.clear(); blockSetValues(0, 1.0);	};
	void blockSetValues( uint32_t Nalleles, double MaxObsPerPos = 1.0 ) { numAlleles = Nalleles; meanWeight = 1.0 / MaxObsPerPos; };
	void addBlock(const Block_t &newBlock) { blocks.push_back(newBlock); };
	inline uint32_t numBlocks() const { return((uint32_t)blocks.size()); };	// 64->32 
	//
	inline uint32_t numMonoPat(void) const  { uint32_t c = 0; for (uint32_t i = 0; i < numBlocks(); i++) c += blocks[i].monoNumPat();  return(c); };
	inline uint32_t numPolyHPat(void) const { uint32_t c = 0; for (uint32_t i = 0; i < numBlocks(); i++) c += blocks[i].polyHNumPat(); return(c); };
	inline uint32_t numPolyLPat(void) const { uint32_t c = 0; for (uint32_t i = 0; i < numBlocks(); i++) c += blocks[i].polyLNumPat(); return(c); };
	// this is tne humber of unique values (may be multiple observations of each), needed for memory allocation
	inline uint32_t numMonoSites(void) const {  uint32_t c = 0; for (uint32_t i = 0; i < numBlocks(); i++) c += blocks[i].monoNumSite();  return(c); };
	inline uint32_t numPolyHSites(void) const { uint32_t c = 0; for (uint32_t i = 0; i < numBlocks(); i++) c += blocks[i].polyHNumSite(); return(c); };
	inline uint32_t numPolyLSites(void) const { uint32_t c = 0; for (uint32_t i = 0; i < numBlocks(); i++) c += blocks[i].polyLNumSite(); return(c); };
	// this is the number of positions, needed for parameter estimates
	inline void polyCountsSites(uint32_t &HF, uint32_t &LF) const {
		HF = LF = 0;
		for (uint32_t i = 0; i < blocks.size(); i++) { 
			HF += blocks[i].polyHNumSites();
			LF += blocks[i].polyLNumSites(); };
		return;
	};
	inline void polyCountsPat(uint32_t &HF, uint32_t &LF) const {
		HF = LF = 0;
		for (uint32_t i = 0; i < blocks.size(); i++) {
			HF += blocks[i].polyHNumPat();
			LF += blocks[i].polyLNumpat();
		};
		return;
	};
};

