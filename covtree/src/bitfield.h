#pragma once
#include <vector>

#include "vecarray.h"
#include "fastfile.h"

#include "FmtSimple.h"
#include "posMapper.h"
// T must be an unsigned intergral type, reccomended uint8_t, uint16_t, uint32_t or uint64_t
// For FitCons2 I reccomend T as uint64_t... sadly, we just missed fitting into 32 bits...

struct fieldsWide {		// 12 bytes
	uint32_t length;	// number of position with these tagBits
	uint64_t tagb;		// tag for these positions, in bit packed format
	static const uint32_t MAX_LENGTH=(uint32_t )-1;
	static const uint8_t  MAX_TAGBITS = 64;
	void clear() { length=0; tagb=0;}
	fieldsWide() { clear();};
};
struct fieldsNarrow { // 8 bytes, requries special handling of addresses to break blocks into smaller units...
	uint64_t length : 16;	// 16 bits of 64
	uint64_t tagb   : 48;	// 48 bits of 64...
	static const uint16_t MAX_LENGTH=(uint16_t)-1;
	static const uint8_t  MAX_TAGBITS = 48;
	void clear() { length = 0; tagb = 0; }
	fieldsNarrow() { clear(); };
};

// Use principally for IO. Defines a conversion between a tag definition and a bitfield representation.
template< typename TBitfieldBlock>
class bitfieldDef {
public:
	typedef uint8_t						slen_t;
	typedef std::vector<slen_t>			vslen_t;
	typedef std::vector<TBitfieldBlock>	vT;
private:
	vslen_t	tagLengths;	// Number of charcters used torepresent each string field... in characters
	vslen_t	bitLengths;	// size of each bit field, in bits...
	vT		masks;		// bit masks for field (all bits 1), dont use the length field, just the tagb field
	vslen_t	offsets;	// shift values to align bit fields into bit record (tagb).
	slen_t	numBits;	// total number of utilized bits in bitfield (bit encoded record, bits)
	slen_t  fieldLen;	// length of string tag record (length of string encoded record, characters)
	static const slen_t MAX_BITS = TBitfieldBlock::MAX_TAGBITS;	// container size, size of maximim definable bit encoded record tfor thsi container (bits)

																// SAFE, high performance, conversion from a known-length ascii string of digits to an unsigned number
																// Checks for null ptr, end of string and nondigit character.
	inline bool qatou(const char *st, slen_t len, slen_t &res) const {
		slen_t v = 0; res = 0;
		// len is always >=1. Gnerally len=1, sometimes len=2....
		if (st == NULL) return false;
		res = (slen_t)(((unsigned char)*st) - '0');    if (res>9) return false;
		while (--len>0) {
			v = (slen_t)(((unsigned char)*(++st)) - '0');  if (v>9) return false;
			res = (res << 3) + (res << 1) + v;
		};
		return true;
	};
public:
	void clear() { tagLengths.clear(); bitLengths.clear(); masks.clear(); offsets.clear(); numBits = 0; fieldLen = 0; return; };
	bitfieldDef() { clear(); };
	~bitfieldDef() { clear(); };

	inline slen_t numFields() const { return( (slen_t) tagLengths.size()); };

	// terribly inefficent, but only done once.. so make it clear.
	// Fields are read left to right, like numerials. So we assemble our field definitions, left to right
	//	lowest indicies are "first" or conceptually "most left" as read, that is most significant in a 32 bit (unsigned) word...
	//  tags are numbered of values (1-N), with 0 meaning "missign data", so we need to accomidate N+1 possibilities
	//	thus, N+1 bits per field. Each possibility is a unique "1" in the bitfield.
	bool addField(slen_t StringTagWidth, slen_t NumTags) {
		slen_t n = numFields(); slen_t newbits = NumTags + 1;
		// make sure we are not trying to add more bits than our base type will accomidate
		if (numBits + (uint64_t)newbits > (uint64_t)MAX_BITS) return false;
		// Shift all existing fields to the left
		for (slen_t i = 0; i < n; i++) {
			masks[i].tagb <<= newbits; offsets[i] += newbits;
		};
		// Add the new field
		TBitfieldBlock newblock;																// add the mask
		for (slen_t i = 0; i<newbits; i++) newblock.tagb = ((newblock.tagb << 1) | 1);
		masks.push_back(newblock);
		tagLengths.push_back(StringTagWidth); bitLengths.push_back(newbits); numBits += newbits;// add tagLen and bitLen values
		offsets.push_back(0);																	// add a new 0 offset item...
		fieldLen += StringTagWidth;																// increment total field length (chars)
		return true;
	};

	//This is fast and SAFE. To make it faster remove checks in first two lines...
	inline bool accumulateMask(uint8_t CovNum, uint8_t CovValNum, TBitfieldBlock &bitFieldBlock) const {
		if (CovNum    >= numFields()) return false;
		if (CovValNum >= bitLengths[CovNum]) return false;
		bitFieldBlock.tagb |= (((uint64_t)1) << (CovValNum + offsets[CovNum]));
		return true;
	}

	//This is SAFE... 
	// Gnerate a bit record from a tagRecord composed of a string of fixed-width tagFields.
	inline bool tagToBits(const char *TagField, TBitfieldBlock &bitFieldBlock) const {
		if ((TagField == NULL) || (strlen(TagField) < fieldLen)) return(false);
		slen_t v = 0, n = numFields(); const char *p = TagField; bitFieldBlock.clear();
		for (slen_t i = 0; i < n; ++i) {
			if (!qatou(p, tagLengths[i], v)) return false;			// 
			if (v >= bitLengths[i]) return false;					// value in required range?
																	// bitField |= ((1 << (v + offsets[i])) & masks[i] );	// safety
			bitFieldBlock.tagb |= (((uint64_t)1) << (v + offsets[i]));
			p += tagLengths[i];
		}
		return true;
	};

	// UNSAFE. Tag Field must contain enough space. BitField must be valid. etc...
	//	Uses low performance IO,  mostly for debugging...
	bool bitsToTag(const TBitfieldBlock &bitFieldBlock, char *TagField) const {
		slen_t n = numFields(); char *p = TagField; unsigned long tagval = 0;
		for (slen_t i = 0; i < n; ++i) {
			tagval = (bitFieldBlock.tagb & masks[i]) >> offsets[i];
			sprintf(p, "%0*lu", (int)tagLengths[i], (int)tagval); p += tagLengths[i];
		}
		return true;
	}

	// A SetMask defines a subset of acceptable patterns. Each mask has at least one bit per field, and each record has 
	//	exactly one bit per field.
	static bool recordInSet(const TBitfieldBlock &bitSet, const TBitfieldBlock &bitRecord) {
		// Handle missing data? Genrally missing data fieldd has first position set, so handle this by using suitable bitset / mask.
		return ((bitRecord.tagb & bitSet.tagb) == bitRecord.tagb);
	};	// is a bit field in a mask?
};


template< typename TBitfieldBlock>
class bitfield {
public:
	typedef bitfieldDef<TBitfieldBlock>		bmap_t;		// map between string of tags, and encoded bitfields
private:
	static const bmap_t *bmap;// = constexpr NULL;
	TBitfieldBlock data;
public:
	void clear() { data.clear(); }
	bitfield() {};
	~bitfield() {};
	inline uint64_t getLength() const { return(data.length); };
	inline uint64_t getTagBits() const { return(data.tagb); };
	static TBitfieldBlock old_data;
	static uint64_t	old_start, old_end;
	// Requires pDed and PositionMapper
	static void setMap(const bmap_t *BMap = NULL) { bmap = BMap; old_data.clear(); old_start = old_end = 0; };
	fastfile::FFStatus readln(fastfile::FILE *f) {
		using fastfile::FFStatus; FFStatus fst;
		static char buf_chr[1024], buf_tag[1024];	// way too big... chr is usually <6 and tags are usually < 10
		data.clear();
		if (old_end != 0) {	// if we are still working on a large locus, jsut return the next "max LEN" segment
			uint64_t size_left = old_end - old_start + 1;
			if (size_left <= data.MAX_LENGTH) {
				data.length = size_left; data.tagb = old_data.tagb;  old_data.clear(); old_start = old_end = 0;
			} else {
				data.length = data.MAX_LENGTH; data.tagb = old_data.tagb; old_start += data.MAX_LENGTH;
			}
			return(FFStatus::SEOL);
		};
		uint64_t pst = 0, pen = 0;
		if (bmap == NULL) return(FFStatus::SERR);
		buf_chr[0] = buf_tag[0] = 0; pst = pen = 0;
		fst = fastfile::getString(f, buf_chr);
		if ((fst == FFStatus::SEOF) && (buf_chr[0] == 0)) return(FFStatus::SEOF);			// no chars read, and EOF, normal end of file...
		if (fst != FFStatus::SEOR) return(FFStatus::SERR);
		if (fastfile::getInt(f, &pst) != FFStatus::SEOR) return(FFStatus::SERR);
		if (fastfile::getInt(f, &pen) != FFStatus::SEOR) return(FFStatus::SERR); pen--;		// Convert from half open to fully closed...
		fst = fastfile::getString(f, buf_tag);
		if (fst == FFStatus::SEOR) fst = fastfile::clearLine(f);
		if (fst == FFStatus::SERR) return(FFStatus::SERR);
		// Got a record, and if here, we must have the end of line or end of file....
		// We assume sorted bedg, that partitions the entire genome, so we jsut need to maintain the width ofthis locus.
		//	NOTE if we move to a 16 bit address, we might need to create multiple 64 k blocks if len is too long....
		if (!bmap->tagToBits(buf_tag, data)) return(FFStatus::SERR);
		if (pen - pst + 1 > data.MAX_LENGTH) {				// big locus? break it down into MAX_LEN chunks...
			old_data = data; old_start = pst; old_end = pen; return(readln(f));
		}
		data.length = pen - pst + 1;			// already converted to closed intervals...
		return(FFStatus::SEOL);
	}
};

template <typename TBitfieldBlock> TBitfieldBlock bitfield<TBitfieldBlock>::old_data;
template <typename TBitfieldBlock> uint64_t bitfield<TBitfieldBlock>::old_start;
template <typename TBitfieldBlock> uint64_t bitfield<TBitfieldBlock>::old_end;
template <typename TBitfieldBlock> bitfieldDef<TBitfieldBlock> const * bitfield<TBitfieldBlock>::bmap;
//static const bmap_t *bmap

// An array of BitfieldBlocks that can be read from an tagged ascii .bed file OR a dense cached binary file using VecArray
template< typename TBitfieldBlock, typename TLinearAddress>
class bitfieldDBElement {
public:
	typedef typename bitfield<TBitfieldBlock>::bmap_t			bmap_t;		// map between string of tags, and encoded bitfields
	typedef stringer string;
private:
	typedef bitfield<TBitfieldBlock>	bitrec_t;	// single address/bitfield record. with .readline() method
	typedef VecArray<bitrec_t>			bitrecs_t;	// structure for accelerated loading and caching of bitrecs.
	bitrecs_t							bitrecs;
	const bmap_t						*bmap = NULL;
	string								elemName = "", fileName = "";
	// Temporary values used during scan... Convert to 64 for larger genomes.... chr is chromosome name (keep it small), buf_tag is string tag for each entry
	char								buf_chr[128], buf_tag[128];
	TLinearAddress						scan_locIdx, scan_locStart, scan_locEnd;
public:
	void clear() { bitrecs.clear(); bmap = NULL; elemName = fileName = "";  buf_chr[0] = buf_tag[0] = 0; scan_locIdx = scan_locStart = scan_locEnd = 0; };
	bitfieldDBElement() { clear(); };
	~bitfieldDBElement() { clear(); };
	bool readDB(const string &fname, const string &ElemName, const bmap_t *BMap, TLinearAddress &PositionsRead) {
		if (BMap== NULL) return(false);
		clear(); bmap = BMap; bitrec_t::setMap( bmap );		// set static pointer to tag and address itnerperters 
		elemName = ElemName; fileName = fname;
		// bool ok = bitrecs_t::FileToVecCached(fname, bitrecs);
		bool ok = bitrecs.FromFile( fname );
		PositionsRead = 0;
		for (typename bitrecs_t::size_v i=0; i<bitrecs.size(); i++) {
			PositionsRead += (TLinearAddress) bitrecs[i].getLength(); };
		return ok;
	}
	bool scanInit() {
		// initialize a monotonicly increasing scan of positions
		scan_locIdx = 0;
		scan_locStart = 0;
		scan_locEnd = (TLinearAddress) (bitrecs[scan_locIdx].getLength() - 1);
		return true;
	}
	inline bool scanPos(const TLinearAddress scanpos, TBitfieldBlock &tagbits) {
		// scanpos must be nonddecreasing between calls to scanInit(). 
		while ((scanpos> scan_locEnd) && (scan_locIdx<bitrecs.size())) {
			scan_locStart = scan_locEnd + 1;															// new locus starts imemdiately after end of current one
			scan_locIdx++;																				// increment locus pointer
			scan_locEnd = (TLinearAddress)  (scan_locStart + bitrecs[scan_locIdx].getLength() - 1);		// last position of new locus.
		}
		if (scanpos>scan_locEnd) return(false);	// check for overrun
		tagbits.length = scan_locEnd - scanpos + 1;
		tagbits.tagb   = bitrecs[scan_locIdx].getTagBits();
		return true;
	}
};


// AddressFull is uint32_t for hg19, but could be enlarged to uint64_t for larger genomes
// TAddressCompact is a locus width and initially is uint32_t. If need to reduce it use uint16_t, but this might mean splitting blocks.
// TBitField is uint64_t, however 
template< typename TBitfieldBlock, typename TLinearAddress >
class bitfieldDB {
public:
	typedef stringer string;
	typedef bitfieldDBElement<TBitfieldBlock, TLinearAddress> dbelem_t;
	typedef typename dbelem_t::bmap_t bitfieldDef_t;
	struct	dbElemSrc { string name; string fname; void clear() {name=fname="";}; dbElemSrc() {}; };
	typedef typename std::vector<dbElemSrc> vdbElemSrc_t;
	struct	scanTarget {							// Each ScanTarget represents a single candidate covariate split
		TBitfieldBlock	patA1, patC1, patA2, patC2;	// Patterns for Annotations & CTS patterns for mask1 (A1&C1) and its conditional compliment (A2&C2)
		uint8_t			hits1, hits2;				// really only need to be uint8 or uint16, large enough count up to number of cell types
		void clear()	{ patA1.clear(); patA2.clear(); patC1.clear(); patC2.clear(); hits1 = hits2 = 0; }
		scanTarget()	{ clear(); }; ~scanTarget() { clear(); };
	};
	typedef std::vector<scanTarget>	vscanTarget_t;
private:
	// database elements
	string					dbase;
	dbElemSrc				annSrc;
	vdbElemSrc_t			ctsSrc;
	// Most of the data is in the following 2 elements... generally many gigabytes of memory...
	dbelem_t				dbAnn;	// Annotations, shared among all
	std::vector<dbelem_t>	dbCts;	// Cell Type Specific bitfields
	// used in scanning
	TLinearAddress			scan_lastStart, scan_lastLength;
	bitfieldDef_t			tagmapAnn, tagmapCTS;
public:
	void clear() { dbase=""; annSrc.clear(); ctsSrc.clear(); dbAnn.clear(); dbCts.clear(); tagmapAnn.clear(); tagmapCTS.clear(); };
	uint32_t NumCells() const { return((uint32_t) dbCts.size()); };
	bitfieldDB()	{ clear(); };
	~bitfieldDB()	{ clear(); };
	const bitfieldDef_t &mapAnn() const { return tagmapAnn; };
	const bitfieldDef_t &mapCTS() const { return tagmapCTS; };
	// given a mapping from tags to bitfields, load each Annotation and CellTypeSpecific (CTS) tagset into memory. We assume that input
	//	file is a partition and uniquely spans the LinearAddress Space with full coverage....
	bool loadDB( const string &DirBase, const dbElemSrc &AnnElem, const bitfieldDef_t &AnnTagmap, const vdbElemSrc_t &CtsElem, const bitfieldDef_t &CtsTagmap, const TLinearAddress NumPos, const int verb=0 ) {
		string fin; TLinearAddress pos_read=0; // dbelem_t tmp_elem;

		if (verb>0) { std::cerr << FmtSimple::strNow() << " CovDBLoad: Starting.\n"; };
		clear();
		dbase=DirBase;
		tagmapAnn = AnnTagmap; tagmapCTS=CtsTagmap;

		// Read Annotations
		fin=(AnnElem.fname.substr(0,1)=="/" ? AnnElem.fname : (DirBase != "" ? string(DirBase + "/" + AnnElem.fname) : AnnElem.fname));
		if (verb>1) { std::cerr << FmtSimple::strNow() << " CovDBLoad: Loading Annotation database from " << fin << ".\n"; };
		if (!dbAnn.readDB( fin, AnnElem.name, &tagmapAnn, pos_read)) return false;
		if (pos_read != NumPos ) return false;

		// Read CTS data
		if (verb>1) { std::cerr << FmtSimple::strNow() << " CovDBLoad: Loading Cell Type Sensative (CTS) database.\n"; };
		if (verb>3) { std::cerr << FmtSimple::strNow() << " CovDBLoad: CTS size " << CtsElem.size() << " \n"; };
		dbCts.clear(); dbCts.resize(CtsElem.size());
		for (uint32_t i = 0; i < CtsElem.size(); i++) {
			fin = (CtsElem[i].fname.substr(0, 1) == "/" ? CtsElem[i].fname : (DirBase != "" ? string(DirBase + "/" + CtsElem[i].fname) : CtsElem[i].fname ) );
			// dbCts.push_back( tmp_elem );
			if (verb>2) { std::cerr << FmtSimple::strNow() << "\tCovDBLoad: database " << fin << "\n"; };
			if ( !dbCts[i].readDB(fin, CtsElem[i].name, &tagmapCTS, pos_read) ) return( false );
			if (pos_read != NumPos) return false;
		}

		if (verb>0) { std::cerr << FmtSimple::strNow() << " CovDBLoad: Database Loadded.\n"; };
		return true;
	};

	bool scanpatInit() {
		dbAnn.scanInit(); for ( auto & dbe : dbCts ) { dbe.scanInit(); };
		scan_lastStart = scan_lastLength = 0;
		return true;
	}

	// Loop through DB, match annotations. PosStart must be nondecreasing, so it is to to caller to add "Length" to next PosStart
	//	and make sure we dont overrun the LinearAddress space. LocusLength is the longest span of positions that has consistent annotations
	//	in each cell type.
	// thsi is arranged so a collecxtion pf patterns representing all possible "covariate splits (about 1000)" is compared at each genomic position.
	//	This avoids thrashing the entire database of positions for each split, as the DB can take 5-100GB of memory.... improving coherency of
	//	memory access is a big win. Modern CPU's have mem throughput of 10-50BG/sec... so just mem access can be slow...
	bool scanpatGetNext(vscanTarget_t &Targets, TLinearAddress &PosStart, TLinearAddress &LocusLength ) {
		LocusLength=0;
		TLinearAddress tmp_len=0; 

		// Get block length
		TBitfieldBlock tmp_block;
		dbAnn.scanPos(PosStart, tmp_block); tmp_len = tmp_block.length;
		if ((LocusLength == 0) || (tmp_len <LocusLength)) LocusLength = tmp_len;
		for (auto & ct : dbCts) {
			ct.scanPos(PosStart, tmp_block); tmp_len = tmp_block.length;
			if ((LocusLength == 0) || (tmp_len <LocusLength)) LocusLength = tmp_len;
		}

		// Count pattern matches
		TBitfieldBlock &tmp_covpat=tmp_block;							// just rename the variable...
		for (auto & targ : Targets) {
			targ.hits1 = 0;
			dbAnn.scanPos(PosStart, tmp_covpat);
			if (bitfieldDef_t::recordInSet(targ.patA1, tmp_covpat)) {	// If annotations dont match, we have 0 hits, don't bother looking at CTS data
				for (auto & ct : dbCts) {								// Count number of cell types with assay patterns that fall within the mask...
					ct.scanPos(PosStart, tmp_covpat); 
					if (bitfieldDef_t::recordInSet(targ.patC1, tmp_covpat)) targ.hits1++;
				}
			}

			targ.hits2 = 0;
			dbAnn.scanPos(PosStart, tmp_covpat);
			if (bitfieldDef_t::recordInSet(targ.patA2, tmp_covpat)) {
				for (auto & ct : dbCts) {
					ct.scanPos(PosStart, tmp_covpat);
					if (bitfieldDef_t::recordInSet(targ.patC2, tmp_covpat)) targ.hits2++;
				}
			}
		}
		return true;	// for each Target pattern mask (subset of covariates) we have coutned the number of cell types with matchng civariates.

		}

};

