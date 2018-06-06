#pragma once

// Types used throught hte project belong here.
// ESPECIALLY types definitions / aliases that represent type-resolved instantiations of templates.

// system includes
#include <cstdint>		// uint8_t ... etc

// butils
#include "stringer.h"

// local includes
#include "bitfield.h"		// fieldsNarrow - use 8 byte/record bitfield packed record definition rather than 12 byte unpacked.
#include "posMapper.h"		// Maps to/from chromosome space, to linear space by presuming sorted, nonoverlapping chromosomes.

typedef uint32_t		LinAddr_t;			// # of positions in DB, fits in uit32 (4gb)
typedef fieldsNarrow	CovBits_t;			// #bitfields, require 8-12 bytes per record, use 8 byte (narrow) version....
typedef bitfieldDB< CovBits_t, LinAddr_t > covDB_t;
typedef stringer string;					// use enhanced string object, publically rerived from std::string, so casting is easy.
typedef posMapper<> posMapperHG19_t;		// the default for this template is hg19 autosomal positions, it may be overridden and loaded from a file, but this is not reccomended.
