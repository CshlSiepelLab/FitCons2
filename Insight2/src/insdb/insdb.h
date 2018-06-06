#pragma once

// These shoudl be first incldued in global namespace, If they have not yet bee, then this IS the 
//	global namespace, so include them now....
#include "insdb_ext.h"
#include "butils/butils.h"

namespace insdb {
	#include "Locus.h"
	#include "Poly.h"
	#include "Block.h"
	#include "Mono.h"
	#include "DB.h"
}
