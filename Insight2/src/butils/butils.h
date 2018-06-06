#pragma once

// if butils is included solo, then this makes sure the needed global declariations are provided
// if butils is encapsulated in another namespace (foo), then butils_ext wehoudl be incldued in THAT
//	namespace's foo_ext.h file...
#include "butils_ext.h"


// wrap all utilities into useful namespace to help isolate it from
//	outside dependencies.
namespace butils {
	#include "util.h"		// utility definitions...
	#include "fastfile.h"	// fast single threaded IO
	//#include "LogProb.h"	// change to lowercase
	//#include "phalf.h"	// half precision floating point probability storage class
	#include "vecarray.h"	// fast disk IO for reading vectors and convertign to arrays.
	#include "RTable.h"
	#include "optimizer.h"	// lbfgs gradient based optimizer. Wrapped in class. Include only.
	#include "ezmatrix.h"
	#include "stringer.h"	// some extensions of std::string, including tokenizer...
	#include "mathplus.h"
}