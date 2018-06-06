#pragma once
#include "inscomp_ext.h"

namespace inscomp {
//	#include "Rtable.h"
	#include "WLocus.h"
	#include "Block.h"
	#include "BlockSet.h"
	typedef double			Real_t;
	typedef unsigned char	AncPriIdx_t;
	typedef uint32_t		PatternCounts_t;
	typedef inscomp::BlockSet<Real_t, AncPriIdx_t, PatternCounts_t>	BlockSet_t;
	typedef inscomp::Block<Real_t, AncPriIdx_t, PatternCounts_t>	Block_t;
	#include "BlockCompiler.h"
	#include "modelInsight.h"
	#include "modelInsightAccelerator.h"
	#include "modelInsightPosterior.h"
	#include "modelInsightSup.h"
	#include "optBeta.h"
	#include "optInsight.h"
}