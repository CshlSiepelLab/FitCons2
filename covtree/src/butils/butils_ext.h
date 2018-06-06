#pragma once

// to enable namespace encapsulation, all global defines from thsi namespace (such as stdio)
// need to be replicated at the namespace / module level. To encapsulate the namespace
//	Say with bar, bar must also have an an external include that proceeds
// bar_ext.h 
//		#include "butils_ext"
// This allows global namespace encapsulation when aggregating, while preserving file
//	autonomy when not aggregating.

#include <cstdio>
#include <climits>
#include <cmath>
#include <iomanip>

#if defined( _WIN32 ) || defined( _WIN64 )
   #include <direct.h>		// Required for _mkdir in windows, not present in unix
   // #include <windows.h>	// FlushFileBuffers
   #define WINBASEAPI __declspec(dllimport)
   typedef int BOOL;
   #define WINAPI __stdcall
   // #define _In_  _SAL2_Source_(_In_, (), _Pre1_impl_(__notnull_impl_notref) _Pre_valid_impl_ _Deref_pre1_impl_(__readaccess_impl_notref))
   typedef void * HANDLE;
   // extern "C" WINBASEAPI BOOL WINAPI FlushFileBuffers(	_In_ HANDLE hFile );
   extern "C" WINBASEAPI BOOL WINAPI FlushFileBuffers(	HANDLE hFile );
   //int FlushFileBuffers( void * Handle );
   #include <errno.h> 
   #include <io.h>			// _get_osfhandle
#else
   #include <cstdlib>
   #include <cstring>
   #include <limits>
#endif


#include <vector>
#include <string>

//RTable
#include <vector>
#include <map>
#include <fstream>

#include <algorithm>
#include <functional>	// templates for common functions

#include <chrono>		// for Sleep functionality
#include <thread>		// for Sleep functionality
