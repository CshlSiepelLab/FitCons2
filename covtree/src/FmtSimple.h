#pragma once

#include <ctime>
#include <cstdio>

#include "stringer.h"

// Simple & slow formatter for value to string conversion
struct FmtSimple {
	typedef stringer string;
	static const uint32_t buf_size = 1024;
	static char *buf() { static char bbuf[buf_size]; return(&(bbuf[0])); };
	static string fmt(const string &s, uint16_t w) { sprintf(buf(), "%*.*s", w, w, s.substr(0, buf_size - 1).c_str()); return(buf()); };
	static string fmt(double d, uint16_t w, uint16_t p, bool sci = false) { 
		string fm = "%" + std::to_string(w) + "." + std::to_string(p) + (sci ? "le" : "lf");	
		sprintf(buf(), fm.c_str(), d); 
		string r(buf());
		return r; };
	// treat double as unsigned integer
	static string fmt2(double d, uint16_t w, bool comma = false, bool zeros = false) { d=std::max(d,0.0); return fmt(uint64_t(d + 0.5), w, comma, zeros); };
	static string fmt(uint64_t v, uint16_t w, bool comma = false, bool zeros = false) {
		string fm = "%";	
#if defined(_WIN32) || defined(_WIN64)
		// Windows comma seperator is broken so ignore it
		fm += (comma ? "" : ""); fm += (zeros ? "0" : ""); fm += std::to_string(w); fm += "I64u";
#else
		fm += (comma ? "'" : ""); fm += (zeros ? "0" : ""); fm += std::to_string(w); fm += "llu";
#endif
		sprintf(buf(), fm.c_str(), (unsigned long long)v);
		string r(buf());
		return r;
	}
	static double toDouble( const string &s ) { double d=0.0; sscanf(s.c_str(),"%lf",&d); return d; };
	static string strNow() { time_t now = time(0);	tm* localtm = localtime(&now);	return( string(asctime(localtm)).trim() ); };
};