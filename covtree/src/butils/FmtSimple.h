#pragma once

#include <ctime>
#include <cstdio>

#include "butils/stringer.h"

#if defined( _WIN32 ) || defined( _WIN64 )
// remove warnings about unsafe sprintf usage...
#pragma warning(disable:4996)
#endif

// Simple & slow formatter for value to string conversion
struct FmtSimple {
	typedef stringer string;
	static const uint32_t buf_size = 1024;
	static char *buf() { static char bbuf[buf_size]; return(&(bbuf[0])); };
	static string fmt(const string &s, long w) { string f="%*.*s"; if (w<0) {w=-w; f = "%-*.*s"; }; sprintf(buf(), f.c_str(), w, w, s.substr(0, buf_size - 1).c_str()); return(buf()); };
	static string fmt(double d, long w, uint16_t p, bool sci = false) { 
		string fm = "%" + std::to_string(w) + "." + std::to_string(p) + (sci ? "le" : "lf");	
		sprintf(buf(), fm.c_str(), d); 
		string r(buf());
		return r; };
	// treat double as unsigned integer
	static string fmt2(double d, long w, bool comma = false, bool zeros = false) { d=std::max(d,0.0); return fmt(uint64_t(d + 0.5), w, comma, zeros); };
	static string fmt(uint64_t v, long  w, bool comma = false, bool zeros = false) {
		string fm = std::to_string(v);	// slow, but machine independent....
		if (comma) fm = addCommas(fm);
		if (zeros && fm.length()< (size_t)w ) { string z; z.resize(w-fm.length(),'0'); fm = z + fm; }
		if (fm.length()<(size_t) (w<0?-w:w)) fm = fmt(fm,w);
		return fm;
	}
	static double toDouble( const string &s ) { double d=0.0; sscanf(s.c_str(),"%lf",&d); return d; };
	static string strNow() { time_t now = time(0);	tm* localtm = localtime(&now);	return( string(asctime(localtm)).trim() ); };
	static string addCommas( string const &sin, long  w = 0 ) { 
		size_t sep = sin.find_first_of(".Ee");
		string wrk = sin; wrk.trim();
		string head = wrk, tail = "", res = "";
		if (sep != stringer::npos ) { head = wrk.substr(0,sep); tail=wrk.substr(sep); }
		while (head.length()>0) {
			res = head.tail(3) + (res.length()>0 ? "," : "")  + res;
			head.eraseTail(3);
		}
		res += tail;						// Add trailing bits, if any (exponent, decimal)
		if (w != 0) res = fmt( res, w );	// Restore formatting, w<0 means left justify, w>0 means right (default, if active). 0 means no justify (inactive default)
		return res;
	}
};