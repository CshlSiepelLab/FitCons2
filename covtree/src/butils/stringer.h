#pragma once

#include <string>
#include <vector>
//#include <utility>
#include <algorithm>    // std::transform

// #include <sstream>
// adds some features that are missing in std::string

class stringer : public std::string {
public:
	typedef std::string string;
	typedef std::vector<stringer> vecS;
	typedef string::size_type size_t;
	stringer(const string &foo) : string(foo) {};
	stringer(const char *foo="") : string(foo) {};
	operator const char *() const { return( this->c_str() ); }
	stringer trim(const string &whitesp = " \t\n\r") { 
		size_t s = find_first_not_of(whitesp);
		if (s == string::npos) { *this = ""; return(stringer("")); };
		size_t e = find_last_not_of(whitesp);
		size_t l = (e < s ? 0 : e - s + 1);
		*this = (*this).substr(s, l);
		return( *this );
	};
	void tolower() { std::transform(begin(), end(), begin(), [](unsigned char c) { return ::tolower(c); }); };
	void toupper() { std::transform(begin(), end(), begin(), [](unsigned char c) { return ::toupper(c); });  };
	string aslower() const { stringer a = *this; a.tolower(); return(a); }; // RVO Retrn Value Optimization makes this OK....
	string asupper() const { stringer a = *this; a.toupper(); return(a); };
	bool startswith(const stringer &match) const { return(substr(0, match.length()) == match); }
	string tr(const string &del ) { if (del.length()>0 ) { erase(std::remove(begin(), end(), del[0]),end()); }; return *this; }
	string tr(const string &src, const string &repl) { if (src.length()>0 && repl.length()>0) { std::replace(begin(), end(), src[0], repl[0]); }; return *this; }
	// tokenize assumes any token splits a field, so two successive tokens indicate a null field between.
	void tokenize(vecS &Toks, const string &delim = " \t\n\r") const {
		Toks.clear();
		size_t last_pos = 0;
		for (size_t pos = find_first_of(delim, last_pos); pos != string::npos; pos = find_first_of(delim, last_pos)) {
			Toks.push_back(substr(last_pos, pos - last_pos)); last_pos = pos+1;
		};
		Toks.push_back(substr(last_pos));
	}
	// tokenizeQuoted() is a bit different from tokenize().
	// tokenize quoted skips consecutive delims and leading delims, no empty fields are inferred.
	//	To explicitely get an empty field, use two sucessive quotation marks.
	void tokenizeQuoted(vecS &Toks, const string &delim = " \t\n\r", const string &quotechars = "\"" ) const {
		//static const string quote="\"";
		Toks.clear();
		size_t pos=0,last_pos = find_first_not_of(delim);
		while ( last_pos != string::npos ) {
			//bool quoted = (substr(last_pos,1) == quote); if (quoted) last_pos++;	
			bool quoted = (substr(last_pos, 1).find_first_of(quotechars) != string::npos); if (quoted) last_pos++;
			pos = find_first_of( (quoted ? quotechars: delim), last_pos);
			if (pos == string::npos) { // string ends before next delimiter......
				Toks.push_back(substr(last_pos)); last_pos=pos; continue; };
			Toks.push_back(   substr(last_pos , pos - last_pos) );
			last_pos= find_first_not_of(delim,pos+1);
		}
		return;
	};
	string tail(size_t maxLen) const {
		if (maxLen >= length()) return *this;
		return substr( length()-maxLen );
	}
	string eraseTail(size_t maxLen) {
		if (maxLen >= length()) { *this = ""; return *this; }
		return( erase(length() - maxLen, maxLen ) );
	}
	string reverse() { 
		std::reverse(begin(), end()); return *this; }

	string asReplaced(const string &Pat, const string &Repl) const {	// return string with literal "pat" replaces with "Repl", can be used to delete if Repl is ""
		string res = ""; size_t oldpos = 0, newpos = 0;
		newpos = find(Pat, oldpos);
		while (newpos!=npos) { 
			res += substr(oldpos,newpos-oldpos) + Repl;
			oldpos=newpos+Pat.length();
			newpos = find(Pat, oldpos);
		}
		res += substr(oldpos); 
		return res; 
	}

};

