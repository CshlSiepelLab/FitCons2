#pragma once

// Deal  with things that shoudl be defined in c++ but may not be...

#if defined( _WIN32 ) || defined( _WIN64 )
// various fixed bit width quantities, not always defiend in usoft C.
//typedef unsigned char			uint8_t;
//typedef unsigned short int 		uint16_t;
//typedef unsigned long int		uint32_t;
//typedef unsigned long long int	uint64_t;
//typedef unsigned long long		int64_t;
#define getc_unlocked(f)		(fgetc(f))
#else

#endif

// useful intergral types 
typedef unsigned char			uchar;
typedef unsigned short int		ushort;
typedef unsigned long int		ulong;
typedef long long int			llong;
typedef unsigned long long int	ullong;


#ifdef FALSE
// small fixed length strings
#define LEN_STRINGT 8		// Tiny:    size of largest integral data type (64 bit ullong)
#define LEN_STRINGS 16		// Small:   typical high precision ascii float + null, or short tag
#define LEN_STRINGM 4096	// Medium:  typical line of text IO
typedef char stringT[ LEN_STRINGT ];
typedef char stringS[ LEN_STRINGS ];
typedef char stringM[ LEN_STRINGM ];
#endif

// fixed length arrays can not be simply assigned. Here is a wrapper
//	that takes no additional memory, but does allow arrays to be treated
//	as structs which can be assignes (id stringT a,b; strcpy(a.c(),"foo"); b=a };
// Genreates a fixed lenth E array of type T
template< typename T, unsigned int N >
struct StringWrapper {
	typedef StringWrapper<T, N> TSW;	// terse reference for argument types
	T		s[N];						// ACTUAL memory allocation, up to padding, struct is just as big as this item.

	// provide access to current length and maximum size for string
	inline ullong  siz()       { return sizeof(T[N]); };
	inline ullong  len()       { return std::strlen(&(s[0])); };	// TODO rewrite, in case T is not char!

	// get a pointer to he underlying string
	      inline T * c()       { return(&(s[0])); };	// emperically for t, &t == t.s == &(t.s[0]), up to type
	const inline T * c() const { return(&(s[0])); };	// emperically for t, &t == t.s == &(t.s[0]).

	// access string elements using the array brackets...
	      inline T & operator[](std::size_t i)       { return(s[i]); };
	const inline T & operator[](std::size_t i) const { return(s[i]); };

	//STL needs these to sort strings, really, only < is needed. To make these safe, convert to strncmp Min(A.siz,B.siz)
	inline bool operator<  (const TSW &B) const { return (strcmp(this->c(), B.c()) <  0); };
	inline bool operator<= (const TSW &B) const { return (strcmp(this->c(), B.c()) <= 0); };
	inline bool operator>  (const TSW &B) const { return (strcmp(this->c(), B.c()) >  0); };
	inline bool operator>= (const TSW &B) const { return (strcmp(this->c(), B.c()) >= 0); };
	inline bool operator== (const TSW &B) const { return (strcmp(this->c(), B.c()) == 0); };

	// Simple & slow & unsafe IO operations. To make safe limit fields size to string siz-1, and insure null termination.
	friend std::ostream& operator<< (std::ostream& os, const TSW & obj) { return os << obj.c(); };
	friend std::istream& operator>> (std::istream& is,       TSW & obj) { return is >> obj.s;   };

	// Simple Assignment
	inline TSW & operator= (const char *a) { strncpy(c(), a, siz() - 1); s[siz() - 1] = '\0'; return(*this); };
};

typedef StringWrapper< char,    8>  stringT; // Tiny:    size of largest integral data type (64 bit ullong)
typedef StringWrapper< char,   16>  stringS; // Small:   typical high precision ascii float + null, or short tag
typedef StringWrapper< char, 4096>  stringM; // Medium:  typical line of text IO

