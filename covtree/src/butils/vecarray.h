#pragma once

// these form part of the interface, and are drawn from global stl, so should be outside the class scope...
#include <vector>
#include <string>

#include "fastfile.h"

// Accelerate tragicly slow stdio processing of numbers by converting object vectors into
//	object arrays and writing the binary form of hte array to a cache file. The process is reversed
//	for quick loads, and performed transeprently. Cache file is created if it is not present, and read
//	from if it is present. Fallback to slow read and cache generation. Cahcread/Write may be supressed w/ user flags
// Original file need not be present if Cache is!
// Objects must be solid, that is, no pointers / references and must have a copy/assignment constructor
//	(usually created by default).
//	

template<class TC> 
class VecArray {
	// limit scope of includes to class, that way we don't need a namespace to keep the default space clean.
public:
	typedef std::vector<long>::size_type	size_v;		// size type for vector indicies
	//typedef std::size_t size_v;
	typedef fastfile::FFStatus FFStatus;
private:
	static const unsigned int CACHE_BUF_SIZE = 1024;	// string representing the number of objects and memory footprint (size+pad) of each object.
	size_v	item_size;									// size of each object, 
	size_v  item_count;									// number of objects in cache
	char	cache_buf[ CACHE_BUF_SIZE ];				// sprintf() of size and count in this buffer for easy dump....
	//
	TC		*base;										// pointer to array form of object vector.

	bool	fileReadable( const std::string &fn )		{	FILE *f=fastfile::ffopen(fn.c_str());if (f==NULL) return( false);fastfile::ffclose(&f); return( true); };
	bool	fileReadCached( const std::string &fn )	{
		bool ret = false;
		if (!fileReadable(fn)) return ret;
		FILE *f=fastfile::ffopen(fn.c_str());
		if (f==NULL) return( ret);
		clear();
		if (fread(cache_buf,sizeof(char),CACHE_BUF_SIZE,f)==CACHE_BUF_SIZE) {	// read string w/ number of objects & object size
			long long unsigned int isv=0, icv=0;
			if (sscanf(cache_buf,"%llu\t%llu",&isv,&icv)==2) {
				item_size =(size_v)isv; item_count=(size_v)icv;
				base=(TC *)calloc(item_count,item_size);						// allocoate array
				if (base != NULL) {
					if (fread(base,item_size,item_count,f)==item_count) ret=true;	// read objects into array, all data, one sys call.
				}
			}
		}
		fastfile::ffclose(&f); return(ret );
	}

	bool	fileWriteCached( const std::string &fn ) {
		bool ret = false;
		FILE *f=fastfile::ffopen(fn.c_str(),false);
		if (f==NULL) return( ret);
		sprintf(cache_buf,"%llu\t%llu",(long long unsigned int)item_size,(long long unsigned int)item_count);				// generate header string
		if (fwrite(cache_buf,sizeof(char),CACHE_BUF_SIZE,f)==CACHE_BUF_SIZE) {	// dump cashed data to disk in 1 syscall...
			if (fwrite(base,item_size,item_count,f)==item_count) ret=true;
		}
		fastfile::ffclose(&f); return(ret );
	}

public:

	void clear( bool first = false ) {
		// if "first" assume pointers may hold garbage and just clear them, otherwise free nonzero pointers..
		item_size = item_count = 0; 
		cache_buf[0]='\0';
		if (first) base = (TC *) NULL;
		if (base!=NULL) free( base );
		base = NULL;
	}

	// for compatibility with vector classes
	inline  size_v size( void ) { return( item_count ); };				// 
	inline	TC & operator[] (size_v index ) { return base[index]; }		// Potentially unsafe!!

	inline TC *arrayBase( void ) { return( base); };					// carefull, this can be null....

	VecArray(void) { clear( true ); }

	virtual ~VecArray(void) { clear(); }

	// Generate a fName.cache file from exusting array object
	// will force overwrite of existing cache file!
	static bool ToFile( const std::string &fName ) { return( fileWriteCached( fName ) ); };

	// Accelerate read of data int array, if possible, otherwise perform slow read on vector file
	//	CacheRead	- use cache if possible
	//	CacheWrite	- write a cache file from input data, if no cache exists now
	//	vin			- if non null, return vector form of data here. otherwise temporaro vector is created then deallocated.
	bool FromFile(const std::string &fName, bool CacheRead = true, bool CacheWrite = true, std::vector<TC> *vin = NULL) {
		using fastfile::SEOF;
		bool	ret=false;
		std::string fnCache=fName + ".cache";
		if (CacheRead && fileReadable(fnCache)) {
			// if a cache exists, use it!
			if (! fileReadCached(fnCache) ) return( false );
			if (vin != NULL) ToVec( *vin );
			ret = true;
		} else {
			// if not read slowly from original data file
			std::vector<TC> v, *vin_t=NULL;
			vin_t=( vin==NULL? &v : vin);
			// Slow read, TC must support a function: FFStatus TC.readln( FILE *f)
			if (FileToVec( fName, *vin_t ) != SEOF ) return( false );
			FromVec( *vin_t ); // convert from Vec to Array
			// write a cache file, unless supressed...
			if (CacheWrite) ret = fileWriteCached( fnCache );
		}	// vector is deallcoated as v leaves scope...
		return( ret );		
	}

	// Allocate and fill an array from a vector
	bool FromVec( std::vector<TC> &v ){
		clear();
		item_size = sizeof(TC); item_count = v.size();
		base=(TC *)calloc(item_size,item_count);				// allcoate array memory
		if (base == NULL) return( false );
		// memcpy might be faster....
		for ( size_v i=0; i< v.size(); i++) { base[i]=v[i]; }	// copy each object from vector to array
		sprintf(cache_buf,"%llu\t%llu\t",(long long unsigned int) item_size, (long long unsigned int )item_count);
		return( true );
	}

	// Export array data to a vector...
	bool ToVec( std::vector<TC> &v ){
		v.resize( item_count );
		for ( size_v i=0; i< item_count; i++) { v[i]=base[i]; }
		return( true );
	}

	// STATIC METHODS for 
	// FileToVec		- simple read a file into a vector, one line at a time, slow / uncached
	// FileToVecCache	- use the object to performe a cached read, if possible. Also can write a cache file
	// VecToFile		- use the object to generate a cache from existing vector, and write it.

	// Read a Fiel into a vector of TC using TC's .readln() method (it must have one).
	static FFStatus FileToVec( const std::string &fName, std::vector<TC> &v ) {
		using fastfile::SEOL;
		using fastfile::SEOF;
		using fastfile::SERR;
		// read a set of records from a file
		FFStatus	s=SEOF;
		TC			item;
		FILE		*f=fastfile::ffopen(fName.c_str());		// open file, read mdoe is default
		if (f==NULL) return SERR;
		do {
			s=item.readln(f);							// read one object, object must supprot readln() function...
			if ( s==SEOL ) v.push_back(item);			// store in vector
		} while ( (s!=SEOF) && (s!=SERR) );
		fastfile::ffclose(&f); return( s );
	};

	// Static form, just get me the vector with cache acceleeration, if possible.
	// CacheWrite - generate a cache file if none exists and slow read suceeded.
	static bool FileToVecCached( const std::string &fName, std::vector<TC> &v, bool CacheWrite = true ) {
		VecArray<TC> va;
		return( va.FromFile( fName, true, CacheWrite, &v ) );	// array data is deallocated by destructor when va leaves scope
	}

	static bool VecToFile( const std::string &fName, std::vector<TC> &v ) {
		VecArray<TC> va;
		va.FromVec( v );
		return( va.ToFile(fName ) );
	};

};
