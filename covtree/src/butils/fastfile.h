#pragma once

#include <cstdint>	// such as uint32_t
#include <cstdio>
#include <cstring>
#include <chrono>	// for Sleep functionality
#include <thread>	// for Sleep functionality

#if defined( _WIN32 ) || defined( _WIN64 )
   #include <direct.h>		// Required for _mkdir in windows, not present in unix
   #include <direct.h>		// Required for _mkdir in windows, not present in unix
   // #include <windows.h>	// FlushFileBuffers
   #define WINBASEAPI __declspec(dllimport)
   typedef int BOOL;
   #define WINAPI __stdcall
   // #define _In_  _SAL2_Source_(_In_, (), _Pre1_impl_(__notnull_impl_notref) _Pre_valid_impl_ _Deref_pre1_impl_(__readaccess_impl_notref))
   typedef void * HANDLE;
   // extern "C" WINBASEAPI BOOL WINAPI FlushFileBuffers(_In_ HANDLE hFile);
   extern "C" WINBASEAPI BOOL WINAPI FlushFileBuffers(HANDLE hFile);
   #include <errno.h> 
   #include <io.h>			// _get_osfhandle
#else
   #include <unistd.h>		// rmdir
#endif
#include <sys/stat.h>	// requires for mkdir() in unix, present in windows
#include <sys/types.h>	// requires for mkdir() in unix, present in windows

#include "util.h"

namespace fastfile {
	typedef std::FILE FILE;
#if defined( _WIN32 ) || defined( _WIN64 )
	#pragma warning(disable: 4996)
	// remove warning about unsafe functions like fscanf.
#endif

	// Generic File status enum for IO
	enum FFStatus { SEOR, SEOL, SEOF, SERR };	// SEOR end of frecord (success), SEOL end of line, SEOF end of file, SERR error during read...
	enum FFExists { FENE, FEFILE, FEDIR };		// ffexists2() FENE Does Not Exist, FEFILE - exists and is a file, FEDIR - exists and is a directory
	inline int		getChar(FILE *f) { return(getc_unlocked(f)); }
	template< typename Tint >	inline FFStatus getUInt(FILE *f, Tint *c);

	// use fopen / fclose to open or close
	// fread / fwrite to get or dump large blocks aof bytes
	inline FILE *	ffopen( const char *fName, bool read = true )	{ return( std::fopen(fName,(read?"rb":"wb")) ); } // binary mode file access, read is default, to write overwrite second arg default w/ "false"
	inline void		ffclose( FILE **f )								{ if (!f) return; if (! *f) return; fclose(*f); *f=NULL; return; }
	inline bool		ffexists(const char *fName)						{ FILE *f = ffopen(fName); bool ok = (f != NULL);  ffclose(&f); return(ok);  }
	inline FFExists ffexists2(const char *fName)					{ // windows and linux...
		struct stat info;
		if (stat(fName, &info) != 0) return( FFExists::FENE );
		if (info.st_mode & S_IFDIR)  return( FFExists::FEDIR );	
		return( FFExists::FEFILE );
	}
	inline int		ffmv(const char *fsrc, const char *fdst)		{ return(rename(fsrc, fdst)); }
	inline int		ffrmf(const char *fsrc )						{ return(remove(fsrc)); }			// remove a FILE, 0 indicates success
	inline int		ffrmd(const char *dName)						{									// remove EMPTY dir
#if defined( _WIN32 ) || defined( _WIN64 )
		int ret = _rmdir(dName);
#else
		int ret = rmdir(dName);
#endif
		return ret;
	}

	inline int		ffmkdir(const char *dName) {	// 0 is success, anything else is a problem... will not make entire path, just last entry... Windows will use / or \ as delimiter
		// fails if directory already exists, so failure code could stand some refinement...
#if defined( _WIN32 ) || defined( _WIN64 )
		int ret= _mkdir( dName );
#else
		int ret= mkdir(dName, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
#endif
		return ret;
	}

	inline int		putChar(FILE *f, int c) {
#if defined( _WIN32 ) || defined( _WIN64 )
		return( _fputc_nolock( c,f) );
#else
		return(fputc_unlocked(c,f));
#endif
	}

	inline int		fflush(FILE *f = NULL) {
#if defined( _WIN32 ) || defined( _WIN64 )
		if (f != NULL) { std::fflush(f); FlushFileBuffers((void *)_get_osfhandle(fileno(f))); }
#else
		if (f != NULL) { std::fflush(f); fsync(fileno(f)); }
		sync();
#endif
		return 0;
	}

	inline int		putString(FILE *f, const char *p) { int s = 0; while ((*p != '\0') && (s != EOF)) s = putChar(f, *(p++)); return s; }
	inline int		putString(FILE *f, const char *p, uint64_t n) { int s = 0; while ((*p != '\0') && (n>0) && (s != EOF)) { s = putChar(f, *(p++)); --n; }; return s; }

	template< typename Tint >
	inline FFStatus putUInt(FILE *f, Tint i ) {
		if (putString(f, std::to_string(i).c_str()) == EOF) return (SERR);
		return SEOF;
	}

	// TODO: someday we may need to inpliment the SIGNED type too, in a different template...
	// TODO: check for overrun.....
	template< typename Tint >
	inline FFStatus getUInt(FILE *f, Tint *c) {
		int i; // *v=(ulong)0; can pass in a value being built, be careful with this!
		while (1) {
			i = getChar(f);
			switch (i) {
			case EOF:	return(SEOF);	break;
			case '\t':	return(SEOR);	break;
			case '\n':  return(SEOL);	break;
			default:
				i -= (int) '0'; if ((i<0) || (i>9)) return(SERR);
				// multiply by 10, then add i.
				*c = ((*c << 3) + (*c << 1) + (Tint)(i));
				break;
			}
		}
	}

	// access methods...
	inline FFStatus	termChar( int c ) { switch (c)					{ case '\t': return(SEOR); case '\n': return(SEOL); case EOF: return(SEOF); default: return (SERR);} }
	inline FFStatus	putInt(FILE *f, uint8_t i)						{ return(putUInt< uint8_t >(f, i)); }
	inline FFStatus	putInt(FILE *f, uint16_t i)						{ return(putUInt< uint16_t >(f, i)); }
	inline FFStatus	putInt(FILE *f, uint32_t i )					{ return(putUInt< uint32_t >( f, i )); }
	inline FFStatus	putInt(FILE *f, uint64_t i )					{ return(putUInt< uint64_t >( f, i )); }
	//inline FFStatus	getInt(FILE *f,  uint64_t *c)					{ return(getUInt< uint64_t >( f, c ) ); }
	inline FFStatus	getInt(FILE *f, uint8_t  *c)					{ return(getUInt< uint8_t >(f, c)); }
	inline FFStatus	getInt(FILE *f, uint16_t *c)					{ return(getUInt< uint16_t >(f, c)); }
	inline FFStatus	getInt(FILE *f, uint32_t *c)					{ return(getUInt< uint32_t >(f, c)); }
	inline FFStatus	getInt(FILE *f, uint64_t *c)					{ return(getUInt< uint64_t >(f, c)); }
#if defined( _WIN32 ) || defined( _WIN64 )
	// ugh needed to appease both g++ and VisC++2015 compilers...
	inline FFStatus	getInt( FILE *f, ulong    *c )					{ return(getUInt< ulong >(f, c) ); }
#endif
	//inline FFStatus	getInt( FILE *f, ullong   *c )					{ return( getUInt< ullong  >( f, c ) ); }
	inline FFStatus	getInt( FILE *f, llong    *c )					{ return( getUInt< llong   >( f, c ) ); }	// Caution: signed return type, but reads unsigned values only
	inline FFStatus getDouble(FILE *f, double *d);
	// FFStatus		getString( FILE *f, char *p );
	// FFStatus		clearLine( FILE *f );

	//template<typename T> 			FFStatus FileToVec( const char *fName, std::vector<T> &v );
	//template< typename Tuint > 	inline FFStatus getUInt( FILE *f, Tint *c ) {

	// TODO: someday we may need to inpliment the SIGNED type too, in a different template...
	// std::to_string could be faster.... sprintf()? ltoa? something else?

	// readan unsigned integer into the provided integral type,


	// unsafe, does not limit string size to buffer size, DOES insure null termination...
	inline FFStatus getString( FILE *f, char *p ) {
		int i;
		do {
			*p='\0'; i=getChar( f );
			switch (i) {
				case EOF:	return( SEOF);		break;
				case '\t':	return( SEOR);		break;
				case '\n':  return( SEOL);		break;
				default:	*p=(char )i; p++;	break;
			}
		} while (1);
	}


	// TODO this can be optimized a lot, and error checkign provided.... see
	// http://tinodidriksen.com/uploads/code/cpp/speed-string-to-double.cpp or
	// http://www.leapsecond.com/tools/fast_atof.c
	inline FFStatus getDouble(FILE *f, double *d) {
		static char buf[256];
		buf[0] = '\0';
		FFStatus s = getString(f, buf);
		if (s == SEOF) return(s);
		if (s == SERR) return(s);
		if (d) *d = (double)atof(buf);
		return(s);
	}


	// clear remainder of line from stream
	inline FFStatus clearLine( FILE *f) {
		int i;
		do {
			i=getChar( f );
			if (i==EOF) return( SEOF);
		} while (i!='\n');
		return( SEOL );
	}

	#include <vector>
	// read a collection of objects, each of which must support a .readln(FILE *F) method, into a vector
	template<typename T> 
	FFStatus FileToVec( const char *fName, std::vector<T> &v ) {
		// read a set of records from a file
		FFStatus	s=SEOF;
		T			item;
		FILE		*f=ffopen(fName);
		if (f==NULL) return SERR;
		do {
			s=item.readln(f);
			if ( s==SEOL ) v.push_back(item);
		} while ( (s!=SEOF) && (s!=SERR) );
		ffclose(&f); return( s );
	}

	void sleepMillis(uint64_t millis) {
		std::this_thread::sleep_for(std::chrono::milliseconds(millis));
	}

}

