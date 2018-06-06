#pragma once

#include <vector>
#include <map>
#include <fstream>
#include <iostream>	//cerr for error reporting, low performance
// This class facilitates the replacement of verbose data types with terse ones.
//	while "double" floats might be needed (8 bytes), if only a small number (say 255) of unique values
//	need to be represented, individual vlaues can can be replaced with an index into a table, saving
//	83% of the memory footprint. If the long value is a string, the savings can be far greater.
//
// Generally Tind is unsigned char (256 vals) or unsigned short (65536 vals), often index 0 is used to represent
//		a "missing" value.
// TVal is usually a double (8 bytes), long long int (8 bytes) or fixed width string.... (char [])

template< class Tval, class Tind>	// Tval is a value to be stored usually string or double, Tind is the type of index, usually uchar or ushort.
class RTable
{
public:
	typedef Tind tIndex;
	//typedef std::cerr cerr;
	//typedef std::endl endl;
	RTable( ) { clear( true ); };
	virtual ~RTable( ) { };
	//
	inline void set( tIndex ind, Tval val ) { vals.at(ind)=val; if (ind>indexTopV) indexTopV=ind; vlookup[val]=ind; };
	const Tval * operator[] ( Tind I ) { return( &(vals[I]) ); };			// get a value from an index... fast (O(1)) 
	const Tval * operator[] (Tind I) const { return(&(vals[I])); };			// get a value from an index... fast (O(1)) 
	inline Tind map( Tval V ) { return( vlookup[V] ); };					// get an index from a value... slow... (lg(n))
	inline unsigned long indexMax( void ) { return( indexMaxV); };			// highest index available to be set
	inline unsigned long indexTop( void ) { return( indexTopV); };			// highest index actually set
	bool readSimple(const char *FName, int V=0) {							// slow, but simple, loading of table OK for small data sets...
		using std::cerr;
		using std::endl;
		Tind i; Tval v; bool ok = false; unsigned long itmp;
		if (V > 0)  cerr << "RTable, opening file " << FName << endl;
		std::fstream fin(FName, std::ios::in);
		clear( );
		while (fin.good()) {
			fin >> itmp; if (!fin.good()) { ok = true; break; };				// normal EOF
			i = (Tind)itmp;
			fin >> v;    if (!fin.good()) { ok = false; break; };				// normal EOF
			set(i, v);
			if (V > 1) cerr << "\t Found " << (unsigned long ) i << " " << v << endl;
		};
		if (V > 1) {
			cerr << "Table Read complete, dumping values max , top =  " << indexMax() << " " << indexTop() << endl;
			for (Tind ind = 0; ind <= indexMax(); ind++) {
				cerr << "\t" << (unsigned long ) ind << " " << *((*this)[ind]) << " " << (ulong) map(*(*this)[ind]) << endl;
				if (ind == indexMax()) break;
			};
			cerr << "Table dumping complete" << endl;
		}
		return(ok);
	};
	void clear(bool first = false) { indexMaxV = std::numeric_limits<tIndex>::max(); vlookup.clear(); vals.clear(); vals.resize(indexMaxV + 1); indexTopV = 0;}
private:
	unsigned long indexMaxV;
	unsigned long indexTopV;
	std::vector<Tval>	vals;		// index to value (fast)
	std::map<Tval,Tind>	vlookup;	// value to index (slow)
};
