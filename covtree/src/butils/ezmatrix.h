#pragma once

#include <cstdio>
#include <vector>
#include <string>

// low performance matrix class for simple operations

class ezmatrix {
private:
	std::vector<std::vector<double>> m;
	double dot(const std::vector<double> &v1, const std::vector<double> &v2) const {
		if (v1.size() != v2.size()) return(0);
		double cum = 0.0;
		for (unsigned int i = 0; i < v1.size(); i++) cum += v1[i] * v2[i];
		return(cum);
	};
	inline static double mylog(double v) { return (v > .5 ? log1p(v - 1.0) : log(v)); };
public:
	inline int nrow() const { return((int)m.size()); };
	inline int ncol() const { return( (int) ( nrow() < 1 ? 0 : m[0].size()) ); };
	inline void clear() { m.clear(); }
	inline void resize(int Size) {	resize(Size, Size);	}
	inline void resize(int Rows, int Colls) {
		if ((Rows < 1) || (Colls < 1)) { clear(); return; }
		m.resize(Rows); for (int i = 0; i < nrow(); i++) m[i].resize(Colls);
	}

	inline double *vp(int i, int j) { return(&(m[i][j])); };
	inline double &v(int i, int j)  { return(m[i][j]); };
	inline const double &v(int i, int j) const { return(m[i][j]); };
	inline bool isSquare() const { return(nrow() == ncol());	}
	inline bool isEmpty() const { return(nrow()==0 || ncol()==0); }

	void delRow(int I) {
		if (I < 0 || I >= nrow()) return;
		for (int i = (I + 1); i < nrow(); i++) m[i - 1] = m[i];
		m.resize(nrow() - 1);
	};

	void delCol(int J) {
		if (J < 0 || J >= ncol()) return;
		int oldcolls = ncol();
		int newcolls = oldcolls-1;
		for (int i = 0; i < nrow(); i++) {
			for (int j = (J + 1); j < oldcolls; j++)  m[i][j - 1] = m[i][j];
			m[i].resize(newcolls);
		}
	};

	// adjoint, that is, the transpose of the cofactor matrix 
	bool adj(ezmatrix &Res) const {
		double fac = 1.0;
		if (!isSquare()) return(false);
		Res.resize(nrow(), ncol());
		// gn3erate cofactor matrix
		ezmatrix tmp;
		for (int i = 0; i < nrow(); i++) {
			for (int j = 0; j < ncol(); j++) {
				tmp = (*this); tmp.delRow(i); tmp.delCol(j);
				Res.v(i,j) = fac * tmp.det();
				fac = -fac;
			}
		}
		Res.t();	// take transpose
		return(true);
	}

	void mult(const double &B) {
		for (int i = 0; i < nrow(); i++)
			for (int j = 0; j < ncol(); j++)
				m[i][j] *= B;
		return;
	}

	void mult(const double &B, ezmatrix &Res) const {
		Res = (*this);
		Res.mult(B);
	};

	bool add(const ezmatrix &B) {
		if (ncol() != B.ncol()) return(false);
		if (nrow() != B.nrow()) return(false);
		for (int i = 0; i < nrow(); i++)
			for (int j = 0; j < ncol(); j++)
				m[i][j] += B.v(i, j);
		return(true);
	};

	bool add(const ezmatrix &B, ezmatrix &Res) const {
		if (ncol() != B.ncol()) return(false);
		if (nrow() != B.nrow()) return(false);
		Res = (*this); Res.add(B);
		return(true);
	};

	bool mult(const ezmatrix &B, ezmatrix &Res) {
		if (ncol() != B.nrow()) return(false);
		Res.resize(nrow(), B.ncol());
		ezmatrix Bt = B;
		Bt.t();
		for (int i = 0; i < nrow(); i++)
			for (int j = 0; j < Res.ncol(); j++)
				Res.v(i, j) = dot(m[i], Bt.m[j]);
		return(true);
	};

	bool makeSym() {
		if (nrow() != ncol()) return( false );
		for (int i = 0; i < nrow(); i++) {
			for (int j = i+1; j < ncol(); j++) {
				double v = ( m[i][j] + m [j][i] ) / 2.0;
				m[i][j] = m[j][i] = v;
			};
		}
		return true;
	};

	bool logElem() {
		bool ok = true;
		for (int i = 0; i < nrow(); i++) {
			for (int j = 0; j < ncol(); j++) {
				double v = m[i][j];
				if (v < 0) ok = false;	// bit of a hack, accept log(0) as -maxval. 
//				m[i][j] = (v==0?-std::numeric_limits<double>::max() : mylog(m[i][j]) );
				m[i][j] = (v==0? 0 : mylog(m[i][j]) );
			};
		}
		return ok;
	};

	// matrix inverse
	bool inv(ezmatrix &Res) const {
		double d = det();
		if (d == 0.0) return(false);
		d = 1.0 / d;
		adj(Res); Res.mult(d);
		return(true);
	};

	// Normalize 1-norm to 1.0, if possible. Particularly useful if matrix contains probabilities...
	double norm1() {
		double tmp = 0.0;
		for (int i = 0; i < nrow(); i++) for (int j = 0; j < ncol(); j++) tmp += fabs(m[i][j]); 
		if (tmp != 0.0) {	
			double tmp2 = 1.0 / tmp;
			for (int i = 0; i < nrow(); i++) for (int j = 0; j < ncol(); j++) m[i][j] *= tmp2;
		}
		return tmp;
	}

	// transpose...
	void t() {
		double tmp = 0;
		for (int i = 0; i < nrow(); i++) {
			for (int j = i+1; j < ncol(); j++) { 
				tmp = m[i][j]; m[i][j] = m[j][i]; m[j][i] = tmp; 
			};
		}
	}

	// determinant, slow but simple. recursive
	double det() const {
		if (!isSquare()) return(0.0); if (nrow() < 1) return(0.0);
		if (nrow() == 1) return(m[0][0]);
		if (nrow() == 2) return(m[0][0] * m[1][1] - m[1][0] * m[0][1]);
		double fac = 1.0, cum=0.0;
		ezmatrix tmp;
		for (int j = 0; j < ncol(); j++) {
			tmp = (*this); tmp.delRow(0);  tmp.delCol(j);
			cum += fac * m[0][j] * tmp.det();
			fac = -fac;
		}
		return(cum);
	}

	std::string toStr(const std::string &fmt="%12.6lf") {
		char buf[1024];
		std::string res;
		for (int i = 0; i < nrow(); i++) {
			for (int j = 0; j < ncol(); j++) {
				if (j != 0) res += "\t";
				sprintf(buf,fmt.c_str(), m[i][j]); res += buf;
			}
			res += "\n";
		}
		return(res);
	}

	ezmatrix(int Rows, int Colls) { resize(Rows, Colls); }
	ezmatrix(int Size) : ezmatrix(Size, Size) {};
	ezmatrix() : ezmatrix(1) {};
	virtual ~ezmatrix() { };
};

