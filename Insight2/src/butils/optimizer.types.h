#pragma once

typedef std::vector<double> vecD;
struct limit { double *vmin, *vmax; };
typedef std::vector<limit>  vecL;

// in the lbfgs namespace
struct options {	// USER options for  lbfgs algorithm...
	int verbose;	// reporting to stderr....
	int m;			// Number of nubmer of bfgs step. Author suggests [3,20]. Yeifi used 30. Memory O(m^2).
	double factr;	// termination condition, relative to machine precision, low 1.0e12, med 1.0e7, v high 1.0e1 . 0.0 -> stop only when iter produces no improvement.
	double pgtol;	// termination condition, terminate on low gradient, 1e-12 typical, 0 to disable.
	options() { clear(); };
	void clear() { verbose = 0; m = 20; factr = 1.0e2;  pgtol = 0.0; };
};

struct wrkBuffer {	// hide and manage details of the WORK buffer used by the algorithm...
public:
	double *parameters, *boundl, *boundu, *gradient, *wa, f, dsave[29];
	int *boundf, *wai, lsave[4], isave[44], n, m;
	char task[60], csave[60];
	wrkBuffer() { clear(true); };
	~wrkBuffer() { clear(false); };
	bool isTask(const char *c) { return (sstrcmp(task, c) == 0); };
	void setTask(const char *c) { strncpy(task, c, strlen(c)); };
	void clear(void) { clear(false); }
	bool init(const vecD &Params, const vecL &Limits, int M) {
		n = (int) Params.size(); m = M;	// loss of precision
		if (!alloc(n, m)) return(false);
		for (int i = 0; i < n; i++) {
			parameters[i] = Params[i]; gradient[i] = 0.0;
			boundf[i] = 0; boundl[i] = boundu[i] = 0.0;
			if (Limits[i].vmin != NULL) { boundf[i] += 1; boundl[i] = *(Limits[i].vmin); };
			if (Limits[i].vmax != NULL) { boundf[i] += 2; boundu[i] = *(Limits[i].vmax); };
			if (boundf[i] >= 2) boundf[i] = (boundf[i] == 2 ? 3 : 2);	// sigh..... author didn't use bit flags
		}
		return true;
	}
private:
	static int sstrcmp(const char *a, const char *b) { return (strncmp(a, b, strlen(b))); };
	void ffree(void **ptr) { if (!ptr) return; if (!*ptr) return; free(*ptr); *ptr = NULL; }
	bool alloc(int SizeN, int SizeM) {
		bool ok = true;
		int dim = (2 * SizeM + 4)*SizeN + 12 * SizeM* SizeM + 12 * SizeM;
		clear();
		parameters	= (double *)calloc(SizeN, sizeof(*parameters));		ok = ok && (parameters != NULL);
		boundl		= (double *)calloc(SizeN, sizeof(*boundl));			ok = ok && (boundl != NULL);
		boundu		= (double *)calloc(SizeN, sizeof(*boundu));			ok = ok && (boundu != NULL);
		gradient	= (double *)calloc(SizeN, sizeof(*gradient));		ok = ok && (gradient != NULL);
		wa			= (double *)calloc(dim, sizeof(*wa));				ok = ok && (wa != NULL);
		wai			= (int *)calloc(3 * SizeN, sizeof(*wai));			ok = ok && (wai != NULL);
		boundf		= (int *)calloc(SizeN, sizeof(*boundf));			ok = ok && (boundf != NULL);
		n = SizeN; m = SizeM; f = 0.0;
		if (!ok) clear();
		return(ok);
	}
	void clear(bool first) {	// call with first=true only to initialize pointers to NULL during construction
		if (first) {
			parameters = boundl = boundu = gradient = wa = NULL;
			boundf = wai = NULL;
		};
		ffree((void **)&parameters); 
		ffree((void **)&boundl); 
		ffree((void **)&boundu); 
		ffree((void **)&gradient); 
		ffree((void **)&wa);
		ffree((void **)&boundf); 
		ffree((void **)&wai);
		//memset(dsave, 0, 29 * sizeof(double)); 
		//memset(lsave, 0, 4 * sizeof(int)); 
		//memset(isave, 0, 44 * sizeof(int));
		//memset(task, 0, 60 * sizeof(char)); 
		//memset(csave, 0, 4 * sizeof(char));
		f = 0.0; n = m = 0;
	};
};