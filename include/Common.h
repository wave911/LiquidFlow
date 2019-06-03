#ifndef COMMON_H_
#define COMMON_H_

#include <ctype.h>
#include <string>
#include "real_type.h"
#include "Constants.h"

#define CFG_PLENGTH 256
#ifdef APPLE
#define CONFIG_FILENAME "/Users/epifanov/Projects/LiquidFlow/bin/cfg.txt"
#define K_MATRIX_FILENAME "/Users/epifanov/Projects/LiquidFlow/bin/kmatrix.bin"
#endif
#ifndef APPLE
#define CONFIG_FILENAME "cfg.txt"
#define K_MATRIX_FILENAME "kmatrix.bin"
#endif

class CConfigFileParser
{
	private:
		std::string m_filename;
	public:
		CConfigFileParser(std::string filename);
		virtual ~CConfigFileParser();

		std::string getParameter(const char *pname);
		
};

void dump2binfile(const real_t *buf, const int count, const char *filename);
void binfile2data(real_t *buf, const int count, const char *filename);
void printMatrix2File(const char *filename, const real_t *m, const real_t *f, const int size);
inline real_t fact(int n) {return (n == 1 || n == 0) ? 1 : fact(n - 1) * n;}

class CGaussRule
{
	typedef real_t (*func_t2)(const int idxN, const int l_row, const int element, const std::vector<real_t> ksi);
	typedef real_t (*func_t3)(const int idxN, const int l_col, const int l_row, const int jdxN, const int element, const std::vector<real_t> ksi);
	private:
		int m_dim,
			m_order;

	public:
		int m_intpoints;
		std::vector<real_t>* m_p;
		real_t* m_wi;
		CGaussRule(const int order, const MeshGeometryType geom_type);
		virtual ~CGaussRule();
		real_t integrate(func_t2 *func, const int idxN, const int l_row, const int element, const real_t elemSize);
		real_t integrate(func_t3 *func, const int idxN, const int jdxN, const int l_col, const int l_row, const int element, const real_t elemSize);
};

#endif //COMMON_H_
