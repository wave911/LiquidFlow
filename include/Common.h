#ifndef COMMON_H_
#define COMMON_H_

#include <ctype.h>
#include <string>
#include "real_type.h"
#include "Constants.h"

#define CFG_PLENGTH 256
#define CONFIG_FILENAME "cfg.txt"
#define K_MATRIX_FILENAME "kmatrix.bin"

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

#endif //COMMON_H_
