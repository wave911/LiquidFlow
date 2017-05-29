#ifndef COMMON_H_
#define COMMON_H_

#include <ctype.h>
#include <string>
#include "real_type.h"

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

enum class MeshGeometryType {G2D = 0, G3D = 1};

#endif //COMMON_H_
