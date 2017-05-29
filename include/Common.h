#ifndef COMMON_H_
#define COMMON_H_

//#include "Mesh.h"
//#include "Point.h"
#include <ctype.h>
#include <string>
//#include <list>
//#include <map>

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

enum class MeshGeometryType {G2D = 0, G3D = 1};

#endif //COMMON_H_
