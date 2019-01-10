/*
 * Common.cpp
 *
 *  Created on: 6.04.2017
 *      Author: Alexander Epifanov
 */
#include "Common.h"
//#include <regex.h>
#include <regex>
#include <iostream>
#include <fstream>
#include "real_type.h"

using namespace std;

CConfigFileParser::CConfigFileParser(std::string filename)
{
	m_filename = filename;
}

CConfigFileParser::~CConfigFileParser()
{

}

std::string CConfigFileParser::getParameter(const char *pname_regexp) {
    ifstream in_stream;
    std::string line;
    std::string res;
	//regex_t rex;
	
	std::string pattern = std::string(pname_regexp);
	std::regex re(pname_regexp);
	std::smatch matches;

    in_stream.open(m_filename);
    if (!in_stream)
    	cout << "ERROR" << endl;

    //regcomp(&rex, pname_regexp, REG_EXTENDED);
    //regmatch_t *matches = new regmatch_t[rex.re_nsub + 1];
    while(in_stream >> line)
    {
//        int rc = 0;
//		if ((rc = regexec(&rex, line.c_str(), rex.re_nsub + 1, matches, 0)) < 1)
//		{
//			int pos = matches[1].rm_so;
//			int len = matches[1].rm_eo - matches[1].rm_so;
//			cout << pos << " " << len << endl;
//			res =std::string(line.substr(pos, len));
//			break;
//		}
    	if (std::regex_search(line, matches, re) && matches.size() > 1) {
    	    res = matches.str(1);
    	    break;
    	}
    }
    in_stream.close();
    //delete (matches);

    return res;
}

void dump2binfile(const real_t *buf, const int count, const char *filename) {
    FILE *file = fopen(filename, "wb");
    fwrite(buf, sizeof(*buf), count, file);
    fflush(file);
    fclose(file);
}

void binfile2data(real_t *buf, const int count, const char *filename) {
    FILE *file = fopen(filename, "rb");
    fread(buf, sizeof(*buf), count, file );
    fclose(file);
}

void printMatrix2File(const char *filename, const real_t *m, const real_t *f, const int size)
{
	FILE *fp;
	fp = fopen(filename, "w");
	if (fp)
	{
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++)
			{
				fprintf(fp, "%f\t", m[j * size + i]);

			}
			if (f != NULL)
				fprintf(fp, "%f\n", f[i]);
		}
		fclose(fp);
	}
};

CGaussRule::CGaussRule(const int order, const MeshGeometryType geom_type) {
	this->m_order = order;
	this->m_wi = nullptr;
	this->m_p = nullptr;

	switch(geom_type) {
		case MeshGeometryType::G2D:
			this->m_dim = 2;
			switch(order) {
				case 1:
					this->m_intpoints = 1;
					m_wi = new real_t[m_intpoints];
					m_p = new vector<real_t>[m_intpoints];
					m_wi[0] = 1.0;
					m_p[0] = {1.0/3.0, 1.0/3.0, 1.0/3.0};
					break;
				case 2:
					this->m_intpoints = 3;
					m_wi = new real_t[m_intpoints];
					m_p = new vector<real_t>[m_intpoints];
					m_wi[0] = 1.0/3.0;
					m_wi[1] = 1.0/3.0;
					m_wi[2] = 1.0/3.0;
					m_p[0] = {2.0/3.0, 1.0/6.0, 1.0/6.0};
					m_p[1] = {1.0/6.0, 1.0/6.0, 2.0/3.0};
					m_p[2] = {1.0/6.0, 2.0/3.0, 1.0/6.0};
					break;
				case 3:
					this->m_intpoints = 4;
					m_wi = new real_t[m_intpoints];
					m_p = new vector<real_t>[m_intpoints];
					m_wi[0] = -0.5625;
					m_wi[1] = 0.5208333333333;
					m_wi[2] = 0.5208333333333;
					m_wi[3] = 0.5208333333333;
					m_p[0] = {1.0/3.0, 1.0/3.0, 1.0/3.0};
					m_p[1] = {0.6, 0.2, 0.2};
					m_p[2] = {0.2, 0.2, 0.6};
					m_p[3] = {0.2, 0.6, 0.2};
					break;
				default:
					break;
			}
			break;
		case MeshGeometryType::G3D:
			this->m_dim = 3;
			switch(order) {
				case 1:
					this->m_intpoints = 1;
					m_wi = new real_t[m_intpoints];
					m_p = new vector<real_t>[m_intpoints];
					m_wi[0] = 1.0;
					m_p[0] = {1.0/4.0, 1.0/4.0, 1.0/4.0, 1.0/4.0};
					break;
				case 2:
					this->m_intpoints = 4;
					m_wi = new real_t[m_intpoints];
					m_p = new vector<real_t>[m_intpoints];
					m_wi[0] = 0.25;
					m_wi[1] = 0.25;
					m_wi[2] = 0.25;
					m_wi[3] = 0.25;
					m_p[0] = {0.585410196624969,
							  0.138196601125011,
							  0.138196601125011,
							  0.138196601125011};
					m_p[1] = {0.138196601125011,
							  0.138196601125011,
							  0.138196601125011,
							  0.585410196624969};
					m_p[2] = {0.138196601125011,
							  0.138196601125011,
							  0.585410196624969,
							  0.138196601125011};
					m_p[3] = {0.138196601125011,
							  0.585410196624969,
							  0.138196601125011,
							  0.138196601125011};
					break;
				case 3:
					this->m_intpoints = 5;
					m_wi = new real_t[m_intpoints];
					m_p = new vector<real_t>[m_intpoints];
					m_wi[0] = -0.8;
					m_wi[1] = 0.45;
					m_wi[2] = 0.45;
					m_wi[3] = 0.45;
					m_wi[4] = 0.45;
					m_p[0] = {0.25, 0.25, 0.25, 0.25};
					m_p[1] = {1.0/2.0, 1.0/6.0, 1.0/6.0, 1.0/6.0};
					m_p[2] = {1.0/6.0, 1.0/6.0, 1.0/6.0, 1.0/2.0};
					m_p[3] = {1.0/6.0, 1.0/6.0, 1.0/2.0, 1.0/6.0};
					m_p[4] = {1.0/6.0, 1.0/2.0, 1.0/6.0, 1.0/6.0};
					break;
				default:
					break;
			}
			break;


		default:
			break;
	}
}

CGaussRule::~CGaussRule() {
//	if (m_wi != nullptr)
//		delete m_wi;
//	if (m_p != nullptr)
//		delete m_wi;
}

real_t CGaussRule::integrate(func_t2* func, const int idxN, const int l_row, const int element, const real_t elemSize) {
	real_t res = 0;

//	for (int i =0; i < this->m_intpoints; ++i) {
//		res += m_wi[i] *  func(idxN, l_row, element, m_p[i]);
//	}

	return res * elemSize;
}

real_t CGaussRule::integrate(func_t3 *func, const int idxN, const int jdxN, const int l_col, const int l_row, const int element, const real_t elemSize) {
	real_t res = 0;

//	for (int i =0; i < this->m_intpoints; ++i) {
//		res += m_wi[i] *  func(idxN, jdxN, l_col, l_row, element, m_p[i]);
//	}

	return res * elemSize;
}

