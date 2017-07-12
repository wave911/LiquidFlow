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

