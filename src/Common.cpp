/*
 * Common.cpp
 *
 *  Created on: 6.04.2017
 *      Author: Alexander Epifanov
 */
#include "Common.h"
#include <regex.h>
#include <iostream>
#include <fstream>

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
	regex_t rex;
	
    in_stream.open(m_filename);

    regcomp(&rex, pname_regexp, REG_EXTENDED);
    regmatch_t *matches = new regmatch_t[rex.re_nsub + 1];
    while(in_stream >> line)
    {
        int rc = 0;
		if ((rc = regexec(&rex, line.c_str(), rex.re_nsub + 1, matches, 0)) < 1)
		{
			int pos = matches[1].rm_so;
			int len = matches[1].rm_eo - matches[1].rm_so;
			res =std::string(line.substr(pos, len));
			break;

		}			
    }
    in_stream.close();
    delete (matches);

    return res;
}