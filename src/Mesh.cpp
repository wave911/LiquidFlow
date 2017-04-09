#include "Mesh.h"
#include <iostream>
#include <fstream>
#include <cctype>
#include <algorithm>
#include "ElementType.h"

using namespace std;

CSalomeMesh::CSalomeMesh() {
}

CSalomeMesh::CSalomeMesh(const std::string filename) : CMesh() {
	m_filename = filename;
	m_pointsNumber = 0;
	m_elementsNumber = 0;
	m_trianglesNumber = 0;
}

CSalomeMesh::~CSalomeMesh() {
	m_points.clear();
	m_borderPoints.clear();
	m_triangles.clear();
	m_mesh.clear();
}

int CSalomeMesh::Init() {
	ifstream in_stream;
	string line;
    in_stream.open(m_filename);

    std::getline(in_stream, line);
    std::vector<std::string> tokens = split(line, " ");
    m_pointsNumber = stoi(tokens[0]);

    int idx = 0;
    while(std::getline(in_stream, line))
    {
    	std::vector<std::string> tokens = split(line, " ");
    	if (idx++ < m_pointsNumber) {
    		CPoint3D point = CPoint3D(stof(tokens[1]), stof(tokens[2]), stof(tokens[3]));
    		m_points.push_back(point);
    	}
    	if (idx >= m_pointsNumber ) {
	    	if (tokens[1].find_first_of("220", 0) == 0) {
	    		addPoints(tokens, m_triangles);
	    	}
	    	if (tokens[1].find_first_of("30", 0) == 0) {
	    		addPoints(tokens, m_mesh);
	    	}    	
    	}
    }
    getBorderPoints();

    cout << m_mesh.size() << endl;
    cout << m_triangles.size() << endl;
    cout << m_borderPoints.size() << endl;
    in_stream.close();
}

int CSalomeMesh::getPointsNumber() {
	return m_points.size();
}

int CSalomeMesh::getElementsNumber() {
	return m_mesh.size();
}

int CSalomeMesh::getTrianglesNumber() {
	return m_triangles.size();
}

std::list<CPoint3D> CSalomeMesh::getPoints() {
	return m_points;
}

std::map<int,std::list<int>> CSalomeMesh::getElements() {
	return m_mesh;
}

void CSalomeMesh::addPoints(std::vector<std::string>& tokens, std::map<int,std::list<int>>& aMap) {
	list<int> temp;
	for (int i = 2; i < tokens.size(); i++) { 
		if (trim(tokens[i]).length() != 0) {
			temp.push_back(stoi(tokens[i]));
		}
	}
	int index = stoi(tokens[0]);
	aMap[index] = temp;
	temp.clear();
}

void CSalomeMesh::getBorderPoints() {
	for (auto const& x : m_triangles) {
		list<int> pts = x.second;
		for (auto const& p: pts) {
			m_borderPoints.insert(p);
		}
	}
}

std::vector<std::string> CSalomeMesh::split(const std::string& text, const std::string& delims) {
    std::vector<std::string> tokens;
    std::size_t start = text.find_first_not_of(delims), end = 0;

    while((end = text.find_first_of(delims, start)) != std::string::npos)
    {
        tokens.push_back(text.substr(start, end - start));
        start = text.find_first_not_of(delims, end);
    }
    if(start != std::string::npos)
        tokens.push_back(text.substr(start));

    return tokens;
}

std::string CSalomeMesh::trim(const std::string &s)
{
   auto wsfront=std::find_if_not(s.begin(),s.end(),[](int c){return std::isspace(c);});
   auto wsback=std::find_if_not(s.rbegin(),s.rend(),[](int c){return std::isspace(c);}).base();
   return (wsback<=wsfront ? std::string() : std::string(wsfront,wsback));
}


