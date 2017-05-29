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
	m_borderElementsNumber = 0;
}

CSalomeMesh::~CSalomeMesh() {
	m_points.clear();
	m_borderPoints.clear();
	m_borderElements.clear();
	m_mesh.clear();
}

int CSalomeMesh::Init(MeshGeometryType meshType) {
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
    		CPoint3D point;
    		if (meshType == MeshGeometryType::G2D) {
    			point = CPoint3D(stof(tokens[1]), stof(tokens[2]), 0);
    		}
    		if (meshType == MeshGeometryType::G3D) {
    			point = CPoint3D(stof(tokens[1]), stof(tokens[2]), stof(tokens[3]));
    		}
    		m_points.push_back(point);
    	}    	
    	if (idx > m_pointsNumber ) {
    		if (meshType == MeshGeometryType::G3D) {
		    	if (tokens[1].find("220", 0) != std::string::npos) {
		    		addPoints(tokens, m_borderElements);
		    	}
		    	if (tokens[1].find("30",0) != std::string::npos) {
		    		addPoints(tokens, m_mesh);
		    	}    	
	    	}
    		if (meshType == MeshGeometryType::G2D) {
		    	if (tokens[1].find("10", 0) != std::string::npos) {
		    		addPoints(tokens, m_borderElements);
		    	}
		    	if (tokens[1].find("20",0) != std::string::npos) {
		    		addPoints(tokens, m_mesh);
		    	}    	
	    	}	    	
    	}
    }
    createBorderPoints();

    cout << m_mesh.size() << endl;
    cout << m_borderElements.size() << endl;
    cout << m_borderPoints.size() << endl;
    cout << m_points.size() << endl;
    in_stream.close();
}

int CSalomeMesh::getPointsNumber() {
	return m_points.size();
}

int CSalomeMesh::getElementsNumber() {
	return m_mesh.size();
}

int CSalomeMesh::getBorderElementsNumber() {
	return m_borderElements.size();
}

std::vector<CPoint3D> CSalomeMesh::getPoints() {
	return m_points;
}

std::map<int,std::vector<int>> CSalomeMesh::getElements() {
	return m_mesh;
}

std::vector<int> CSalomeMesh::getElementByIndex(const int idx) {
	return m_mesh[idx];
}

CPoint3D CSalomeMesh::getPointByIndex(const int idx) {
	return m_points[idx];
}

bool CSalomeMesh::isBorderPoint(const int idx) {
	if (m_borderPoints.count(idx) != 0)
		return true;
	else
		return false;
}

std::set<int> CSalomeMesh::getBorderPoints() {
	return m_borderPoints;
}

void CSalomeMesh::addPoints(std::vector<std::string>& tokens, std::map<int,std::vector<int>>& aMap) {
	vector<int> temp;
	for (int i = 2; i < tokens.size(); i++) { 
		if (trim(tokens[i]).length() != 0) {
			temp.push_back(stoi(tokens[i]));
		}
	}
	int index = aMap.size();
	aMap[index] = temp;
	temp.clear();
}

void CSalomeMesh::createBorderPoints() {
	for (auto const& x : m_borderElements) {
		vector<int> pts = x.second;
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


