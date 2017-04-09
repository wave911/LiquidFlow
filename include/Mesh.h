#ifndef MESH_H_
#define MESH_H_

#include "real_type.h"
#include "Point.h"
#include <string.h>
#include <list>
#include <map>
#include <vector>
#include <set>

class CMesh
{
	public:
		CMesh() {};
		virtual ~CMesh() {};

		virtual int Init() = 0;
		virtual int getPointsNumber() = 0;
		virtual int getElementsNumber() = 0;
		virtual int getTrianglesNumber() = 0;
		virtual std::list<CPoint3D> getPoints() = 0;
		virtual std::map<int,std::list<int>> getElements() = 0;
};

class CSalomeMesh : public CMesh {
	private:
		std::list<CPoint3D> m_points;
		std::set<int> m_borderPoints;
		std::map<int,std::list<int>> m_triangles;
		std::map<int,std::list<int>> m_mesh;
		int m_pointsNumber;
		int m_elementsNumber;
		int m_trianglesNumber;

		std::string m_filename;
	public:
		CSalomeMesh();
		CSalomeMesh(const std::string filename);
		virtual ~CSalomeMesh();

		virtual int Init();
		virtual int getPointsNumber();
		virtual int getElementsNumber();
		virtual int getTrianglesNumber();
		virtual std::list<CPoint3D> getPoints();
		virtual std::map<int,std::list<int>> getElements();

	private:
		std::vector<std::string> split(const std::string& text, const std::string& delims);
		void addPoints(std::vector<std::string>& tokens, std::map<int,std::list<int>>& aMap);
		std::string trim(const std::string &s);
		void getBorderPoints();

};

#endif /* MESH_H_ */