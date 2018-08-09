#ifndef MESH_H_
#define MESH_H_

#include "real_type.h"
#include "Point.h"
#include "Common.h"
#include <string.h>
#include <vector>
#include <map>
#include <vector>
#include <set>

class CMesh
{
	public:
		CMesh() {};
		virtual ~CMesh() {};

		virtual int Init(MeshGeometryType meshType) = 0;
		virtual int getPointsNumber() = 0;
		virtual int getElementsNumber() = 0;
		virtual int getBorderElementsNumber() = 0;
		virtual std::vector<CPoint3D> getPoints() = 0;
		virtual std::set<int> getBorderPoints() = 0;
		virtual std::map<int,std::vector<int> > getElements() = 0;
		virtual std::vector<int> getElementByIndex(const int idx) = 0;
		virtual CPoint3D getPointByIndex(const int idx) = 0;
		virtual bool isBorderPoint(const int idx) = 0;
		virtual int getPointsNumberPerElement() = 0;
};

class CSalomeMesh : public CMesh 
{
	private:
		std::vector<CPoint3D> m_points;
		std::set<int> m_borderPoints;
		std::map<int,std::vector<int> > m_borderElements;
		std::map<int,std::vector<int> > m_mesh;
		int m_pointsNumber;
		int m_elementsNumber;
		int m_borderElementsNumber;

		std::string m_filename;
	public:
		CSalomeMesh();
		CSalomeMesh(const std::string filename);
		virtual ~CSalomeMesh();

		virtual int Init(MeshGeometryType meshType);
		virtual int getPointsNumber();
		virtual int getElementsNumber();
		virtual int getBorderElementsNumber();
		virtual std::vector<CPoint3D> getPoints();
		virtual std::set<int> getBorderPoints();
		virtual std::map<int,std::vector<int> > getElements();
		virtual std::vector<int> getElementByIndex(const int idx);
		virtual CPoint3D getPointByIndex(const int idx);
		virtual bool isBorderPoint(const int idx);
		virtual int getPointsNumberPerElement();

	private:
		std::vector<std::string> split(const std::string& text, const std::string& delims);
		void addPoints(std::vector<std::string>& tokens, std::map<int,std::vector<int> >& aMap);
		std::string trim(const std::string &s);
		void createBorderPoints();

};

class CFreeFemMesh: public CMesh {
private:
		std::vector<CPoint3D> m_points;
		std::set<int> m_borderPoints;
		std::map<int,std::vector<int> > m_borderElements;
		std::map<int,std::vector<int> > m_mesh;
		int m_pointsNumber;
		int m_elementsNumber;
		int m_borderElementsNumber;

		std::string m_filename;
	public:
		CFreeFemMesh();
		CFreeFemMesh(const std::string filename);
		virtual ~CFreeFemMesh();

		virtual int Init(MeshGeometryType meshType);
		virtual int getPointsNumber();
		virtual int getElementsNumber();
		virtual int getBorderElementsNumber();
		virtual std::vector<CPoint3D> getPoints();
		virtual std::set<int> getBorderPoints();
		virtual std::map<int,std::vector<int> > getElements();
		virtual std::vector<int> getElementByIndex(const int idx);
		virtual CPoint3D getPointByIndex(const int idx);
		virtual bool isBorderPoint(const int idx);
		virtual int getPointsNumberPerElement();

	private:
		std::vector<std::string> split(const std::string& text, const std::string& delims);
		void addPoints(std::vector<std::string>& tokens, std::map<int,std::vector<int> >& aMap);
		std::string trim(const std::string &s);
		void createBorderPoints();
};
#endif /* MESH_H_ */
