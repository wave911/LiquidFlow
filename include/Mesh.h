#ifndef MESH_H_
#define MESH_H_

#include "common_headers/real_type.h"
#include "common_headers/all_includes.h"
#include "Point.h"
#include <string.h>

class CExtGenMesh
{
	protected:
		CPoint3D			*m_point;
		CPointProperties	*m_ptprop;	
		int		m_elnumber,
				m_ptnumber,
				m_meshgen;
		int*	m_elptnum;
		char*   m_filename;

	private:
		int m_order;
		int m_trianglenum;
		int *p_trian_pointnum;

	public:
		CExtGenMesh();
		CExtGenMesh(const char *filename, const int meshgen, int order);
		virtual ~CExtGenMesh();
		virtual int SalomeDat2DParser();
		virtual int SalomeDat3DParser();
		
		CPoint3D*	GetPoints();
		CPointProperties* GetPtProp();
		int*		GetElpoints();
		int*		GetTrianglePoints();
		int			GetElnumber();
		int			GetPtnumber();
		int			GetTriangleNumber();	
		int			GetPtsPerElement();			
		int			PrintElements(const char *filename);

	protected:
		virtual void SetPtNumber(int ptnumber);
		virtual void SetElPtNumber(int elnumber);
		virtual void SetElPoints(int *elpoints);
		virtual void SetTriangleNumber(int trianglenumber);
		virtual void SetTrianglePoints(int *trianglepoints);
		virtual void SetPoints(CPoint3D *pt, CPointProperties *prp);
};

// int Init(CMeshBaseD *mb);
// int MakeNet_(CMeshBaseD *mb);
// int InitPoints(CMeshBaseD *mb);
// CPipeMeshSelfGen3D* new_CPipeMeshSelfGen3D(real_t R, real_t H, int n, int nz, int schemeorder);
// CCircleMeshSelfGen2D* new_CCircleMeshSelfGen2D(const real_t R, const int n, const int schemeorder);
// CSquareTriangleMesh2D* new_CSquareTriangleMesh2D(const real_t W, const real_t H, const int wn, const int hn);
// CExtGenMesh2D * new_CExtGenMesh2D (char *filename, const int filenamelen, const int meshgen);
// CMeshLine1D* new_CMeshLine1D(const real_t sL, const real_t eL, const real_t h);

#endif /* MESH_H_ */