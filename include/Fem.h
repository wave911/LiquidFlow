#ifndef FEM_H_
#define FEM_H_

#include "Mesh.h"
#include "Problem.h"

#ifdef LIBBLASLAPACK
#endif

class CFem 
{
	protected:
		CMesh *m_mesh;
		CProblem *m_pr;
		MeshGeometryType m_mesh_geometry_type;
		real_t *m_K,
			   *m_C,
			   *m_F,
			   *m_U_temp;
		int m_equations_number;
	public:
		CFem(const int equations_number, const MeshGeometryType mgt);
		virtual ~CFem() {};
		virtual void init(CProblem *pr);
		virtual void assembleKMatrix();
		virtual void assembleRightVector(const int timestep);
		virtual void setBorderConditions(const int timestep);
		virtual void perform(const int timesteps);
		virtual real_t getVolume(const int element) = 0;
		virtual std::vector<real_t> getLocalCoordinates(const int element,
													  const CPoint3D p) = 0;

		virtual real_t getN(const int idxN, const std::vector<real_t> ksi) = 0;
		virtual real_t getdNdKsi(const int idxN, const int idxKsi, const std::vector<real_t> ksi) = 0;

		void dgemv(char *trans, int m, int n, real_t alpha, real_t *a,
				   int lda, real_t *x, int incx, real_t beta, real_t *y, int incy);
		void dgesv(int n, real_t *M, real_t *B );
		void dgemm(char *TransA, char *TransB, int m, int n, int k, real_t alpha, real_t *A, int lda,
					real_t *B, int ldb, real_t beta, real_t *C, int ldc);

		virtual real_t getKK(const int idxN, const int jdxN, const int l_col, const int l_row, const int element, const std::vector<real_t> ksi) = 0;
		virtual real_t getCC(const int idxN, const int jdxN, const int l_col, const int l_row, const int element, const std::vector<real_t> ksi) = 0;
		virtual real_t getFF(const int idxN, const int l_row, const int element, const std::vector<real_t> ksi) = 0;

		virtual real_t getdNdX(const int idxN, const int element, const std::vector<real_t> ksi) = 0;
		virtual real_t getdNdY(const int idxN, const int element, const std::vector<real_t> ksi) = 0;
		virtual real_t getdNdZ(const int idxN, const int element) {return 0;};
};

class CFemLocalLinear2D : public CFem
{
	public:
		CFemLocalLinear2D(CMesh *mesh, const int equations_number, const MeshGeometryType mgt);
		virtual ~CFemLocalLinear2D();
		//virtual void init(CProblem *pr);
//		virtual void assembleKMatrix();
//		virtual void assembleRightVector(const int timestep);
//		virtual void setBorderConditions(const int timestep);
//		virtual void perform(const int timesteps);
	protected:
		virtual real_t getVolume(const int element);
		virtual std::vector<real_t> getLocalCoordinates(const int element,
													  const CPoint3D p);
		virtual real_t getN(const int idxN, const std::vector<real_t> ksi);
		virtual real_t getdNdKsi(const int idxN, const int idxKsi, const std::vector<real_t> ksi);
		virtual real_t getKK(const int idxN, const int jdxN, const int l_col, const int l_row, const int element, const std::vector<real_t> ksi);
		virtual real_t getCC(const int idxN, const int jdxN, const int l_col, const int l_row, const int element, const std::vector<real_t> ksi);
		virtual real_t getFF(const int idxN, const int l_row, const int element, const std::vector<real_t> ksi);
		virtual real_t getdNdX(const int idxN, const int element, const std::vector<real_t> ksi);
		virtual real_t getdNdY(const int idxN, const int element, const std::vector<real_t> ksi);
		virtual real_t getdKsidX(const int idx, const int element);
		virtual real_t getdKsidY(const int idx, const int element);
		virtual real_t getdUdX(const int element, const int dim, const std::vector<real_t> ksi);
		virtual real_t getdUdY(const int element, const int dim, const std::vector<real_t> ksi);
};

class CFemLocalQuad2D : public CFemLocalLinear2D
{
	public:
		CFemLocalQuad2D(CMesh *mesh, const int equations_number, const MeshGeometryType mgt);
		virtual ~CFemLocalQuad2D();
		//virtual void init(CProblem *pr);
	protected:
		virtual std::vector<real_t> getLocalCoordinates(const int element,
													  const CPoint3D p);
		virtual real_t getN(const int idxN, std::vector<real_t> ksi);
//		virtual real_t dNdKsi(const int idxN, const int idxKsi, std::vector<real_t> ksi);
		virtual real_t getdNdX(const int idxN, const int element, std::vector<real_t> ksi);
		virtual real_t getdNdY(const int idxN, const int element, std::vector<real_t> ksi);
};

class CFemLocalLinear3D : public CFem
{
	public:
		CFemLocalLinear3D(CMesh *mesh, const int equations_number, const MeshGeometryType mgt);
		virtual ~CFemLocalLinear3D();
		//virtual void init(CProblem *pr);
//		virtual void assembleKMatrix();
//		virtual void assembleRightVector(const int timestep);
//		virtual void setBorderConditions(const int timestep);
//		virtual void perform(const int timesteps);
	protected:
		virtual std::vector<real_t> getLocalCoordinates(const int element,
													  const CPoint3D p);
		virtual real_t getVolume(const int element);
		virtual real_t getN(const int idxN, std::vector<real_t> ksi);
		virtual real_t getdNdKsi(const int idxN, const int idxKsi, std::vector<real_t> ksi);

		virtual real_t getKK(const int idxN, const int jdxN, const int l_col, const int l_row, const int element, const std::vector<real_t> ksi);
		virtual real_t getFF(const int idxN, const int l_row, const int element, const std::vector<real_t> ksi);
		virtual real_t getCC(const int idxN, const int jdxN, const int l_col, const int l_row, const int element, const std::vector<real_t> ksi);
		real_t getVolume(const int element, const int subVolume);

		virtual real_t getdNdX(const int idxN, const int element, std::vector<real_t> ksi);
		virtual real_t getdNdY(const int idxN, const int element, std::vector<real_t> ksi);
		virtual real_t getdNdZ(const int idxN, const int element, std::vector<real_t> ksi);
		virtual real_t getdUdX(const int element, const int dim, std::vector<real_t> ksi);
		virtual real_t getdUdY(const int element, const int dim, std::vector<real_t> ksi);
		virtual real_t getdUdZ(const int element, const int dim, std::vector<real_t> ksi);
		virtual real_t getdKsidX(const int idx, const int element);
		virtual real_t getdKsidY(const int idx, const int element);
		virtual real_t getdKsidZ(const int idx, const int element);
};

class CFemLocalQuad3D : public CFemLocalLinear3D
{
	public:
		CFemLocalQuad3D(CMesh *mesh, const int equations_number, const MeshGeometryType mgt);
		virtual ~CFemLocalQuad3D();
	protected:
		virtual real_t getN(const int idxN, std::vector<real_t> ksi);
		virtual real_t getdNdKsi(const int idxN, const int idxKsi, std::vector<real_t> ksi);
		std::vector<real_t> getdKsidXYZ(const int dim, const int element, std::vector<real_t> ksi);
		virtual real_t getdNdX(const int idxN, const int element, std::vector<real_t> ksi);
		virtual real_t getdNdY(const int idxN, const int element, std::vector<real_t> ksi);
		virtual real_t getdNdZ(const int idxN, const int element, std::vector<real_t> ksi);
};


#endif /* FEM_H_ */
