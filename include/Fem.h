#ifndef FEM_H_
#define FEM_H_

#include "Mesh.h"
#include "Problem.h"

#ifdef LIBBLASLAPACK
#endif

class CFem 
{
	public:
		CFem() {};
		virtual ~CFem() {};
		virtual void init(CProblem *pr) = 0;
		virtual std::vector<real_t> getLocalCoordinates(const int element, 
													  const CPoint3D p) = 0;		
		virtual real_t getdNdX(const int idxN, const int element) = 0;	
		virtual real_t getdNdY(const int idxN, const int element) = 0;
		virtual void assembleKMatrix() = 0;
		virtual void assembleRightVector(const int timestep) = 0;
		virtual void setBorderConditions(const int timestep) = 0;
		virtual void perform(const int timesteps) = 0;
		void dgemv(char *trans, int m, int n, real_t alpha, real_t *a, 
				   int lda, real_t *x, int incx, real_t beta, real_t *y, int incy);	
		void dgesv(int n, real_t *M, real_t *B );
};

class CFemLocalLinear2D : public CFem
{
	private:
		CMesh *m_mesh;
		CProblem *m_pr;
		real_t *m_K,
			   *m_C,
			   *m_F,
			   *m_U_temp;

	public:
		CFemLocalLinear2D(CMesh *mesh);
		virtual ~CFemLocalLinear2D();
		virtual void init(CProblem *pr);
		virtual void assembleKMatrix();
		virtual void assembleRightVector(const int timestep);
		virtual void setBorderConditions(const int timestep);
		virtual void perform(const int timesteps);
		real_t getSquare(const int element);
	protected:
		virtual std::vector<real_t> getLocalCoordinates(const int element,
													  const CPoint3D p);
		virtual real_t getdNdX(const int idxN, const int element);
		virtual real_t getdNdY(const int idxN, const int element);
		virtual real_t getdKsidX(const int idx, const int element);
		virtual real_t getdKsidY(const int idx, const int element);
		virtual real_t getdNdKsi(const int idxN, const int idxKsi);
		virtual real_t getdUdX(const int element, const int dim);
		virtual real_t getdUdY(const int element, const int dim);
};

class CFemLocalQuad2D : public CFemLocalLinear2D
{
	private:
		CMesh *m_mesh;
		CProblem *m_pr;
		real_t *m_K,
			   *m_C,
			   *m_F,
			   *m_U_temp;
	public:
		CFemLocalQuad2D(CMesh *mesh);
		virtual ~CFemLocalQuad2D();
		void test();
	protected:
		virtual std::vector<real_t> getLocalCoordinates(const int element,
													  const CPoint3D p);
		virtual real_t getdNdX(const int idxN, const int element);
		virtual real_t getdNdY(const int idxN, const int element);
		virtual real_t getdNdKsi(const int idxN, const int idxKsi,
											const std::vector<real_t> ksi);
};

class CFemLocalLinear3D : public CFem
{
	private:
		CMesh *m_mesh;

	public:
		CFemLocalLinear3D(CMesh *mesh);
		virtual ~CFemLocalLinear3D();
		virtual void init(CProblem *pr);
		virtual std::vector<real_t> getLocalCoordinates(const int element, 
													  const CPoint3D p);		
		virtual real_t getdNdX(const int idxN, const int element);
		virtual real_t getdNdY(const int idxN, const int element);
		virtual void assembleKMatrix();
		virtual void assembleRightVector(const int timestep);
		virtual void setBorderConditions(const int timestep);
		virtual void perform(const int timesteps);
	protected:
		real_t getVolume(const int element, const int subVolume);
		virtual real_t getdKsidX(const int idx, const int element);
		virtual real_t getdKsidY(const int idx, const int element);
		virtual real_t getdKsidZ(const int idx, const int element);
		virtual real_t getdNdKsi(const int idxN, const int idxKsi);
};

#endif /* FEM_H_ */
