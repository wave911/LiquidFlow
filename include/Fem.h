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
		virtual real_t integrate(int elem_idx, const int ksi1, const int ksi2, const int ksi3) = 0;
};

class CFemLocalLinear2D : public CFem
{
	protected:
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
		virtual real_t integrateiNjN(const int iN, const int jN, const int element);
		virtual real_t integrateiNjdN(const int iN, const int jN, const int element);
		virtual real_t integrateidNjdN(const int iN, const int jN, const int element);
};

class CFemLocalQuad2D : public CFemLocalLinear2D
{
//	private:
//		CMesh *m_mesh;
//		CProblem *m_pr;
//		real_t *m_K,
//			   *m_C,
//			   *m_F,
//			   *m_U_temp;
	public:
		CFemLocalQuad2D(CMesh *mesh);
		virtual ~CFemLocalQuad2D();
		//virtual void init(CProblem *pr);
	protected:
		virtual std::vector<real_t> getLocalCoordinates(const int element,
													  const CPoint3D p);
		virtual real_t getdNdX(const int idxN, const int element);
		virtual real_t getdNdY(const int idxN, const int element);
		virtual real_t getdNdKsi(const int idxN, const int idxKsi,
											const std::vector<real_t> ksi);
		virtual void assembleKMatrix();
		virtual void assembleRightVector(const int timestep);
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
		real_t integrate(int elem_idx, const int ksi1, const int ksi2, const int ksi3, const int ksi4);
};

#endif /* FEM_H_ */
