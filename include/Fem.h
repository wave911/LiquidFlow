#ifndef FEM_H_
#define FEM_H_

#include "Mesh.h"

#ifdef LIBBLASLAPACK
#endif

class CFem 
{
	public:
		CFem() {};
		virtual ~CFem() {};
		virtual void init();
		virtual std::vector<real_t> getLocalCoordinates(const int element, 
													  const CPoint3D p) = 0;		
		virtual real_t getdNdX(const int idxN, const int element) = 0;	
		virtual real_t getdNdY(const int idxN, const int element) = 0;
		virtual void assembleKMatrix();
		virtual void assembleRightVector();
		void dgemv(char *trans, int m, int n, real_t alpha, real_t *a, 
				   int lda, real_t *x, int incx, real_t beta, real_t *y, int incy);	
};

class CFemLocalLinear2D : public CFem
{
	private:
		CMesh *m_mesh;
		CProblem *m_pr;
		real_t *m_K,
			   *m_C,
			   *m_F;	

	public:
		CFemLocalLinear2D(CMesh *mesh);
		virtual ~CFemLocalLinear2D();
		virtual void init(CProblem *pr);
		virtual std::vector<real_t> getLocalCoordinates(const int element, 
													  const CPoint3D p);		
		virtual real_t getdNdX(const int idxN, const int element);
		virtual real_t getdNdY(const int idxN, const int element);
		virtual void assembleKMatrix();
		virtual void assembleRightVector();
		real_t getSquare(const int element);
	protected:
		virtual real_t getdKsidX(const int idx, const int element);
		virtual real_t getdKsidY(const int idx, const int element);
		virtual real_t getdNdKsi(const int idxN, const int idxKsi);
		virtual real_t getdUdX(const int element, const int dim);
		virtual real_t getdUdY(const int element, const int dim);
};

class CFemLocalLinear3D : public CFem
{
	private:
		CMesh *m_mesh;

	public:
		CFemLocalLinear3D(CMesh *mesh);
		virtual ~CFemLocalLinear3D();
		virtual void init();
		virtual std::vector<real_t> getLocalCoordinates(const int element, 
													  const CPoint3D p);		
		virtual real_t getdNdX(const int idxN, const int element);
		virtual real_t getdNdY(const int idxN, const int element);
	
	protected:
		real_t getVolume(const int element, const int subVolume);
		virtual real_t getdKsidX(const int idx, const int element);
		virtual real_t getdKsidY(const int idx, const int element);
		virtual real_t getdKsidZ(const int idx, const int element);
		virtual real_t getdNdKsi(const int idxN, const int idxKsi);
};

#endif /* FEM_H_ */