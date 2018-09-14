#ifndef PROBLEM_H_
#define PROBLEM_H_

#include "Mesh.h"

class CProblem {
	protected:
		real_t m_Re;
		real_t m_t;
	public:
		CProblem() {m_Re = 0; m_t = 0;};
		virtual ~CProblem() {};
		virtual void init() {};
		real_t getRe();
		real_t getTau();
		virtual real_t setRe(const real_t Re);
		virtual real_t setTau(const real_t t);
		virtual void setU(const int idx, const short dim, const real_t value) {};
		virtual void setU(real_t *value) {};
		virtual real_t *getU() {};
		virtual real_t getU(const int idx, const short dim) {};
		virtual real_t getBorderCondition(const int idx, const int dim, const real_t time) {};		
};

class CProblem2DCircle : public CProblem {
	protected:
		CMesh *m_mesh;
		real_t *m_U;
	public:
		CProblem2DCircle(CMesh *mesh);
		virtual ~CProblem2DCircle();
		virtual void init();
		virtual void setU(const int idx, const short dim, const real_t value);
		virtual void setU(real_t *value);
		virtual real_t *getU();
		virtual real_t getU(const int idx, const short dim);
		virtual real_t getBorderCondition(const int idx, const int dim, const real_t time);
};

class CProblem2DMixer : public CProblem {
	protected:
		CMesh *m_mesh;
		real_t *m_U;
	public:
		CProblem2DMixer(CMesh *mesh);
		virtual ~CProblem2DMixer();
		virtual void init();
		virtual void setU(const int idx, const short dim, const real_t value);
		virtual void setU(real_t *value);
		virtual real_t *getU();
		virtual real_t getU(const int idx, const short dim);
		virtual real_t getBorderCondition(const int idx, const int dim, const real_t time);
		virtual real_t eta(const real_t time);
		virtual real_t fi(const real_t time);
		virtual real_t dtEta(const real_t time);
		virtual real_t dtFi(const real_t time);
};

class CProblem3DPipe : public CProblem {
	protected:
		CMesh *m_mesh;
		real_t *m_U;
		real_t constA;
	public:
		CProblem3DPipe(CMesh *mesh);
		CProblem3DPipe(CMesh *mesh, const real_t A);
		virtual ~CProblem3DPipe();
		virtual void init();
		virtual void setU(const int idx, const short dim, const real_t value);
		virtual void setU(real_t *value);
		virtual real_t *getU();
		virtual real_t getU(const int idx, const short dim);
		virtual real_t getBorderCondition(const int idx, const int dim, const real_t time);
};

class CProblem3DCubeTest4 : public CProblem3DPipe {
	protected:
		CMesh *m_mesh;
		real_t *m_U;
		real_t constA;
		real_t constDelta;
	public:
		CProblem3DCubeTest4(CMesh *mesh);
		CProblem3DCubeTest4(CMesh *mesh, const real_t A, const real_t constDelta);
		virtual ~CProblem3DCubeTest4();
		virtual void init();
		virtual real_t getBorderCondition(const int idx, const int dim, const real_t time);
	private:
		virtual real_t getV1(const int idx, const real_t time);
		virtual real_t getV2(const int idx, const real_t time);
		virtual real_t getV3(const int idx, const real_t time);
		virtual real_t getP(const int idx, const real_t time);
};

#endif /* PROBLEM_H_ */
