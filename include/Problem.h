#ifndef PROBLEM_H_
#define PROBLEM_H_

#include "Mesh.h"

class CProblem {
	protected:
		real_t m_Re;
		real_t m_t;
	public:
		CProblem() {};
		virtual ~CProblem() {};
		virtual void init() {m_Re = 0; m_t = 0;};
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
	private:
		CMesh *m_mesh;
		real_t *m_U;

	public:
		CProblem2DCircle(CMesh *mesh);
		virtual ~CProblem2DCircle();
		virtual void setU(const int idx, const short dim, const real_t value);
		virtual void setU(real_t *value);
		virtual real_t *getU();
		virtual real_t getU(const int idx, const short dim);
		virtual real_t getBorderCondition(const int idx, const int dim, const real_t time);
};

#endif /* PROBLEM_H_ */