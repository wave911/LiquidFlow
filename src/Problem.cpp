#include "Problem.h"

real_t CProblem::getRe() {
	return m_Re;
}

real_t CProblem::getTau() {
	return m_t;
}

real_t CProblem::setRe(const real_t Re) {
	m_Re = Re;
}

real_t CProblem::setTau(const real_t t) {
	m_t = t;
}

CProblem2DCircle::CProblem2DCircle(CMesh *mesh) : CProblem() {
	m_mesh = mesh;
	m_U = NULL;
}

CProblem2DCircle::~CProblem2DCircle() {
	delete [] m_U;
}

void CProblem2DCircle::init() {
	int n = 3;
	int count = m_mesh->getPointsNumber();
	m_U = new real_t[count * n];

	for (int i = 0; i < count; i++) {
		CPoint3D p = m_mesh->getPointByIndex(const int idx);
		m_U[i * n + 0] = -p.m_y;
		m_U[i * n + 1] = p.m_x;
		m_U[i * n + 2] = 0;
	}
}

void CProblem2DCircle::setU(const int idx, const short dim, const real_t value) {
	int n = 3;
	m_U[idx * n + dim] = value;
}

void CProblem2DCircle::setU(real_t *value) {
	*m_U = *value;
}

real_t* CProblem2DCircle::getU() {
	return m_U;
}

real_t CProblem2DCircle::getU(const int idx, const short dim) {
	int n = 3;
	return m_U[idx * n + dim];
}
