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
		CPoint3D p = m_mesh->getPointByIndex(i);
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
	int n = 3;
	memcpy(m_U, value, m_mesh->getPointsNumber() * n * sizeof(real_t));
}

real_t* CProblem2DCircle::getU() {
	return m_U;
}

real_t CProblem2DCircle::getU(const int idx, const short dim) {
	int n = 3;
	return m_U[idx * n + dim];
}

real_t CProblem2DCircle::getBorderCondition(const int idx, const int dim, const real_t time) {
	CPoint3D p = m_mesh->getPointByIndex(idx);
	switch(dim) {
		case 0:
			return -p.m_y;
		case 1:
			return p.m_x;
		case 2:
			return -0.5 * (p.m_x * p.m_x + p.m_y * p.m_y);
		default:
			return 0;
	}
}

CProblem3DPipe::CProblem3DPipe(CMesh *mesh) : CProblem() {
	m_mesh = mesh;
	m_U = NULL;
	constA = 1;
}

CProblem3DPipe::CProblem3DPipe(CMesh *mesh, const real_t A) : CProblem() {
	m_mesh = mesh;
	m_U = NULL;
	constA = A;
}

CProblem3DPipe::~CProblem3DPipe() {
	delete [] m_U;
}

void CProblem3DPipe::init() {
	int n = 4;
	int count = m_mesh->getPointsNumber();
	m_U = new real_t[count * n];

	for (int i = 0; i < count; i++) {
		CPoint3D p = m_mesh->getPointByIndex(i);
		m_U[i * n + 0] = -p.m_y;
		m_U[i * n + 1] = p.m_x;
		m_U[i * n + 2] = constA * p.m_z;
		m_U[i * n + 3] = 0;
	}
}

void CProblem3DPipe::setU(const int idx, const short dim, const real_t value) {
	int n = 4;
	m_U[idx * n + dim] = value;
}

void CProblem3DPipe::setU(real_t *value) {
	int n = 4;
	memcpy(m_U, value, m_mesh->getPointsNumber() * n * sizeof(real_t));
}

real_t* CProblem3DPipe::getU() {
	return m_U;
}

real_t CProblem3DPipe::getU(const int idx, const short dim) {
	int n = 4;
	return m_U[idx * n + dim];
}

real_t CProblem3DPipe::getBorderCondition(const int idx, const int dim, const real_t time) {
	CPoint3D p = m_mesh->getPointByIndex(idx);
	switch(dim) {
		case 0:
			return -p.m_y;
		case 1:
			return p.m_x;
		case 2:
			return constA * time;
		case 3:
			return -0.5 * (1 - (p.m_x * p.m_x + p.m_y * p.m_y)) - constA * p.m_z;
		default:
			return 0;
	}
}

