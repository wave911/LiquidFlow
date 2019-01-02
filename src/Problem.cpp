#include "Problem.h"
#include <math.h>

#include <iostream>
using namespace std;


double sgn(real_t x) {
  if (x > 0.0) return 1.0;
  if (x < 0.0) return -1.0;
  return x;
}

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
	//delete [] m_U;
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

CProblem2DMixer::CProblem2DMixer(CMesh *mesh) : CProblem() {
	m_mesh = mesh;
	m_U = NULL;
}

CProblem2DMixer::~CProblem2DMixer() {
	//delete [] m_U;
}

void CProblem2DMixer::init() {
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

void CProblem2DMixer::setU(const int idx, const short dim, const real_t value) {
	int n = 3;
	m_U[idx * n + dim] = value;
}

void CProblem2DMixer::setU(real_t *value) {
	int n = 3;
	memcpy(m_U, value, m_mesh->getPointsNumber() * n * sizeof(real_t));
}

real_t* CProblem2DMixer::getU() {
	return m_U;
}

real_t CProblem2DMixer::getU(const int idx, const short dim) {
	int n = 3;
	return m_U[idx * n + dim];
}

real_t CProblem2DMixer::getBorderCondition(const int idx, const int dim, const real_t time) {
	CPoint3D p = m_mesh->getPointByIndex(idx);
	switch(dim) {
		case 0:
			return -p.m_y+ fi(time);
		case 1:
			return p.m_x + eta(time);
		case 2:
			return 0.5 * (pow(-p.m_y + fi(time), 2) + pow(p.m_x + eta(time), 2)) - \
					dtFi(time) * (p.m_x + eta(time)) - dtEta(time) * (-p.m_y + fi(time));
		default:
			return 0;
	}
}

real_t CProblem2DMixer::eta(const real_t time) {
	return cosf(time);
}

real_t CProblem2DMixer::fi(const real_t time) {
	return sinf(time);
}

real_t CProblem2DMixer::dtEta(const real_t time) {
	return -sinf(time);
}

real_t CProblem2DMixer::dtFi(const real_t time) {
	return cosf(time);
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
	//delete [] m_U;
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

CProblem3DCubeTest4::CProblem3DCubeTest4(CMesh *mesh) : CProblem3DPipe(mesh) {
	m_mesh = mesh;
	m_U = NULL;
	constA = 1;
	constDelta = pow(10, -6);
}

CProblem3DCubeTest4::CProblem3DCubeTest4(CMesh *mesh, const real_t A, const real_t delta) : CProblem3DPipe(mesh, A) {
	m_mesh = mesh;
	m_U = NULL;
	constA = A;
	constDelta = delta;
}

CProblem3DCubeTest4::~CProblem3DCubeTest4() {
//	delete [] m_U;
}

void CProblem3DCubeTest4::init() {
	int n = 4;
	int count = m_mesh->getPointsNumber();
	m_U = new real_t[count * n];
	std::cout << "CProblem3DCubeTest4::init() " << count << std::endl;
	for (int i = 0; i < count; i++) {
		CPoint3D p = m_mesh->getPointByIndex(i);
		m_U[i * n + 0] = -p.m_y;
		m_U[i * n + 1] = p.m_x;
		m_U[i * n + 2] = constA * p.m_z;
		m_U[i * n + 3] = 0;
	}
}

real_t CProblem3DCubeTest4::getBorderCondition(const int idx, const int dim, const real_t time) {
	CPoint3D p = m_mesh->getPointByIndex(idx);
	switch(dim) {
		case 0:
			return getV1(idx, time);
		case 1:
			return getV2(idx, time);
		case 2:
			return getV3(idx, time);
		case 3:
			return getP(idx, time);
		default:
			return 0;
	}
}

real_t CProblem3DCubeTest4::getV1(const int idx, const real_t t) {
	CPoint3D p = m_mesh->getPointByIndex(idx);
	real_t s1, s2, s3;

	s1 = sin(p.m_z - exp(t) * sin(3 * t)) + (p.m_x + p.m_z - exp(t) * (sin(3 * t) + pow(t - 0.5, 2) * sin(1/(abs(t - 0.5) + constA) ) )) * cos(p.m_y - exp(t) * sin(2 * t));
	s2 = sin(p.m_x - exp(t) * pow(t - 0.5, 2) * sin(1/(abs(t - 0.5) + constA))) + cos(p.m_x - exp(t) * pow(t - 0.5, 2) * sin(1/(abs(t - 0.5) + constA))) + cos(p.m_z - exp(t) * sin(3 * t));
	s3 = exp(t) * ((t - 0.5) * (t + 1.5) * sin(1/(abs(t - 0.5) + constA)) - sgn(t - 0.5) * pow((t - 0.5)/(abs(t - 0.5) + constA),2) * cos(1/(abs(t - 0.5) + constA)) );

	return exp(-t) * (s1 + s2) + s3;
}

real_t CProblem3DCubeTest4::getV2(const int idx, const real_t t) {
	CPoint3D p = m_mesh->getPointByIndex(idx);
	real_t s1, s2, s3;
	s1 = -sin(p.m_z - exp(t) * sin(3 * t)) - (p.m_y + p.m_z - exp(t) * (sin(2 * t) + sin(3 * t)) ) * cos(p.m_x - exp(t) * pow(t - 0.5, 2) * sin(1/(abs(t - 0.5) + constA)));
	s2 = -sin(p.m_y - exp(t) * sin(2 * t)) + cos(p.m_y - exp(t) * sin(2 * t)) + cos(p.m_y - exp(t) * sin(3 * t));
	s3 = exp(t) * (sin(2 * t) + 2 * cos(2 * t));
	return exp(-t) * (s1 + s2) + s3;
}

real_t CProblem3DCubeTest4::getV3(const int idx, const real_t t) {
	CPoint3D p = m_mesh->getPointByIndex(idx);
	real_t s1, s2, s3, s4;
	s1 = (p.m_y + p.m_z - exp(t) * (sin(2 * t) + sin(3 * t) )) * sin(p.m_x - exp(t) * pow(t - 0.5, 2) * sin(1/(abs(t - 0.5) + constA)));
	s2 = (p.m_x+ p.m_z - exp(t) * (sin(3 * t) + pow(t - 0.5, 2) * sin(1/(abs(t - 0.5) + constA)) )) * sin(p.m_y - exp(t) * sin(2 * t));
	s3 = exp(t) * (sin(3 * t) + 3 * cos(3 * t));
	return exp(-t) * (s1 + s2) + s3;
}

real_t CProblem3DCubeTest4::getP(const int idx, const real_t t) {
	CPoint3D p = m_mesh->getPointByIndex(idx);
	real_t s0, s1, s2, s3, s4, s5, s6, s7;

	s0 = -0.5 * (pow(getV1(idx, t), 2) + pow(getV2(idx, t), 2) + pow(getV3(idx, t), 2));
	s1 = (p.m_x + p.m_z - exp(t) * (sin(3 * t) + pow(t - 0.5, 2) * sin(1/(abs(t - 0.5) + constA)))) * cos(p.m_y - exp(t) * sin(2 * t)) * \
				((t - 0.5) * (t + 1.5) * sin(1/(abs(t - 0.5) + constA)) - sgn(t - 0.5) * pow((t - 0.5)/(abs(t - 0.5) + constA), 2) * cos(1/(abs(t - 0.5) + constA)));
	s2 = (sin(3 * t) + 3 * cos(3 * t)) * sin(p.m_y - exp(t) * sin(2 * t)) * (p.m_x + p.m_z - exp(t) * (sin(3 * t) + pow(t - 0.5, 2) * sin(1/(abs(t - 0.5) + constA))));
	s3 = p.m_x * (-exp(t) * sin(1/(abs(t - 0.5) + constA)) * (t * t + 3 * t + 1/4.0) + 2 * constDelta * sgn(t - 0.5) * (1/pow((abs(t - 0.5) + constA), 2)) * \
			cos(1/(abs(t - 0.5) + constA)) * (t * t - 1/4.0 + constA * ( (t - 0.5)/(abs(t - 0.5) + constA) )) + exp(t) * pow((t - 0.5)/(abs(t - 0.5) + constA),2) * \
			(2 * constDelta * (t - 0.5) * cos(1/(abs(t - 0.5) + constA)) + sin(1/(abs(t - 0.5) + constA)) * (1/pow(abs(t - 0.5) + constA ,2))));
	s4 = sin(p.m_x - exp(t) * pow(t - 0.5, 2) * sin(1/(abs(t - 0.5) + constA))) * ( (t - 0.5) * (t + 1.5) * sin(1/(abs(t - 0.5) + constA)) - \
			sgn(t - 0.5) * pow((t - 0.5)/(abs(t - 0.5) + constA), 2) * cos(1/(abs(t - 0.5) + constA)) + (sin(3 * t) + 3 * cos(3 * t)) * \
			(p.m_y + p.m_z - exp(t) * (sin(2 * t) + sin(3 * t))));
	s5 = cos(p.m_x - exp(t) * pow(t - 0.5, 2) * sin(1/(abs(t - 0.5) + constA))) * ( (t - 0.5) * (t + 1.5) * sin(1/(abs(t - 0.5) + constA)) - \
			sgn(t - 0.5) * pow((t - 0.5)/(abs(t - 0.5) + constA), 2) * cos(1/(abs(t - 0.5) + constA)) - (sin(2 * t) + 2 * cos(2 * t)) * \
			(p.m_y + p.m_z - exp(t) * (sin(2 * t) + sin(3 * t))) ) - p.m_y * exp(t) * (-3 * sin(2 * t) + 4 * cos(2 * t));
	s6 = (sin(2 * t) + 2 * cos(2 * t)) * (cos(p.m_y - exp(t) * sin(2 * t)) - sin(p.m_y - exp(t) * sin(2 * t))) - \
			p.m_z * exp(t) * (6 * cos(3 * t) - 8 * sin(3 * t));
	s7 = (sin(p.m_z - exp(t) * sin(3 * t)) + cos(p.m_z - exp(t) * sin(3 * t))) * ( (t - 0.5) * (t + 1.5) * sin(1/(abs(t - 0.5) + constA)) - \
			sgn(t - 0.5) * pow((t - 0.5)/(abs(t - 0.5) + constA), 2) * cos(1/(abs(t - 0.5) + constA)) ) - \
			(sin(2 * t) + 2 * cos(2 * t)) * (sin (p.m_z - exp(t) * sin(3 * t)) - cos(p.m_z - exp(t) * sin(3 * t)) );
	return s0 + s1 + s2 + s3 + s4 + s5 + s6 + s7;
}


