#include "Fem.h"
#include <iostream>
#include "f2c.h"
#include "dgemv.h"
#include "dgesv.h"

using namespace std;

CFemLocalLinear3D::CFemLocalLinear3D(CMesh *mesh) {
	m_mesh = mesh;
}

CFemLocalLinear3D::~CFemLocalLinear3D() {
}

std::vector<real_t> CFemLocalLinear3D::getLocalCoordinates(const int element,
											  const CPoint3D p) {
	real_t *matrix = new real_t[16];
	real_t *vect = new real_t[4];
	real_t *ksi = new real_t[4];
	vector<int> points = m_mesh->getElementByIndex(element);
	CPoint3D p1 = m_mesh->getPointByIndex(points[0]);
	CPoint3D p2 = m_mesh->getPointByIndex(points[1]);
	CPoint3D p3 = m_mesh->getPointByIndex(points[2]);
	CPoint3D p4 = m_mesh->getPointByIndex(points[3]);

	real_t volume = getVolume(element, 0);
	real_t vo1 = getVolume(element, 1);
	real_t vo2 = getVolume(element, 2);
	real_t vo3 = getVolume(element, 3);
	real_t vo4 = getVolume(element, 4);

	matrix[0] = 6 * vo1;
	matrix[1] = 6 * vo2;
	matrix[2] = 6 * vo3;
	matrix[3] = 6 * vo4;

	cout << vo1 + vo2 + vo3 + vo4 << " " << volume << endl;

	matrix[4] = (p4.m_y - p2.m_y) * (p3.m_z - p2.m_z) - (p3.m_y - p2.m_y) * (p4.m_z - p2.m_z);
	matrix[5] = (p3.m_y - p1.m_y) * (p4.m_z - p3.m_z) - (p3.m_y - p4.m_y) * (p1.m_z - p3.m_z);
	matrix[6] = (p2.m_y - p4.m_y) * (p1.m_z - p4.m_z) - (p1.m_y - p4.m_y) * (p2.m_z - p4.m_z);
	matrix[7] = (p1.m_y - p3.m_y) * (p2.m_z - p1.m_z) - (p1.m_y - p2.m_y) * (p3.m_z - p1.m_z);

	matrix[8] = (p3.m_x - p2.m_x) * (p4.m_z - p2.m_z) - (p4.m_x - p2.m_x) * (p3.m_z - p2.m_z);
	matrix[9] = (p4.m_x - p3.m_x) * (p3.m_z - p1.m_z) - (p1.m_x - p3.m_x) * (p3.m_z - p4.m_z);
	matrix[10] = (p1.m_x - p4.m_x) * (p2.m_z - p4.m_z) - (p2.m_x - p4.m_x) * (p1.m_z - p4.m_z);
	matrix[11] = (p2.m_x - p1.m_x) * (p1.m_z - p3.m_z) - (p3.m_x - p1.m_x) * (p1.m_z - p2.m_z);

	matrix[12] = (p4.m_x - p2.m_x) * (p3.m_y - p2.m_y) - (p3.m_x - p2.m_x) * (p4.m_y - p2.m_y);
	matrix[13] = (p3.m_x - p1.m_x) * (p4.m_y - p3.m_y) - (p3.m_x - p4.m_x) * (p1.m_y - p3.m_y);
	matrix[14] = (p2.m_x - p4.m_x) * (p1.m_y - p4.m_y) - (p1.m_x - p4.m_x) * (p2.m_y - p4.m_y);
	matrix[15] = (p1.m_x - p3.m_x) * (p2.m_y - p1.m_y) - (p1.m_x - p2.m_x) * (p3.m_y - p1.m_y);

	vect[0] = 1;
	vect[1] = p.m_x;
	vect[2] = p.m_y;
	vect[3] = p.m_z;

	char *ch = "N";
	int m_m = 4,
		m_n = 4,
		lda = 4,
		incx = 1,
		incy = 1;
	real_t tau = 1/(6 * volume),
		   beta = 0;

	dgemv(ch, m_m, m_n, tau, &matrix[0], lda, &vect[0],
				 incx, beta, ksi, incy);

	std::vector<real_t> res = std::vector<real_t>(ksi, ksi + 4);
	delete(ksi);
	delete(matrix);
	delete(vect);

	return res;
}

real_t CFemLocalLinear3D::getdNdX(const int idxN, const int element) {

}

real_t CFemLocalLinear3D::getdNdY(const int idxN, const int element) {

}

real_t CFemLocalLinear3D::getVolume(const int element, const int subVolume) {
	vector<int> points = m_mesh->getElementByIndex(element);
	CPoint3D p1 = m_mesh->getPointByIndex(points[0]);
	CPoint3D p2 = m_mesh->getPointByIndex(points[1]);
	CPoint3D p3 = m_mesh->getPointByIndex(points[2]);
	CPoint3D p4 = m_mesh->getPointByIndex(points[3]);
	real_t vo1 = (p2.m_x * (p3.m_y * p4.m_z - p4.m_y * p3.m_z) + \
				  p3.m_x * (p4.m_y * p2.m_z - p2.m_y * p4.m_z) + \
				  p4.m_x * (p2.m_y * p3.m_z - p3.m_y * p2.m_z))/6;

	real_t vo2 = (p1.m_x * (p4.m_y * p3.m_z - p3.m_y * p4.m_z) + \
				  p3.m_x * (p1.m_y * p4.m_z - p4.m_y * p1.m_z) + \
				  p4.m_x * (p3.m_y * p1.m_z - p1.m_y * p3.m_z))/6;

	real_t vo3 = (p1.m_x * (p2.m_y * p4.m_z - p4.m_y * p2.m_z) + \
				  p2.m_x * (p4.m_y * p1.m_z - p1.m_y * p4.m_z) + \
				  p4.m_x * (p1.m_y * p2.m_z - p2.m_y * p1.m_z))/6;

	real_t vo4 = (p1.m_x * (p3.m_y * p2.m_z - p2.m_y * p3.m_z) + \
				  p2.m_x * (p1.m_y * p3.m_z - p3.m_y * p1.m_z) + \
				  p3.m_x * (p2.m_y * p1.m_z - p1.m_y * p2.m_z))/6;

	switch(subVolume) {
		case 0:
			return vo1 + vo2 + vo3 + vo4;
		case 1:
			return vo1;
		case 2:
			return vo2;
		case 3:
			return vo3;
		case 4:
			return vo4;
	}
}

real_t CFemLocalLinear3D::getdKsidX(const int idx, const int element) {

}

real_t CFemLocalLinear3D::getdKsidY(const int idx, const int element) {

}

real_t CFemLocalLinear3D::getdKsidZ(const int idx, const int element) {

}

real_t CFemLocalLinear3D::getdNdKsi(const int idxN, const int idxKsi) {

}

void CFemLocalLinear3D::assembleKMatrix() {

}

void CFemLocalLinear3D::assembleRightVector(const int timestep) {

}

void CFemLocalLinear3D::setBorderConditions(const int timestep) {

}

void CFemLocalLinear3D::perform(const int timesteps) {

}

void CFemLocalLinear3D::init(CProblem *pr) {

}
