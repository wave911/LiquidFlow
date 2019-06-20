#include "Fem.h"
#include <iostream>
#include "f2c.h"
#include "dgemv.h"
#include "dgesv.h"

using namespace std;

CFemLocalQuad2D::CFemLocalQuad2D(CMesh *mesh, const int equations_number, const MeshGeometryType mgt) : CFemLocalLinear2D(mesh, equations_number, mgt) {
	m_mesh = mesh;
	m_mesh_geometry_type = mgt;
//	m_K = nullptr;
//	m_C = nullptr;
//	m_F = nullptr;
//	m_U_temp = nullptr;
//	m_pr = nullptr;
}

CFemLocalQuad2D::~CFemLocalQuad2D() {
//	if (this->m_K != nullptr)
//		delete [] this->m_K;
//	if (this->m_C != nullptr)
//		delete [] this->m_C;
//	if (this->m_F != nullptr)
//		delete [] this->m_F;
//	if (this->m_U_temp != nullptr)
//		delete [] this->m_U_temp;
}

std::vector<real_t> CFemLocalQuad2D::getLocalCoordinates(const int element,
											  const CPoint3D p) {
	real_t *matrix = new real_t[9];
	real_t *vect = new real_t[3];
	real_t *ksi = new real_t[3];
	vector<int> points = m_mesh->getElementByIndex(element);
	CPoint3D p1 = m_mesh->getPointByIndex(points[0]);
	CPoint3D p2 = m_mesh->getPointByIndex(points[1]);
	CPoint3D p3 = m_mesh->getPointByIndex(points[2]);

	real_t square = getVolume(element);
	matrix[0] = p2.m_x * p3.m_y - p3.m_x * p2.m_y;
	matrix[1] = p3.m_x * p1.m_y - p1.m_x * p3.m_y;
	matrix[2] = p1.m_x * p2.m_y - p2.m_x * p1.m_y;
	matrix[3] = p2.m_y - p3.m_y;
	matrix[4] = p3.m_y - p1.m_y;
	matrix[5] = p1.m_y - p2.m_y;
	matrix[6] = p3.m_x - p2.m_x;
	matrix[7] = p1.m_x - p3.m_x;
	matrix[8] = p2.m_x - p1.m_x;

	vect[0] = 1;
	vect[1] = p.m_x;
	vect[2] = p.m_y;

	char *ch = "N";
	int m_m = 3,
		m_n = 3,
		lda = 3,
		incx = 1,
		incy = 1;
	real_t tau = 1/(2 * square),
		   beta = 0;

	dgemv(ch, m_m, m_n, tau, &matrix[0], lda, &vect[0],
				 incx, beta, ksi, incy);

	std::vector<real_t> res = std::vector<real_t>(ksi, ksi + 3);
	delete(ksi);
	delete(matrix);
	delete(vect);

	return res;
}

real_t CFemLocalQuad2D::getN(const int idxN, const std::vector<real_t> ksi) {
	if (ksi.size() == 0) {
		return 0;
	}
	real_t res = 0;

	switch(idxN) {
	case 0:
		res = ksi[0] * (2 * ksi[0] - 1);
		break;
	case 1:
		res = ksi[1] * (2 * ksi[1] - 1);
		break;
	case 2:
		res = ksi[2] * (2 * ksi[2] - 1);
		break;
	case 3:
		res = 4 * ksi[0] * ksi[1];
		break;
	case 4:
		res = 4 * ksi[1] * ksi[2];
		break;
	case 5:
		res = 4 * ksi[2] * ksi[0];
		break;
	default:
		break;
	}
	return res;
}

//real_t CFemLocalQuad2D::getdNdX(const int idxN, const int element, std::vector<real_t> ksi) {
//	real_t sum = 0;
////	std::vector<int> e = m_mesh->getElementByIndex(element);
////	CPoint3D p = m_mesh->getPointByIndex(e[idxN]);
////	std::vector<real_t> ksi = this->getLocalCoordinates(element, p);
//
//	for (int i = 0; i < ksi.size(); i++) {
//		sum += getdNdKsi(idxN, i, ksi) * getdKsidX(i, element);
//	}
//	return sum;
//}

real_t CFemLocalQuad2D::getdNdX(const int idxN, const int element, std::vector<real_t> ksi) {
	vector<int> points = m_mesh->getElementByIndex(element);
	CPoint3D p1 = m_mesh->getPointByIndex(points[0]);
	CPoint3D p2 = m_mesh->getPointByIndex(points[1]);
	CPoint3D p3 = m_mesh->getPointByIndex(points[2]);
	//std::vector<real_t> ksi = this->getLocalCoordinates(element, m_mesh->getPointByIndex(points[idxN]));

	real_t J2 = (p2.m_x - p1.m_x) * (p3.m_y - p1.m_y) - (p1.m_y - p2.m_y) * (p1.m_x - p3.m_x),
		 Jy23 = p2.m_y - p3.m_y,
		 Jy31 = p3.m_y - p1.m_y,
		 Jy12 = p1.m_y - p2.m_y;

	switch(idxN) {
		case 0:
			return (4 * ksi[0] - 1) * Jy23/J2;
		case 1:
			return (4 * ksi[1] - 1) * Jy31/J2;
		case 2:
			return (4 * ksi[2] - 1) * Jy12/J2;
		case 3:
			return 4 * (ksi[1] * Jy23 + ksi[0] * Jy31) / J2;
		case 4:
			return 4 * (ksi[2] * Jy31 + ksi[1] * Jy12) / J2;
		case 5:
			return 4 * (ksi[0] * Jy12 + ksi[2] * Jy23) / J2;
	}
}

//real_t CFemLocalQuad2D::getdNdY(const int idxN, const int element, std::vector<real_t> ksi) {
//	real_t sum = 0;
////	std::vector<int> e = m_mesh->getElementByIndex(element);
////	CPoint3D p = m_mesh->getPointByIndex(e[idxN]);
////	std::vector<real_t> ksi = this->getLocalCoordinates(element, p);
//
//	for (int i = 0; i < ksi.size(); i++) {
//		sum += getdNdKsi(idxN, i, ksi) * getdKsidY(i, element);
//	}
//	return sum;
//}

real_t CFemLocalQuad2D::getdNdY(const int idxN, const int element, std::vector<real_t> ksi) {
	vector<int> points = m_mesh->getElementByIndex(element);
	CPoint3D p1 = m_mesh->getPointByIndex(points[0]);
	CPoint3D p2 = m_mesh->getPointByIndex(points[1]);
	CPoint3D p3 = m_mesh->getPointByIndex(points[2]);
	//std::vector<real_t> ksi = this->getLocalCoordinates(element, m_mesh->getPointByIndex(points[idxN]));

	real_t J2 = (p2.m_x - p1.m_x) * (p3.m_y - p1.m_y) - (p1.m_y - p2.m_y) * (p1.m_x - p3.m_x),
		 Jx32 = p3.m_x - p2.m_x,
		 Jx13 = p1.m_x - p3.m_x,
		 Jx21 = p2.m_x - p1.m_x;

	switch(idxN) {
		case 0:
			return (4 * ksi[0] - 1) * Jx32/J2;
		case 1:
			return (4 * ksi[1] - 1) * Jx13/J2;
		case 2:
			return (4 * ksi[2] - 1) * Jx21/J2;
		case 3:
			return 4 * (ksi[1] * Jx32 + ksi[0] * Jx13) / J2;
		case 4:
			return 4 * (ksi[2] * Jx13 + ksi[1] * Jx21) / J2;
		case 5:
			return 4 * (ksi[0] * Jx21 + ksi[2] * Jx32) / J2;
	}
}
