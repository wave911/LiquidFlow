#include "Fem.h"
#include <iostream>
#include "f2c.h"
#include "dgemv.h"
#include "dgesv.h"

using namespace std;

CFemLocalLinear3D::CFemLocalLinear3D(CMesh *mesh, const int equations_number, const MeshGeometryType mgt) :
		CFem(equations_number, mgt) {
	m_mesh = mesh;
	m_mesh_geometry_type = mgt;
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

real_t CFemLocalLinear3D::getN(const int idxN, const std::vector<real_t> ksi) {
	if (ksi.size() == 0) {
		return 0;
	}
	real_t res = 0;

	switch(idxN) {
	case 0:
		res = ksi[0];
		break;
	case 1:
		res = ksi[1];
		break;
	case 2:
		res = ksi[2];
		break;
	case 3:
		res = ksi[3];
		break;
	default:
		break;
	}
	return res;
}

real_t CFemLocalLinear3D::getdNdX(const int idxN, const int element, std::vector<real_t> ksi) {
	real_t sum = 0;
	for (int i = 0; i < ksi.size(); i++) {
		sum += getdNdKsi(idxN, i, ksi) * getdKsidX(i, element);
	}
	return sum;
}

real_t CFemLocalLinear3D::getdNdY(const int idxN, const int element, std::vector<real_t> ksi) {
	real_t sum = 0;
	for (int i = 0; i < ksi.size(); i++) {
		sum += getdNdKsi(idxN, i, ksi) * getdKsidY(i, element);
	}
	return sum;
}

real_t CFemLocalLinear3D::getdNdZ(const int idxN, const int element, std::vector<real_t> ksi) {
	real_t sum = 0;
	for (int i = 0; i < ksi.size(); i++) {
		sum += getdNdKsi(idxN, i, ksi) * getdKsidZ(i, element);
	}
	return sum;
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

real_t CFemLocalLinear3D::getVolume(const int element) {
	return getVolume(element, 0);
}

real_t CFemLocalLinear3D::getdKsidX(const int idx, const int element) {
	vector<int> points = m_mesh->getElementByIndex(element);
	CPoint3D p1 = m_mesh->getPointByIndex(points[0]);
	CPoint3D p2 = m_mesh->getPointByIndex(points[1]);
	CPoint3D p3 = m_mesh->getPointByIndex(points[2]);
	CPoint3D p4 = m_mesh->getPointByIndex(points[3]);
	real_t volume = this->getVolume(element, 0);

	real_t yy1, yy2, zz1, zz2;
	switch(idx) {
		case 0: //i = 1 j = 2 k = 3  <<
			yy1 = p4.m_y - p2.m_y;
			zz1 = p3.m_z - p2.m_z;
			yy2 = p3.m_y - p2.m_y;
			zz2 = p4.m_z - p2.m_z;
			return (yy1 * zz1 - yy2 * zz2)/(6 * volume);

		case 1: //i = 2 j = 3 k = 1  <<
			yy1 = p3.m_y - p1.m_y;
			zz1 = p4.m_z - p3.m_z;
			yy2 = p3.m_y - p4.m_y;
			zz2 = p1.m_z - p3.m_z;
			return (yy1 * zz1 - yy2 * zz2)/(6 * volume);

		case 2: //i = 3 j = 1 k = 2  <<
			yy1 = p2.m_y - p4.m_y;
			zz1 = p1.m_z - p4.m_z;
			yy2 = p1.m_y - p4.m_y;
			zz2 = p2.m_z - p4.m_z;
			return (yy1 * zz1 - yy2 * zz2)/(6 * volume);

		case 3: //i = 3 j = 1 k = 2  <<
			yy1 = p1.m_y - p3.m_y;
			zz1 = p2.m_z - p1.m_z;
			yy2 = p1.m_y - p2.m_y;
			zz2 = p3.m_z - p1.m_z;
			return (yy1 * zz1 - yy2 * zz2)/(6 * volume);
	}
}

real_t CFemLocalLinear3D::getdKsidY(const int idx, const int element) {
	vector<int> points = m_mesh->getElementByIndex(element);
	CPoint3D p1 = m_mesh->getPointByIndex(points[0]);
	CPoint3D p2 = m_mesh->getPointByIndex(points[1]);
	CPoint3D p3 = m_mesh->getPointByIndex(points[2]);
	CPoint3D p4 = m_mesh->getPointByIndex(points[3]);
	real_t volume = this->getVolume(element, 0);

	real_t xx1, xx2, zz1, zz2;
	switch(idx) {
		case 0: //i = 1 j = 2 k = 3  <<
			zz1 = p4.m_z - p2.m_z;
			xx1 = p3.m_x - p2.m_x;
			zz2 = p3.m_z - p2.m_z;
			xx2 = p4.m_x - p2.m_x;
			return (xx1 * zz1 - zz2 * xx2)/(6 * volume);
		case 1: //i = 2 j = 3 k = 1  <<
			zz1 = p3.m_z - p1.m_z;
			xx1 = p4.m_x - p3.m_x;
			zz2 = p3.m_z - p4.m_z;
			xx2 = p1.m_x - p3.m_x;
			return (zz1 * xx1 - zz2 * xx2)/(6 * volume);
		case 2: //i = 3 j = 1 k = 2  <<
			zz1 = p2.m_z - p4.m_z;
			xx1 = p1.m_x - p4.m_x;
			zz2 = p1.m_z - p4.m_z;
			xx2 = p2.m_x - p4.m_x;
			return (zz1 * xx1 - zz2 * xx2)/(6 * volume);
		case 3: //i = 3 j = 1 k = 2  <<
			zz1 = p1.m_z - p3.m_z;
			xx1 = p2.m_x - p1.m_x;
			zz2 = p1.m_z - p2.m_z;
			xx2 = p3.m_x - p1.m_x;
			return (zz1 * xx1 - zz2 * xx2)/(6 * volume);
	}
}

real_t CFemLocalLinear3D::getdKsidZ(const int idx, const int element) {
	vector<int> points = m_mesh->getElementByIndex(element);
	CPoint3D p1 = m_mesh->getPointByIndex(points[0]);
	CPoint3D p2 = m_mesh->getPointByIndex(points[1]);
	CPoint3D p3 = m_mesh->getPointByIndex(points[2]);
	CPoint3D p4 = m_mesh->getPointByIndex(points[3]);
	real_t volume = this->getVolume(element, 0);

	real_t xx1, xx2, yy1, yy2;
	switch(idx) {
		case 0: //i = 1 j = 2 k = 3  <<
			xx1 = p4.m_x - p2.m_x;
			yy1 = p3.m_y - p2.m_y;
			xx2 = p3.m_x - p2.m_x;
			yy2 = p4.m_y - p2.m_y;
			return (xx1 * yy1 - xx2 * yy2)/(6 * volume);
		case 1: //i = 2 j = 3 k = 1  <<
			xx1 = p3.m_x - p1.m_x;
			yy1 = p4.m_y - p3.m_y;
			xx2 = p3.m_x - p4.m_x;
			yy2 = p1.m_y - p3.m_y;
			return (xx1 * yy1 - xx2 * yy2)/(6 * volume);
		case 2: //i = 3 j = 1 k = 2  <<
			xx1 = p2.m_x - p4.m_x;
			yy1 = p1.m_y - p4.m_y;
			xx2 = p1.m_x - p4.m_x;
			yy2 = p2.m_y - p4.m_y;
			return (xx1 * yy1 - xx2 * yy2)/(6 * volume);
		case 3: //i = 3 j = 1 k = 2  <<
			xx1 = p1.m_x - p3.m_x;
			yy1 = p2.m_y - p1.m_y;
			xx2 = p1.m_x - p2.m_x;
			yy2 = p3.m_y - p1.m_y;
			return (xx1 * yy1 - xx2 * yy2)/(6 * volume);
	}
}

real_t CFemLocalLinear3D::getdNdKsi(const int idxN, const int idxKsi, const std::vector<real_t> ksi) {
	if (ksi.size() == 0) {
		return 0;
	}
	if (idxN == idxKsi)
		return 1;
	else
		return 0;
}

real_t CFemLocalLinear3D::getdUdX(const int element_idx, const int dim, std::vector<real_t> ksi) {
	std::vector<int> element = m_mesh->getElementByIndex(element_idx);
	real_t res = 0;
	for (int i = 0; i < element.size(); i++) {
		real_t U = m_pr->getU(element[i], dim);
		res += U * getdNdX(i, element_idx, ksi);
	}
	return res;
}

real_t CFemLocalLinear3D::getdUdY(const int element_idx, const int dim, std::vector<real_t> ksi) {
	std::vector<int> element = m_mesh->getElementByIndex(element_idx);
	real_t res = 0;
	for (int i = 0; i < element.size(); i++) {
		real_t U = m_pr->getU(element[i], dim);
		res += U * getdNdY(i, element_idx, ksi);
	}

	return res;
}

real_t CFemLocalLinear3D::getdUdZ(const int element_idx, const int dim, std::vector<real_t> ksi) {
	std::vector<int> element = m_mesh->getElementByIndex(element_idx);
	real_t res = 0;
	for (int i = 0; i < element.size(); i++) {
		real_t U = m_pr->getU(element[i], dim);
		res += U * getdNdZ(i, element_idx,ksi);
	}

	return res;
}

real_t CFemLocalLinear3D::getKK(const int idxN, const int jdxN, const int l_col, const int l_row, const int element, const std::vector<real_t> ksi) {
	//int n = 4;
	const int n = m_equations_number;
	real_t res = 0;
	if (l_col == l_row) {
		res = getdNdX(idxN,element, ksi) * getdNdX(jdxN,element, ksi) + getdNdY(idxN,element, ksi) * getdNdY(jdxN,element, ksi) + getdNdZ(idxN,element, ksi) * getdNdZ(jdxN,element, ksi);
		if (l_row < n - 1)
			res = res/m_pr->getRe();
	}
	else {
		if ((3 == l_col) && (0 == l_row))
			res = getdNdX(idxN, element, ksi);
		if ((3 == l_col) && (1 == l_row))
			res = getdNdY(idxN, element, ksi);
		if ((3 == l_col) && (2 == l_row))
			res = getdNdZ(idxN, element, ksi);
	}
	return res;
}

real_t CFemLocalLinear3D::getCC(const int idxN, const int jdxN, const int l_col, const int l_row, const int element, const std::vector<real_t> ksi) {
	//int n = 4;
	const int n = m_equations_number;
	real_t res = 0;
	if (l_row < n - 1) {
		if (idxN == jdxN)
			res = this->getVolume(element, 0)/10;
		else
			res = this->getVolume(element, 0)/20;
	}
	return res;
}

//void CFemLocalLinear3D::assembleKMatrix() {
//	//const int n = 4;
//	const int n = m_equations_number;
//	real_t cc = 0,
//		   kk = 0;
//	CGaussRule *gr = new CGaussRule(3, m_mesh_geometry_type);
//	int elementsNum = m_mesh->getElementsNumber();
//	if (elementsNum > 0) {
//		for (int i = 0; i < elementsNum; i++) {
//			std::vector<int> element = m_mesh->getElementByIndex(i);
//			for (int g_col = 0; g_col < element.size(); g_col++) {
//				for (int g_row = 0; g_row < element.size(); g_row++) {
//					for (int l_col = 0; l_col < n; l_col++) {
//						for (int l_row = 0; l_row < n; l_row++) {
//							int idx = (element[g_col] * n + l_col) * n * m_mesh->getPointsNumber() + element[g_row] * n + l_row;
//
//							if (l_col == l_row) {
//								if (l_row < n - 1) {
//									if (g_col == g_row)
//										cc = this->getVolume(i, 0)/10;
//									else
//										cc = this->getVolume(i, 0)/20;
//								}
//								for (int l = 0; l < gr->m_intpoints; l++) {
//									//m_K[idx] += getVolume(i, 0) * (getKK(g_col, g_row, l_col, l_row, i, gr->m_p[l])) * gr->m_wi[l];
//									kk += getVolume(i, 0) * (getKK(g_col, g_row, l_col, l_row, i, gr->m_p[l])) * gr->m_wi[l];
//								}
//							}
//							else {
//								for (int l = 0; l < gr->m_intpoints; l++) {
//									//m_K[idx] += getVolume(i, 0) * (getKK(g_col, g_row, l_col, l_row, i, gr->m_p[l]) * this->getN(g_row, gr->m_p[l]) * gr->m_wi[l]);
//									kk += getVolume(i, 0) * (getKK(g_col, g_row, l_col, l_row, i, gr->m_p[l]) * this->getN(g_row, gr->m_p[l]) * gr->m_wi[l]);
//								}
//							}
//
//							m_C[idx] += cc;
//							m_K[idx] += kk + cc/m_pr->getTau();
//							kk = 0;
//							cc = 0;
//						}
//					}
//				}
//			}
//		}
//	}
//	delete gr;
//}

real_t CFemLocalLinear3D::getFF(const int idxN, const int l_row, const int element, const std::vector<real_t> ksi) {
	real_t res = 0;
	real_t U1 = m_pr->getU(idxN, 0);
	real_t U2 = m_pr->getU(idxN, 1);
	real_t U3 = m_pr->getU(idxN, 2);
	switch(l_row) {
		case 0:
			res = U1 * getdUdX(element, 0, ksi) + U2 * getdUdY(element, 0, ksi) + U3 * getdUdZ(element, 0, ksi);
			break;
		case 1:
			res = U1 * getdUdX(element, 1, ksi) + U2 * getdUdY(element, 1, ksi) + U3 * getdUdZ(element, 1, ksi);
			break;
		case 2:
			res = U1 * getdUdX(element, 2, ksi) + U2 * getdUdY(element, 2, ksi) + U3 * getdUdZ(element, 2, ksi);
			break;
		case 3:
			res = -2 * (getdUdY(element, 0, ksi) * getdUdX(element, 1, ksi) + \
					getdUdZ(element, 0, ksi) * getdUdX(element, 2, ksi) + \
					getdUdZ(element, 1, ksi) * getdUdY(element, 2, ksi) );
			break;
	}
	return res;
};

//void CFemLocalLinear3D::assembleRightVector(const int timestep) {
//	//const int n = 4;
//	const int n = m_equations_number;
//	int elnumber = m_mesh->getElementsNumber();
//	int ptnumber = m_mesh->getPointsNumber();
//	CGaussRule *gr = new CGaussRule(3, m_mesh_geometry_type);
//
//	for (int i = 0; i < elnumber; i++) {
//		std::vector<int> element = m_mesh->getElementByIndex(i);
//		for (int j = 0; j < element.size(); j++) {
//			for (int k = 0; k < n; k++) {
//				//integration over RHS
//				for (int l = 0; l < gr->m_intpoints; l++) {
//					m_F[element[j] * n + k] -= getVolume(i, 0) * ( getFF(element[j], k, i, gr->m_p[l]) * getN(j, gr->m_p[l]) * gr->m_wi[l] );
//				}
//				if (k < n - 1) {
//					m_U_temp[element[j] * n + k] = m_pr->getU(element[j], k);
//				}
//			}
//		}
//	}
//
//	char *ch = "N";
//	int m_m = ptnumber * n,
//		m_n = ptnumber * n,
//		lda = ptnumber * n,
//		incx = 1,
//		incy = 1;
//	real_t tau = 1/m_pr->getTau(),
//		   beta = 1;
//
//	dgemv(ch, m_m, m_n, tau, &m_C[0], lda, &m_U_temp[0], incx, beta, &m_F[0], incy);
//
//	delete gr;
//}

//void CFemLocalLinear3D::setBorderConditions(const int timestep) {
//	//const int n = 4;
//	const int n = m_equations_number;
//	int ptnumber = m_mesh->getPointsNumber();
//	std::set<int> borderPoints = m_mesh->getBorderPoints();
//
//	for (auto const& i : borderPoints) {
//		for (int k = 0; k < n; k++) {
//			m_F[i * n + k] = m_pr->getBorderCondition(i, k, (timestep) * m_pr->getTau());
//		}
//
//		for (int k = 0; k < ptnumber; k++) {
//			for (int ii = 0; ii < n; ii++) {
//				for (int jj = 0; jj < n; jj++) {
//					int pivot = (ptnumber * n * n) * k + n * i;
//					int g_idx = pivot + jj * m_mesh->getPointsNumber() * n + ii;
//					if (k == i) {
//						if (ii == jj) {
//							m_K[g_idx] = 1;
//						}
//						else
//							m_K[g_idx] = 0;
//					}
//					else {
//						m_K[g_idx] = 0;
//					}
//				}
//			}
//		}
//	}
//}

//void CFemLocalLinear3D::perform(const int timesteps) {
//	//int n = 4;
//	const int n = m_equations_number;
//	int count = m_mesh->getPointsNumber();
//
//	this->assembleKMatrix();
//
//	//printMatrix2File("k_matrix.txt", m_K, m_F, count * n);
//	dump2binfile(m_K, count * n * count * n, K_MATRIX_FILENAME);
//	for (int step = 1; step < timesteps; step++) {
//		this->assembleRightVector(step);
//		this->setBorderConditions(step);
//		//printMatrix2File("k_matrix.txt", m_K, m_F, count * n);
//		dgesv(count * n, m_K, m_F);
//		m_pr->setU(m_F);
//		memset(m_F, 0, count * n * sizeof(real_t));
//		binfile2data(m_K, count * n * count * n, K_MATRIX_FILENAME);
//		cout << "step " << step << " / " << timesteps << " is done" << endl;
//	}
//	cout << "number of points = " << count << endl;
//	for (int i = 0; i < count; i++) {
//		cout << i * n + 0 << "=" << m_pr->getU(i, 0) << " " << m_pr->getBorderCondition(i, 0, (timesteps - 1) * m_pr->getTau()) << endl;
//		cout << i * n + 1 << "=" << m_pr->getU(i, 1) << " " << m_pr->getBorderCondition(i, 1, (timesteps - 1) * m_pr->getTau()) <<endl;
//		cout << i * n + 2 << "=" << m_pr->getU(i, 2) << " " << m_pr->getBorderCondition(i, 2, (timesteps - 1) * m_pr->getTau()) <<endl;
//		cout << i * n + 3 << "=" << m_pr->getU(i, 3) << " " << m_pr->getBorderCondition(i, 3, (timesteps - 1) * m_pr->getTau()) <<endl;
//	}
//}

//void CFemLocalLinear3D::init(CProblem *pr) {
//	//const int n = 4; //number of equations to solve
//	const int n = m_equations_number;
//	int count = m_mesh->getPointsNumber();
//	m_K = new real_t[count * n * count * n];
//	m_C = new real_t[count * n * count * n];
//	m_F = new real_t[count * n];
//	m_U_temp = new real_t[count * n];
//
//	memset(m_F, 0, sizeof(real_t) * count * n);
//	memset(m_U_temp, 0, sizeof(real_t) * count * n);
//	memset(m_K, 0, sizeof(real_t) * count * n * count * n);
//	memset(m_C, 0, sizeof(real_t) * count * n * count * n);
//
//	m_pr = pr;
//	m_pr->init();
//	cout << " TAU1 " << m_pr->getTau() << endl;
//}
