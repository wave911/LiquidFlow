//#include "Fem.h"
//#include <iostream>
//#include "f2c.h"
//#include "dgemv.h"
//#include "dgesv.h"
//
//using namespace std;
//
//CFemLocalLinear3D::CFemLocalLinear3D(CMesh *mesh) {
//	m_mesh = mesh;
//}
//
//CFemLocalLinear3D::~CFemLocalLinear3D() {
//}
//
//std::vector<real_t> CFemLocalLinear3D::getLocalCoordinates(const int element,
//											  const CPoint3D p) {
//	real_t *matrix = new real_t[16];
//	real_t *vect = new real_t[4];
//	real_t *ksi = new real_t[4];
//	vector<int> points = m_mesh->getElementByIndex(element);
//	CPoint3D p1 = m_mesh->getPointByIndex(points[0]);
//	CPoint3D p2 = m_mesh->getPointByIndex(points[1]);
//	CPoint3D p3 = m_mesh->getPointByIndex(points[2]);
//	CPoint3D p4 = m_mesh->getPointByIndex(points[3]);
//
//	real_t volume = getVolume(element, 0);
//	real_t vo1 = getVolume(element, 1);
//	real_t vo2 = getVolume(element, 2);
//	real_t vo3 = getVolume(element, 3);
//	real_t vo4 = getVolume(element, 4);
//
//	matrix[0] = 6 * vo1;
//	matrix[1] = 6 * vo2;
//	matrix[2] = 6 * vo3;
//	matrix[3] = 6 * vo4;
//
//	matrix[4] = (p4.m_y - p2.m_y) * (p3.m_z - p2.m_z) - (p3.m_y - p2.m_y) * (p4.m_z - p2.m_z);
//	matrix[5] = (p3.m_y - p1.m_y) * (p4.m_z - p3.m_z) - (p3.m_y - p4.m_y) * (p1.m_z - p3.m_z);
//	matrix[6] = (p2.m_y - p4.m_y) * (p1.m_z - p4.m_z) - (p1.m_y - p4.m_y) * (p2.m_z - p4.m_z);
//	matrix[7] = (p1.m_y - p3.m_y) * (p2.m_z - p1.m_z) - (p1.m_y - p2.m_y) * (p3.m_z - p1.m_z);
//
//	matrix[8] = (p3.m_x - p2.m_x) * (p4.m_z - p2.m_z) - (p4.m_x - p2.m_x) * (p3.m_z - p2.m_z);
//	matrix[9] = (p4.m_x - p3.m_x) * (p3.m_z - p1.m_z) - (p1.m_x - p3.m_x) * (p3.m_z - p4.m_z);
//	matrix[10] = (p1.m_x - p4.m_x) * (p2.m_z - p4.m_z) - (p2.m_x - p4.m_x) * (p1.m_z - p4.m_z);
//	matrix[11] = (p2.m_x - p1.m_x) * (p1.m_z - p3.m_z) - (p3.m_x - p1.m_x) * (p1.m_z - p2.m_z);
//
//	matrix[12] = (p4.m_x - p2.m_x) * (p3.m_y - p2.m_y) - (p3.m_x - p2.m_x) * (p4.m_y - p2.m_y);
//	matrix[13] = (p3.m_x - p1.m_x) * (p4.m_y - p3.m_y) - (p3.m_x - p4.m_x) * (p1.m_y - p3.m_y);
//	matrix[14] = (p2.m_x - p4.m_x) * (p1.m_y - p4.m_y) - (p1.m_x - p4.m_x) * (p2.m_y - p4.m_y);
//	matrix[15] = (p1.m_x - p3.m_x) * (p2.m_y - p1.m_y) - (p1.m_x - p2.m_x) * (p3.m_y - p1.m_y);
//
//	vect[0] = 1;
//	vect[1] = p.m_x;
//	vect[2] = p.m_y;
//	vect[3] = p.m_z;
//
//	char *ch = "N";
//	int m_m = 4,
//		m_n = 4,
//		lda = 4,
//		incx = 1,
//		incy = 1;
//	real_t tau = 1/(6 * volume),
//		   beta = 0;
//
//	dgemv(ch, m_m, m_n, tau, &matrix[0], lda, &vect[0],
//				 incx, beta, ksi, incy);
//
//	std::vector<real_t> res = std::vector<real_t>(ksi, ksi + 4);
//	delete(ksi);
//	delete(matrix);
//	delete(vect);
//
//	return res;
//}
//
//real_t CFemLocalLinear3D::getdNdX(const int idxN, const int element) {
//	real_t sum = 0;
//	for (int i = 0; i < 4; i++) {
//		sum += getdNdKsi(idxN, i) * getdKsidX(i, element);
//	}
//	return sum;
//}
//
//real_t CFemLocalLinear3D::getdNdY(const int idxN, const int element) {
//	real_t sum = 0;
//	for (int i = 0; i < 4; i++) {
//		sum += getdNdKsi(idxN, i) * getdKsidY(i, element);
//	}
//	return sum;
//}
//
//real_t CFemLocalLinear3D::getdNdZ(const int idxN, const int element) {
//	real_t sum = 0;
//	for (int i = 0; i < 4; i++) {
//		sum += getdNdKsi(idxN, i) * getdKsidZ(i, element);
//	}
//	return sum;
//}
//
//real_t CFemLocalLinear3D::getVolume(const int element, const int subVolume) {
//	vector<int> points = m_mesh->getElementByIndex(element);
//	CPoint3D p1 = m_mesh->getPointByIndex(points[0]);
//	CPoint3D p2 = m_mesh->getPointByIndex(points[1]);
//	CPoint3D p3 = m_mesh->getPointByIndex(points[2]);
//	CPoint3D p4 = m_mesh->getPointByIndex(points[3]);
//	real_t vo1 = (p2.m_x * (p3.m_y * p4.m_z - p4.m_y * p3.m_z) + \
//				  p3.m_x * (p4.m_y * p2.m_z - p2.m_y * p4.m_z) + \
//				  p4.m_x * (p2.m_y * p3.m_z - p3.m_y * p2.m_z))/6;
//
//	real_t vo2 = (p1.m_x * (p4.m_y * p3.m_z - p3.m_y * p4.m_z) + \
//				  p3.m_x * (p1.m_y * p4.m_z - p4.m_y * p1.m_z) + \
//				  p4.m_x * (p3.m_y * p1.m_z - p1.m_y * p3.m_z))/6;
//
//	real_t vo3 = (p1.m_x * (p2.m_y * p4.m_z - p4.m_y * p2.m_z) + \
//				  p2.m_x * (p4.m_y * p1.m_z - p1.m_y * p4.m_z) + \
//				  p4.m_x * (p1.m_y * p2.m_z - p2.m_y * p1.m_z))/6;
//
//	real_t vo4 = (p1.m_x * (p3.m_y * p2.m_z - p2.m_y * p3.m_z) + \
//				  p2.m_x * (p1.m_y * p3.m_z - p3.m_y * p1.m_z) + \
//				  p3.m_x * (p2.m_y * p1.m_z - p1.m_y * p2.m_z))/6;
//
//	switch(subVolume) {
//		case 0:
//			return vo1 + vo2 + vo3 + vo4;
//		case 1:
//			return vo1;
//		case 2:
//			return vo2;
//		case 3:
//			return vo3;
//		case 4:
//			return vo4;
//	}
//}
//
//real_t CFemLocalLinear3D::getdKsidX(const int idx, const int element) {
//	vector<int> points = m_mesh->getElementByIndex(element);
//	CPoint3D p1 = m_mesh->getPointByIndex(points[0]);
//	CPoint3D p2 = m_mesh->getPointByIndex(points[1]);
//	CPoint3D p3 = m_mesh->getPointByIndex(points[2]);
//	CPoint3D p4 = m_mesh->getPointByIndex(points[3]);
//	real_t volume = this->getVolume(element, 0);
//
//	real_t yy1, yy2, zz1, zz2;
//	switch(idx) {
//		case 0: //i = 1 j = 2 k = 3  <<
//			yy1 = p4.m_y - p2.m_y;
//			zz1 = p3.m_z - p2.m_z;
//			yy2 = p3.m_y - p2.m_y;
//			zz2 = p4.m_z - p2.m_z;
//			return (yy1 * zz1 - yy2 * zz2)/(6 * volume);
//
//		case 1: //i = 2 j = 3 k = 1  <<
//			yy1 = p3.m_y - p1.m_y;
//			zz1 = p4.m_z - p3.m_z;
//			yy2 = p3.m_y - p4.m_y;
//			zz2 = p1.m_z - p3.m_z;
//			return (yy1 * zz1 - yy2 * zz2)/(6 * volume);
//
//		case 2: //i = 3 j = 1 k = 2  <<
//			yy1 = p2.m_y - p4.m_y;
//			zz1 = p1.m_z - p4.m_z;
//			yy2 = p1.m_y - p4.m_y;
//			zz2 = p2.m_z - p4.m_z;
//			return (yy1 * zz1 - yy2 * zz2)/(6 * volume);
//
//		case 3: //i = 3 j = 1 k = 2  <<
//			yy1 = p1.m_y - p3.m_y;
//			zz1 = p2.m_z - p1.m_z;
//			yy2 = p1.m_y - p2.m_y;
//			zz2 = p3.m_z - p1.m_z;
//			return (yy1 * zz1 - yy2 * zz2)/(6 * volume);
//	}
//}
//
//real_t CFemLocalLinear3D::getdKsidY(const int idx, const int element) {
//	vector<int> points = m_mesh->getElementByIndex(element);
//	CPoint3D p1 = m_mesh->getPointByIndex(points[0]);
//	CPoint3D p2 = m_mesh->getPointByIndex(points[1]);
//	CPoint3D p3 = m_mesh->getPointByIndex(points[2]);
//	CPoint3D p4 = m_mesh->getPointByIndex(points[3]);
//	real_t volume = this->getVolume(element, 0);
//
//	real_t xx1, xx2, zz1, zz2;
//	switch(idx) {
//		case 0: //i = 1 j = 2 k = 3  <<
//			zz1 = p4.m_z - p2.m_z;
//			xx1 = p3.m_x - p2.m_x;
//			zz2 = p3.m_z - p2.m_z;
//			xx2 = p4.m_x - p2.m_x;
//			return (xx1 * zz1 - zz2 * xx2)/(6 * volume);
//		case 1: //i = 2 j = 3 k = 1  <<
//			zz1 = p3.m_z - p1.m_z;
//			xx1 = p4.m_x - p3.m_x;
//			zz2 = p3.m_z - p4.m_z;
//			xx2 = p1.m_x - p3.m_x;
//			return (zz1 * xx1 - zz2 * xx2)/(6 * volume);
//		case 2: //i = 3 j = 1 k = 2  <<
//			zz1 = p2.m_z - p4.m_z;
//			xx1 = p1.m_x - p4.m_x;
//			zz2 = p1.m_z - p4.m_z;
//			xx2 = p2.m_x - p4.m_x;
//			return (zz1 * xx1 - zz2 * xx2)/(6 * volume);
//		case 3: //i = 3 j = 1 k = 2  <<
//			zz1 = p1.m_z - p3.m_z;
//			xx1 = p2.m_x - p1.m_x;
//			zz2 = p1.m_z - p2.m_z;
//			xx2 = p3.m_x - p1.m_x;
//			return (zz1 * xx1 - zz2 * xx2)/(6 * volume);
//	}
//}
//
//real_t CFemLocalLinear3D::getdKsidZ(const int idx, const int element) {
//	vector<int> points = m_mesh->getElementByIndex(element);
//	CPoint3D p1 = m_mesh->getPointByIndex(points[0]);
//	CPoint3D p2 = m_mesh->getPointByIndex(points[1]);
//	CPoint3D p3 = m_mesh->getPointByIndex(points[2]);
//	CPoint3D p4 = m_mesh->getPointByIndex(points[3]);
//	real_t volume = this->getVolume(element, 0);
//
//	real_t xx1, xx2, yy1, yy2;
//	switch(idx) {
//		case 0: //i = 1 j = 2 k = 3  <<
//			xx1 = p4.m_x - p2.m_x;
//			yy1 = p3.m_y - p2.m_y;
//			xx2 = p3.m_x - p2.m_x;
//			yy2 = p4.m_y - p2.m_y;
//			return (xx1 * yy1 - xx2 * yy2)/(6 * volume);
//		case 1: //i = 2 j = 3 k = 1  <<
//			xx1 = p3.m_x - p1.m_x;
//			yy1 = p4.m_y - p3.m_y;
//			xx2 = p3.m_x - p4.m_x;
//			yy2 = p1.m_y - p3.m_y;
//			return (xx1 * yy1 - xx2 * yy2)/(6 * volume);
//		case 2: //i = 3 j = 1 k = 2  <<
//			xx1 = p2.m_x - p4.m_x;
//			yy1 = p1.m_y - p4.m_y;
//			xx2 = p1.m_x - p4.m_x;
//			yy2 = p2.m_y - p4.m_y;
//			return (xx1 * yy1 - xx2 * yy2)/(6 * volume);
//		case 3: //i = 3 j = 1 k = 2  <<
//			xx1 = p1.m_x - p3.m_x;
//			yy1 = p2.m_y - p1.m_y;
//			xx2 = p1.m_x - p2.m_x;
//			yy2 = p3.m_y - p1.m_y;
//			return (xx1 * yy1 - xx2 * yy2)/(6 * volume);
//	}
//}
//
//real_t CFemLocalLinear3D::getdNdKsi(const int idxN, const int idxKsi) {
//	if (idxN == idxKsi)
//		return 1;
//	else
//		return 0;
//}
//
//real_t CFemLocalLinear3D::getdUdX(const int element_idx, const int dim) {
//	std::vector<int> element = m_mesh->getElementByIndex(element_idx);
//	real_t res = 0;
//	for (int i = 0; i < element.size(); i++) {
//		real_t U = m_pr->getU(element[i], dim);
//		res += U * getdNdX(i, element_idx);
//	}
//
//	return res;
//}
//
//real_t CFemLocalLinear3D::getdUdY(const int element_idx, const int dim) {
//	std::vector<int> element = m_mesh->getElementByIndex(element_idx);
//	real_t res = 0;
//	for (int i = 0; i < element.size(); i++) {
//		real_t U = m_pr->getU(element[i], dim);
//		res += U * getdNdY(i, element_idx);
//	}
//
//	return res;
//}
//
//real_t CFemLocalLinear3D::getdUdZ(const int element_idx, const int dim) {
//	std::vector<int> element = m_mesh->getElementByIndex(element_idx);
//	real_t res = 0;
//	for (int i = 0; i < element.size(); i++) {
//		real_t U = m_pr->getU(element[i], dim);
//		res += U * getdNdZ(i, element_idx);
//	}
//
//	return res;
//}
//
//void CFemLocalLinear3D::assembleKMatrix() {
//	const int n = 4;
//	real_t cc = 0,
//		   kk = 0;
//
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
//									cc = integrateiNjN(g_col, g_row, i);
//								}
//								kk = getdNdX(g_col,i) * getdNdX(g_row,i) + getdNdY(g_col,i) * getdNdY(g_row,i) + getdNdZ(g_col,i) * getdNdZ(g_row,i);
//								kk = kk * integrateidNjdN(g_row, g_col, i);
//								if (l_row < n - 1)
//									kk = kk/m_pr->getRe();
//							}
//							else {
//								if ((n - 1 == l_col) && (0 == l_row))
//									kk = getdNdX(g_col, i) * integrateiNjdN(g_row, g_col, i);
//								if ((n - 1 == l_col) && (1 == l_row))
//									kk = getdNdY(g_col, i) * integrateiNjdN(g_row, g_col, i);
//								if ((n - 1 == l_col) && (2 == l_row))
//									kk = getdNdZ(g_col, i) * integrateiNjdN(g_row, g_col, i);
//								cc = 0;
//							}
//
//							m_K[idx] += kk + cc/m_pr->getTau();
//							m_C[idx] += cc;
//							kk = 0;
//							cc = 0;
//						}
//					}
//				}
//			}
//		}
//	}
//}
//
//void CFemLocalLinear3D::assembleRightVector(const int timestep) {
//	const int n = 4;
//	int elnumber = m_mesh->getElementsNumber();
//	int ptnumber = m_mesh->getPointsNumber();
//
//	for (int i = 0; i < elnumber; i++) {
//		std::vector<int> element = m_mesh->getElementByIndex(i);
//		for (int j = 0; j < element.size(); j++) {
//			real_t U1 = m_pr->getU(element[j], 0);
//			real_t U2 = m_pr->getU(element[j], 1);
//			real_t U3 = m_pr->getU(element[j], 2);
//
//			m_F[element[j] * n + 0] += (U1 * getdUdX(i, 0) + U2 * getdUdY(i, 0) + U3 * getdUdZ(i, 0)) * integrateiNjdN(0, j, i);
//			m_F[element[j] * n + 1] += (U1 * getdUdX(i, 1) + U2 * getdUdY(i, 1) + U3 * getdUdZ(i, 1)) * integrateiNjdN(1, j, i);
//			m_F[element[j] * n + 2] += (U1 * getdUdX(i, 2) + U2 * getdUdY(i, 2) + U3 * getdUdZ(i, 2)) * integrateiNjdN(2, j, i);
//			m_F[element[j] * n + 3] += -2 * (getdUdY(i, 0) * getdUdX(i, 1) + getdUdZ(i, 0) * getdUdX(i, 2) + getdUdZ(i, 1) * getdUdY(i, 2) )  * integrateiNjdN(3, j, i);
//
//			m_U_temp[element[j] * n + 0] = U1;
//			m_U_temp[element[j] * n + 1] = U2;
//			m_U_temp[element[j] * n + 2] = U3;
//			m_U_temp[element[j] * n + 3] = 0;
//		}
//	}
//
//	char *ch = "N";
//	int m_m = m_mesh->getPointsNumber() * n,
//		m_n = m_mesh->getPointsNumber() * n,
//		lda = m_mesh->getPointsNumber() * n,
//		incx = 1,
//		incy = 1;
//	real_t tau = 1/m_pr->getTau(),
//		   beta = 1;
//
//	dgemv(ch, m_m, m_n, tau, &m_C[0], lda, &m_U_temp[0], incx, beta, &m_F[0], incy);
//}
//
//void CFemLocalLinear3D::setBorderConditions(const int timestep) {
//	const int n = 4;
//	int ptnumber = m_mesh->getPointsNumber();
//	std::set<int> borderPoints = m_mesh->getBorderPoints();
//
//	for (auto const& i : borderPoints) {
//		m_F[i * n + 0] = m_pr->getBorderCondition(i, 0, (timestep) * m_pr->getTau());
//		m_F[i * n + 1] = m_pr->getBorderCondition(i, 1, (timestep) * m_pr->getTau());
//		m_F[i * n + 2] = m_pr->getBorderCondition(i, 2, (timestep) * m_pr->getTau());
//		m_F[i * n + 3] = m_pr->getBorderCondition(i, 3, (timestep) * m_pr->getTau());
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
//
//void CFemLocalLinear3D::perform(const int timesteps) {
//	int n = 4;
//	int count = m_mesh->getPointsNumber();
//
//	this->assembleKMatrix();
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
//	//}
//	cout << "number of points = " << count << endl;
//	for (int i = 0; i < count; i++) {
//		cout << i * n + 0 << "=" << m_pr->getU(i, 0) << " " << m_pr->getBorderCondition(i, 0, (timesteps - 1) * m_pr->getTau()) << endl;
//		cout << i * n + 1 << "=" << m_pr->getU(i, 1) << " " << m_pr->getBorderCondition(i, 1, (timesteps - 1) * m_pr->getTau()) <<endl;
//		cout << i * n + 2 << "=" << m_pr->getU(i, 2) << " " << m_pr->getBorderCondition(i, 2, (timesteps - 1) * m_pr->getTau()) <<endl;
//		cout << i * n + 3 << "=" << m_pr->getU(i, 3) << " " << m_pr->getBorderCondition(i, 3, (timesteps - 1) * m_pr->getTau()) <<endl;
//
////		cout << i * n + 0 << "=" << m_pr->getU(i, 0) - m_pr->getBorderCondition(i, 0, (timesteps - 1) * m_pr->getTau()) << endl;
////		cout << i * n + 1 << "=" << m_pr->getU(i, 1) - m_pr->getBorderCondition(i, 1, (timesteps - 1) * m_pr->getTau()) <<endl;
////		cout << i * n + 2 << "=" << m_pr->getU(i, 2) - m_pr->getBorderCondition(i, 2, (timesteps - 1) * m_pr->getTau()) <<endl;
////		cout << i * n + 3 << "=" << m_pr->getU(i, 3) - m_pr->getBorderCondition(i, 3, (timesteps - 1) * m_pr->getTau()) <<endl;
//	}
//	//
//	}
//}
//
//void CFemLocalLinear3D::init(CProblem *pr) {
//	const int n = 4; //number of equations to solve
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
//
//real_t CFemLocalLinear3D::integrateiNjN(const int iN, const int jN, const int element) {
//	if (iN == jN) {
//		return getVolume(element, 0)/10;
//	}
//	else {
//		return getVolume(element, 0)/20;
//	}
//}
//real_t CFemLocalLinear3D::integrateiNjdN(const int iN, const int jN, const int element) {
//	return getVolume(element, 0)/4;
//}
//real_t CFemLocalLinear3D::integrateidNjdN(const int iN, const int jN, const int element) {
//	return getVolume(element, 0);
//}
