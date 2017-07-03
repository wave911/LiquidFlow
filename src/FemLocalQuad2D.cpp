#include "Fem.h"
#include <iostream>
#include "f2c.h"
#include "dgemv.h"
#include "dgesv.h"

using namespace std;

CFemLocalQuad2D::CFemLocalQuad2D(CMesh *mesh) : CFemLocalLinear2D(mesh) {
	m_mesh = mesh;
	m_K = nullptr;
	m_C = nullptr;
	m_F = nullptr;
	m_U_temp = nullptr;
	m_pr = nullptr;
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

	real_t square = getSquare(element);
	matrix[0] = p2.m_x * p3.m_y - p3.m_x * p2.m_y;//2 * square/3;
	matrix[1] = p3.m_x * p1.m_y - p1.m_x * p3.m_y;//2 * square/3;
	matrix[2] = p1.m_x * p2.m_y - p2.m_x * p1.m_y;//2 * square/3;
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

real_t CFemLocalQuad2D::getdNdX(const int idxN, const int element) {
	real_t sum = 0;
	std::vector<int> e = m_mesh->getElementByIndex(element);
	CPoint3D p = m_mesh->getPointByIndex(e[idxN]);
	std::vector<real_t> ksi = this->getLocalCoordinates(element, p);

	for (int i = 0; i < 3; i++) {
		sum += getdNdKsi(idxN, i, ksi) * getdKsidX(i, element);
	}
	return sum;
}

real_t CFemLocalQuad2D::getdNdY(const int idxN, const int element) {
	real_t sum = 0;
	std::vector<int> e = m_mesh->getElementByIndex(element);
	CPoint3D p = m_mesh->getPointByIndex(e[idxN]);
	std::vector<real_t> ksi = this->getLocalCoordinates(element, p);

	for (int i = 0; i < 3; i++) {
		sum += getdNdKsi(idxN, i, ksi) * getdKsidY(i, element);
	}
	return sum;
}

real_t CFemLocalQuad2D::getdNdKsi(const int idxN, const int idxKsi, const std::vector<real_t> ksi) {
	switch(idxN) {
		case 0:
			if (idxKsi == 0)
				return 4 * ksi[0] - 1;
			else
				return 0;
			break;
		case 1:
			if (idxKsi == 1)
				return 4 * ksi[1] - 1;
			else
				return 0;
			break;
		case 2:
			if (idxKsi == 2)
				return 4 * ksi[2] - 1;
			else
				return 0;
			break;
		case 3:
			if (idxKsi == 0)
				return 4 * ksi[1];
			else if (idxKsi == 1)
				return 4 * ksi[0];
			else if (idxKsi == 2)
				return 0;
			break;
		case 4:
			if (idxKsi == 0)
				return 0;
			else if (idxKsi == 1)
				return 4 * ksi[2];
			else if (idxKsi == 2)
				return 4 * ksi[1];
			break;
		case 5:
			if (idxKsi == 0)
				return 4 * ksi[2];
			else if (idxKsi == 1)
				return 0;
			else if (idxKsi == 2)
				return 4 * ksi[0];
			break;
	}
}

void CFemLocalQuad2D::assembleKMatrix() {
	const int n = 3;
	real_t cc = 0,
		   kk = 0;

	int elementsNum = m_mesh->getElementsNumber();
	if (elementsNum > 0) {
		for (int i = 0; i < elementsNum; i++) {
			std::vector<int> element = m_mesh->getElementByIndex(i);
			for (int g_col = 0; g_col < element.size(); g_col++) {
				for (int g_row = 0; g_row < element.size(); g_row++) {
					for (int l_col = 0; l_col < n; l_col++) {
						for (int l_row = 0; l_row < n; l_row++) {
							int idx = (element[g_col] * n + l_col) * n * m_mesh->getPointsNumber() + element[g_row] * n + l_row;

							if (l_col == l_row) {
								if (l_row < n - 1) {
									if (g_col == g_row)
										cc = getSquare(i)/6;
									else
										cc = getSquare(i)/12;
								}
								kk = getdNdX(g_col,i) * getdNdX(g_row,i) + getdNdY(g_col,i) * getdNdY(g_row,i);
								kk = kk * getSquare(i);
								if (l_row < n - 1)
									kk = kk/m_pr->getRe();
							}
							else {
								if ((2 == l_col) && (0 == l_row))
									kk = getdNdX(g_col, i) * getSquare(i)/6;
								if ((2 == l_col) && (1 == l_row))
									kk = getdNdY(g_col, i) * getSquare(i)/6;
								cc = 0;
							}

							m_K[idx] += kk + cc/m_pr->getTau();
							m_C[idx] += cc;
							kk = 0;
							cc = 0;
						}
					}
				}
			}
		}
	}
}

void CFemLocalQuad2D::assembleRightVector(const int timestep) {
	const int n = 3;
	int elnumber = m_mesh->getElementsNumber();
	int ptnumber = m_mesh->getPointsNumber();

	for (int i = 0; i < elnumber; i++) {
		std::vector<int> element = m_mesh->getElementByIndex(i);
		for (int j = 0; j < element.size(); j++) {
			real_t U1 = m_pr->getU(element[j], 0);
			real_t U2 = m_pr->getU(element[j], 1);

			m_F[element[j] * n + 0] = getdUdX(i, 0);
			m_F[element[j] * n + 0] += (U1 * getdUdX(i, 0) + U2 * getdUdY(i, 0)) * (this->getSquare(i)/6);
			m_F[element[j] * n + 1] += (U1 * getdUdX(i, 1) + U2 * getdUdY(i, 1)) * (this->getSquare(i)/6);
			m_F[element[j] * n + 2] += -(2 * getdUdY(i, 0) * getdUdX(i, 1)) * (this->getSquare(i)/3);

			m_U_temp[element[j] * n + 0] = m_pr->getU(element[j], 0);
			m_U_temp[element[j] * n + 1] = m_pr->getU(element[j], 1);
			m_U_temp[element[j] * n + 2] = 0;
		}
	}

	char *ch = "N";
	int m_m = m_mesh->getPointsNumber() * n,
		m_n = m_mesh->getPointsNumber() * n,
		lda = m_mesh->getPointsNumber() * n,
		incx = 1,
		incy = 1;
	real_t tau = 1/m_pr->getTau(),
		   beta = 1;

	dgemv(ch, m_m, m_n, tau, &m_C[0], lda, &m_U_temp[0], incx, beta, &m_F[0], incy);
}
