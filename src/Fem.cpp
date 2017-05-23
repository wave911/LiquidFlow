#include "Fem.h"
#include <iostream>
#include "f2c.h"
#include "dgemv.h"
using namespace std;

void CFem::dgemv(char *trans, int m, int n, real_t alpha, real_t *a, 
				   int lda, real_t *x, int incx, real_t beta, real_t *y, int incy) {
#ifdef LIBBLASLAPACK
	dgemv_(trans, (integer*)&m, (integer*)&n, (doublereal*)&alpha, 
				(doublereal*)a, (integer*)&lda, (doublereal*)x, 
				(integer*)&incx, (doublereal*)&beta, (doublereal*)y, 
				(integer*)&incy);
#endif
}

CFemLocalLinear2D::CFemLocalLinear2D(CMesh *mesh) {
	m_mesh = mesh;
}

CFemLocalLinear2D::~CFemLocalLinear2D() {
	delete [] m_K;
	delete [] m_C;
	delete [] m_F;
}

void CFemLocalLinear2D::init(CProblem *pr) {
	const int n = 3; //number of equations to solve
	int count = m_mesh->getPointsNumber();
	m_K = new real_t[count * n * count * n];
	m_C = new real_t[count * n * count * n];
	m_F = new real_t[count * n];
	m_pr = pr;
}

std::vector<real_t> CFemLocalLinear2D::getLocalCoordinates(const int element, 
											  const CPoint3D p) {
	real_t *matrix = new real_t[9];
	real_t *vect = new real_t[3];
	real_t *ksi = new real_t[3];
	vector<int> points = m_mesh->getElementByIndex(element);
	CPoint3D p1 = m_mesh->getPointByIndex(points[0]);
	CPoint3D p2 = m_mesh->getPointByIndex(points[1]);
	CPoint3D p3 = m_mesh->getPointByIndex(points[2]);

	cout << "triangle" << endl;
	cout << p1.m_x << " " << p1.m_y << " " << p1.m_z << endl;
	cout << p2.m_x << " " << p2.m_y << " " << p2.m_z << endl;
	cout << p3.m_x << " " << p3.m_y << " " << p3.m_z << endl;

	real_t square = getSquare(element);
	cout << "square " << square << endl;  
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

real_t CFemLocalLinear2D::getdNdX(const int idxN, const int element) {
	real_t sum = 0;
	for (int i = 0; int < 3; i++) {
		sum += getdNdKsi(idxN, i) * getdKsidX(i, element);
	}
	return sum;
}

real_t CFemLocalLinear2D::getdNdY(const int idxN, const int element) {
	real_t sum = 0;
	for (int i = 0; int < 3; i++) {
		sum += getdNdKsi(idxN, i) * getdKsidY(i, element);
	}
	return sum;
}

real_t CFemLocalLinear2D::getSquare(const int element) {
	vector<int> points = m_mesh->getElementByIndex(element);
	CPoint3D p1 = m_mesh->getPointByIndex(points[0]);
	CPoint3D p2 = m_mesh->getPointByIndex(points[1]);
	CPoint3D p3 = m_mesh->getPointByIndex(points[2]);

	real_t square = ((p2.m_x * p3.m_y - p3.m_x * p2.m_y) + \
					 (p3.m_x * p1.m_y - p1.m_x * p3.m_y) + \
					 (p1.m_x * p2.m_y - p2.m_x * p1.m_y))/2;
	return square;
}

real_t CFemLocalLinear2D::getdKsidX(const int idx, const int element) {
	vector<int> points = m_mesh->getElementByIndex(element);
	CPoint3D p1 = m_mesh->getPointByIndex(points[0]);
	CPoint3D p2 = m_mesh->getPointByIndex(points[1]);
	CPoint3D p3 = m_mesh->getPointByIndex(points[2]);
	real_t square = getSquare(element);

	switch(idx) {
		case 0: //i = 1 j = 2 k = 3  <<
			return (p2.m_y - p3.m_y)/(2 * square);
		case 1: //i = 2 j = 3 k = 1  <<
			return (p3.m_y - p1.m_y)/(2 * square);
		case 2: //i = 3 j = 1 k = 2  <<
			return (p1.m_y - p2.m_y)/(2 * square);						
	}
}

real_t CFemLocalLinear2D::getdKsidY(const int idx, const int element) {
	vector<int> points = m_mesh->getElementByIndex(element);
	CPoint3D p1 = m_mesh->getPointByIndex(points[0]);
	CPoint3D p2 = m_mesh->getPointByIndex(points[1]);
	CPoint3D p3 = m_mesh->getPointByIndex(points[2]);
	real_t square = getSquare(element);

	switch(idx) {
		case 0: //i = 1 j = 2 k = 3  <<
			return (p3.m_x - p2.m_x)/(2 * square);
		case 1: //i = 2 j = 3 k = 1  <<
			return (p1.m_x - p3.m_x)/(2 * square);
		case 2: //i = 3 j = 1 k = 2  <<
			return (p2.m_x - p1.m_x)/(2 * square);						
	}
}

real_t CFemLocalLinear2D::getdNdKsi(const int idxN, const int idxKsi) {
	if (idxN == idxKsi)
		return 1;
	else
		return 0;
}

real_t CFemLocalLinear2D::getdUdX(const int element, const int dim) {
	std::vector<int> element = m_mesh->getElementByIndex(i);
	real_t res = 0;
	for (int i = 0; i < element.size(); i++) {
		real_t U = m_pr->getU(element[i], dim);
		res += U * getdNdX(i, element);
	}

	return res;
}

real_t CFemLocalLinear2D::getdUdY(const int element, const int dim) {
	std::vector<int> element = m_mesh->getElementByIndex(i);
	real_t res = 0;
	for (int i = 0; i < element.size(); i++) {
		real_t U = m_pr->getU(element[i], dim);
		res += U * getdNdY(i, element);
	}

	return res;
}

void CFemLocalLinear2D::assembleKMatrix() {
	const int n = 3;
	real_t cc = 0,
		   kk = 0;

	int elementsNum = m_mesh->getElementsNumber();
	if (elementsNum > 0) {
		for (int i = 0; i < elementsNum; i++) {
			std::vector<int> element = m_mesh->getElementByIndex(i);
			for (int g_col = 0; g_col < element.size(); g_col++) {
				for (int g_row = 0; g_row < element.size(); g_row++) {
					int pivot = (m_mesh->getPointsNumber() * n * n) * element[g_col] + n * element[g_row];
					for (int l_col = 0; l_col < n; l_col++) {
						for (int l_row = 0; l_row < n; l_row++) {
							if (l_col == g_row) {
								if (g_col == g_row)
									cc = getSquare(i)/6;
								else
									cc = getSquare(i)/12;						
								kk = getdNdX(g_col,i) * getdNdX(g_row,i) + getdNdY(g_col,i) * getdNdY(g_row,i);
								kk = kk * getSquare(i)/12;
								if (g_row < n - 1)
									kk = kk/m_pr->getRe();

							} 
							else {
								if ((2 == l_col) && (0 == l_row))
									kk = getdNdX(g_col, i) * getSquare(i)/3;
								if ((2 == l_col) && (1 == l_row))
									kk = getdNdY(g_col, i) * getSquare(i)/3;
							}

							int g_idx = pivot + l_col * m_mesh->getPointsNumber() * n + l_row;
							m_K[g_idx] += (kk + cc)/m_pr->getTau();
							m_C[g_idx] += cc;
						}
					}
				}
			}
		}
	}
}

void CFemLocalLinear2D::assembleRightVector() {
	const int n = 3;
	real_t ff1 = 0, ff2 = 0, ff3 = 0;
	int elnumber = m_mesh->getElementsNumber();
	int ptnumber = m_mesh->getPointsNumber();

	for (int i = 0; i < elnumber; i++) {
		std::list<int> element = m_mesh->getElementByIndex(i);
		for (int j = 0; j < element.size(); j++) {
			if (!m_mesh->isBorderPoint(element[j])) {
				real_t U1 = m_pr->getU(i, 0);
				real_t U2 = m_pr->getU(i, 1);
				m_F[element[j] * ptnumber + 0] += (U1 * getdUdX(i, 0) + U2 * getdUdY(i, 0))/(this->getSquare(i) * 3);
				m_F[element[j] * ptnumber + 1] += (U1 * getdUdX(i, 1) + U2 * getdUdY(i, 1))/(this->getSquare(i) * 3);
				m_F[element[j] * ptnumber + 2] += 0;
			}
			else {
				m_F[element[j] * ptnumber + 0] = 1;
				m_F[element[j] * ptnumber + 1] = 1;
				m_F[element[j] * ptnumber + 2] = 1;

				for (int k = 0; k < ptnumber; k++) {
					for (int ii = 0; ii < n; ii++) {
						for (int jj = 0; jj < n; jj++) {
							int pivot = (ptnumber * n * n) * k + n * element[j];
							int g_idx = pivot + jj * m_mesh->getPointsNumber() * n + ii;
							if (k == element[j]) {
								if (ii == jj) {
									m_K[g_idx] = 1;
								}
								else
									m_K[g_idx] = 0;
							}
							else {
								m_K[g_idx] = 0;
							}
						}
					}
				}

			}
		}
	}

}

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

	cout << "tetrahedra" << endl;
	cout << p1.m_x << " " << p1.m_y << " " << p1.m_z << endl;
	cout << p2.m_x << " " << p2.m_y << " " << p2.m_z << endl;
	cout << p3.m_x << " " << p3.m_y << " " << p3.m_z << endl;
	cout << p4.m_x << " " << p4.m_y << " " << p4.m_z << endl;

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
