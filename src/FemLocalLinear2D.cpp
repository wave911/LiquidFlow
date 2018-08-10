#include "Fem.h"
#include <iostream>
#include "f2c.h"
#include "dgemv.h"
#include "dgesv.h"

using namespace std;

CFemLocalLinear2D::CFemLocalLinear2D(CMesh *mesh) {
	m_mesh = mesh;
	m_K = nullptr;
	m_C = nullptr;
	m_F = nullptr;
	m_U_temp = nullptr;
	m_pr = nullptr;

}

CFemLocalLinear2D::~CFemLocalLinear2D() {
	if (this->m_K != nullptr)
		delete [] this->m_K;
	if (this->m_C != nullptr)
		delete [] this->m_C;
	if (this->m_F != nullptr)
		delete [] this->m_F;
	if (this->m_U_temp != nullptr)
		delete [] this->m_U_temp;
}

void CFemLocalLinear2D::init(CProblem *pr) {
	const int n = 3; //number of equations to solve
	int count = m_mesh->getPointsNumber();
	m_K = new real_t[count * n * count * n];
	m_C = new real_t[count * n * count * n];
	m_F = new real_t[count * n];
	m_U_temp = new real_t[count * n];

	memset(m_F, 0, sizeof(real_t) * count * n);
	memset(m_U_temp, 0, sizeof(real_t) * count * n);
	memset(m_K, 0, sizeof(real_t) * count * n * count * n);
	memset(m_C, 0, sizeof(real_t) * count * n * count * n);

	m_pr = pr;
	m_pr->init();
	cout << " TAU1 " << m_pr->getTau() << endl;
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

real_t CFemLocalLinear2D::getdNdX(const int idxN, const int element) {
	real_t sum = 0;
	for (int i = 0; i < 3; i++) {
		sum += getdNdKsi(idxN, i) * getdKsidX(i, element);
	}
	return sum;
}

real_t CFemLocalLinear2D::getdNdY(const int idxN, const int element) {
	real_t sum = 0;
	for (int i = 0; i < 3; i++) {
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

real_t CFemLocalLinear2D::getdUdX(const int element_idx, const int dim) {
	std::vector<int> element = m_mesh->getElementByIndex(element_idx);
	real_t res = 0;
	for (int i = 0; i < element.size(); i++) {
		real_t U = m_pr->getU(element[i], dim);
		res += U * getdNdX(i, element_idx);
	}

	return res;
}

real_t CFemLocalLinear2D::getdUdY(const int element_idx, const int dim) {
	std::vector<int> element = m_mesh->getElementByIndex(element_idx);
	real_t res = 0;
	for (int i = 0; i < element.size(); i++) {
		real_t U = m_pr->getU(element[i], dim);
		res += U * getdNdY(i, element_idx);
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
					for (int l_col = 0; l_col < n; l_col++) {
						for (int l_row = 0; l_row < n; l_row++) {
							int idx = (element[g_col] * n + l_col) * n * m_mesh->getPointsNumber() + element[g_row] * n + l_row;

							if (l_col == l_row) {
								if (l_row < n - 1) {
									cc = integrateiNjN(g_col, g_row, i);
//									if (g_col == g_row)
//										cc = getSquare(i)/6;
//									else
//										cc = getSquare(i)/12;
								}
								kk = getdNdX(g_col,i) * getdNdX(g_row,i) + getdNdY(g_col,i) * getdNdY(g_row,i);
								kk = kk * integrateidNjdN(g_row, g_col, i); //getSquare(i);
								if (l_row < n - 1)
									kk = kk/m_pr->getRe();
							}
							else {
								if ((2 == l_col) && (0 == l_row))
									//kk = getdNdX(g_col, i) * getSquare(i)/3;
									kk = getdNdX(g_col, i) * integrateiNjdN(g_row, g_col, i);
								if ((2 == l_col) && (1 == l_row))
									//kk = getdNdY(g_col, i) * getSquare(i)/3;
									kk = getdNdY(g_col, i) * integrateiNjdN(g_row, g_col, i);
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

void CFemLocalLinear2D::assembleRightVector(const int timestep) {
	const int n = 3;
	int elnumber = m_mesh->getElementsNumber();
	int ptnumber = m_mesh->getPointsNumber();

	for (int i = 0; i < elnumber; i++) {
		std::vector<int> element = m_mesh->getElementByIndex(i);
		for (int j = 0; j < element.size(); j++) {
			real_t U1 = m_pr->getU(element[j], 0);
			real_t U2 = m_pr->getU(element[j], 1);

			m_F[element[j] * n + 0] += (U1 * getdUdX(i, 0) + U2 * getdUdY(i, 0)) * integrateiNjdN(0, j, i);//(this->getSquare(i)/3);
			m_F[element[j] * n + 1] += (U1 * getdUdX(i, 1) + U2 * getdUdY(i, 1)) * integrateiNjdN(1, j, i);//(this->getSquare(i)/3);
			m_F[element[j] * n + 2] += -(2 * getdUdY(i, 0) * getdUdX(i, 1)) * integrateiNjdN(2, j, i);//(this->getSquare(i)/3);

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

void CFemLocalLinear2D::setBorderConditions(const int timestep) {
	const int n = 3;
	int ptnumber = m_mesh->getPointsNumber();
	std::set<int> borderPoints = m_mesh->getBorderPoints();

	for (auto const& i : borderPoints) {
		m_F[i * n + 0] = m_pr->getBorderCondition(i, 0, 0);
		m_F[i * n + 1] = m_pr->getBorderCondition(i, 1, 0);
		m_F[i * n + 2] = m_pr->getBorderCondition(i, 2, 0);

		for (int k = 0; k < ptnumber; k++) {
			for (int ii = 0; ii < n; ii++) {
				for (int jj = 0; jj < n; jj++) {
					int pivot = (ptnumber * n * n) * k + n * i;
					int g_idx = pivot + jj * m_mesh->getPointsNumber() * n + ii;
					if (k == i) {
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

void CFemLocalLinear2D::perform(const int timesteps) {
	int n = 3;
	int count = m_mesh->getPointsNumber();

	this->assembleKMatrix();
	//printMatrix2File("k_matrix.txt", m_K, m_F, count * n);
	dump2binfile(m_K, count * n * count * n, K_MATRIX_FILENAME);
	for (int step = 1; step < timesteps; step++) {
	 	this->assembleRightVector(step);
	 	this->setBorderConditions(step);
	 	//printMatrix2File("k_matrix.txt", m_K, m_F, count * n);
	 	dgesv(count * n, m_K, m_F);
	 	m_pr->setU(m_F);
	 	memset(m_F, 0, count * n * sizeof(real_t));
	 	binfile2data(m_K, count * n * count * n, K_MATRIX_FILENAME);
	}
	cout << "number of points = " << count << endl;
	for (int i = 0; i < count; i++) {
		//cout << i * n + 0 << "=" << m_pr->getU(i, 0) << " " << m_pr->getBorderCondition(i, 0, 0) <<endl;
		//cout << i * n + 1 << "=" << m_pr->getU(i, 1) << " " << m_pr->getBorderCondition(i, 1, 0) <<endl;
		cout << abs(abs(m_pr->getU(i, 2)) - abs(m_pr->getBorderCondition(i, 2, 0)) )<<endl;
	}
}

real_t CFemLocalLinear2D::integrateiNjN(const int iN, const int jN, const int elementIdx) {

	if (iN == jN) {
		return getSquare(elementIdx)/6;
	}
	else {
		return getSquare(elementIdx)/12;
	}
	return 0;
}

real_t CFemLocalLinear2D::integrateiNjdN(const int iN, const int jN, const int elementIdx) {

	return getSquare(elementIdx)/3;
}

real_t CFemLocalLinear2D::integrateidNjdN(const int iN, const int jN, const int elementIdx) {

	return getSquare(elementIdx);
}
