#include "Fem.h"
#include <iostream>
#include "f2c.h"
#include "dgemv.h"
#include "dgesv.h"
#include "Common.h"
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


real_t CFemLocalLinear2D::getN(const int idxN, const std::vector<real_t> ksi) {
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
	default:
		break;
	}
	return res;
}

real_t CFemLocalLinear2D::getdNdKsi(const int idxN, const int idxKsi, const std::vector<real_t> ksi) {
	if (ksi.size() == 0) {
		return 0;
	}
	if (idxN == idxKsi)
		return 1;
	else
		return 0;
}

real_t CFemLocalLinear2D::getdNdX(const int idxN, const int element, const std::vector<real_t> ksi) {
	real_t sum = 0;
	for (int i = 0; i < ksi.size(); i++) {

		sum += getdNdKsi(idxN, i, ksi) * getdKsidX(i, element);
	}
	return sum;
}

real_t CFemLocalLinear2D::getdNdY(const int idxN, const int element, const std::vector<real_t> ksi) {
	real_t sum = 0;
	for (int i = 0; i < 3; i++) {
		sum += getdNdKsi(idxN, i, ksi) * getdKsidY(i, element);
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

real_t CFemLocalLinear2D::getdUdX(const int element_idx, const int dim, const std::vector<real_t> ksi) {
	std::vector<int> element = m_mesh->getElementByIndex(element_idx);
	real_t res = 0;
	for (int i = 0; i < element.size(); i++) {
		real_t U = m_pr->getU(element[i], dim);
		res += U * getdNdX(i, element_idx, ksi);
	}

	return res;
}

real_t CFemLocalLinear2D::getdUdY(const int element_idx, const int dim, const std::vector<real_t> ksi) {
	std::vector<int> element = m_mesh->getElementByIndex(element_idx);
	real_t res = 0;
	for (int i = 0; i < element.size(); i++) {
		real_t U = m_pr->getU(element[i], dim);
		res += U * getdNdY(i, element_idx, ksi);
	}

	return res;
}

real_t CFemLocalLinear2D::getKK(const int idxN, const int jdxN, const int l_col, const int l_row, const int element, const std::vector<real_t> ksi) {
	real_t res = 0;
	if (l_col == l_row) {
		res = getdNdX(idxN,element, ksi) * getdNdX(jdxN,element, ksi) + getdNdY(idxN,element, ksi) * getdNdY(jdxN,element, ksi);
	}
	else {
		if ((2 == l_col) && (0 == l_row))
			res = getdNdX(idxN, element, ksi);
		if ((2 == l_col) && (1 == l_row))
			res = getdNdY(idxN, element, ksi);
	}
	return res;
}

void CFemLocalLinear2D::assembleKMatrix() {
	const int n = 3;
	real_t cc = 0,
		   kk = 0,
		   alfa = 0;
	CGaussRule *gr = new CGaussRule(3, MeshGeometryType::G2D);
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
								for (int l = 0; l < gr->m_intpoints; l++) {
									kk += getSquare(i) * (getKK(g_col, g_row, l_col, l_row, i, gr->m_p[l])) * gr->m_wi[l];
								}
								if (l_row < n - 1)
									kk = kk/m_pr->getRe();
							}
							else {
								for (int l = 0; l < gr->m_intpoints; l++) {
									if ((2 == l_col) && (0 == l_row))
										kk += getSquare(i) * (getKK(g_col, g_row, l_col, l_row, i, gr->m_p[l]) * this->getN(g_row, gr->m_p[l]) * gr->m_wi[l]);
									if ((2 == l_col) && (1 == l_row))
										kk += getSquare(i) * (getKK(g_col, g_row, l_col, l_row, i, gr->m_p[l]) * this->getN(g_row, gr->m_p[l]) * gr->m_wi[l]);
								}
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
	delete gr;
}

real_t CFemLocalLinear2D::getFF(const int idxN, const int l_row, const int element, const std::vector<real_t> ksi) {
	real_t res = 0;
	real_t U1 = m_pr->getU(idxN, 0);
	real_t U2 = m_pr->getU(idxN, 1);
	switch(l_row) {
		case 0:
			res = U1 * getdUdX(element, 0, ksi) + U2 * getdUdY(element, 0, ksi);
			break;
		case 1:
			res = U1 * getdUdX(element, 1, ksi) + U2 * getdUdY(element, 1, ksi);
			break;
		case 2:
			res = -2 * getdUdY(element, 0, ksi) * getdUdX(element, 1, ksi);
			break;
	}
	return res;
}

void CFemLocalLinear2D::assembleRightVector(const int timestep) {
	const int n = 3;
	int elnumber = m_mesh->getElementsNumber();
	int ptnumber = m_mesh->getPointsNumber();
	CGaussRule *gr = new CGaussRule(3, MeshGeometryType::G2D);
	real_t ff = 0;
	for (int i = 0; i < elnumber; i++) {
		std::vector<int> element = m_mesh->getElementByIndex(i);
		for (int j = 0; j < element.size(); j++) {
			real_t U1 = m_pr->getU(element[j], 0);
			real_t U2 = m_pr->getU(element[j], 1);
			//integration over RHS
			for (int l = 0; l < gr->m_intpoints; l++) {
				m_F[element[j] * n + 0] += getSquare(i) * ( getFF(element[j], 0, i, gr->m_p[l]) * getN(j, gr->m_p[l]) * gr->m_wi[l]);
				m_F[element[j] * n + 1] += getSquare(i) * ( getFF(element[j], 1, i, gr->m_p[l]) * getN(j, gr->m_p[l]) * gr->m_wi[l]);
				m_F[element[j] * n + 2] += getSquare(i) * ( getFF(element[j], 2, i, gr->m_p[l]) * getN(j, gr->m_p[l]) * gr->m_wi[l] );
			}

			m_U_temp[element[j] * n + 0] = m_pr->getU(element[j], 0);
			m_U_temp[element[j] * n + 1] = m_pr->getU(element[j], 1);
			m_U_temp[element[j] * n + 2] = 0;

			ff = 0;
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

	delete gr;
}

void CFemLocalLinear2D::setBorderConditions(const int timestep) {
	const int n = 3;
	int ptnumber = m_mesh->getPointsNumber();
	std::set<int> borderPoints = m_mesh->getBorderPoints();

	for (auto const& i : borderPoints) {
		m_F[i * n + 0] = m_pr->getBorderCondition(i, 0, (timestep) * m_pr->getTau());
		m_F[i * n + 1] = m_pr->getBorderCondition(i, 1, (timestep) * m_pr->getTau());
		m_F[i * n + 2] = m_pr->getBorderCondition(i, 2, (timestep) * m_pr->getTau());

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
	//printMatrix2File("k_matrix_init.txt", m_K, m_F, count * n);
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
	//for (int i = 0; i < m_mesh->getElementsNumber(); i++) {
		//cout << i * n + 0 << "=" << m_pr->getU(i, 0) << " " << m_pr->getBorderCondition(i, 0, 0) << " " << m_F[i * n + 0] << endl;
		cout << i * n + 2 << "=" << m_pr->getU(i, 2) << " " << m_pr->getBorderCondition(i, 2, 0 ) << " " << m_F[i * n + 2] <<endl;
		//cout << abs(abs(m_pr->getU(i, 2)) - abs(m_pr->getBorderCondition(i, 2, (timesteps - 1) * m_pr->getTau())) )<<endl;
		//out << "points " << m_mesh->getElementByIndex(i)[0] << " " << m_mesh->getElementByIndex(i)[1] << " " << m_mesh->getElementByIndex(i)[2] << endl;
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
