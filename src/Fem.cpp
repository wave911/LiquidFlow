#include "Fem.h"
#include <iostream>
#ifdef LIBBLASLAPACK
#include "f2c.h"
#include "dgemv.h"
#include "dgesv.h"
#include "dgemm.h"
#endif
#ifdef APPLE
#include <Accelerate/Accelerate.h>
#endif

using namespace std;

CFem::CFem(const int equations_number, const MeshGeometryType mgt) {
	m_pr = nullptr;
	m_U_temp = nullptr;
	m_mesh = nullptr;
	m_K = nullptr;
	m_F = nullptr;
	m_C = nullptr;
	m_equations_number = equations_number;
	m_mesh_geometry_type = mgt;
};

void CFem::init(CProblem *pr) {
	//const int n = 4; //number of equations to solve
	const int n = m_equations_number;
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

void CFem::assembleKMatrix() {
	const int n = m_equations_number;
	real_t cc = 0,
		   kk = 0;
	CGaussRule *gr = new CGaussRule(3, m_mesh_geometry_type);
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
								cc = getCC(g_col, g_row, l_col, l_row, i, gr->m_p[0]);
								for (int l = 0; l < gr->m_intpoints; l++) {
									kk += getVolume(i) * (getKK(g_col, g_row, l_col, l_row, i, gr->m_p[l])) * gr->m_wi[l];
								}
							}
							else {
								for (int l = 0; l < gr->m_intpoints; l++) {
									kk += getVolume(i) * (getKK(g_col, g_row, l_col, l_row, i, gr->m_p[l]) * this->getN(g_row, gr->m_p[l]) * gr->m_wi[l]);
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

void CFem::assembleRightVector(const int timestep) {
	const int n = m_equations_number;
	int elnumber = m_mesh->getElementsNumber();
	int ptnumber = m_mesh->getPointsNumber();
	CGaussRule *gr = new CGaussRule(3, m_mesh_geometry_type);
	for (int i = 0; i < elnumber; i++) {
		std::vector<int> element = m_mesh->getElementByIndex(i);
		for (int j = 0; j < element.size(); j++) {
			real_t U1 = m_pr->getU(element[j], 0);
			real_t U2 = m_pr->getU(element[j], 1);
			//integration over RHS
			for (int k = 0; k < n; k++) {
				for (int l = 0; l < gr->m_intpoints; l++) {
					m_F[element[j] * n + k] -= getVolume(i) * ( getFF(element[j], k, i, gr->m_p[l]) * getN(j, gr->m_p[l]) * gr->m_wi[l]);
				}
				if (k < n - 1) {
					m_U_temp[element[j] * n + k] = m_pr->getU(element[j], k);
				}
			}
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

void CFem::setBorderConditions(const int timestep) {
	const int n = m_equations_number;
	int ptnumber = m_mesh->getPointsNumber();
	std::set<int> borderPoints = m_mesh->getBorderPoints();

	for (auto const& i : borderPoints) {
		for (int k = 0; k < n; k++) {
			m_F[i * n + k] = m_pr->getBorderCondition(i, k, (timestep) * m_pr->getTau());
		}
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

void CFem::perform(const int timesteps) {
	const int n = m_equations_number;
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
	vector<vector<real_t>> error_u(n, vector<real_t>(0));
	cout << "number of points = " << count << endl;
	for (int i = 0; i < count; i++) {
		for (int j = 0; j < n; j++) {
			error_u[j].push_back( abs(m_pr->getU(i, j) - m_pr->getBorderCondition(i, j, (timesteps - 1) * m_pr->getTau())) );
			//cout << m_pr->getU(i, j) << " " << m_pr->getBorderCondition(i, j, (timesteps - 1) * m_pr->getTau()) << endl;
		}
	}
	for (int i = 0; i < n; i++) {
		cout << "error U" <<  i << " max " << *max_element(error_u[i].begin(), error_u[i].end()) << endl;
	}
}

void CFem::dgemv(char *trans, int m, int n, real_t alpha, real_t *a, 
				   int lda, real_t *x, int incx, real_t beta, real_t *y, int incy) {
#ifdef LIBBLASLAPACK
	dgemv_(trans, (integer*)&m, (integer*)&n, (doublereal*)&alpha,
				(doublereal*)a, (integer*)&lda, (doublereal*)x,
				(integer*)&incx, (doublereal*)&beta, (doublereal*)y,
				(integer*)&incy);
#endif
#ifdef APPLE
	cblas_dgemv(CBLAS_ORDER::CblasColMajor, CBLAS_TRANSPOSE::CblasNoTrans,
			    m, n, alpha, a, lda, x, incx, beta, y, incy);
#endif
}

void CFem::dgesv(int n, real_t *M, real_t *B ) {
	int status = 0;
	int NRHS = 1;
	int* IPIV = new int[n];
	int N = n;
    int lda = n;
    int ldb = n;
#ifdef LIBBLASLAPACK
	dgesv_((integer*)&N, (integer*)&NRHS, M, (integer*)&lda, (integer*)IPIV, B, (integer*)&ldb, (integer*)&status);
#endif
#ifdef APPLE
	dgesv_((__CLPK_integer*)&N, (__CLPK_integer*)&NRHS, M, (__CLPK_integer*)&lda, (__CLPK_integer*)IPIV, B, (__CLPK_integer*)&ldb, (__CLPK_integer*)&status);
#endif
	delete [] IPIV;
}

void CFem::dgemm(char *TransA, char *TransB, int m, int n, int k, real_t alpha, real_t *A, int lda,
				real_t *B, int ldb, real_t beta, real_t *C, int ldc) {
#ifdef LIBBLASLAPACK
	dgemm_(TransA, TransB, (integer*)&m, (integer*)&n, (integer*)&k, (doublereal*)&alpha,
			  (doublereal*)&A, (integer*)&lda, (doublereal*)&B, (integer*)&ldb, (doublereal*)&beta,
			  (doublereal*)&C, (integer*)&ldc);
#endif
#ifdef APPLE
	cblas_dgemm(CBLAS_ORDER::CblasColMajor, CBLAS_TRANSPOSE::CblasNoTrans,
			CBLAS_TRANSPOSE::CblasNoTrans, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
#endif
}




