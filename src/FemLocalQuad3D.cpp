#include "Fem.h"
#include <iostream>
#include "f2c.h"
#include "dgemv.h"
#include "dgesv.h"

using namespace std;

CFemLocalQuad3D::CFemLocalQuad3D(CMesh *mesh) : CFemLocalLinear3D(mesh) {
	m_mesh = mesh;
}

CFemLocalQuad3D::~CFemLocalQuad3D() {
}

real_t CFemLocalQuad3D::getN(const int idxN, std::vector<real_t> ksi) {
	if (ksi.size() == 0) {
		return 0;
	}
	real_t res = 0;
	switch(idxN) {
	case 0:
		res = (2 * ksi[0] - 1) * ksi[0];
		break;
	case 1:
		res = (2 * ksi[1] - 1) * ksi[1];
		break;
	case 2:
		res = (2 * ksi[2] - 1) * ksi[2];
		break;
	case 3:
		res = (2 * ksi[3] - 1) * ksi[3];
		break;
	case 4:
		res = 4 * ksi[1] * ksi[2];
		break;
	case 5:
		res = 4 * ksi[0] * ksi[2];
		break;
	case 6:
		res = 4 * ksi[0] * ksi[1];
		break;
	case 7:
		res = 4 * ksi[0] * ksi[3];
		break;
	case 8:
		res = 4 * ksi[1] * ksi[3];
		break;
	case 9:
		res = 4 * ksi[2] * ksi[3];
		break;
	default:
		break;
	}
	return res;
}

real_t CFemLocalQuad3D::getdNdKsi(const int idxN, const int idxKsi, const std::vector<real_t> ksi) {
	if (ksi.size() == 0) {
		return 0;
	}
	real_t res = 0;

	switch(idxN) {
		case 0:
			switch(idxKsi) {
				case 0:
					res = 4 * ksi[0] - 1;
					break;
				case 1:
					res = 0;
					break;
				case 2:
					res = 0;
					break;
				case 3:
					res = 0;
					break;
			}
			break;
		case 1:
			switch(idxKsi) {
				case 0:
					res = 0;
					break;
				case 1:
					res = 4 * ksi[1] - 1;
					break;
				case 2:
					res = 0;
					break;
				case 3:
					res = 0;
					break;
			}
			break;
		case 2:
			switch(idxKsi) {
				case 0:
					res = 0;
					break;
				case 1:
					res = 0;
					break;
				case 2:
					res = 4 * ksi[2] - 1;
					break;
				case 3:
					res = 0;
					break;
			}
			break;
		case 3:
			switch(idxKsi) {
				case 0:
					res = 0;
					break;
				case 1:
					res = 0;
					break;
				case 2:
					res = 0;
					break;
				case 3:
					res = 4 * ksi[3] - 1;
					break;
			}
			break;
		case 4:
			switch(idxKsi) {
				case 0:
					res = 4 * ksi[1];
					break;
				case 1:
					res = 4 * ksi[0];
					break;
				case 2:
					res = 0;
					break;
				case 3:
					res = 0;
					break;
			}
			break;
		case 5:
			switch(idxKsi) {
				case 0:
					res = 0;
					break;
				case 1:
					res = 4 * ksi[2];
					break;
				case 2:
					res = 4 * ksi[1];
					break;
				case 3:
					res = 0;
					break;
			}
			break;
		case 6:
			switch(idxKsi) {
				case 0:
					res = 4 * ksi[2];
					break;
				case 1:
					res = 0;
					break;
				case 2:
					res = 4 * ksi[0];
					break;
				case 3:
					res = 0;
					break;
			}
			break;
		case 7:
			switch(idxKsi) {
				case 0:
					res = 4 * ksi[3];
					break;
				case 1:
					res = 0;
					break;
				case 2:
					res = 0;
					break;
				case 3:
					res = 4 * ksi[0];
					break;
			}
			break;
		case 8:
			switch(idxKsi) {
				case 0:
					res = 0;
					break;
				case 1:
					res = 4 * ksi[3];
					break;
				case 2:
					res = 0;
					break;
				case 3:
					res = 4 * ksi[1];
					break;
			}
			break;
		case 9:
			switch(idxKsi) {
				case 0:
					res = 0;
					break;
				case 1:
					res = 0;
					break;
				case 2:
					res = 4 * ksi[3];
					break;
				case 3:
					res = 4 * ksi[2];
					break;
			}
			break;
		}
	return res;
}

real_t CFemLocalQuad3D::getdNdX(const int idxN, const int element, std::vector<real_t> ksi) {
	std::vector<real_t> dKsidX = getdKsidXYZ(0, element, ksi);
	real_t N = 0;
	real_t sum = 0;
	for (int i = 0; i < ksi.size(); i++) {
		sum += getdNdKsi(idxN, i, ksi) * dKsidX[i];
	}

	return sum;
}

real_t CFemLocalQuad3D::getdNdY(const int idxN, const int element, std::vector<real_t> ksi) {
	std::vector<real_t> dKsidY = getdKsidXYZ(1, element, ksi);
	real_t sum = 0;
	for (int i = 0; i < ksi.size(); i++) {
		sum += getdNdKsi(idxN, i, ksi) * dKsidY[i];
	}

	return sum;
}

real_t CFemLocalQuad3D::getdNdZ(const int idxN, const int element, std::vector<real_t> ksi) {
	std::vector<real_t> dKsidZ = getdKsidXYZ(2, element, ksi);
	real_t sum = 0;
	for (int i = 0; i < ksi.size(); i++) {
		sum += getdNdKsi(idxN, i, ksi) * dKsidZ[i];

	}

	return sum;
}

std::vector<real_t> CFemLocalQuad3D::getdKsidXYZ(const int dim, const int element, std::vector<real_t> ksi) {
	vector<int> points = m_mesh->getElementByIndex(element);
	CPoint3D p1 = m_mesh->getPointByIndex(points[0]);
	CPoint3D p2 = m_mesh->getPointByIndex(points[1]);
	CPoint3D p3 = m_mesh->getPointByIndex(points[2]);
	CPoint3D p4 = m_mesh->getPointByIndex(points[3]);
	CPoint3D p5 = m_mesh->getPointByIndex(points[4]);
	CPoint3D p6 = m_mesh->getPointByIndex(points[5]);
	CPoint3D p7 = m_mesh->getPointByIndex(points[6]);
	CPoint3D p8 = m_mesh->getPointByIndex(points[7]);
	CPoint3D p9 = m_mesh->getPointByIndex(points[8]);
	CPoint3D p10 = m_mesh->getPointByIndex(points[9]);

	real_t Jx1, Jx2, Jx3, Jx4,
		   Jy1, Jy2, Jy3, Jy4,
		   Jz1, Jz2, Jz3, Jz4;

	Jx1 = p1.m_x * (4 * ksi[0] - 1) + 4 * p5.m_x * ksi[1] + 4 * p7.m_x * ksi[2] + 4 * p8.m_x * ksi[3];
	Jx2 = p2.m_x * (4 * ksi[1] - 1) + 4 * p6.m_x * ksi[2] + 4 * p5.m_x * ksi[0] + 4 * p9.m_x * ksi[3];
	Jx3 = p3.m_x * (4 * ksi[2] - 1) + 4 * p7.m_x * ksi[0] + 4 * p6.m_x * ksi[1] + 4 * p10.m_x * ksi[3];
	Jx4 = p4.m_x * (4 * ksi[3] - 1) + 4 * p8.m_x * ksi[0] + 4 * p9.m_x * ksi[1] + 4 * p10.m_x * ksi[2];

	Jy1 = p1.m_y * (4 * ksi[0] - 1) + 4 * p5.m_y * ksi[1] + 4 * p7.m_y * ksi[2] + 4 * p8.m_y * ksi[3];
	Jy2 = p2.m_y * (4 * ksi[1] - 1) + 4 * p6.m_y * ksi[2] + 4 * p5.m_y * ksi[0] + 4 * p9.m_y * ksi[3];
	Jy3 = p3.m_y * (4 * ksi[2] - 1) + 4 * p7.m_y * ksi[0] + 4 * p6.m_y * ksi[1] + 4 * p10.m_y * ksi[3];
	Jy4 = p4.m_y * (4 * ksi[3] - 1) + 4 * p8.m_y * ksi[0] + 4 * p9.m_y * ksi[1] + 4 * p10.m_y * ksi[2];

	Jz1 = p1.m_z * (4 * ksi[0] - 1) + 4 * p5.m_z * ksi[1] + 4 * p7.m_z * ksi[2] + 4 * p8.m_z * ksi[3];
	Jz2 = p2.m_z * (4 * ksi[1] - 1) + 4 * p6.m_z * ksi[2] + 4 * p5.m_z * ksi[0] + 4 * p9.m_z * ksi[3];
	Jz3 = p3.m_z * (4 * ksi[2] - 1) + 4 * p7.m_z * ksi[0] + 4 * p6.m_z * ksi[1] + 4 * p10.m_z * ksi[3];
	Jz4 = p4.m_z * (4 * ksi[3] - 1) + 4 * p8.m_z * ksi[0] + 4 * p9.m_z * ksi[1] + 4 * p10.m_z * ksi[2];

	real_t *matrix = new real_t[9];

	matrix[0] = Jx2;
	matrix[1] = Jy2;
	matrix[2] = Jz2;
	matrix[3] = Jx3;
	matrix[4] = Jy3;
	matrix[5] = Jz3;
	matrix[6] = Jx4;
	matrix[7] = Jy4;
	matrix[8] = Jz4;

	real_t a11 = matrix[0] * matrix[4] * matrix[8] + matrix[1] * matrix[5] * matrix[6] + matrix[3] * matrix[7] * matrix[2] - \
				 matrix[2] * matrix[4] * matrix[6] - matrix[5] * matrix[7] * matrix[0] - matrix[3] * matrix[1] * matrix[8];

	matrix[0] = Jx1;
	matrix[1] = Jy1;
	matrix[2] = Jz1;
	matrix[3] = Jx3;
	matrix[4] = Jy3;
	matrix[5] = Jz3;
	matrix[6] = Jx4;
	matrix[7] = Jy4;
	matrix[8] = Jz4;

	real_t a12 = matrix[0] * matrix[4] * matrix[8] + matrix[1] * matrix[5] * matrix[6] + matrix[3] * matrix[7] * matrix[2] - \
				 matrix[2] * matrix[4] * matrix[6] - matrix[5] * matrix[7] * matrix[0] - matrix[3] * matrix[1] * matrix[8];
	a12 = -a12;

	matrix[0] = Jx1;
	matrix[1] = Jy1;
	matrix[2] = Jz1;
	matrix[3] = Jx2;
	matrix[4] = Jy2;
	matrix[5] = Jz2;
	matrix[6] = Jx4;
	matrix[7] = Jy4;
	matrix[8] = Jz4;

	real_t a13 = matrix[0] * matrix[4] * matrix[8] + matrix[1] * matrix[5] * matrix[6] + matrix[3] * matrix[7] * matrix[2] - \
				 matrix[2] * matrix[4] * matrix[6] - matrix[5] * matrix[7] * matrix[0] - matrix[3] * matrix[1] * matrix[8];

	matrix[0] = Jx1;
	matrix[1] = Jy1;
	matrix[2] = Jz1;
	matrix[3] = Jx2;
	matrix[4] = Jy2;
	matrix[5] = Jz2;
	matrix[6] = Jx3;
	matrix[7] = Jy3;
	matrix[8] = Jz3;

	real_t a14 = matrix[0] * matrix[4] * matrix[8] + matrix[1] * matrix[5] * matrix[6] + matrix[3] * matrix[7] * matrix[2] - \
				 matrix[2] * matrix[4] * matrix[6] - matrix[5] * matrix[7] * matrix[0] - matrix[3] * matrix[1] * matrix[8];
	a14 = -a14;
	real_t detA = a11 + a12 + a13 + a14;

	matrix[0] = 1;
	matrix[1] = Jy2;
	matrix[2] = Jz2;
	matrix[3] = 1;
	matrix[4] = Jy3;
	matrix[5] = Jz3;
	matrix[6] = 1;
	matrix[7] = Jy4;
	matrix[8] = Jz4;

	real_t a21 = matrix[0] * matrix[4] * matrix[8] + matrix[1] * matrix[5] * matrix[6] + matrix[3] * matrix[7] * matrix[2] - \
				 matrix[2] * matrix[4] * matrix[6] - matrix[5] * matrix[7] * matrix[0] - matrix[3] * matrix[1] * matrix[8];
	a21 = -a21;

	matrix[0] = 1;
	matrix[1] = Jy1;
	matrix[2] = Jz1;
	matrix[3] = 1;
	matrix[4] = Jy3;
	matrix[5] = Jz3;
	matrix[6] = 1;
	matrix[7] = Jy4;
	matrix[8] = Jz4;

	real_t a22 = matrix[0] * matrix[4] * matrix[8] + matrix[1] * matrix[5] * matrix[6] + matrix[3] * matrix[7] * matrix[2] - \
				 matrix[2] * matrix[4] * matrix[6] - matrix[5] * matrix[7] * matrix[0] - matrix[3] * matrix[1] * matrix[8];

	matrix[0] = 1;
	matrix[1] = Jy1;
	matrix[2] = Jz1;
	matrix[3] = 1;
	matrix[4] = Jy2;
	matrix[5] = Jz2;
	matrix[6] = 1;
	matrix[7] = Jy4;
	matrix[8] = Jz4;

	real_t a23 = matrix[0] * matrix[4] * matrix[8] + matrix[1] * matrix[5] * matrix[6] + matrix[3] * matrix[7] * matrix[2] - \
				 matrix[2] * matrix[4] * matrix[6] - matrix[5] * matrix[7] * matrix[0] - matrix[3] * matrix[1] * matrix[8];
	a23 = -a23;

	matrix[0] = 1;
	matrix[1] = Jy1;
	matrix[2] = Jz1;
	matrix[3] = 1;
	matrix[4] = Jy2;
	matrix[5] = Jz2;
	matrix[6] = 1;
	matrix[7] = Jy3;
	matrix[8] = Jz3;

	real_t a24 = matrix[0] * matrix[4] * matrix[8] + matrix[1] * matrix[5] * matrix[6] + matrix[3] * matrix[7] * matrix[2] - \
				 matrix[2] * matrix[4] * matrix[6] - matrix[5] * matrix[7] * matrix[0] - matrix[3] * matrix[1] * matrix[8];

	matrix[0] = 1;
	matrix[1] = Jx2;
	matrix[2] = Jz2;
	matrix[3] = 1;
	matrix[4] = Jx3;
	matrix[5] = Jz3;
	matrix[6] = 1;
	matrix[7] = Jx4;
	matrix[8] = Jz4;

	real_t a31 = matrix[0] * matrix[4] * matrix[8] + matrix[1] * matrix[5] * matrix[6] + matrix[3] * matrix[7] * matrix[2] - \
				 matrix[2] * matrix[4] * matrix[6] - matrix[5] * matrix[7] * matrix[0] - matrix[3] * matrix[1] * matrix[8];

	matrix[0] = 1;
	matrix[1] = Jx1;
	matrix[2] = Jz1;
	matrix[3] = 1;
	matrix[4] = Jx3;
	matrix[5] = Jz3;
	matrix[6] = 1;
	matrix[7] = Jx4;
	matrix[8] = Jz4;

	real_t a32 = matrix[0] * matrix[4] * matrix[8] + matrix[1] * matrix[5] * matrix[6] + matrix[3] * matrix[7] * matrix[2] - \
				 matrix[2] * matrix[4] * matrix[6] - matrix[5] * matrix[7] * matrix[0] - matrix[3] * matrix[1] * matrix[8];
	a32 = -a32;
	matrix[0] = 1;
	matrix[1] = Jx1;
	matrix[2] = Jz1;
	matrix[3] = 1;
	matrix[4] = Jx2;
	matrix[5] = Jz2;
	matrix[6] = 1;
	matrix[7] = Jx4;
	matrix[8] = Jz4;

	real_t a33 = matrix[0] * matrix[4] * matrix[8] + matrix[1] * matrix[5] * matrix[6] + matrix[3] * matrix[7] * matrix[2] - \
				 matrix[2] * matrix[4] * matrix[6] - matrix[5] * matrix[7] * matrix[0] - matrix[3] * matrix[1] * matrix[8];

	matrix[0] = 1;
	matrix[1] = Jx1;
	matrix[2] = Jz1;
	matrix[3] = 1;
	matrix[4] = Jx2;
	matrix[5] = Jz2;
	matrix[6] = 1;
	matrix[7] = Jx3;
	matrix[8] = Jz3;

	real_t a34 = matrix[0] * matrix[4] * matrix[8] + matrix[1] * matrix[5] * matrix[6] + matrix[3] * matrix[7] * matrix[2] - \
				 matrix[2] * matrix[4] * matrix[6] - matrix[5] * matrix[7] * matrix[0] - matrix[3] * matrix[1] * matrix[8];
	a34 = -a34;

	matrix[0] = 1;
	matrix[1] = Jx2;
	matrix[2] = Jy2;
	matrix[3] = 1;
	matrix[4] = Jx3;
	matrix[5] = Jy3;
	matrix[6] = 1;
	matrix[7] = Jx4;
	matrix[8] = Jy4;

	real_t a41 = matrix[0] * matrix[4] * matrix[8] + matrix[1] * matrix[5] * matrix[6] + matrix[3] * matrix[7] * matrix[2] - \
				 matrix[2] * matrix[4] * matrix[6] - matrix[5] * matrix[7] * matrix[0] - matrix[3] * matrix[1] * matrix[8];
	a41 = -a41;

	matrix[0] = 1;
	matrix[1] = Jx1;
	matrix[2] = Jy1;
	matrix[3] = 1;
	matrix[4] = Jx3;
	matrix[5] = Jy3;
	matrix[6] = 1;
	matrix[7] = Jx4;
	matrix[8] = Jy4;

	real_t a42 = matrix[0] * matrix[4] * matrix[8] + matrix[1] * matrix[5] * matrix[6] + matrix[3] * matrix[7] * matrix[2] - \
				 matrix[2] * matrix[4] * matrix[6] - matrix[5] * matrix[7] * matrix[0] - matrix[3] * matrix[1] * matrix[8];

	matrix[0] = 1;
	matrix[1] = Jx1;
	matrix[2] = Jy1;
	matrix[3] = 1;
	matrix[4] = Jx2;
	matrix[5] = Jy2;
	matrix[6] = 1;
	matrix[7] = Jx4;
	matrix[8] = Jy4;

	real_t a43 = matrix[0] * matrix[4] * matrix[8] + matrix[1] * matrix[5] * matrix[6] + matrix[3] * matrix[7] * matrix[2] - \
				 matrix[2] * matrix[4] * matrix[6] - matrix[5] * matrix[7] * matrix[0] - matrix[3] * matrix[1] * matrix[8];
	a43 = -a43;

	matrix[0] = 1;
	matrix[1] = Jx1;
	matrix[2] = Jy1;
	matrix[3] = 1;
	matrix[4] = Jx2;
	matrix[5] = Jy2;
	matrix[6] = 1;
	matrix[7] = Jx3;
	matrix[8] = Jy3;

	real_t a44 = matrix[0] * matrix[4] * matrix[8] + matrix[1] * matrix[5] * matrix[6] + matrix[3] * matrix[7] * matrix[2] - \
				 matrix[2] * matrix[4] * matrix[6] - matrix[5] * matrix[7] * matrix[0] - matrix[3] * matrix[1] * matrix[8];

	real_t *m1 = new real_t[16];
	m1[0] = a11/detA;
	m1[1] = a12/detA;
	m1[2] = a13/detA;
	m1[3] = a14/detA;
	m1[4] = a21/detA;
	m1[5] = a22/detA;
	m1[6] = a23/detA;
	m1[7] = a24/detA;
	m1[8] = a31/detA;
	m1[9] = a32/detA;
	m1[10] = a33/detA;
	m1[11] = a34/detA;
	m1[12] = a41/detA;
	m1[13] = a42/detA;
	m1[14] = a43/detA;
	m1[15] = a44/detA;

	real_t *m2 = new real_t[12];
	m2[0] = 0;
	m2[1] = 1;
	m2[2] = 0;
	m2[3] = 0;
	m2[4] = 0;
	m2[5] = 0;
	m2[6] = 1;
	m2[7] = 0;
	m2[8] = 0;
	m2[9] = 0;
	m2[10] = 0;
	m2[11] = 1;



	real_t *c = new real_t[12];

	char *ch = "N";
	int m_m = 4, //num of rows of m1
		m_n = 3, //num of cols of m2
		m_k = 4, //num of rows of m2
		lda = 4, //max(M,1)
		ldb = 4, // max(K,1)
		ldc = 4; // max(M,1)
	real_t alpha = 1,///detA,
		   beta = 0;
	dgemm(ch, ch, m_m, m_n, m_k, alpha, &m1[0], lda, &m2[0], ldb, beta, &c[0], ldc);

	std::vector<real_t> dKsi;
	int n_start = 0, n_end = 0;
	switch(dim) {
		case 0:
			n_start = 0;
			n_end = 4;
			break;
		case 1:
			n_start = 4;
			n_end = 8;
			break;
		case 2:
			n_start = 8;
			n_end = 12;
			break;
	}
	for (int i = n_start; i < n_end; i++) {
		dKsi.push_back(c[i]);
	}

	delete [] m1;
	delete [] m2;
	delete [] c;
	delete [] matrix;
	return dKsi;
}
