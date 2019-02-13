#include "Fem.h"
#include <iostream>
#include "f2c.h"
#include "dgemv.h"
#include "dgesv.h"

using namespace std;

CFemLocalQuad3D::CFemLocalQuad3D(CMesh *mesh) {
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
		res = 4 * ksi[0] * ksi[1];
		break;
	case 2:
		res = 4 * ksi[0] * ksi[2];
		break;
	case 3:
		res = 4 * ksi[0] * ksi[3];
		break;
	case 4:
		res = (2 * ksi[1] - 1) * ksi[1];
		break;
	case 5:
		res = 4 * ksi[1] * ksi[2];
		break;
	case 6:
		res = 4 * ksi[1] * ksi[3];
		break;
	case 7:
		res = (2 * ksi[2] - 1) * ksi[2];
		break;
	case 8:
		res = 4 * ksi[2] * ksi[3];
		break;
	case 9:
		res = (2 * ksi[3] - 1) * ksi[3];
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
		case 2:
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
		case 3:
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
		case 4:
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
		case 7:
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
		case 8:
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
		case 9:
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
		}
	return res;
}
