import numpy as np
import mesh_parser as mp
import gauss_integration as gi

class FemLocalLinear2D(object):

    def __init__(self, mesh, U, equations_number):
        self.mesh = mesh
        self.U = U
        self.N_count = len(mesh[mp.MESH_KEY.T2D][0]) #3
        self.equations_number = equations_number
        pass

    def get_volume(self, element_idx):
        elem_points = self.mesh[mp.MESH_KEY.T2D][element_idx]

        p1 = self.get_point_by_index(elem_points[0])
        p2 = self.get_point_by_index(elem_points[1])
        p3 = self.get_point_by_index(elem_points[2])

        p1_x, p2_x, p3_x = p1[0], p2[0], p3[0]
        p1_y, p2_y, p3_y = p1[1], p2[1], p3[1]

        square = ((p2_x * p3_y - p3_x * p2_y) + (p3_x * p1_y - p1_x * p3_y) + (p1_x * p2_y - p2_x * p1_y)) / 2
        return square

    # def get_local_coordinates(self, element_idx, point):
    #     elem_points = self.mesh[mp.MESH_KEY.T2D][element_idx]
    #
    #     p1 = self.get_point_by_index(elem_points[0])
    #     p2 = self.get_point_by_index(elem_points[1])
    #     p3 = self.get_point_by_index(elem_points[2])
    #     p = self.get_point_by_index(point)
    #
    #     p1_x, p2_x, p3_x = p1[0], p2[0], p3[0]
    #     p1_y, p2_y, p3_y = p1[1], p2[1], p3[1]
    #
    #     matrix = [p2_x * p3_y - p3_x * p2_y,
    #               p3_x * p1_y - p1_x * p3_y,
    #               p1_x * p2_y - p2_x * p1_y,
    #               p2_y - p3_y,
    #               p3_y - p1_y,
    #               p1_y - p2_y,
    #               p3_x - p2_x,
    #               p1_x - p3_x,
    #               p2_x - p1_x ]
    #
    #     vector = [1, p[0], p[1]]
    #     vector = np.array(vector)
    #
    #     matrix = np.array(matrix)
    #     matrix = matrix.reshape(3, 3).swapaxes(0, 1)
    #     matrix = matrix / (2 * self.get_volume(element_idx))
    #
    #     ksi = np.matmul(matrix, vector)
    #
    #     return ksi

    def get_kk(self, idxN, jdxN, l_col, l_row, element_idx, ksi):
        n = self.equations_number
        res = 0
        if l_col == l_row:
            res = self.get_dN_dx(element_idx, ksi)[idxN] * self.get_dN_dx(element_idx, ksi)[jdxN] + \
                self.get_dN_dy(element_idx, ksi)[idxN] * self.get_dN_dy(element_idx, ksi)[jdxN]
            if l_row < n - 1:
                res = res/1 #m_pr->getRe()
        else:
            if (2 == l_col) and (0 == l_row):
                res = self.get_dN_dx(element_idx, ksi)[idxN]
            if (2 == l_col) and (1 == l_row):
                res = self.get_dN_dy(element_idx, ksi)[idxN]
        return res

    def get_cc(self, idxN, jdxN, l_col, l_row, element_idx, ksi):
        n = self.equations_number
        res = 0
        if l_row < n - 1:
            if idxN == jdxN:
                res = self.get_volume(element_idx)/6
            else:
                res = self.get_volume(element_idx)/12
        return res

    def get_ff(self, idx, l_row, element_idx, ksi):
        res = 0
        U1 = self.U[idx][0]
        U2 = self.U[idx][1]
        if l_row == 0:
            res = U1 * self.get_du_dx(element_idx, ksi)[0] + U2 * self.get_du_dy(element_idx, ksi)[0]
        if l_row == 1:
            res = U1 * self.get_du_dx(element_idx, ksi)[1] + U2 * self.get_du_dy(element_idx, ksi)[1]
        if l_row == 2:
            res = -2 * self.get_du_dy(element_idx, ksi)[0] * self.get_du_dx(element_idx, ksi)[1]

        return res

    def get_du_dx(self, element_idx, ksi):
        elem_points = self.mesh[mp.MESH_KEY.T2D][element_idx]

        dn_dx = self.get_dN_dx(element_idx, ksi)
        du_dx = [0.0, 0.0]
        for i in range(len(elem_points)):
            du_dx[0] += self.U[elem_points[i]][0] * dn_dx[i]
            du_dx[1] += self.U[elem_points[i]][1] * dn_dx[i]
        return du_dx

    def get_du_dy(self, element_idx, ksi):
        elem_points = self.mesh[mp.MESH_KEY.T2D][element_idx]

        dn_dy = self.get_dN_dy(element_idx, ksi)
        du_dy = [0.0, 0.0]
        for i in range(len(elem_points)):
            du_dy[0] += self.U[elem_points[i]][0] * dn_dy[i]
            du_dy[1] += self.U[elem_points[i]][1] * dn_dy[i]
        return du_dy

    def get_N(self, ksi):
        N = [ksi[0], ksi[1], ksi[2]]
        return N

    def __get_dN_dKsi__(self, idx_N, ksi):
        dN = {}
        dN[0] = [1, 0, 0]
        dN[1] = [0, 1, 0]
        dN[2] = [0, 0, 1]
        return dN[idx_N]

    def get_dN_dx(self, element_idx, ksi):

        square = self.get_volume(element_idx)
        dn_dx = []

        elem_points = self.mesh[mp.MESH_KEY.T2D][element_idx]
        p1 = self.get_point_by_index(elem_points[0])
        p2 = self.get_point_by_index(elem_points[1])
        p3 = self.get_point_by_index(elem_points[2])
        p1_y, p2_y, p3_y = p1[1], p2[1], p3[1]

        y23 = p2_y - p3_y
        y31 = p3_y - p1_y
        y12 = p1_y - p2_y

        for i in range(0, self.N_count):
            dn_dksi = self.__get_dN_dKsi__(i, ksi)
            dn_dx.append((dn_dksi[0] * y23 + dn_dksi[1] * y31 + dn_dksi[2] * y12)/(2 * square))

        # in fact returns an array of 3 elements [1,2,3]
        return dn_dx

    def get_dN_dy(self, element_idx, ksi):

        square = self.get_volume(element_idx)
        dn_dy = []

        elem_points = self.mesh[mp.MESH_KEY.T2D][element_idx]
        p1 = self.get_point_by_index(elem_points[0])
        p2 = self.get_point_by_index(elem_points[1])
        p3 = self.get_point_by_index(elem_points[2])
        p1_x, p2_x, p3_x = p1[0], p2[0], p3[0]

        x32 = p3_x - p2_x
        x13 = p1_x - p3_x
        x21 = p2_x - p1_x

        for i in range(0, self.N_count):
            dn_dksi = self.__get_dN_dKsi__(i, ksi)
            dn_dy.append((dn_dksi[0] * x32 + dn_dksi[1] * x13 + dn_dksi[2] * x21)/(2 * square))

        # in fact returns an array of 3 elements [1,2,3
        return dn_dy

    def get_border_condition(self, point_idx, dim, time = 0.0):
        p = self.get_point_by_index(point_idx)
        bc = [-p[1],
              p[0],
              -0.5 * (1 - (p[0] * p[0] + p[1] * p[1]))]
        return bc[dim]

    def get_point_by_index(self, point_idx):
        return self.mesh[mp.MESH_KEY.Points][point_idx]

    def get_volume_elements(self):
        return self.mesh[mp.MESH_KEY.T2D]

    def get_border_elements(self):
        return self.mesh[mp.MESH_KEY.T1D]

    def get_gauss_rule(self, order):
        return gi.GaussRule(mp.MESH_TYPE.T2D, order)

    def get_Uc(self):
        Uc = [[p[0], p[1], 0] for p in self.U]
        return [j for i in Uc for j in i]