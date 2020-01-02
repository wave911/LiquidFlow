import numpy as np
import mesh_parser as mp
import gauss_integration as gi

class FemLocalLinear3D(object):

    def __init__(self, mesh, U, equations_number):
        self.mesh = mesh
        self.U = U
        self.N_count = len(mesh[mp.MESH_KEY.T3D][0]) #4
        self.equations_number = equations_number
        pass

    def __get_volume__(self, element_idx, vol_id):
        elem_points = self.mesh[mp.MESH_KEY.T3D][element_idx]

        p1 = self.get_point_by_index(elem_points[0])
        p2 = self.get_point_by_index(elem_points[1])
        p3 = self.get_point_by_index(elem_points[2])
        p4 = self.get_point_by_index(elem_points[3])

        p1_x, p2_x, p3_x, p4_x = p1[0], p2[0], p3[0], p4[0]
        p1_y, p2_y, p3_y, p4_y = p1[1], p2[1], p3[1], p4[1]
        p1_z, p2_z, p3_z, p4_z = p1[2], p2[2], p3[2], p4[2]

        vo1 = (p2_x * (p3_y * p4_z - p4_y * p3_z) + p3_x * (p4_y * p2_z - p2_y * p4_z) +
               p4_x * (p2_y * p3_z - p3_y * p2_z))/6
        vo2 = (p1_x * (p4_y * p3_z - p3_y * p4_z) + p3_x * (p1_y * p4_z - p4_y * p1_z) +
               p4_x * (p3_y * p1_z - p1_y * p3_z))/6
        vo3 = (p1_x * (p2_y * p4_z - p4_y * p2_z) + p2_x * (p4_y * p1_z - p1_y * p4_z) +
               p4_x * (p1_y * p2_z - p2_y * p1_z))/6
        vo4 = (p1_x * (p3_y * p2_z - p2_y * p3_z) + p2_x * (p1_y * p3_z - p3_y * p1_z) +
               p3_x * (p2_y * p1_z - p1_y * p2_z))/6
        volume = [vo1 + vo2 + vo3 + vo4,
                  vo1,
                  vo2,
                  vo3,
                  vo4]

        return volume[vol_id]

    def get_volume(self, element_idx):
        return self.__get_volume__(element_idx, 0)

    def get_N(self, ksi):
        N = [ksi[0], ksi[1], ksi[2], ksi[3]]
        return N

    def __get_dN_dKsi__(self, idx_N, ksi):
        dN = {}
        dN[0] = [1, 0, 0, 0]
        dN[1] = [0, 1, 0, 0]
        dN[2] = [0, 0, 1, 0]
        dN[3] = [0, 0, 0, 1]
        return dN[idx_N]

    def __get_dKsi_dx__(self, element_idx, ksi):
        volume = self.get_volume(element_idx)
        elem_points = self.mesh[mp.MESH_KEY.T3D][element_idx]

        p1 = self.get_point_by_index(elem_points[0])
        p2 = self.get_point_by_index(elem_points[1])
        p3 = self.get_point_by_index(elem_points[2])
        p4 = self.get_point_by_index(elem_points[3])

        p1_y, p2_y, p3_y, p4_y = p1[1], p2[1], p3[1], p4[1]
        p1_z, p2_z, p3_z, p4_z = p1[2], p2[2], p3[2], p4[2]

        dKsi_dx = [
            ((p4_y - p2_y) * (p3_z - p2_z) - (p3_y - p2_y) * (p4_z - p2_z))/(6 * volume),
            ((p3_y - p1_y) * (p4_z - p3_z) - (p3_y - p4_y) * (p1_z - p3_z))/(6 * volume),
            ((p2_y - p4_y) * (p1_z - p4_z) - (p1_y - p4_y) * (p2_z - p4_z))/(6 * volume),
            ((p1_y - p3_y) * (p2_z - p1_z) - (p1_y - p2_y) * (p3_z - p1_z))/(6 * volume)
        ]
        return dKsi_dx

    def __get_dKsi_dy__(self, element_idx, ksi):
        volume = self.get_volume(element_idx)
        elem_points = self.mesh[mp.MESH_KEY.T3D][element_idx]

        p1 = self.get_point_by_index(elem_points[0])
        p2 = self.get_point_by_index(elem_points[1])
        p3 = self.get_point_by_index(elem_points[2])
        p4 = self.get_point_by_index(elem_points[3])

        p1_x, p2_x, p3_x, p4_x = p1[0], p2[0], p3[0], p4[0]
        p1_y, p2_y, p3_y, p4_y = p1[1], p2[1], p3[1], p4[1]
        p1_z, p2_z, p3_z, p4_z = p1[2], p2[2], p3[2], p4[2]

        dKsi_dy = [
            ((p4_z - p2_z) * (p3_x - p2_x) - (p3_z - p2_z) * (p4_x - p2_x))/(6 * volume),
            ((p3_z - p1_z) * (p4_x - p3_x) - (p3_z - p4_z) * (p1_x - p3_x))/(6 * volume),
            ((p2_z - p4_z) * (p1_x - p4_x) - (p1_z - p4_z) * (p2_x - p4_x))/(6 * volume),
            ((p1_z - p3_z) * (p2_x - p1_x) - (p1_z - p2_z) * (p3_x - p1_x))/(6 * volume)
        ]
        return dKsi_dy

    def __get_dKsi_dz__(self, element_idx, ksi):
        volume = self.get_volume(element_idx)
        elem_points = self.mesh[mp.MESH_KEY.T3D][element_idx]

        p1 = self.get_point_by_index(elem_points[0])
        p2 = self.get_point_by_index(elem_points[1])
        p3 = self.get_point_by_index(elem_points[2])
        p4 = self.get_point_by_index(elem_points[3])

        p1_x, p2_x, p3_x, p4_x = p1[0], p2[0], p3[0], p4[0]
        p1_y, p2_y, p3_y, p4_y = p1[1], p2[1], p3[1], p4[1]
        p1_z, p2_z, p3_z, p4_z = p1[2], p2[2], p3[2], p4[2]

        dKsi_dz = [
            ((p4_x - p2_x) * (p3_y - p2_y) - (p3_x - p2_x) * (p4_y - p2_y))/(6 * volume),
            ((p3_x - p1_x) * (p4_y - p3_y) - (p3_x - p4_x) * (p1_y - p3_y))/(6 * volume),
            ((p2_x - p4_x) * (p1_y - p4_y) - (p1_x - p4_x) * (p2_y - p4_y))/(6 * volume),
            ((p1_x - p3_x) * (p2_y - p1_y) - (p1_x - p2_x) * (p3_y - p1_y))/(6 * volume)
        ]
        return dKsi_dz

    def get_dN_dx(self, element_idx, ksi):
        dn_dx = []
        dKsi_dx = self.__get_dKsi_dx__(element_idx, ksi)
        dn_dx = [self.__get_dN_dKsi__(i, ksi)[0] * dKsi_dx[0] +
                 self.__get_dN_dKsi__(i, ksi)[1] * dKsi_dx[1] +
                 self.__get_dN_dKsi__(i, ksi)[2] * dKsi_dx[2] +
                 self.__get_dN_dKsi__(i, ksi)[3] * dKsi_dx[3] for i in range(0, self.N_count)]
        # for i in range(0, self.N_count):
        #     dn_dksi = self.__get_dN_dKsi__(i, ksi)
        #     dn_dx.append(dn_dksi[0] * dKsi_dx[0] +
        #                  dn_dksi[1] * dKsi_dx[1] +
        #                  dn_dksi[2] * dKsi_dx[2] +
        #                  dn_dksi[3] * dKsi_dx[3])

        # in fact returns an array of 3 elements [1,2,3,4]
        return dn_dx

    def get_dN_dy(self, element_idx, ksi):
        dn_dy = []
        dKsi_dy = self.__get_dKsi_dy__(element_idx, ksi)
        dn_dy = [self.__get_dN_dKsi__(i, ksi)[0] * dKsi_dy[0] +
                 self.__get_dN_dKsi__(i, ksi)[1] * dKsi_dy[1] +
                 self.__get_dN_dKsi__(i, ksi)[2] * dKsi_dy[2] +
                 self.__get_dN_dKsi__(i, ksi)[3] * dKsi_dy[3] for i in range(0, self.N_count)]
        # for i in range(0, self.N_count):
        #     dn_dksi = self.__get_dN_dKsi__(i, ksi)
        #     dn_dy.append(dn_dksi[0] * dKsi_dy[0] +
        #                  dn_dksi[1] * dKsi_dy[1] +
        #                  dn_dksi[2] * dKsi_dy[2] +
        #                  dn_dksi[3] * dKsi_dy[3])

        # in fact returns an array of 3 elements [1,2,3,4]
        return dn_dy

    def get_dN_dz(self, element_idx, ksi):
        dn_dz = []
        dKsi_dz = self.__get_dKsi_dz__(element_idx, ksi)
        dn_dz = [self.__get_dN_dKsi__(i, ksi)[0] * dKsi_dz[0] +
                 self.__get_dN_dKsi__(i, ksi)[1] * dKsi_dz[1] +
                 self.__get_dN_dKsi__(i, ksi)[2] * dKsi_dz[2] +
                 self.__get_dN_dKsi__(i, ksi)[3] * dKsi_dz[3] for i in range(0, self.N_count)]
        # for i in range(0, self.N_count):
        #     dn_dksi = self.__get_dN_dKsi__(i, ksi)
        #     dn_dz.append(dn_dksi[0] * dKsi_dz[0] +
        #                  dn_dksi[1] * dKsi_dz[1] +
        #                  dn_dksi[2] * dKsi_dz[2] +
        #                  dn_dksi[3] * dKsi_dz[3])

        # in fact returns an array of 3 elements [1,2,3,4]
        return dn_dz

    def get_du_dx(self, element_idx, ksi):
        elem_points = self.mesh[mp.MESH_KEY.T3D][element_idx]

        dn_dx = self.get_dN_dx(element_idx, ksi)
        du_dx = [0.0, 0.0, 0.0]
        for i in range(len(elem_points)):
            du_dx[0] += self.U[elem_points[i]][0] * dn_dx[i]
            du_dx[1] += self.U[elem_points[i]][1] * dn_dx[i]
            du_dx[2] += self.U[elem_points[i]][2] * dn_dx[i]
        return du_dx

    def get_du_dy(self, element_idx, ksi):
        elem_points = self.mesh[mp.MESH_KEY.T3D][element_idx]

        dn_dy = self.get_dN_dy(element_idx, ksi)
        du_dy = [0.0, 0.0, 0.0]
        for i in range(len(elem_points)):
            du_dy[0] += self.U[elem_points[i]][0] * dn_dy[i]
            du_dy[1] += self.U[elem_points[i]][1] * dn_dy[i]
            du_dy[2] += self.U[elem_points[i]][2] * dn_dy[i]
        return du_dy

    def get_du_dz(self, element_idx, ksi):
        elem_points = self.mesh[mp.MESH_KEY.T3D][element_idx]

        dn_dz = self.get_dN_dz(element_idx, ksi)
        du_dz = [0.0, 0.0, 0.0]
        for i in range(len(elem_points)):
            du_dz[0] += self.U[elem_points[i]][0] * dn_dz[i]
            du_dz[1] += self.U[elem_points[i]][1] * dn_dz[i]
            du_dz[2] += self.U[elem_points[i]][2] * dn_dz[i]
        return du_dz

    def get_point_by_index(self, point_idx):
        return self.mesh[mp.MESH_KEY.Points][point_idx]

    def get_volume_elements(self):
        return self.mesh[mp.MESH_KEY.T3D]

    def get_border_elements(self):
        return self.mesh[mp.MESH_KEY.T2D]

    def get_gauss_rule(self, order):
        return gi.GaussRule(mp.MESH_TYPE.T3D, order)

    def get_Uc(self):
        Uc = [[p[0], p[1], p[2], 0] for p in self.U]
        return [j for i in Uc for j in i]

    def get_kk(self, idxN, jdxN, l_col, l_row, element_idx, ksi):
        n = self.equations_number
        res = 0
        if l_col == l_row:
            res = self.get_dN_dx(element_idx, ksi)[idxN] * self.get_dN_dx(element_idx, ksi)[jdxN] + \
                  self.get_dN_dy(element_idx, ksi)[idxN] * self.get_dN_dy(element_idx, ksi)[jdxN] + \
                  self.get_dN_dz(element_idx, ksi)[idxN] * self.get_dN_dz(element_idx, ksi)[jdxN]

            if l_row < n - 1:
                res = res / 1  # m_pr->getRe()
        else:
            if (3 == l_col) and (0 == l_row):
                res = self.get_dN_dx(element_idx, ksi)[idxN]
            if (3 == l_col) and (1 == l_row):
                res = self.get_dN_dy(element_idx, ksi)[idxN]
            if (3 == l_col) and (2 == l_row):
                res = self.get_dN_dz(element_idx, ksi)[idxN]
        return res

    def get_cc(self, idxN, jdxN, l_col, l_row, element_idx, ksi):
        n = self.equations_number
        res = 0
        if l_row < n - 1:
            if idxN == jdxN:
                res = self.get_volume(element_idx)/10
            else:
                res = self.get_volume(element_idx)/20
        return res

    def get_ff(self, idx, l_row, element_idx, ksi):
        U1 = self.U[idx][0]
        U2 = self.U[idx][1]
        U3 = self.U[idx][2]
        res = [
                U1 * self.get_du_dx(element_idx, ksi)[0] + U2 * self.get_du_dy(element_idx, ksi)[0] +
                U3 * self.get_du_dz(element_idx, ksi)[0],
                U1 * self.get_du_dx(element_idx, ksi)[1] + U2 * self.get_du_dy(element_idx, ksi)[1] +
                U3 * self.get_du_dz(element_idx, ksi)[1],
                U1 * self.get_du_dx(element_idx, ksi)[2] + U2 * self.get_du_dy(element_idx, ksi)[2] +
                U3 * self.get_du_dz(element_idx, ksi)[2],
                -2 * (self.get_du_dy(element_idx, ksi)[0] * self.get_du_dx(element_idx, ksi)[1] +
                      self.get_du_dz(element_idx, ksi)[0] * self.get_du_dx(element_idx, ksi)[2] +
                      self.get_du_dz(element_idx, ksi)[1] * self.get_du_dy(element_idx, ksi)[2])
            ]

        return res[l_row]

    def get_border_condition(self, point_idx, dim, time = 0.0):
        p = self.get_point_by_index(point_idx)
        bc = [-p[1],
              p[0],
              1 * time,
              -0.5 * (1 - (p[0] * p[0] + p[1] * p[1])) - 1 * p[2]
              ]
        return bc[dim]