import numpy as np
import mesh_parser as mp
import fem_local_linear_3d as fll3d


class FemLocalQuad3D(fll3d.FemLocalLinear3D):
    def __init__(self, mesh, U, equations_number):
        super().__init__(mesh, U, equations_number)
        pass

    def get_N(self, ksi):
        N = [ksi[0] * (2 * ksi[0] - 1),
             ksi[1] * (2 * ksi[1] - 1),
             ksi[2] * (2 * ksi[2] - 1),
             ksi[3] * (2 * ksi[3] - 1),
             4 * ksi[0] * ksi[1],
             4 * ksi[1] * ksi[2],
             4 * ksi[2] * ksi[0],
             4 * ksi[0] * ksi[3],
             4 * ksi[1] * ksi[3],
             4 * ksi[2] * ksi[3]]
        return N

    def __get_dN_dKsi__(self, idx_N, ksi):
        dN = {}
        dN[0] = [4 * ksi[0] - 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        dN[1] = [0, 4 * ksi[1] - 1, 0, 0, 0, 0, 0, 0, 0, 0]
        dN[2] = [0, 0, 4 * ksi[2] - 1, 0, 0, 0, 0, 0, 0, 0]
        dN[3] = [0, 0, 0, 4 * ksi[3] - 1, 0, 0, 0, 0, 0, 0]
        dN[4] = [4 * ksi[1], 4 * ksi[0], 0, 0, 0, 0, 0, 0, 0, 0]
        dN[5] = [0, 4 * ksi[2], 4 * ksi[1], 0, 0, 0, 0, 0, 0, 0]
        dN[6] = [4 * ksi[2], 0, 4 * ksi[0], 0, 0, 0, 0, 0, 0, 0]
        dN[7] = [4 * ksi[3], 0, 0, 4 * ksi[0], 0, 0, 0, 0, 0, 0]
        dN[8] = [0, 4 * ksi[3], 0, 4 * ksi[1], 0, 0, 0, 0, 0, 0]
        dN[9] = [0, 0, 4 * ksi[3], 4 * ksi[1], 0, 0, 0, 0, 0, 0]
        return dN[idx_N]

    def __get_dKsi_dXYZ__(self, dim, element_idx, ksi):
        elem_points = self.mesh[mp.MESH_KEY.T3D][element_idx]

        p1 = self.get_point_by_index(elem_points[0])
        p2 = self.get_point_by_index(elem_points[1])
        p3 = self.get_point_by_index(elem_points[2])
        p4 = self.get_point_by_index(elem_points[3])
        p5 = self.get_point_by_index(elem_points[4])
        p6 = self.get_point_by_index(elem_points[5])
        p7 = self.get_point_by_index(elem_points[6])
        p8 = self.get_point_by_index(elem_points[7])
        p9 = self.get_point_by_index(elem_points[8])
        p10 = self.get_point_by_index(elem_points[9])

        p1_x, p2_x, p3_x, p4_x, p5_x, p6_x, p7_x, p8_x, p9_x, p10_x = p1[0], p2[0], p3[0], p4[0], p5[0], \
                                                                      p6[0], p7[0], p8[0], p9[0], p10[0]
        p1_y, p2_y, p3_y, p4_y, p5_y, p6_y, p7_y, p8_y, p9_y, p10_y = p1[1], p2[1], p3[1], p4[1], p5[1], \
                                                                      p6[1], p7[1], p8[1], p9[1], p10[1]
        p1_z, p2_z, p3_z, p4_z, p5_z, p6_z, p7_z, p8_z, p9_z, p10_z = p1[2], p2[2], p3[2], p4[2], p5[2], \
                                                                      p6[2], p7[2], p8[2], p9[2], p10[2]

        Jx1 = p1_x * (4 * ksi[0] - 1) + 4 * p5_x * ksi[1] + 4 * p7_x * ksi[2] + 4 * p8_x * ksi[3]
        Jx2 = p2_x * (4 * ksi[1] - 1) + 4 * p6_x * ksi[2] + 4 * p5_x * ksi[0] + 4 * p9_x * ksi[3]
        Jx3 = p3_x * (4 * ksi[2] - 1) + 4 * p7_x * ksi[0] + 4 * p6_x * ksi[1] + 4 * p10_x * ksi[3]
        Jx4 = p4_x * (4 * ksi[3] - 1) + 4 * p8_x * ksi[0] + 4 * p9_x * ksi[1] + 4 * p10_x * ksi[2]

        Jy1 = p1_y * (4 * ksi[0] - 1) + 4 * p5_y * ksi[1] + 4 * p7_y * ksi[2] + 4 * p8_y * ksi[3]
        Jy2 = p2_y * (4 * ksi[1] - 1) + 4 * p6_y * ksi[2] + 4 * p5_y * ksi[0] + 4 * p9_y * ksi[3]
        Jy3 = p3_y * (4 * ksi[2] - 1) + 4 * p7_y * ksi[0] + 4 * p6_y * ksi[1] + 4 * p10_y * ksi[3]
        Jy4 = p4_y * (4 * ksi[3] - 1) + 4 * p8_y * ksi[0] + 4 * p9_y * ksi[1] + 4 * p10_y * ksi[2]

        Jz1 = p1_z * (4 * ksi[0] - 1) + 4 * p5_z * ksi[1] + 4 * p7_z * ksi[2] + 4 * p8_z * ksi[3]
        Jz2 = p2_z * (4 * ksi[1] - 1) + 4 * p6_z * ksi[2] + 4 * p5_z * ksi[0] + 4 * p9_z * ksi[3]
        Jz3 = p3_z * (4 * ksi[2] - 1) + 4 * p7_z * ksi[0] + 4 * p6_z * ksi[1] + 4 * p10_z * ksi[3]
        Jz4 = p4_z * (4 * ksi[3] - 1) + 4 * p8_z * ksi[0] + 4 * p9_z * ksi[1] + 4 * p10_z * ksi[2]

        J = [[1, 1, 1, 1],
             [Jx1, Jx2, Jx3, Jx4],
             [Jy1, Jy2, Jy3, Jy4],
             [Jz1, Jz2, Jz3, Jz4]]

        J = np.array(J)
        Zeros = [[0, 0, 0],
                 [1, 0, 0],
                 [0, 1, 0],
                 [0, 0, 1]]
        Zeros = np.array(Zeros)
        Jinv = np.linalg.inv(J)

        res = np.matmul(Jinv, Zeros)

        c = [0.0, 0.0, 0.0, 0.0]
        c[0] = res[0][dim]  # c[0]
        c[1] = res[1][dim]  # c[1]
        c[2] = res[2][dim]  # c[2]
        c[3] = res[3][dim]  # c[3]

        return c

    def __get_dKsi_dx__(self, element_idx, ksi):
        dKsi_dx = self.__get_dKsi_dXYZ__(0, element_idx, ksi)
        return dKsi_dx

    def __get_dKsi_dy__(self, element_idx, ksi):
        dKsi_dy = self.__get_dKsi_dXYZ__(1, element_idx, ksi)
        return dKsi_dy

    def __get_dKsi_dz__(self, element_idx, ksi):
        dKsi_dz = self.__get_dKsi_dXYZ__(2, element_idx, ksi)
        return dKsi_dz

    def get_dN_dx(self, element_idx, ksi):
        vol = 1  # self.get_volume(element_idx)
        c = self.__get_dKsi_dx__(element_idx, ksi)
        a1 = c[0]
        a2 = c[1]
        a3 = c[2]
        a4 = c[3]
        dn_dx = [
            (4 * ksi[0] - 1) * a1 / vol,
            (4 * ksi[1] - 1) * a2 / vol,
            (4 * ksi[2] - 1) * a3 / vol,
            (4 * ksi[3] - 1) * a4 / vol,
            4 * (ksi[1] * a1 + ksi[0] * a2) / vol,
            4 * (ksi[2] * a2 + ksi[1] * a3) / vol,
            4 * (ksi[2] * a1 + ksi[0] * a3) / vol,
            4 * (ksi[3] * a1 + ksi[0] * a4) / vol,
            4 * (ksi[3] * a2 + ksi[1] * a4) / vol,
            4 * (ksi[3] * a3 + ksi[2] * a4) / vol
        ]

        # in fact returns an array of 3 elements [1,2,3,4]
        return dn_dx

    def get_dN_dy(self, element_idx, ksi):
        vol = 1  # self.get_volume(element_idx)
        c = self.__get_dKsi_dy__(element_idx, ksi)
        b1 = c[0]
        b2 = c[1]
        b3 = c[2]
        b4 = c[3]
        dn_dy = [
            (4 * ksi[0] - 1) * b1 / vol,
            (4 * ksi[1] - 1) * b2 / vol,
            (4 * ksi[2] - 1) * b3 / vol,
            (4 * ksi[3] - 1) * b4 / vol,
            4 * (ksi[1] * b1 + ksi[0] * b2) / vol,
            4 * (ksi[2] * b2 + ksi[1] * b3) / vol,
            4 * (ksi[2] * b1 + ksi[0] * b3) / vol,
            4 * (ksi[3] * b1 + ksi[0] * b4) / vol,
            4 * (ksi[3] * b2 + ksi[1] * b4) / vol,
            4 * (ksi[3] * b3 + ksi[2] * b4) / vol
        ]

        # in fact returns an array of 3 elements [1,2,3,4]
        return dn_dy

    def get_dN_dz(self, element_idx, ksi):
        vol = 1  # self.get_volume(element_idx)
        c = self.__get_dKsi_dz__(element_idx, ksi)
        c1 = c[0]
        c2 = c[1]
        c3 = c[2]
        c4 = c[3]
        dn_dz = [
            (4 * ksi[0] - 1) * c1 / vol,
            (4 * ksi[1] - 1) * c2 / vol,
            (4 * ksi[2] - 1) * c3 / vol,
            (4 * ksi[3] - 1) * c4 / vol,
            4 * (ksi[1] * c1 + ksi[0] * c2) / vol,
            4 * (ksi[2] * c2 + ksi[1] * c3) / vol,
            4 * (ksi[2] * c1 + ksi[0] * c3) / vol,
            4 * (ksi[3] * c1 + ksi[0] * c4) / vol,
            4 * (ksi[3] * c2 + ksi[1] * c4) / vol,
            4 * (ksi[3] * c3 + ksi[2] * c4) / vol
        ]
        # in fact returns an array of 3 elements [1,2,3,4]
        return dn_dz