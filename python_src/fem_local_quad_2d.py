import numpy as np
import mesh_parser as mp
import fem_local_linear_2d as f2d

class FemLocalQuad2D(f2d.FemLocalLinear2D):

    def __init__(self, mesh, U, equations_number):
        super().__init__(mesh, U, equations_number)
        pass

    def get_N(self, ksi):
        N = [ksi[0] * (2 * ksi[0] - 1),
             ksi[1] * (2 * ksi[1] - 1),
             ksi[2] * (2 * ksi[2] - 1),
             4 * ksi[0] * ksi[1],
             4 * ksi[1] * ksi[2],
             4 * ksi[2] * ksi[0]]
        return N

    def get_dN_dx(self, element_idx, ksi):
        elem_points = self.mesh[mp.MESH_KEY.T2D][element_idx]
        p1 = self.get_point_by_index(elem_points[0])
        p2 = self.get_point_by_index(elem_points[1])
        p3 = self.get_point_by_index(elem_points[2])
        p1_x, p2_x, p3_x = p1[0], p2[0], p3[0]
        p1_y, p2_y, p3_y = p1[1], p2[1], p3[1]

        J2 = (p2_x - p1_x) * (p3_y - p1_y) - (p1_y - p2_y) * (p1_x - p3_x)
        Jy23 = p2_y - p3_y
        Jy31 = p3_y - p1_y
        Jy12 = p1_y - p2_y

        dn_dx = [(4 * ksi[0] - 1) * Jy23/J2,
                 (4 * ksi[1] - 1) * Jy31 / J2,
                 (4 * ksi[2] - 1) * Jy12 / J2,
                 4 * (ksi[1] * Jy23 + ksi[0] * Jy31) / J2,
                 4 * (ksi[2] * Jy31 + ksi[1] * Jy12) / J2,
                 4 * (ksi[0] * Jy12 + ksi[2] * Jy23) / J2 ]
        # in fact returns an array of 3 elements [1,2,3...]
        return dn_dx

    def get_dN_dy(self, element_idx, ksi):
        elem_points = self.mesh[mp.MESH_KEY.T2D][element_idx]
        p1 = self.get_point_by_index(elem_points[0])
        p2 = self.get_point_by_index(elem_points[1])
        p3 = self.get_point_by_index(elem_points[2])
        p1_x, p2_x, p3_x = p1[0], p2[0], p3[0]
        p1_y, p2_y, p3_y = p1[1], p2[1], p3[1]

        J2 = (p2_x - p1_x) * (p3_y - p1_y) - (p1_y - p2_y) * (p1_x - p3_x)
        Jx32 = p3_x - p2_x
        Jx13 = p1_x - p3_x
        Jx21 = p2_x - p1_x

        dn_dy = [(4 * ksi[0] - 1) * Jx32/J2,
                 (4 * ksi[1] - 1) * Jx13 / J2,
                 (4 * ksi[2] - 1) * Jx21 / J2,
                 4 * (ksi[1] * Jx32 + ksi[0] * Jx13) / J2,
                 4 * (ksi[2] * Jx13 + ksi[1] * Jx21) / J2,
                 4 * (ksi[0] * Jx21 + ksi[2] * Jx32) / J2 ]

        # in fact returns an array of 3 elements [1,2,3...]
        return dn_dy
