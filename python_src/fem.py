import numpy as np
import mesh_parser as mp


class FEM(object):

    def __init__(self, mesh, fem_imp, U, tau):
        self.mesh = mesh
        self.fem_imp = fem_imp
        self.point_count = self.mesh[mp.MESH_KEY.Points_Number]
        self.n = self.fem_imp.equations_number
        self.K = [0] * (self.n * self.point_count) ** 2
        self.C = [0] * (self.n * self.point_count) ** 2
        self.F = [0] * self.n * self.point_count

        self.U = U
        self.tau = tau
        pass

    def process_element_matrix(self, i):
        cc = 0
        kk = 0
        gr = self.fem_imp.get_gauss_rule(3)
        print("process {0} element".format(i))
        element = self.fem_imp.get_volume_elements()[i]
        elem_size = len(element)
        for g_col in range(0, elem_size):
            for g_row in range(0, elem_size):
                for l_col in range(0, self.n):
                    for l_row in range(0, self.n):
                        if l_col == l_row:
                            cc += self.fem_imp.get_cc(g_col, g_row, l_col, l_row, i, gr.p[0])
                            for l in range(0, gr.int_points):
                                kk += self.fem_imp.get_volume(i) * (self.fem_imp.get_kk(g_col, g_row,
                                                                                         l_col, l_row,
                                                                                         i,
                                                                                         gr.p[l])) * gr.wi[l]
                        else:
                            for l in range(0, gr.int_points):
                                kk += self.fem_imp.get_volume(i) * \
                                       (self.fem_imp.get_kk(g_col, g_row, l_col, l_row, i, gr.p[l]) *
                                        self.fem_imp.get_N(gr.p[l])[g_row] * gr.wi[l])
                            cc = 0

                        idx = (element[g_col] * self.n + l_col) * self.n * self.point_count + \
                                                                    element[g_row] * self.n + l_row
                        self.K[idx] += kk + cc/self.tau
                        self.C[idx] += cc
                        kk = 0
                        cc = 0
        pass

    def assemble_K_matrix(self):
        elements_count = len(self.fem_imp.get_volume_elements())
        if (elements_count <= 0):
            pass
        #     print("started element {0} of {1}".format(i, elements_count))
        #     print("finished element {0} of {1}".format(i, elements_count))
        _ = [self.process_element_matrix(i) for i in range(0, elements_count)]
        pass

    def process_element_vector(self, i):
        #print("started element {0}".format(i))
        gr = self.fem_imp.get_gauss_rule(3)
        element = self.fem_imp.get_volume_elements()[i]
        elem_size = len(element)
        for j in range(0, elem_size):
            for k in range(0, self.n):
                for l in range(0, gr.int_points):
                    self.F[element[j] * self.n + k] -= self.fem_imp.get_N(gr.p[l])[j] * \
                                                       self.fem_imp.get_volume(i) * \
                                                       self.fem_imp.get_ff(element[j], k, i, gr.p[l]) * gr.wi[l]

        #print("finished element {0}".format(i))
        pass

    def assemble_F_vector(self):
        elements_count = len(self.fem_imp.get_volume_elements())
        if (elements_count <= 0):
            pass

        size = self.mesh[mp.MESH_KEY.Points_Number] * self.n
        C1 = np.array(self.C)
        Uc1 = np.array(self.fem_imp.get_Uc())
        self.F = np.matmul(C1.reshape(size, size).swapaxes(0, 1)/self.tau, Uc1)

        _ = [self.process_element_vector(i) for i in range(0, elements_count)]
        pass

    def process_point_border_condition(self, i, time, points_count):
        #print("started border point {0}".format(i))
        for k in range(0, self.n):
            self.F[i * self.n + k] = self.fem_imp.get_border_condition(i, k, time)

        for k in range(0, points_count):
            for ii in range(0, self.n):
                for jj in range(0, self.n):
                    pivot = (points_count * self.n * self.n) * k + self.n * i
                    g_idx = pivot + jj * points_count * self.n + ii
                    if k == i:
                        if ii == jj:
                            self.K[g_idx] = 1
                        else:
                            self.K[g_idx] = 0
                    else:
                        self.K[g_idx] = 0
        #print("finished border point {0}".format(i))
        pass

    def set_border_conditions(self, time):
        points_count = self.mesh[mp.MESH_KEY.Points_Number]
        bolder_points = self.fem_imp.get_border_elements()
        bolder_points = set([j for i in bolder_points for j in i])
        for i in bolder_points:
            self.process_point_border_condition(i, time, points_count)
        pass

    def run_calculations(self, timesteps):

        print("started preparing K matrix")
        self.assemble_K_matrix()
        self.K = np.array(self.K)
        np.save("k_matrix.out", self.K)
        np.save("c_matrix.out", self.C)
        print("finished preparing K matrix")
        # self.K = np.load("k_matrix.out.npy")
        # self.C = np.load("c_matrix.out.npy")
        K1 = self.K
        for step in range(1, timesteps):
            print("start {0} of {1}".format(step, timesteps))
            self.assemble_F_vector()
            self.set_border_conditions(self.tau * step)
            size = self.mesh[mp.MESH_KEY.Points_Number] * self.n
            self.U = np.linalg.solve(np.array(self.K).reshape(size, size).swapaxes(0, 1), np.array(self.F))
            self.U = self.U.reshape(int(len(self.U)/self.n), self.n)
            self.F = [0] * self.n * self.point_count
            self.K = K1
#             #self.K = np.load("/content/drive/My Drive/FEM/k_matrix.out.npy")
            print("finished {0} of {1}".format(step, timesteps))
        pass