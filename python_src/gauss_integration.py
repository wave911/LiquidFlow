import mesh_parser as mp


class GaussRule(object):

    def __init__(self, mesh_type, order):
        self.mesh_type = mesh_type
        self.order = order
        self.int_points = 0
        self.p = {}
        self.__init_points__()

    def __init_points__(self):
        if self.mesh_type == mp.MESH_TYPE.T2D:
            if self.order == 1:
                self.int_points = 1
                self.wi = [1.0]
                self.p[0] = [1.0/3.0, 1.0/3.0, 1.0/3.0]
            if self.order == 2:
                self.int_points = 3
                self.wi = [1.0/3.0, 1.0/3.0, 1.0/3.0]
                self.p[0] = [2.0/3.0, 1.0/6.0, 1.0/6.0]
                self.p[1] = [1.0/6.0, 1.0/6.0, 2.0/3.0]
                self.p[2] = [1.0/6.0, 2.0/3.0, 1.0/6.0]
            if self.order == 3:
                self.int_points = 4
                self.wi = [-0.5625, 0.5208333333333, 0.5208333333333, 0.5208333333333]
                self.p[0] = [1.0/3.0, 1.0/3.0, 1.0/3.0]
                self.p[1] = [0.6, 0.2, 0.2]
                self.p[2] = [0.2, 0.2, 0.6]
                self.p[3] = [0.2, 0.6, 0.2]

        if self.mesh_type == mp.MESH_TYPE.T3D:
            if self.order == 1:
                self.int_points = 1
                self.wi = [1.0]
                self.p[0] = [1.0/4.0, 1.0/4.0, 1.0/4.0, 1.0/4.0]
            if self.order == 2:
                self.int_points = 4
                self.wi = [0.25, 0.25, 0.25, 0.25]
                self.p[0] = [0.585410196624969, 0.138196601125011, 0.138196601125011, 0.138196601125011]
                self.p[1] = [0.138196601125011, 0.138196601125011, 0.138196601125011, 0.585410196624969]
                self.p[2] = [0.138196601125011, 0.138196601125011, 0.585410196624969, 0.138196601125011]
                self.p[3] = [0.138196601125011, 0.585410196624969, 0.138196601125011, 0.138196601125011]
            if self.order == 3:
                self.int_points = 5
                self.wi = [-0.8, 0.45, 0.45, 0.45, 0.45]
                self.p[0] = [0.25, 0.25, 0.25, 0.25]
                self.p[1] = [1.0/2.0, 1.0/6.0, 1.0/6.0, 1.0/6.0]
                self.p[2] = [1.0/6.0, 1.0/6.0, 1.0/6.0, 1.0/2.0]
                self.p[3] = [1.0/6.0, 1.0/6.0, 1.0/2.0, 1.0/6.0]
                self.p[4] = [1.0/6.0, 1.0/2.0, 1.0/6.0, 1.0/6.0]
        pass