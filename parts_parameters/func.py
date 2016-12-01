import os
import numpy as np

number_type = 'Number'
mm_unit = 'MilliMeter'
rad_unit = 'Radian'
m_unit = 'Meter'
deg_unit = 'Degrees'
integer_type = 'Integer'
nd_unit = 'nd'


class NXExpression:
    def __init__(self, exp_type, exp_name, exp_value, exp_unit):
        self.type = exp_type
        self.name = exp_name
        self.value = exp_value
        self.unit = exp_unit


def dict_combination(*args):
    result = dict()
    for i1 in args:
        for i2 in i1:
            result[i2] = i1[i2]
    return result


def number_into_string(exp_value, exp_type):
    if exp_type != 'Point':
        return str(exp_value)
    else:
        return 'Point(' + str(exp_value[0]) + ',' + str(exp_value[1]) + ',' + str(exp_value[2]) + ')'


def add_expressions_into_file(file, expressions_dict):
    for i in expressions_dict:
        file.write(expressions_dict[i].type + ' ' + expressions_dict[i].name + ' ' +
                   number_into_string(expressions_dict[i].value, expressions_dict[i].type) + ' ' +
                   expressions_dict[i].unit + '\n')


def create_expressions_file(exp_filename, part_filename, exp_dict):
    current_dir = os.path.dirname(__file__)
    abs_exp_filename = os.path.join(current_dir, 'output', 'current', exp_filename)
    abs_path_filename = os.path.join(current_dir, 'output', 'model_paths', exp_filename)
    path_file = open(abs_path_filename, 'w')
    path_file.write(os.path.abspath(part_filename) + '\n')
    path_file.close()
    exp_file = open(abs_exp_filename, 'w')
    add_expressions_into_file(exp_file, exp_dict)
    exp_file.close()
    old_versions_dir = os.path.join(os.path.dirname(__file__), 'output', 'old')
    filename_arr = os.listdir(old_versions_dir)
    if filename_arr.count(exp_filename) == 0:
        file = open(os.path.join(old_versions_dir, exp_filename), 'w')
        file.close()


class LockTeethCoordinates:
    def __init__(self, y0, z0, s, r1, r2, h0, phi, gamma, beta, count):
        self.y0 = y0
        self.z0 = z0
        self.s = s
        self.r1 = r1
        self.r2 = r2
        self.h0 = h0
        self.phi = phi
        self.gamma = gamma
        self.beta = beta
        self.count = count
        self.l = s / 2 * np.tan(np.pi - beta - gamma) / (np.tan(np.pi - beta - gamma) + np.tan(beta)) / np.cos(beta) - \
                 r1 / np.tan(gamma / 2)
        self.angle1 = np.pi / 2 + self.phi / 2 - self.beta
        self.angle2 = self.gamma - self.angle1
        self.y1 = self.y0
        self.z1 = self.z0 - self.h0
        self.y2 = self.y1 + self.r2 * (1 - np.sin(self.angle1))
        self.z2 = self.z1 - self.r2 * np.cos(self.angle1)
        self.y3 = self.y2 + self.l * np.cos(self.angle1)
        self.z3 = self.z2 - self.l * np.sin(self.angle1)
        self.y4 = self.y3 + self.l * np.cos(self.angle1)
        self.z4 = self.z3 - self.l * np.sin(self.angle1)
        self.y5 = self.y4 + self.r1 / np.tan(self.gamma / 2) * (np.cos(self.angle1) - np.cos(self.angle2))
        self.z5 = self.z4 - self.r1 / np.tan(self.gamma / 2) * (np.sin(self.angle1) + np.sin(self.angle2))
        self.y6 = self.y3 - 0.5 * self.s * np.sin(self.phi / 2)
        self.z6 = self.z3 - 0.5 * self.s * np.cos(self.phi / 2)
        self.y7 = self.y6 - (self.y5 - self.y6)
        self.z7 = self.z6 - (self.z5 - self.z6)
        self.y_last = self.y3 - self.s / 2 * np.sin(self.phi / 2) * (2 * self.count - 1)
        self.z_last = self.z3 - self.s / 2 * np.cos(self.phi / 2) * (2 * self.count - 1)
        self.y_last_next = self.y2 - self.s * self.count * np.sin(self.phi / 2)
        self.z_last_next = self.z2 - self.s * self.count * np.cos(self.phi / 2)



