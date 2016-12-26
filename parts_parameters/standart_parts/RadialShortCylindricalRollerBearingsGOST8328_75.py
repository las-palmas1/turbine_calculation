from parts_parameters.func import *
import math


d_arr = [40, 45, 50, 55, 60, 65, 70, 75, 80]


def func_from_d(d, *args):
    for i in range(len(d_arr)):
        if d == d_arr[i]:
            return args[i]


def D_av_ser(d):
    return func_from_d(d, 90, 100, 110, 120, 130, 140, 150, 160, 170)


def B_av_ser(d):
    return func_from_d(d, 23, 25, 27, 29, 31, 33, 35, 37, 39)


def r_av_ser(d):
    return func_from_d(d, 2.5, 2.5, 3, 3, 3.5, 3.5, 3.5, 3.5, 3.5)


def r1_av_ser(d):
    return func_from_d(d, 2.5, 2.5, 3, 3, 3.5, 3.5, 3.5, 3.5, 3.5)


def D_light_ser(d):
    return func_from_d(d, 80, 85, 90, 100, 110, 120, 125, 130, 140)


def B_light_ser(d):
    return func_from_d(d, 18, 19, 20, 21, 22, 23, 24, 25, 26)


def r_light_ser(d):
    return func_from_d(d, 2, 2, 2, 2.5, 2.5, 2.5, 2.5, 2.5, 3)


def r1_light_ser(d):
    return func_from_d(d, 2, 2, 2, 2.5, 2.5, 2.5, 2.5, 2.5, 3)


def S1(d, D):
    return 0.1 * (D - d)


def S(d, D):
    return 0.16 * (D - d)


def D_we(d, D):
    return 0.28 * (D - d)


def L_we(d, D):
    return D_we(d, D)


def D_pw(d, D):
    return 0.5 * (d + D)


class RollerBearing:
    def __init__(self, d, series='light'):
        deg = math.pi / 180
        self.d_in = NXExpression(number_type, 'd_in', d, mm_unit)
        if series == 'light':
            self.D = NXExpression(number_type, 'D', D_light_ser(d), mm_unit)
            self.B = NXExpression(number_type, 'B', B_light_ser(d), mm_unit)
            self.r = NXExpression(number_type, 'r', r_light_ser(d), mm_unit)
            self.r1 = NXExpression(number_type, 'r1', r1_light_ser(d), mm_unit)
        elif series == 'average':
            self.D = NXExpression(number_type, 'D', D_av_ser(d), mm_unit)
            self.B = NXExpression(number_type, 'B', B_av_ser(d), mm_unit)
            self.r = NXExpression(number_type, 'r', r_av_ser(d), mm_unit)
            self.r1 = NXExpression(number_type, 'r1', r1_av_ser(d), mm_unit)
        self.S = NXExpression(number_type, 'S', S(d, self.D.value), mm_unit)
        self.S1 = NXExpression(number_type, 'S1', S1(d, self.D.value), mm_unit)
        self.D_we = NXExpression(number_type, 'D_we', D_we(d, self.D.value), mm_unit)
        self.L_we = NXExpression(number_type, 'L_we', L_we(d, self.D.value), mm_unit)
        self.D_pw = NXExpression(number_type, 'D_pw', D_pw(d, self.D.value), mm_unit)
        self.z = NXExpression(integer_type, 'z', int(math.pi * self.D_pw.value / self.D_we.value) - 4, nd_unit)
        self.rad = NXExpression(number_type, 'rad', 0.5, mm_unit)
        self.delta = NXExpression(number_type, 'delta', self.rad.value, mm_unit)
        self.h_sep = NXExpression(number_type, 'h_sep', 0.4 * self.D_we.value, mm_unit)
        self.thickness_sep = NXExpression(number_type, 'thickness_sep', 2, mm_unit)
