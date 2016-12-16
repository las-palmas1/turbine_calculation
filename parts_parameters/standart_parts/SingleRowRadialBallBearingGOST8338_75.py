from parts_parameters.func import *
import math


d_arr = [45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100]


def func_from_d(d, *args):
    for i in range(len(d_arr)):
        if d == d_arr[i]:
            return args[i]


def D_av_ser(d):
    return func_from_d(d, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 215)


def B_av_ser(d):
    return func_from_d(d, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 47)


def r_av_ser(d):
    return func_from_d(d, 2.5, 3, 3, 3.5, 3.5, 3.5, 3.5, 3.5, 4, 4, 4, 4)


def D_light_ser(d):
    return func_from_d(d, 85, 90, 100, 110, 120, 125, 130, 140, 150, 160, 170, 180)


def B_light_ser(d):
    return func_from_d(d, 19, 20, 21, 22, 23, 24, 25, 26, 28, 30, 32, 34)


def r_light_ser(d):
    return func_from_d(d, 2, 2, 2.5, 2.5, 2.5, 2.5, 2.5, 3, 3, 3, 3.5, 3.5)


def D_w(d, D):
    return 0.32 * (D - d)


def S(d, D):
    return 0.15 * (D - d)


def D_pw(d, D):
    return 0.5 * (D + d)


class BallBearing:
    def __init__(self, d, series='light'):
        number_type = 'Number'
        mm_unit = 'MilliMeter'
        rad_unit = 'Radian'
        m_unit = 'Meter'
        deg_unit = 'Degrees'
        integer_type = 'Integer'
        nd_unit = 'nd'
        deg = math.pi / 180
        if series == 'light':
            self.D = NXExpression(number_type, 'D', D_light_ser(d), mm_unit)
            self.B = NXExpression(number_type, 'B', B_light_ser(d), mm_unit)
            self.r = NXExpression(number_type, 'r', r_light_ser(d), mm_unit)
        elif series == 'average':
            self.D = NXExpression(number_type, 'D', D_av_ser(d), mm_unit)
            self.B = NXExpression(number_type, 'B', B_av_ser(d), mm_unit)
            self.r = NXExpression(number_type, 'r', r_av_ser(d), mm_unit)
        self.D_w = NXExpression(number_type, 'D_w', D_w(d, self.D.value), mm_unit)
        self.S = NXExpression(number_type, 'S', S(d, self.D.value), mm_unit)
        self.D_pw = NXExpression(number_type, 'D_pw', D_pw(d, self.D.value), mm_unit)
        self.d_in = NXExpression(number_type, 'd_in', d, mm_unit)
        self.z = NXExpression(integer_type, 'z', int(math.pi * self.D_pw.value / self.D_w.value) - 3, nd_unit)
        self.h_sep = NXExpression(number_type, 'h_sep',
                                  0.8 * (0.5 * self.D.value - 2 * self.S.value - self.d_in.value / 2), mm_unit)
        self.thickness_sep = NXExpression(number_type, 'thickness_sep', 0.7 * 0.5 * (self.B.value - self.D_w.value),
                                          mm_unit)

