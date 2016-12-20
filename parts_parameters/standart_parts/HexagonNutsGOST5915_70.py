from parts_parameters.func import *

d_arr = [6, 8, 10, 12, 16, 20, 24]


def S(d):
    return func_from_par(d, d_arr, 10, 13, 17, 19, 24, 30, 36)


def m(d):
    return func_from_par(d, d_arr, 5, 6.5, 8, 10, 13, 16, 19)


class HexagonNutGOST:
    def __init__(self, d, d_thread_in):
        number_type = 'Number'
        mm_unit = 'MilliMeter'
        rad_unit = 'Radian'
        m_unit = 'Meter'
        deg_unit = 'Degrees'
        integer_type = 'Integer'
        nd_unit = 'nd'
        self.d = NXExpression(number_type, 'd', d, mm_unit)
        self.S = NXExpression(number_type, 'S', S(d), mm_unit)
        self.m = NXExpression(number_type, 'm', m(d), mm_unit)
        self.alpha = NXExpression(number_type, 'alpha', 30, deg_unit)
        self.d_thread_in = NXExpression(number_type, 'd_thread_in', d_thread_in, mm_unit)
        self.d_w = NXExpression(number_type, 'd_w', 0.9 * self.S.value, mm_unit)

