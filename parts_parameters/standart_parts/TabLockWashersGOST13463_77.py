from parts_parameters.func import *
from parts_parameters.standart_parts.HexagonBoltsGOST7798_70 import S

d_arr = [6, 8, 10, 12, 16, 20, 24]


def d1(d):
    return func_from_par(d, d_arr, 6.4, 8.4, 10.5, 13, 17, 21, 25)


def d2(d):
    return func_from_par(d, d_arr, 10, 13, 17, 19, 24, 30, 36)


def B(d):
    return func_from_par(d, d_arr, 6, 8, 10, 12, 15, 18, 20)


def B1(d):
    return func_from_par(d, d_arr, 7.5, 9, 10, 12, 15, 18, 20)


def L(d):
    return func_from_par(d, d_arr, 18, 20, 22, 28, 32, 36, 42)


def L1(d):
    return func_from_par(d, d_arr, 9, 11, 13, 15, 20, 24, 28)


def s(d):
    return func_from_par(d, d_arr, 0.8, 1, 1, 1, 1, 1, 1)


def r(d):
    return func_from_par(d, d_arr, 0.5, 0.5, 1.2, 1.2, 1.2, 1.2, 1.6)


def r1(d):
    return func_from_par(d, d_arr, 0.5, 1, 1, 2, 2, 2, 3)


def r2(d):
    return func_from_par(d, d_arr, 0.8, 1.2, 1.2, 1.6, 1.6, 2, 2)


class TabLockWasherGOST:
    def __init__(self, d, l1, l2):
        number_type = 'Number'
        mm_unit = 'MilliMeter'
        rad_unit = 'Radian'
        m_unit = 'Meter'
        deg_unit = 'Degrees'
        integer_type = 'Integer'
        nd_unit = 'nd'
        self.d = NXExpression(number_type, 'd', d, mm_unit)
        self.l1 = NXExpression(number_type, 'l1', l1, mm_unit)
        self.l2 = NXExpression(number_type, 'l2', l2, mm_unit)
        self.d1 = NXExpression(number_type, 'd1', d1(d), mm_unit)
        self.d2 = NXExpression(number_type, 'd2', d2(d), mm_unit)
        self.B = NXExpression(number_type, 'B', B(d), mm_unit)
        self.B1 = NXExpression(number_type, 'B1', B1(d), mm_unit)
        self.L = NXExpression(number_type, 'L', L(d), mm_unit)
        self.L11 = NXExpression(number_type, 'L11', L1(d), mm_unit)
        self.s = NXExpression(number_type, 's', s(d), mm_unit)
        self.S_bolt = NXExpression(number_type, 'S_bolt', S(d), mm_unit)
        self.r = NXExpression(number_type, 'r', r(d), mm_unit)
        self.r1 = NXExpression(number_type, 'r1', r1(d), mm_unit)
        self.r2 = NXExpression(number_type, 'r2', r2(d), mm_unit)
        self.alpha = NXExpression(number_type, 'alpha', 60, deg_unit)
        self.r3 = NXExpression(number_type, 'r3', 0.6 * self.s.value, mm_unit)






