from parts_parameters.func import *


d_arr = [27, 30, 33, 36, 39, 42, 45, 48, 52, 56, 60, 64, 80, 85, 90]


def d1(d):
    return func_from_par(d, d_arr, 27.5, 30.5, 33.5, 36.5, 39.5, 42.5, 45.5, 48.5, 52.5, 57, 61, 65, 81, 86, 91)


def d2(d):
    return func_from_par(d, d_arr, 47, 50, 54, 58, 62, 67, 72, 77, 82, 87, 92, 98, 117, 122, 127)


def d3(d):
    return func_from_par(d, d_arr, 36, 39, 42, 45, 48, 52, 56, 60, 65, 70, 75, 80, 100, 105, 110)


def b(d):
    return func_from_par(d, d_arr, 4.8, 4.8, 5.8, 5.8, 5.8, 5.8, 5.8, 7.8, 7.8, 7.8, 7.8, 7.8, 9.5, 9.5, 11.5)


def l(d):
    return func_from_par(d, d_arr, 24, 27, 30, 33, 36, 39, 42, 45, 49, 53, 57, 61, 76, 81, 86)


def h(d):
    return func_from_par(d, d_arr, 6, 6, 6, 6, 6, 6, 6, 6, 7.5, 7.5, 7.5, 7.5, 8.5, 8.5, 9.5)


def r(d):
    return func_from_par(d, d_arr, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 1)


def s(d):
    return func_from_par(d, d_arr, 1, 1, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 2)


def r1(d):
    return 1.5 * s(d)


class MultilegLockWasherGOST:
    def __init__(self, d, nut_h, nut_D):
        number_type = 'Number'
        mm_unit = 'MilliMeter'
        rad_unit = 'Radian'
        m_unit = 'Meter'
        deg_unit = 'Degrees'
        integer_type = 'Integer'
        nd_unit = 'nd'
        deg = np.pi / 180
        self.alpha1 = NXExpression(number_type, 'alpha1', 30, deg_unit)
        self.alpha2 = NXExpression(number_type, 'alpha2', 15, deg_unit)
        self.alpha3 = NXExpression(number_type, 'alpha3', 25, deg_unit)
        self.d1 = NXExpression(number_type, 'd1', d1(d), mm_unit)
        self.d2 = NXExpression(number_type, 'd2', d2(d), mm_unit)
        self.d3 = NXExpression(number_type, 'd3', d3(d), mm_unit)
        self.b = NXExpression(number_type, 'b', b(d), mm_unit)
        self.l = NXExpression(number_type, 'l', l(d), mm_unit)
        self.h = NXExpression(number_type, 'h', h(d), mm_unit)
        self.r = NXExpression(number_type, 'r', r(d), mm_unit)
        self.s = NXExpression(number_type, 's', s(d), mm_unit)
        self.r1 = NXExpression(number_type, 'r1', r1(d), mm_unit)
        self.nut_h = NXExpression(number_type, 'nut_h', nut_h, mm_unit)
        self.nut_D = NXExpression(number_type, 'nut_D', nut_D, mm_unit)




