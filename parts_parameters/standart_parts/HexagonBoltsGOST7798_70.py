from parts_parameters.func import *
d_arr = [6, 8, 10, 12, 16, 20, 24]


def S(d):
    return func_from_par(d, d_arr, 10, 13, 17, 19, 24, 30, 36)


def k(d):
    return func_from_par(d, d_arr, 4, 5.3, 6.4, 7.5, 10, 12.5, 15)


class HexagonBoltGOST:
    def __init__(self, d, l, b, chamfer):
        # b - длина резьбы
        number_type = 'Number'
        mm_unit = 'MilliMeter'
        rad_unit = 'Radian'
        m_unit = 'Meter'
        deg_unit = 'Degrees'
        integer_type = 'Integer'
        nd_unit = 'nd'
        self.d = NXExpression(number_type, 'd', d, mm_unit)
        self.l = NXExpression(number_type, 'l', l, mm_unit)
        self.b = NXExpression(number_type, 'b', b, mm_unit)
        self.S = NXExpression(number_type, 'S', S(d), mm_unit)
        self.D = NXExpression(number_type, 'D1', 0.8 * self.S.value, mm_unit)
        self.k = NXExpression(number_type, 'k', k(d), mm_unit)
        self.chamfer = NXExpression(number_type, 'chamfer', chamfer, mm_unit)
        self.alpha = NXExpression(number_type, 'alpha', 30, deg_unit)
        self.d_thread_in = NXExpression(number_type, 'd_thread_in', self.d.value -  chamfer, mm_unit)


