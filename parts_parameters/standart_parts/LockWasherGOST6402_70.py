from parts_parameters.func import *

d_thread_arr = [2.5, 3, 4, 5, 6, 8, 10]


def d(d_thread):
    return func_from_par(d_thread, d_thread_arr, 2.6, 3.1, 4.1, 5.1, 6.1, 8.2, 10.2)


def b(d_thread):
    return func_from_par(d_thread, d_thread_arr, 0.6, 0.8, 1, 1.2, 1.4, 2, 2.5)


def s(d_thread):
    return func_from_par(d_thread, d_thread_arr, 0.6, 0.8, 1, 1.2, 1.4, 2, 2.5)


class LockWasherGOST:
    def __init__(self, d_thread):
        number_type = 'Number'
        mm_unit = 'MilliMeter'
        rad_unit = 'Radian'
        m_unit = 'Meter'
        deg_unit = 'Degrees'
        integer_type = 'Integer'
        nd_unit = 'nd'
        self.d = NXExpression(number_type, 'd', d(d_thread), mm_unit)
        self.b = NXExpression(number_type, 'b', b(d_thread), mm_unit)
        self.s = NXExpression(number_type, 's', s(d_thread), mm_unit)
        self.m = NXExpression(number_type, 'm', 0.65 * self.s.value, mm_unit)
        self.angle = NXExpression(number_type, 'angle', 70, deg_unit)
