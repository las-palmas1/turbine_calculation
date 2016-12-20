from parts_parameters.func import *
import math


d_arr = [27, 30, 33, 36, 39, 42, 45, 48, 52, 56, 60, 64, 80, 85, 90]


def func_from_d(d, *args):
    for i in range(len(d_arr)):
        if d == d_arr[i]:
            return args[i]


def D(d):
    return func_from_d(d, 45, 48, 52, 55, 60, 65, 70, 75, 80, 85, 90, 95, 115, 120, 125)


def m(d):
    return func_from_d(d, 10, 10, 10, 10, 10, 10, 10, 12, 12, 12, 12, 12, 15, 15, 18)


def m1(d):
    return func_from_d(d, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 10, 10, 10)


def D2(d):
    return func_from_d(d, 35, 38, 40, 42, 48, 52, 55, 58, 61, 65, 70, 75, 90, 98, 102)


def b(d):
    return func_from_d(d, 6, 6, 8, 8, 8, 8, 8, 8, 10, 10, 10, 10, 10, 10, 12)


def h(d):
    return func_from_d(d, 2.5, 2.5, 3, 3, 3, 3, 3, 3.5, 3.5, 4, 4, 4, 4, 4, 4)


def c(d):
    return func_from_d(d, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6)


def n(d):
    return func_from_d(d, 4, 4, 4, 4, 4, 4, 6, 6, 6, 6, 6, 6, 6, 6, 6)


class SlottedRoundNutGOST:
    def __init__(self, d, chamfer, series='wide'):
        deg = math.pi / 180
        self.D = NXExpression(number_type, 'D', D(d), mm_unit)
        if series == 'wide':
            self.m = NXExpression(number_type, 'm', m(d), mm_unit)
        elif series == 'narrow':
            self.m = NXExpression(number_type, 'm', m1(d), mm_unit)
        self.D2 = NXExpression(number_type, 'D2', D2(d), mm_unit)
        self.b = NXExpression(number_type, 'b', b(d), mm_unit)
        self.h = NXExpression(number_type, 'h', h(d), mm_unit)
        self.c = NXExpression(number_type, 'c', c(d), mm_unit)
        self.n = NXExpression(integer_type, 'n', n(d), mm_unit)
        self.chamfer = NXExpression(number_type, 'chamfer', chamfer, mm_unit)
        self.d_in = NXExpression(number_type, 'd_in', d, mm_unit)






