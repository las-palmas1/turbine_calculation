import pickle as pk
from parts_parameters.func import *

file = open(os.path.join(os.path.dirname(os.path.dirname(__file__)), 'profiling', 'profiling_results'), 'rb')
stages = pk.load(file)
file.close()

rk_blades = []

for n1, i1 in enumerate(stages):
    rk_blades.append(dict())
    rk_blades[n1]['D1_in'] = NXExpression(number_type, 'D1_in', i1['rk']['D1_in'], mm_unit)
    rk_blades[n1]['D2_in'] = NXExpression(number_type, 'D2_in', i1['rk']['D2_in'], mm_unit)
    rk_blades[n1]['D1_out'] = NXExpression(number_type, 'D1_out', i1['rk']['D1_out'], mm_unit)
    rk_blades[n1]['D2_out'] = NXExpression(number_type, 'D2_out', i1['rk']['D2_out'], mm_unit)
    rk_blades[n1]['z'] = NXExpression(integer_type, 'z', i1['rk']['z'], nd_unit)
    rk_blades[n1]['b_a'] = NXExpression(number_type, 'b_a', i1['rk']['b_a'], mm_unit)
    for n2, i2 in enumerate(i1['rk']['sections']):
        for n3, i3 in enumerate(i2.x_s):
            s_point_name = 's%s_%s' % (n3, n2)
            k_point_name = 'k%s_%s' % (n3, n2)
            rk_blades[n1][s_point_name] = NXExpression('Point', s_point_name, [i2.x_s[n3], i2.y_s[n3], i2.r], mm_unit)
            rk_blades[n1][k_point_name] = NXExpression('Point', k_point_name, [i2.x_k[n3], i2.y_k[n3], i2.r], mm_unit)


sa_blades = []

for n1, i1 in enumerate(stages):
    sa_blades.append(dict())
    sa_blades[n1]['D1_in'] = NXExpression(number_type, 'D1_in', i1['sa']['D0_in'], mm_unit)
    sa_blades[n1]['D2_in'] = NXExpression(number_type, 'D2_in', i1['sa']['D05_in'], mm_unit)
    sa_blades[n1]['D1_out'] = NXExpression(number_type, 'D1_out', i1['sa']['D0_out'], mm_unit)
    sa_blades[n1]['D2_out'] = NXExpression(number_type, 'D2_out', i1['sa']['D05_out'], mm_unit)
    sa_blades[n1]['z'] = NXExpression(integer_type, 'z', i1['rk']['z'], nd_unit)
    sa_blades[n1]['b_a'] = NXExpression(number_type, 'b_a', i1['rk']['b_a'], mm_unit)
    for n2, i2 in enumerate(i1['sa']['sections']):
        for n3, i3 in enumerate(i2.x_s):
            s_point_name = 's%s_%s' % (n3, n2)
            k_point_name = 'k%s_%s' % (n3, n2)
            sa_blades[n1][s_point_name] = NXExpression('Point', s_point_name, [i2.x_s[n3], i2.y_s[n3], i2.r], mm_unit)
            sa_blades[n1][k_point_name] = NXExpression('Point', k_point_name, [i2.x_k[n3], i2.y_k[n3], i2.r], mm_unit)
