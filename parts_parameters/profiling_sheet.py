from parts_parameters.func import *


file = open(os.path.join(os.path.dirname(os.path.dirname(__file__)), 'profiling', 'profiling_results'), 'rb')
stages = pk.load(file)
file.close()


rk_blade = dict()
sa_blade = dict()


rk_blade['z_rk'] = NXExpression(integer_type, 'z_rk', stages[0]['rk']['z'], nd_unit)
rk_blade['b_a_rk'] = NXExpression(number_type, 'b_a_rk', stages[0]['rk']['b_a']*1e3, mm_unit)
rk_blade['r_in'] = NXExpression(number_type, 'r_in', get_blade_section(stages, 0, 0, 'rk').r * 1e3, mm_unit)
rk_blade['r_av'] = NXExpression(number_type, 'r_av', get_blade_section(stages, 0, 2, 'rk').r * 1e3, mm_unit)
rk_blade['r_out'] = NXExpression(number_type, 'r_out', get_blade_section(stages, 0, 4, 'rk').r * 1e3, mm_unit)
rk_blade['r1_rk'] = NXExpression(number_type, 'r1_rk', get_blade_section(stages, 0, 0, 'rk').r1 * 1e3, mm_unit)
for n2, i2 in enumerate(stages[0]['rk']['sections']):
    if n2 == 0 or n2 == 2 or n2 == 4:
        for n3, i3 in enumerate(i2.x_s):
            s_point_name = 'rk_s%s_%s' % (n3, n2)
            k_point_name = 'rk_k%s_%s' % (n3, n2)
            rk_blade[s_point_name] = NXExpression('Point', s_point_name,
                                                  [i2.x_s[n3]*1e3, i2.y_s[n3]*1e3 + n2 * 50, 0], nd_unit)
            rk_blade[k_point_name] = NXExpression('Point', k_point_name,
                                                  [i2.x_k[n3]*1e3, i2.y_k[n3]*1e3 + n2 * 50, 0], nd_unit)

sa_blade['z_sa'] = NXExpression(integer_type, 'z_sa', stages[0]['sa']['z'], nd_unit)
sa_blade['b_a_sa'] = NXExpression(number_type, 'b_a_sa', stages[0]['sa']['b_a']*1e3, mm_unit)
sa_blade['r1_sa'] = NXExpression(number_type, 'r1_sa', get_blade_section(stages, 0, 0, 'sa').r1 * 1e3, mm_unit)
for n2, i2 in enumerate(stages[0]['sa']['sections']):
    if n2 == 0 or n2 == 2 or n2 == 4:
        for n3, i3 in enumerate(i2.x_s):
            s_point_name = 'sa_s%s_%s' % (n3, n2)
            k_point_name = 'sa_k%s_%s' % (n3, n2)
            sa_blade[s_point_name] = NXExpression('Point', s_point_name,
                                                  [i2.x_s[n3]*1e3 - 70, i2.y_s[n3]*1e3 + n2 * 50, 0], nd_unit)
            sa_blade[k_point_name] = NXExpression('Point', k_point_name,
                                                  [i2.x_k[n3]*1e3 - 70, i2.y_k[n3]*1e3 + n2 * 50, 0], nd_unit)


class VelocityTriangles:
    def __init__(self, in_scale=4, av_scale=4, out_scale=4):
        st1 = stages[0]
        self.c1_in = NXExpression(number_type, 'c1_in', st1['c1_in'] / in_scale, mm_unit)
        self.alpha1_in = NXExpression(number_type, 'alpha1_in', np.degrees(st1['alpha1_in']), deg_unit)
        self.w1_in = NXExpression(number_type, 'w1_in', st1['w1_in'] / in_scale, mm_unit)
        self.beta1_in = NXExpression(number_type, 'beta1_in', np.degrees(st1['beta1_in']), deg_unit)
        self.c2_in = NXExpression(number_type, 'c2_in', st1['c2_in'] / in_scale, mm_unit)
        self.alpha2_in = NXExpression(number_type, 'alpha2_in', np.degrees(st1['alpha2_in']), deg_unit)
        self.w2_in = NXExpression(number_type, 'w2_in', st1['w2_in'] / in_scale, mm_unit)
        self.beta2_in = NXExpression(number_type, 'beta2_in', np.degrees(st1['beta2_in']), deg_unit)

        self.c1_av = NXExpression(number_type, 'c1_av', st1['c1_av'] / av_scale, mm_unit)
        self.alpha1_av = NXExpression(number_type, 'alpha1_av', np.degrees(st1['alpha1_av']), deg_unit)
        self.w1_av = NXExpression(number_type, 'w1_av', st1['w1_av'] / av_scale, mm_unit)
        self.beta1_av = NXExpression(number_type, 'beta1_av', np.degrees(st1['beta1_av']), deg_unit)
        self.c2_av = NXExpression(number_type, 'c2_av', st1['c2_av'] / av_scale, mm_unit)
        self.alpha2_av = NXExpression(number_type, 'alpha2_av', np.degrees(st1['alpha2_av']), deg_unit)
        self.w2_av = NXExpression(number_type, 'w2_av', st1['w2_av'] / av_scale, mm_unit)
        self.beta2_av = NXExpression(number_type, 'beta2_av', np.degrees(st1['beta2_av']), deg_unit)

        self.c1_out = NXExpression(number_type, 'c1_out', st1['c1_out'] / out_scale, mm_unit)
        self.alpha1_out = NXExpression(number_type, 'alpha1_out', np.degrees(st1['alpha1_out']), deg_unit)
        self.w1_out = NXExpression(number_type, 'w1_out', st1['w1_out'] / out_scale, mm_unit)
        self.beta1_out = NXExpression(number_type, 'beta1_out', np.degrees(st1['beta1_out']), deg_unit)
        self.c2_out = NXExpression(number_type, 'c2_out', st1['c2_out'] / out_scale, mm_unit)
        self.alpha2_out = NXExpression(number_type, 'alpha2_out', np.degrees(st1['alpha2_out']), deg_unit)
        self.w2_out = NXExpression(number_type, 'w2_out', st1['w2_out'] / out_scale, mm_unit)
        self.beta2_out = NXExpression(number_type, 'beta2_out', np.degrees(st1['beta2_out']), deg_unit)