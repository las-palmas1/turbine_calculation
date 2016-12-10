import pickle as pk
from parts_parameters.func import *
from profiling.profiling import BladeSection

file = open(os.path.join(os.path.dirname(os.path.dirname(__file__)), 'profiling', 'profiling_results'), 'rb')
stages = pk.load(file)
file.close()

rk_blades = []

for n1, i1 in enumerate(stages):
    rk_blades.append(dict())
    rk_blades[n1]['D1_in'] = NXExpression(number_type, 'D1_in', i1['rk']['D1_in'], m_unit)
    rk_blades[n1]['D2_in'] = NXExpression(number_type, 'D2_in', i1['rk']['D2_in'], m_unit)
    rk_blades[n1]['D1_out'] = NXExpression(number_type, 'D1_out', i1['rk']['D1_out'], m_unit)
    rk_blades[n1]['D2_out'] = NXExpression(number_type, 'D2_out', i1['rk']['D2_out'], m_unit)
    rk_blades[n1]['z'] = NXExpression(integer_type, 'z', i1['rk']['z'], nd_unit)
    rk_blades[n1]['b_a'] = NXExpression(number_type, 'b_a', i1['rk']['b_a'], m_unit)
    rk_blades[n1]['delta_r'] = NXExpression(number_type, 'delta_r', i1['rk']['delta_r']*1e3, mm_unit)
    for n2, i2 in enumerate(i1['rk']['sections']):
        for n3, i3 in enumerate(i2.x_s):
            s_point_name = 's%s_%s' % (n3, n2)
            k_point_name = 'k%s_%s' % (n3, n2)
            rk_blades[n1][s_point_name] = NXExpression('Point', s_point_name, [i2.x_s[n3]*1e3, i2.y_s[n3]*1e3,
                                                                               i2.r*1e3], nd_unit)
            rk_blades[n1][k_point_name] = NXExpression('Point', k_point_name, [i2.x_k[n3]*1e3, i2.y_k[n3]*1e3,
                                                                               i2.r*1e3], nd_unit)


sa_blades = []

for n1, i1 in enumerate(stages):
    sa_blades.append(dict())
    sa_blades[n1]['D1_in'] = NXExpression(number_type, 'D1_in', i1['sa']['D0_in'], m_unit)
    sa_blades[n1]['D2_in'] = NXExpression(number_type, 'D2_in', i1['sa']['D05_in'], m_unit)
    sa_blades[n1]['D1_out'] = NXExpression(number_type, 'D1_out', i1['sa']['D0_out'], m_unit)
    sa_blades[n1]['D2_out'] = NXExpression(number_type, 'D2_out', i1['sa']['D05_out'], m_unit)
    sa_blades[n1]['z'] = NXExpression(integer_type, 'z', i1['sa']['z'], nd_unit)
    sa_blades[n1]['b_a'] = NXExpression(number_type, 'b_a', i1['sa']['b_a'], m_unit)
    for n2, i2 in enumerate(i1['sa']['sections']):
        for n3, i3 in enumerate(i2.x_s):
            s_point_name = 's%s_%s' % (n3, n2)
            k_point_name = 'k%s_%s' % (n3, n2)
            sa_blades[n1][s_point_name] = NXExpression('Point', s_point_name, [i2.x_s[n3]*1e3, i2.y_s[n3]*1e3,
                                                                               i2.r*1e3], nd_unit)
            sa_blades[n1][k_point_name] = NXExpression('Point', k_point_name, [i2.x_k[n3]*1e3, i2.y_k[n3]*1e3,
                                                                               i2.r*1e3], nd_unit)


def angle_rotate(blade_section: BladeSection, b1, **kwargs):
    if 'alpha' not in kwargs:
        rot_angle = np.pi / 2 - blade_section.alpha
        rot_matrix = np.array([[np.cos(rot_angle), np.sin(rot_angle)],
                              [-np.sin(rot_angle), np.cos(rot_angle)]])
        y_k_new = np.dot(rot_matrix, [blade_section.x_k, blade_section.y_k])[1]
        y_s_new = np.dot(rot_matrix, [blade_section.x_s, blade_section.y_s])[1]
        chord_blade = (max(y_k_new) - min(y_s_new)) / np.sin(blade_section.alpha)
        chord1 = blade_section.y_k[0] - b1 / np.tan(blade_section.alpha)
        result = (0.5 * chord_blade - chord1) / blade_section.r
    else:
        rot_angle = np.pi / 2 - kwargs['alpha']
        rot_matrix = np.array([[np.cos(rot_angle), np.sin(rot_angle)],
                               [-np.sin(rot_angle), np.cos(rot_angle)]])
        y_k_new = np.dot(rot_matrix, [blade_section.x_k, blade_section.y_k])[1]
        y_s_new = np.dot(rot_matrix, [blade_section.x_s, blade_section.y_s])[1]
        chord_blade = (max(y_k_new) - min(y_s_new)) / np.sin(kwargs['alpha'])
        chord1 = max(blade_section.y_k[0] - b1 / np.tan(kwargs['alpha']),
                     blade_section.y_k[len(blade_section.y_k) - 1] - (b1 + blade_section.b_a) / np.tan(kwargs['alpha']))
        result = (0.5 * chord_blade - chord1) / blade_section.r
    return result


class StageTail:
    def __init__(self, n, b1_rel=0.6, b2_rel=0.5, br_rel=0.2, b_a_tail_rel=0.8, w2_rel=0.6, teeth_count=2,
                 s=5, r1=0.6, r2=1.1, delta_D=3, phi=40, gamma=65, beta=70):
        deg = np.pi / 180
        self.delta_a_sa = NXExpression(number_type, 'delta_a_sa', stages[n]['rk']['delta_a_sa']*1e3, mm_unit)
        self.delta_a_rk = NXExpression(number_type, 'delta_a_rk', stages[n]['rk']['delta_a_rk']*1e3, mm_unit)
        b1_rel = b1_rel
        b2_rel = b2_rel
        self.b_a1 = NXExpression(number_type, 'b_a1', stages[n]['rk']['b_a']*1e3 + self.delta_a_sa.value * b1_rel +
                                 self.delta_a_rk.value * b2_rel, mm_unit)
        self.b1 = NXExpression(number_type, 'b1', self.delta_a_sa.value * b1_rel, mm_unit)
        self.b2 = NXExpression(number_type, 'b2', self.delta_a_rk.value * b2_rel, mm_unit)
        self.h = NXExpression(number_type, 'h', 2, mm_unit)
        self.r3 = NXExpression(number_type, 'r3', 4, mm_unit)
        br_rel = br_rel
        self.br = NXExpression(number_type, 'br', self.b1.value * br_rel, mm_unit)
        b_a_tail_rel = b_a_tail_rel
        self.b_a_tail = NXExpression(number_type, 'b_a_tail', self.b_a1.value * b_a_tail_rel, mm_unit)
        self.b3 = NXExpression(number_type, 'b3', (self.b_a1.value - self.b_a_tail.value) / 2, mm_unit)
        gamma_in = np.arctan(0.5 * (stages[n]['rk']['D1_in'] - stages[n]['rk']['D2_in']) / stages[n]['rk']['b_a'])
        self.w1 = NXExpression(number_type, 'w1', 1, mm_unit)
        self.D1_tail = NXExpression(number_type, 'D1_tail', stages[n]['rk']['D1_in']*1e3 +
                                    2 * (self.b1.value - self.b3.value) * np.tan(gamma_in) -
                                    self.w1.value / np.cos(gamma_in), mm_unit)
        self.D2_tail = NXExpression(number_type, 'D2_tail', self.D1_tail.value -
                                    2 * self.b_a_tail.value * np.tan(gamma_in), mm_unit)
        self.psi = NXExpression(number_type, 'psi', 360 / stages[n]['rk']['z'], deg_unit)
        self.l1 = NXExpression(number_type, 'l1', 1.5, mm_unit)
        self.l3 = NXExpression(number_type, 'l3', self.b_a1.value - self.b_a_tail.value - self.b3.value, mm_unit)
        self.D1 = NXExpression(number_type, 'D1', self.D2_tail.value - 5, mm_unit)
        self.r4 = NXExpression(number_type, 'r4', 0.5, mm_unit)
        self.r5 = NXExpression(number_type, 'r5', 0.5, mm_unit)
        self.delta_D = NXExpression(number_type, 'delta_D', delta_D, mm_unit)
        self.alpha = NXExpression(number_type, 'alpha', stages[n]['rk']['sections'][0].alpha / deg, deg_unit)
        w2_rel = w2_rel
        self.w2 = NXExpression(number_type, 'w2', w2_rel * 0.5 * self.D1_tail.value * np.radians(self.psi.value),
                               mm_unit)
        self.z01 = NXExpression(number_type, 'z01', np.sqrt((0.5 * self.D1_tail.value)**2 - (0.5 * self.w2.value)**2),
                                mm_unit)
        self.z02 = NXExpression(number_type, 'z02', np.sqrt((0.5 * self.D2_tail.value) ** 2 - (0.5 * self.w2.value) ** 2),
                                mm_unit)
        self.z0 = NXExpression(number_type, 'z0', min(self.z01.value, self.z02.value) - self.delta_D.value, mm_unit)
        self.y0 = NXExpression(number_type, 'y0', 0.5 * self.w2.value, mm_unit)
        self.s = NXExpression(number_type, 's', s, mm_unit)
        self.teeth_count = NXExpression(integer_type, 'teeth_count', teeth_count, nd_unit)
        self.r1 = NXExpression(number_type, 'r1', r1, mm_unit)
        self.r2 = NXExpression(number_type, 'r2', r2, mm_unit)
        self.phi = NXExpression(number_type, 'phi', phi, deg_unit)
        self.gamma = NXExpression(number_type, 'gamma', gamma, deg_unit)
        self.beta = NXExpression(number_type, 'beta', beta, deg_unit)
        self.h0 = NXExpression(number_type, 'h0', 1, mm_unit)
        lock_coord = BladeLockTeethCoordinates(self.y0.value, self.z0.value, self.s.value, self.r1.value, self.r2.value,
                                               self.h0.value, np.radians(self.phi.value), np.radians(self.gamma.value),
                                               np.radians(self.beta.value), self.teeth_count.value)
        self.y1 = NXExpression(number_type, 'y1', lock_coord.y1, mm_unit)
        self.z1 = NXExpression(number_type, 'z1', lock_coord.z1, mm_unit)
        self.y2 = NXExpression(number_type, 'y2', lock_coord.y2, mm_unit)
        self.z2 = NXExpression(number_type, 'z2', lock_coord.z2, mm_unit)
        self.y3 = NXExpression(number_type, 'y3', lock_coord.y3, mm_unit)
        self.z3 = NXExpression(number_type, 'z3', lock_coord.z3, mm_unit)
        self.y4 = NXExpression(number_type, 'y4', lock_coord.y4, mm_unit)
        self.z4 = NXExpression(number_type, 'z4', lock_coord.z4, mm_unit)
        self.y5 = NXExpression(number_type, 'y5', lock_coord.y5, mm_unit)
        self.z5 = NXExpression(number_type, 'z5', lock_coord.z5, mm_unit)
        self.y6 = NXExpression(number_type, 'y6', lock_coord.y6, mm_unit)
        self.z6 = NXExpression(number_type, 'z6', lock_coord.z6, mm_unit)
        self.y7 = NXExpression(number_type, 'y7', lock_coord.y7, mm_unit)
        self.z7 = NXExpression(number_type, 'z7', lock_coord.z7, mm_unit)
        self.y_last = NXExpression(number_type, 'y_last', lock_coord.y_last, mm_unit)
        self.z_last = NXExpression(number_type, 'z_last', lock_coord.z_last, mm_unit)
        self.y_last_next = NXExpression(number_type, 'y_last_next', lock_coord.y_last_next, mm_unit)
        self.z_last_next = NXExpression(number_type, 'z_last_next', lock_coord.z_last_next, mm_unit)
        ang_rot = angle_rotate(stages[n]['rk']['sections'][0], self.b1.value / 1e3)
        self.angle_rotation = NXExpression(number_type, 'angle_rotation', ang_rot / deg, deg_unit)
        self.r6 = NXExpression(number_type, 'r6', 1.2, mm_unit)
        self.r_blade = NXExpression(number_type, 'r_blade', 1, mm_unit)
        self.b1_out = NXExpression(number_type, 'b1_out', 1.2, mm_unit)
        self.alpha_out = NXExpression(number_type, 'alpha_out',
                                      self.alpha.value,
                                      deg_unit)
        ang_rot_out = angle_rotate(stages[n]['rk']['sections'][len(stages[n]['rk']['sections']) - 1],
                                   self.b1_out.value / 1e3, alpha=np.radians(self.alpha_out.value))
        self.ang_rotation_out = NXExpression(number_type, 'angle_rotation_out', ang_rot_out / deg, deg_unit)
        self.d1_tail_ledge = NXExpression(number_type, 'd1_tail_ledge', 2 * self.z_last.value + 6, mm_unit)
        self.d2_tail_ledge = NXExpression(number_type, 'd2_tail_ledge', self.d1_tail_ledge.value + 6, mm_unit)
        self.angle_tail_ledge = NXExpression(number_type, 'angle_tail_ledge', 30, deg_unit)
        self.l1_tail_ledge = NXExpression(number_type, 'l1_tail_ledge', 4, mm_unit)
        self.l2_tail_ledge = NXExpression(number_type, 'l2_tail_ledge', self.l1_tail_ledge.value / 2, mm_unit)
        self.t1_tail_ledge = NXExpression(number_type, 't1_tail_ledge', 0.5 * self.d1_tail_ledge.value -
                                          self.z_last.value, mm_unit)


first_stage_tail = StageTail(0, teeth_count=2, s=4, delta_D=4, w2_rel=0.5)
second_stage_tail = StageTail(1, teeth_count=2, s=4, delta_D=4, w2_rel=0.5)
stage_tails = [first_stage_tail.__dict__, second_stage_tail.__dict__]
