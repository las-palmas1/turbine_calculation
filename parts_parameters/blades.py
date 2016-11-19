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


def angle_rotate(blade_section: BladeSection, b1):
    rot_angle = np.pi / 2 - blade_section.alpha
    rot_matrix = np.array([[np.cos(rot_angle), np.sin(rot_angle)],
                          [-np.sin(rot_angle), np.cos(rot_angle)]])
    y_k_new = np.dot(rot_matrix, [blade_section.x_k, blade_section.y_k])[1]
    y_s_new = np.dot(rot_matrix, [blade_section.x_s, blade_section.y_s])[1]
    chord_blade = (max(y_k_new) - min(y_s_new)) / np.sin(blade_section.alpha)
    chord1 = blade_section.y_k[0] - b1 / np.tan(blade_section.alpha)
    result = (0.5 * chord_blade - chord1) / blade_section.r
    return result


class FirstStageTail:
    def __init__(self):
        deg = np.pi / 180
        self.delta_a_sa = NXExpression(number_type, 'delta_a_sa', stages[0]['rk']['delta_a_sa']*1e3, mm_unit)
        self.delta_a_rk = NXExpression(number_type, 'delta_a_rk', stages[0]['rk']['delta_a_rk']*1e3, mm_unit)
        self.D_tail_in = NXExpression(number_type, 'D_tail_in', min(stages[0]['rk']['D1_in'], stages[0]['rk']['D2_in'])*1e3 - 8,
                                 mm_unit)
        b1_rel = 0.3
        b2_rel = 0.3
        self.b_a1 = NXExpression(number_type, 'b_a1', stages[0]['rk']['b_a']*1e3 + self.delta_a_sa.value * b1_rel +
                                 self.delta_a_rk.value * b2_rel, mm_unit)
        self.b1 = NXExpression(number_type, 'b1', self.delta_a_sa.value * b1_rel, mm_unit)
        b_a_tail_rel = 0.8
        self.b_a_tail = NXExpression(number_type, 'b_a_tail', self.b_a1.value * b_a_tail_rel, mm_unit)
        gamma_in = np.arctan((stages[0]['rk']['D1_in'] - stages[0]['rk']['D2_in']) / stages[0]['rk']['b_a'])
        self.D1_tail_out = NXExpression(number_type, 'D1_tail_out', stages[0]['rk']['D1_in']*1e3 + self.b1.value * np.tan(gamma_in),
                                   mm_unit)
        self.D2_tail_out = NXExpression(number_type, 'D2_tail_out', self.D1_tail_out.value - self.b_a1.value * np.tan(gamma_in),
                                   mm_unit)
        self.alpha = NXExpression(number_type, 'alpha', stages[0]['rk']['sections'][0].alpha / deg, deg_unit)
        self.psi = NXExpression(number_type, 'psi', 360 / stages[0]['rk']['z'], deg_unit)
        self.theta = NXExpression(number_type, 'theta', 0.9 * self.psi.value, deg_unit)
        self.s = NXExpression(number_type, 's', 10, mm_unit)
        self.teeth_count = NXExpression(integer_type, 'teeth_count', 3, nd_unit)
        self.r1 = NXExpression(number_type, 'r1', 0.9, mm_unit)
        self.phi = NXExpression(number_type, 'phi', 40, deg_unit)
        self.gamma = NXExpression(number_type, 'gamma', 30, deg_unit)
        self.beta = NXExpression(number_type, 'beta', 30, deg_unit)
        self.y0 = NXExpression(number_type, 'y0', self.D_tail_in.value / 2 * np.sin(np.radians(self.theta.value) / 2),
                               mm_unit)
        self.z0 = NXExpression(number_type, 'z0', self.D_tail_in.value / 2 * np.cos(np.radians(self.theta.value) / 2),
                               mm_unit)
        angle1 = np.pi / 2 - self.phi.value * deg / 2 - self.beta.value * deg
        self.z1 = NXExpression(number_type, 'z1', self.z0.value - 2, mm_unit)
        self.y1 = NXExpression(number_type, 'y1', self.y0.value - (self.z0.value - self.z1.value) / np.tan(angle1),
                               mm_unit)
        lock_teeth = LockTeethCoordinates(self.y1.value, self.z1.value, self.s.value, self.r1.value,
                                          np.radians(self.phi.value), np.radians(self.gamma.value),
                                          np.radians(self.beta.value), self.teeth_count.value)
        self.y2 = NXExpression(number_type, 'y2', lock_teeth.y2, mm_unit)
        self.z2 = NXExpression(number_type, 'z2', lock_teeth.z2, mm_unit)
        self.y3 = NXExpression(number_type, 'y3', lock_teeth.y3, mm_unit)
        self.z3 = NXExpression(number_type, 'z3', lock_teeth.z3, mm_unit)
        self.y4 = NXExpression(number_type, 'y4', lock_teeth.y4, mm_unit)
        self.z4 = NXExpression(number_type, 'z4', lock_teeth.z4, mm_unit)
        self.y5 = NXExpression(number_type, 'y5', lock_teeth.y5, mm_unit)
        self.z5 = NXExpression(number_type, 'z5', lock_teeth.z5, mm_unit)
        self.y6 = NXExpression(number_type, 'y6', lock_teeth.y6, mm_unit)
        self.z6 = NXExpression(number_type, 'z6', lock_teeth.z6, mm_unit)
        self.y7 = NXExpression(number_type, 'y7', lock_teeth.y7, mm_unit)
        self.z7 = NXExpression(number_type, 'z7', lock_teeth.z7, mm_unit)
        self.y_last = NXExpression(number_type, 'y_last', lock_teeth.y_last, mm_unit)
        self.z_last = NXExpression(number_type, 'z_last', lock_teeth.z_last, mm_unit)
        ang_rot = angle_rotate(stages[0]['rk']['sections'][0], self.b1.value / 1e3)
        self.angle_rotation = NXExpression(number_type, 'angle_rotation', ang_rot / deg, deg_unit)


class SecondStageTail:
    def __init__(self):
        deg = np.pi / 180
        self.delta_a_sa = NXExpression(number_type, 'delta_a_sa', stages[1]['rk']['delta_a_sa']*1e3, mm_unit)
        self.delta_a_rk = NXExpression(number_type, 'delta_a_rk', stages[1]['rk']['delta_a_rk']*1e3, mm_unit)
        self.D_tail_in = NXExpression(number_type, 'D_tail_in', min(stages[1]['rk']['D1_in'], stages[1]['rk']['D2_in']) * 1e3 - 8,
                                      mm_unit)
        b1_rel = 0.3
        b2_rel = 0.3
        self.b_a1 = NXExpression(number_type, 'b_a1', stages[1]['rk']['b_a']*1e3 + self.delta_a_sa.value * b1_rel +
                                 self.delta_a_rk.value * b2_rel, mm_unit)
        self.b1 = NXExpression(number_type, 'b1', self.delta_a_sa.value * b1_rel, mm_unit)
        b_a_tail_rel = 0.8
        self.b_a_tail = NXExpression(number_type, 'b_a_tail', self.b_a1.value * b_a_tail_rel, mm_unit)
        gamma_in = np.arctan((stages[1]['rk']['D1_in'] - stages[1]['rk']['D2_in']) / stages[1]['rk']['b_a'])
        self.D1_tail_out = NXExpression(number_type, 'D1_tail_out', stages[1]['rk']['D1_in'] * 1e3 + self.b1.value * np.tan(gamma_in),
                                        mm_unit)
        self.D2_tail_out = NXExpression(number_type, 'D2_tail_out', self.D1_tail_out.value - self.b_a1.value * np.tan(gamma_in),
                                        mm_unit)
        self.alpha = NXExpression(number_type, 'alpha', stages[1]['rk']['sections'][0].alpha / deg, deg_unit)
        self.psi = NXExpression(number_type, 'psi', 360 / stages[1]['rk']['z'], deg_unit)
        self.theta = NXExpression(number_type, 'theta', 0.9 * self.psi.value, deg_unit)
        self.s = NXExpression(number_type, 's', 10, mm_unit)
        self.teeth_count = NXExpression(integer_type, 'teeth_count', 3, nd_unit)
        self.r1 = NXExpression(number_type, 'r1', 0.9, mm_unit)
        self.phi = NXExpression(number_type, 'phi', 40, deg_unit)
        self.gamma = NXExpression(number_type, 'gamma', 30, deg_unit)
        self.beta = NXExpression(number_type, 'beta', 30, deg_unit)
        self.y0 = NXExpression(number_type, 'y0', self.D_tail_in.value / 2 * np.sin(np.radians(self.theta.value) / 2),
                               mm_unit)
        self.z0 = NXExpression(number_type, 'z0', self.D_tail_in.value / 2 * np.cos(np.radians(self.theta.value) / 2),
                               mm_unit)
        angle1 = np.pi / 2 - self.phi.value * deg / 2 - self.beta.value * deg
        self.z1 = NXExpression(number_type, 'z1', self.z0.value - 2, mm_unit)
        self.y1 = NXExpression(number_type, 'y1', self.y0.value - (self.z0.value - self.z1.value) / np.tan(angle1),
                               mm_unit)
        lock_teeth = LockTeethCoordinates(self.y1.value, self.z1.value, self.s.value, self.r1.value,
                                          np.radians(self.phi.value), np.radians(self.gamma.value),
                                          np.radians(self.beta.value), self.teeth_count.value)
        self.y2 = NXExpression(number_type, 'y2', lock_teeth.y2, mm_unit)
        self.z2 = NXExpression(number_type, 'z2', lock_teeth.z2, mm_unit)
        self.y3 = NXExpression(number_type, 'y3', lock_teeth.y3, mm_unit)
        self.z3 = NXExpression(number_type, 'z3', lock_teeth.z3, mm_unit)
        self.y4 = NXExpression(number_type, 'y4', lock_teeth.y4, mm_unit)
        self.z4 = NXExpression(number_type, 'z4', lock_teeth.z4, mm_unit)
        self.y5 = NXExpression(number_type, 'y5', lock_teeth.y5, mm_unit)
        self.z5 = NXExpression(number_type, 'z5', lock_teeth.z5, mm_unit)
        self.y6 = NXExpression(number_type, 'y6', lock_teeth.y6, mm_unit)
        self.z6 = NXExpression(number_type, 'z6', lock_teeth.z6, mm_unit)
        self.y7 = NXExpression(number_type, 'y7', lock_teeth.y7, mm_unit)
        self.z7 = NXExpression(number_type, 'z7', lock_teeth.z7, mm_unit)
        self.y_last = NXExpression(number_type, 'y_last', lock_teeth.y_last, mm_unit)
        self.z_last = NXExpression(number_type, 'z_last', lock_teeth.z_last, mm_unit)
        ang_rot = angle_rotate(stages[1]['rk']['sections'][0], self.b1.value / 1e3)
        self.angle_rotation = NXExpression(number_type, 'angle_rotation', ang_rot / deg, deg_unit)