from parts_parameters.blades import first_stage_tail, second_stage_tail, rk_blades, sa_blades, stages
from parts_parameters.func import *


turbine = get_average_streamline_calculation_results(os.path.join(os.path.dirname(os.path.dirname(__file__)),
                                                                  'average_streamline_calculation',
                                                                  'average_streamline_calculation_results'))


class FirstStageDisk:
    def __init__(self):
        blade_teeth = BladeLockTeethCoordinates(first_stage_tail.y0.value, first_stage_tail.z0.value,
                                                first_stage_tail.s.value, first_stage_tail.r1.value,
                                                first_stage_tail.r2.value, first_stage_tail.h0.value,
                                                np.radians(first_stage_tail.phi.value),
                                                np.radians(first_stage_tail.gamma.value),
                                                np.radians(first_stage_tail.beta.value),
                                                first_stage_tail.teeth_count.value)
        self.z_blade = NXExpression(integer_type, 'z_blade', rk_blades[0]['z'].value, nd_unit)
        self.alpha = first_stage_tail.alpha
        self.s = first_stage_tail.s
        self.r1 = first_stage_tail.r1
        self.phi = first_stage_tail.phi
        self.gamma = first_stage_tail.gamma
        self.beta = first_stage_tail.beta
        self.teeth_count = first_stage_tail.teeth_count
        k = 0.05
        self.h1 = NXExpression(number_type, 'h1', 2, mm_unit)
        self.b_a_tail = first_stage_tail.b_a_tail
        disk_teeth = DiskLockTeethCoordinates(blade_teeth, k, self.h1.value)
        self.angle1 = NXExpression(number_type, 'angle1', np.degrees(blade_teeth.angle1), deg_unit)
        self.ym2 = NXExpression(number_type, 'ym2', disk_teeth.ym2, mm_unit)
        self.zm2 = NXExpression(number_type, 'zm2', disk_teeth.zm2, mm_unit)
        self.ym1 = NXExpression(number_type, 'ym1', disk_teeth.ym1, mm_unit)
        self.zm1 = NXExpression(number_type, 'zm1', disk_teeth.zm1, mm_unit)
        self.y0 = NXExpression(number_type, 'y0', disk_teeth.y0, mm_unit)
        self.z0 = NXExpression(number_type, 'z0', disk_teeth.z0, mm_unit)
        self.y1 = NXExpression(number_type, 'y1', disk_teeth.y1, mm_unit)
        self.z1 = NXExpression(number_type, 'z1', disk_teeth.z1, mm_unit)
        self.y2 = NXExpression(number_type, 'y2', disk_teeth.y2, mm_unit)
        self.z2 = NXExpression(number_type, 'z2', disk_teeth.z2, mm_unit)
        self.y3 = NXExpression(number_type, 'y3', disk_teeth.y3, mm_unit)
        self.z3 = NXExpression(number_type, 'z3', disk_teeth.z3, mm_unit)
        self.y4 = NXExpression(number_type, 'y4', disk_teeth.y4, mm_unit)
        self.z4 = NXExpression(number_type, 'z4', disk_teeth.z4, mm_unit)
        self.y5 = NXExpression(number_type, 'y5', disk_teeth.y5, mm_unit)
        self.z5 = NXExpression(number_type, 'z5', disk_teeth.z5, mm_unit)
        self.y6 = NXExpression(number_type, 'y6', disk_teeth.y6, mm_unit)
        self.z6 = NXExpression(number_type, 'z6', disk_teeth.z6, mm_unit)
        self.y7 = NXExpression(number_type, 'y7', disk_teeth.y7, mm_unit)
        self.z7 = NXExpression(number_type, 'z7', disk_teeth.z7, mm_unit)
        self.y_last = NXExpression(number_type, 'y_last', disk_teeth.y_last, mm_unit)
        self.z_last = NXExpression(number_type, 'z_last', disk_teeth.z_last, mm_unit)
        self.y_last_next = NXExpression(number_type, 'y_last_next', disk_teeth.y_last_next, mm_unit)
        self.z_last_next = NXExpression(number_type, 'z_last_next', disk_teeth.z_last_next, mm_unit)
        self.z_down = NXExpression(number_type, 'z_down', disk_teeth.z_down, mm_unit)
        self.D_out = NXExpression(number_type, 'D_out', 2 * disk_teeth.z0, mm_unit)
        self.D_tail_in = NXExpression(number_type, 'D_tail_in', 2 * disk_teeth.z_down, mm_unit)
        self.angle_disk1 = NXExpression(number_type, 'angle_disk1', 88, deg_unit)
        self.angle_disk2 = NXExpression(number_type, 'angle_disk2', 88, deg_unit)
        self.b1 = NXExpression(number_type, 'b1', 0.55 * self.b_a_tail.value, mm_unit)
        self.h_head = NXExpression(number_type, 'h_head', 0.5 * (self.D_out.value - self.D_tail_in.value) + 5, mm_unit)
        self.h3 = NXExpression(number_type, 'h3', 4, mm_unit)
        self.Dr_bolt = NXExpression(number_type, 'Dr_bolt', 125, mm_unit)
        self.d_bolt = NXExpression(number_type, 'd_bolt', 8.5, mm_unit)
        self.t1 = NXExpression(number_type, 't1', 2.1, mm_unit)
        self.l_l = NXExpression(number_type, 'l_l', 20, mm_unit)
        self.l_r = NXExpression(number_type, 'l_r', 15, mm_unit)
        self.b3 = NXExpression(number_type, 'b3', 1.15 * self.b_a_tail.value, mm_unit)
        self.b5 = NXExpression(number_type, 'b5', 0.9 * self.b3.value, mm_unit)
        self.z_hole1 = NXExpression(integer_type, 'z_hole1', 12, nd_unit)
        self.z_bolt = NXExpression(integer_type, 'z_bolt', 12, nd_unit)
        self.r2 = NXExpression(number_type, 'r2', 0.8, mm_unit)
        self.r3 = NXExpression(number_type, 'r3', 0.8, mm_unit)
        self.z_tooth_hirth = NXExpression(integer_type, 'z_tooth_hirth', 96, nd_unit)
        self.gamma_hirth = NXExpression(number_type, 'gamma_hirth', 360 / (self.z_tooth_hirth.value * 4), deg_unit)
        self.r_hirth = NXExpression(number_type, 'r_hirth', self.Dr_bolt.value / 2 + self.d_bolt.value / 2 +
                                    self.t1.value + 50, mm_unit)
        self.a_hirth = NXExpression(number_type, 'a_hirth', self.r_hirth.value *
                                    np.sin(np.radians(self.gamma_hirth.value)), mm_unit)
        self.l_hirth = NXExpression(number_type, 'l_hirth', self.r_hirth.value *
                                    np.cos(np.radians(self.gamma_hirth.value)), mm_unit)
        self.alpha_hirth = NXExpression(number_type, 'alpha_hirth', 60, deg_unit)
        self.t_hirth = NXExpression(number_type, 't_hirth', self.a_hirth.value /
                                    np.tan(0.5 * np.radians(self.alpha_hirth.value)), mm_unit)
        self.beta_hirth = NXExpression(number_type, 'beta_hirth',
                                       np.degrees(np.arcsin(self.t_hirth.value / self.l_hirth.value)), deg_unit)
        self.alpha_butt_hirth = NXExpression(number_type, 'alpha_butt_hirth',
                                             2 * np.degrees(np.arctan(np.tan(np.radians(self.alpha_hirth.value) / 2) *
                                                            np.cos(np.radians(self.beta_hirth.value)))), deg_unit)
        self.l1_tail_ledge = first_stage_tail.l1_tail_ledge
        self.l2_tail_ledge = first_stage_tail.l2_tail_ledge
        self.d1_tail_ledge = first_stage_tail.d1_tail_ledge
        self.d2_tail_ledge = first_stage_tail.d2_tail_ledge
        self.angle_tail_ledge = first_stage_tail.angle_tail_ledge
        self.t1_tail_ledge = first_stage_tail.t1_tail_ledge
        self.d_hole3 = NXExpression(number_type, 'd_hole3', round(self.t1_tail_ledge.value / 2), mm_unit)
        self.T_rk_blade_in = NXExpression(number_type, 'T_rk_blade_in', stages[0]['rk']['T_rk_blade_in'] - 273,
                                          'Celsius')
        self.angle_velocity = NXExpression(number_type, 'angle_velocity', turbine.n * np.pi / 30, 'DegreesPerSecond')


class SecondStageDisk:
    def __init__(self):
        blade_teeth = BladeLockTeethCoordinates(second_stage_tail.y0.value, second_stage_tail.z0.value,
                                                second_stage_tail.s.value, second_stage_tail.r1.value,
                                                second_stage_tail.r2.value, second_stage_tail.h0.value,
                                                np.radians(second_stage_tail.phi.value),
                                                np.radians(second_stage_tail.gamma.value),
                                                np.radians(second_stage_tail.beta.value),
                                                second_stage_tail.teeth_count.value)
        self.z_blade = NXExpression(integer_type, 'z_blade', rk_blades[1]['z'].value, nd_unit)
        self.alpha = second_stage_tail.alpha
        self.s = second_stage_tail.s
        self.r1 = second_stage_tail.r1
        self.phi = second_stage_tail.phi
        self.gamma = second_stage_tail.gamma
        self.beta = second_stage_tail.beta
        self.teeth_count = second_stage_tail.teeth_count
        k = 0.05
        self.h1 = NXExpression(number_type, 'h1', 2, mm_unit)
        self.b_a_tail = second_stage_tail.b_a_tail
        self.disks_tail_interval = NXExpression(number_type, 'disks_tail_interval', stages[0]['rk']['delta_a_rk']*1e3 +
                                                stages[1]['sa']['b_a']*1e3 + stages[1]['rk']['delta_a_sa']*1e3 +
                                                (first_stage_tail.b_a1.value - first_stage_tail.b3.value -
                                                 first_stage_tail.b_a_tail.value - first_stage_tail.b2.value) +
                                                (second_stage_tail.b3.value - second_stage_tail.b1.value), mm_unit)
        disk_teeth = DiskLockTeethCoordinates(blade_teeth, k, self.h1.value)
        self.angle1 = NXExpression(number_type, 'angle1', np.degrees(blade_teeth.angle1), deg_unit)
        self.ym2 = NXExpression(number_type, 'ym2', disk_teeth.ym2, mm_unit)
        self.zm2 = NXExpression(number_type, 'zm2', disk_teeth.zm2, mm_unit)
        self.ym1 = NXExpression(number_type, 'ym1', disk_teeth.ym1, mm_unit)
        self.zm1 = NXExpression(number_type, 'zm1', disk_teeth.zm1, mm_unit)
        self.y0 = NXExpression(number_type, 'y0', disk_teeth.y0, mm_unit)
        self.z0 = NXExpression(number_type, 'z0', disk_teeth.z0, mm_unit)
        self.y1 = NXExpression(number_type, 'y1', disk_teeth.y1, mm_unit)
        self.z1 = NXExpression(number_type, 'z1', disk_teeth.z1, mm_unit)
        self.y2 = NXExpression(number_type, 'y2', disk_teeth.y2, mm_unit)
        self.z2 = NXExpression(number_type, 'z2', disk_teeth.z2, mm_unit)
        self.y3 = NXExpression(number_type, 'y3', disk_teeth.y3, mm_unit)
        self.z3 = NXExpression(number_type, 'z3', disk_teeth.z3, mm_unit)
        self.y4 = NXExpression(number_type, 'y4', disk_teeth.y4, mm_unit)
        self.z4 = NXExpression(number_type, 'z4', disk_teeth.z4, mm_unit)
        self.y5 = NXExpression(number_type, 'y5', disk_teeth.y5, mm_unit)
        self.z5 = NXExpression(number_type, 'z5', disk_teeth.z5, mm_unit)
        self.y6 = NXExpression(number_type, 'y6', disk_teeth.y6, mm_unit)
        self.z6 = NXExpression(number_type, 'z6', disk_teeth.z6, mm_unit)
        self.y7 = NXExpression(number_type, 'y7', disk_teeth.y7, mm_unit)
        self.z7 = NXExpression(number_type, 'z7', disk_teeth.z7, mm_unit)
        self.y_last = NXExpression(number_type, 'y_last', disk_teeth.y_last, mm_unit)
        self.z_last = NXExpression(number_type, 'z_last', disk_teeth.z_last, mm_unit)
        self.y_last_next = NXExpression(number_type, 'y_last_next', disk_teeth.y_last_next, mm_unit)
        self.z_last_next = NXExpression(number_type, 'z_last_next', disk_teeth.z_last_next, mm_unit)
        self.z_down = NXExpression(number_type, 'z_down', disk_teeth.z_down, mm_unit)
        self.D_out = NXExpression(number_type, 'D_out', 2 * disk_teeth.z0, mm_unit)
        self.D_tail_in = NXExpression(number_type, 'D_tail_in', 2 * disk_teeth.z_down, mm_unit)
        self.angle_disk1 = NXExpression(number_type, 'angle_disk1', 88, deg_unit)
        self.angle_disk2 = NXExpression(number_type, 'angle_disk2', 88, deg_unit)
        self.b1 = NXExpression(number_type, 'b1', 0.55 * self.b_a_tail.value, mm_unit)
        self.h_head = NXExpression(number_type, 'h_head', 0.5 * (self.D_out.value - self.D_tail_in.value) + 5, mm_unit)
        self.h3 = NXExpression(number_type, 'h3', 4, mm_unit)
        self.l_l = NXExpression(number_type, 'l_l', self.disks_tail_interval.value, mm_unit)
        self.l_r = NXExpression(number_type, 'l_r', 6, mm_unit)
        self.b4 = NXExpression(number_type, 'b4', 0.95 * self.b_a_tail.value, mm_unit)
        self.z_bolt = FirstStageDisk().z_bolt
        self.r2 = NXExpression(number_type, 'r2', 0.8, mm_unit)
        self.r3 = NXExpression(number_type, 'r3', 0.8, mm_unit)
        self.z_tooth_hirth = FirstStageDisk().z_tooth_hirth
        self.gamma_hirth = FirstStageDisk().gamma_hirth
        self.r_hirth = FirstStageDisk().r_hirth
        self.alpha_hirth = FirstStageDisk().alpha_hirth
        self.a_hirth = FirstStageDisk().a_hirth
        self.l_hirth = FirstStageDisk().l_hirth
        self.alpha_butt_hirth = FirstStageDisk().alpha_butt_hirth
        self.t_hirth = FirstStageDisk().t_hirth
        self.beta_hirth = FirstStageDisk().beta_hirth
        self.l1_tail_ledge = second_stage_tail.l1_tail_ledge
        self.l2_tail_ledge = second_stage_tail.l2_tail_ledge
        self.d1_tail_ledge = second_stage_tail.d1_tail_ledge
        self.d2_tail1_ledge = second_stage_tail.d2_tail_ledge
        self.angle_tail_ledge = second_stage_tail.angle_tail_ledge
        self.t1_tail_ledge = second_stage_tail.t1_tail_ledge
        self.d_hole3 = NXExpression(number_type, 'd_hole3', round(self.t1_tail_ledge.value / 2), mm_unit)
        self.angle_velocity = FirstStageDisk().angle_velocity
        self.T_rk_blade_in = NXExpression(number_type, 'T_rk_blade_in', stages[1]['rk']['T_rk_blade_in'] - 273,
                                          'Celsius')


class FirstStageDiskComputationModel:
    def __init__(self):
        blade_teeth = BladeLockTeethCoordinates(first_stage_tail.y0.value, first_stage_tail.z0.value,
                                                first_stage_tail.s.value, first_stage_tail.r1.value,
                                                first_stage_tail.r2.value, first_stage_tail.h0.value,
                                                np.radians(first_stage_tail.phi.value),
                                                np.radians(first_stage_tail.gamma.value),
                                                np.radians(first_stage_tail.beta.value),
                                                first_stage_tail.teeth_count.value)
        disk1 = FirstStageDisk()
        self.d1_z_blade = NXExpression(integer_type, 'd1_z_blade', rk_blades[0]['z'].value, nd_unit)
        k = 0.05
        self.d1_h1 = NXExpression(number_type, 'h1', disk1.h1.value, mm_unit)
        self.d1_b_a_tail = NXExpression(number_type, 'd1_b_a_tail', first_stage_tail.b_a_tail.value, mm_unit)
        disk_teeth = DiskLockTeethCoordinates(blade_teeth, k, disk1.h1.value)
        self.d1_angle1 = NXExpression(number_type, 'd1_angle1', np.degrees(blade_teeth.angle1), deg_unit)
        self.d1_ym2 = NXExpression(number_type, 'd1_ym2', disk_teeth.ym2, mm_unit)
        self.d1_zm2 = NXExpression(number_type, 'd1_zm2', disk_teeth.zm2, mm_unit)
        self.d1_ym1 = NXExpression(number_type, 'd1_ym1', disk_teeth.ym1, mm_unit)
        self.d1_zm1 = NXExpression(number_type, 'd1_zm1', disk_teeth.zm1, mm_unit)
        self.d1_y0 = NXExpression(number_type, 'd1_y0', disk_teeth.y0, mm_unit)
        self.d1_z0 = NXExpression(number_type, 'd1_z0', disk_teeth.z0, mm_unit)
        self.d1_y1 = NXExpression(number_type, 'd1_y1', disk_teeth.y1, mm_unit)
        self.d1_z1 = NXExpression(number_type, 'd1_z1', disk_teeth.z1, mm_unit)
        self.d1_y2 = NXExpression(number_type, 'd1_y2', disk_teeth.y2, mm_unit)
        self.d1_z2 = NXExpression(number_type, 'd1_z2', disk_teeth.z2, mm_unit)
        self.d1_y3 = NXExpression(number_type, 'd1_y3', disk_teeth.y3, mm_unit)
        self.d1_z3 = NXExpression(number_type, 'd1_z3', disk_teeth.z3, mm_unit)
        self.d1_y4 = NXExpression(number_type, 'd1_y4', disk_teeth.y4, mm_unit)
        self.d1_z4 = NXExpression(number_type, 'd1_z4', disk_teeth.z4, mm_unit)
        self.d1_y5 = NXExpression(number_type, 'd1_y5', disk_teeth.y5, mm_unit)
        self.d1_z5 = NXExpression(number_type, 'd1_z5', disk_teeth.z5, mm_unit)
        self.d1_y6 = NXExpression(number_type, 'd1_y6', disk_teeth.y6, mm_unit)
        self.d1_z6 = NXExpression(number_type, 'd1_z6', disk_teeth.z6, mm_unit)
        self.d1_y7 = NXExpression(number_type, 'd1_y7', disk_teeth.y7, mm_unit)
        self.d1_z7 = NXExpression(number_type, 'd1_z7', disk_teeth.z7, mm_unit)
        self.d1_y_last = NXExpression(number_type, 'd1_y_last', disk_teeth.y_last, mm_unit)
        self.d1_z_last = NXExpression(number_type, 'd1_z_last', disk_teeth.z_last, mm_unit)
        self.d1_y_last_next = NXExpression(number_type, 'd1_y_last_next', disk_teeth.y_last_next, mm_unit)
        self.d1_z_last_next = NXExpression(number_type, 'd1_z_last_next', disk_teeth.z_last_next, mm_unit)
        self.d1_z_down = NXExpression(number_type, 'd1_z_down', disk_teeth.z_down, mm_unit)
        self.d1_D_out = NXExpression(number_type, 'd1_D_out', 2 * disk_teeth.z0, mm_unit)
        self.d1_D_tail_in = NXExpression(number_type, 'd1_D_tail_in', 2 * disk_teeth.z_down, mm_unit)
        self.d1_angle_disk1 = NXExpression(number_type, 'd1_angle_disk1', disk1.angle_disk1.value, deg_unit)
        self.d1_angle_disk2 = NXExpression(number_type, 'd1_angle_disk2', disk1.angle_disk2.value, deg_unit)
        self.d1_b1 = NXExpression(number_type, 'd1_b1', disk1.b1.value, mm_unit)
        self.d1_h_cut = NXExpression(number_type, 'd1_h_cut', 0.5 * (self.d1_D_out.value - self.d1_D_tail_in.value) + 1, mm_unit)
        self.d1_h4 = NXExpression(number_type, 'd1_h4', 4, mm_unit)
        self.d1_h3 = NXExpression(number_type, 'd1_h3', disk1.h3.value, mm_unit)
        self.Dr_bolt = NXExpression(number_type, 'Dr_bolt', disk1.Dr_bolt.value, mm_unit)
        self.d_bolt = NXExpression(number_type, 'd_bolt', disk1.d_bolt.value, mm_unit)
        self.t1 = NXExpression(number_type, 'd1_t1', disk1.t1.value, mm_unit)
        self.d1_l_l = NXExpression(number_type, 'd1_l_l', disk1.l_l.value, mm_unit)
        self.d1_l_r = NXExpression(number_type, 'd1_l_r', disk1.l_r.value, mm_unit)
        self.d1_b3 = NXExpression(number_type, 'd1_b3', disk1.b3.value, mm_unit)
        self.d1_b5 = NXExpression(number_type, 'd1_b5', disk1.b5.value, mm_unit)
        self.d1_z_hole1 = NXExpression(integer_type, 'd1_z_hole1', disk1.z_hole1.value, nd_unit)
        self.z_bolt = NXExpression(integer_type, 'z_bolt', disk1.z_bolt.value, nd_unit)
        self.d1_r2 = NXExpression(number_type, 'd1_r2', disk1.r2.value, mm_unit)
        self.d1_r3 = NXExpression(number_type, 'd1_r3', disk1.r3.value, mm_unit)
        self.z_tooth_hirth = NXExpression(integer_type, 'z_tooth_hirth', disk1.z_tooth_hirth.value, nd_unit)
        self.gamma_hirth = NXExpression(number_type, 'gamma_hirth', 360 / (self.z_tooth_hirth.value * 4), deg_unit)
        self.r_hirth = NXExpression(number_type, 'r_hirth', self.Dr_bolt.value / 2 + self.d_bolt.value / 2 +
                                    self.t1.value + 50, mm_unit)
        self.a_hirth = NXExpression(number_type, 'a_hirth', self.r_hirth.value *
                                    np.sin(np.radians(self.gamma_hirth.value)), mm_unit)
        self.l_hirth = NXExpression(number_type, 'l_hirth', self.r_hirth.value *
                                    np.cos(np.radians(self.gamma_hirth.value)), mm_unit)
        self.alpha_hirth = NXExpression(number_type, 'alpha_hirth', disk1.alpha_hirth.value, deg_unit)
        self.t_hirth = NXExpression(number_type, 't_hirth', self.a_hirth.value /
                                    np.tan(0.5 * np.radians(self.alpha_hirth.value)), mm_unit)
        self.beta_hirth = NXExpression(number_type, 'beta_hirth',
                                       np.degrees(np.arcsin(self.t_hirth.value / self.l_hirth.value)), deg_unit)
        self.alpha_butt_hirth = NXExpression(number_type, 'alpha_butt_hirth',
                                             2 * np.degrees(np.arctan(np.tan(np.radians(self.alpha_hirth.value) / 2) *
                                                            np.cos(np.radians(self.beta_hirth.value)))), deg_unit)
        self.d1_T_rk_blade_in = NXExpression(number_type, 'd1_T_rk_blade_in', stages[0]['rk']['T_rk_blade_in'] - 273,
                                          'Celsius')
        self.angle_velocity = NXExpression(number_type, 'angle_velocity', turbine.n * np.pi / 30, 'DegreesPerSecond')
        self.d1_pres = NXExpression(number_type, 'd1_press', (1 / 12 * self.angle_velocity.value ** 2 * 8000 /
                                    (self.d1_D_out.value/1e3 - 2 * self.d1_h_cut.value/1e3) *
                                    (self.d1_D_out.value**3 - (self.d1_D_out.value - 2 * self.d1_h_cut.value)**3)/1e9)/1e6,
                                    'PressureNewtonPerSquareMilliMeter')


class SecondStageDiskComputationModel:
    def __init__(self):
        blade_teeth = BladeLockTeethCoordinates(second_stage_tail.y0.value, second_stage_tail.z0.value,
                                                second_stage_tail.s.value, second_stage_tail.r1.value,
                                                second_stage_tail.r2.value, second_stage_tail.h0.value,
                                                np.radians(second_stage_tail.phi.value),
                                                np.radians(second_stage_tail.gamma.value),
                                                np.radians(second_stage_tail.beta.value),
                                                second_stage_tail.teeth_count.value)
        disk1 = FirstStageDisk()
        disk2 = SecondStageDisk()
        self.d2_z_blade = NXExpression(integer_type, 'd2_z_blade', rk_blades[1]['z'].value, nd_unit)
        k = 0.05
        self.d2_h1 = NXExpression(number_type, 'd2_h1', disk2.h1.value, mm_unit)
        self.d2_b_a_tail = NXExpression(number_type, 'd2_b_a_tail', second_stage_tail.b_a_tail.value, mm_unit)
        self.disks_tail_interval = NXExpression(number_type, 'disks_tail_interval',
                                                stages[0]['rk']['delta_a_rk'] * 1e3 +
                                                stages[1]['sa']['b_a'] * 1e3 + stages[1]['rk']['delta_a_sa'] * 1e3 +
                                                (first_stage_tail.b_a1.value - first_stage_tail.b3.value -
                                                 first_stage_tail.b_a_tail.value - first_stage_tail.b2.value) +
                                                (second_stage_tail.b3.value - second_stage_tail.b1.value), mm_unit)
        disk_teeth = DiskLockTeethCoordinates(blade_teeth, k, self.d2_h1.value)
        self.d2_angle1 = NXExpression(number_type, 'd2_angle1', np.degrees(blade_teeth.angle1), deg_unit)
        self.d2_ym2 = NXExpression(number_type, 'd2_ym2', disk_teeth.ym2, mm_unit)
        self.d2_zm2 = NXExpression(number_type, 'd2_zm2', disk_teeth.zm2, mm_unit)
        self.d2_ym1 = NXExpression(number_type, 'd2_ym1', disk_teeth.ym1, mm_unit)
        self.d2_zm1 = NXExpression(number_type, 'd2_zm1', disk_teeth.zm1, mm_unit)
        self.d2_y0 = NXExpression(number_type, 'd2_y0', disk_teeth.y0, mm_unit)
        self.d2_z0 = NXExpression(number_type, 'd2_z0', disk_teeth.z0, mm_unit)
        self.d2_y1 = NXExpression(number_type, 'd2_y1', disk_teeth.y1, mm_unit)
        self.d2_z1 = NXExpression(number_type, 'd2_z1', disk_teeth.z1, mm_unit)
        self.d2_y2 = NXExpression(number_type, 'd2_y2', disk_teeth.y2, mm_unit)
        self.d2_z2 = NXExpression(number_type, 'd2_z2', disk_teeth.z2, mm_unit)
        self.d2_y3 = NXExpression(number_type, 'd2_y3', disk_teeth.y3, mm_unit)
        self.d2_z3 = NXExpression(number_type, 'd2_z3', disk_teeth.z3, mm_unit)
        self.d2_y4 = NXExpression(number_type, 'd2_y4', disk_teeth.y4, mm_unit)
        self.d2_z4 = NXExpression(number_type, 'd2_z4', disk_teeth.z4, mm_unit)
        self.d2_y5 = NXExpression(number_type, 'd2_y5', disk_teeth.y5, mm_unit)
        self.d2_z5 = NXExpression(number_type, 'd2_z5', disk_teeth.z5, mm_unit)
        self.d2_y6 = NXExpression(number_type, 'd2_y6', disk_teeth.y6, mm_unit)
        self.d2_z6 = NXExpression(number_type, 'd2_z6', disk_teeth.z6, mm_unit)
        self.d2_y7 = NXExpression(number_type, 'd2_y7', disk_teeth.y7, mm_unit)
        self.d2_z7 = NXExpression(number_type, 'd2_z7', disk_teeth.z7, mm_unit)
        self.d2_y_last = NXExpression(number_type, 'd2_y_last', disk_teeth.y_last, mm_unit)
        self.d2_z_last = NXExpression(number_type, 'd2_z_last', disk_teeth.z_last, mm_unit)
        self.d2_y_last_next = NXExpression(number_type, 'd2_y_last_next', disk_teeth.y_last_next, mm_unit)
        self.d2_z_last_next = NXExpression(number_type, 'd2_z_last_next', disk_teeth.z_last_next, mm_unit)
        self.d2_z_down = NXExpression(number_type, 'd2_z_down', disk_teeth.z_down, mm_unit)
        self.d2_D_out = NXExpression(number_type, 'd2_D_out', 2 * disk_teeth.z0, mm_unit)
        self.d2_D_tail_in = NXExpression(number_type, 'd2_D_tail_in', 2 * disk_teeth.z_down, mm_unit)
        self.d2_angle_disk1 = NXExpression(number_type, 'd2_angle_disk1', disk2.angle_disk1.value, deg_unit)
        self.d2_angle_disk2 = NXExpression(number_type, 'd2_angle_disk2', disk2.angle_disk2.value, deg_unit)
        self.d2_b1 = NXExpression(number_type, 'd2_b1', disk2.b1.value, mm_unit)
        self.d2_h_cut = NXExpression(number_type, 'd2_h_cut', 0.5 * (self.d2_D_out.value - self.d2_D_tail_in.value) + 1, mm_unit)
        self.d2_h4 = NXExpression(number_type, 'd2_h4', 4, mm_unit)
        self.d2_h3 = NXExpression(number_type, 'd2_h3', disk2.h3.value, mm_unit)
        self.d2_l_l = NXExpression(number_type, 'd2_l_l', self.disks_tail_interval.value, mm_unit)
        self.d2_l_r = NXExpression(number_type, 'd2_l_r', disk2.l_r.value, mm_unit)
        self.d2_b4 = NXExpression(number_type, 'd2_b4', disk2.b4.value, mm_unit)
        self.d2_r2 = NXExpression(number_type, 'd2_r2', disk2.r2.value, mm_unit)
        self.d2_r3 = NXExpression(number_type, 'd2_r3', disk2.r3.value, mm_unit)
        self.d2_T_rk_blade_in = NXExpression(number_type, 'd2_T_rk_blade_in', stages[1]['rk']['T_rk_blade_in'] - 273,
                                          'Celsius')
        self.d2_pres = NXExpression(number_type, 'd2_press', (1 / 12 * disk2.angle_velocity.value ** 2 * 8000 /
                                    (self.d2_D_out.value/1e3 - 2 * self.d2_h_cut.value/1e3) *
                                    (self.d2_D_out.value ** 3 - (self.d2_D_out.value - 2 * self.d2_h_cut.value) ** 3)/1e9)/1e6,
                                    'PressureNewtonPerSquareMilliMeter')
