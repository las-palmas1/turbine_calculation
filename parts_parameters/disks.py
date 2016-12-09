from parts_parameters.blades import first_stage_tail, second_stage_tail, rk_blades, sa_blades, stages
from parts_parameters.func import *


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
        self.angle_disk1 = NXExpression(number_type, 'angle_disk1', 87, deg_unit)
        self.angle_disk2 = NXExpression(number_type, 'angle_disk2', 86, deg_unit)
        self.b1 = NXExpression(number_type, 'b1', 0.5 * self.b_a_tail.value, mm_unit)
        self.h_head = NXExpression(number_type, 'h_head', 0.5 * (self.D_out.value - self.D_tail_in.value) + 10, mm_unit)
        self.h3 = NXExpression(number_type, 'h3', 4, mm_unit)
        self.Dr_bolt = NXExpression(number_type, 'Dr_bolt', 0.5 * self.D_out.value, mm_unit)
        self.d_bolt = NXExpression(number_type, 'd_bolt', 8.5, mm_unit)
        self.t1 = NXExpression(number_type, 't1', 2, mm_unit)
        self.l_l = NXExpression(number_type, 'l_l', 20, mm_unit)
        self.l_r = NXExpression(number_type, 'l_r', 15, mm_unit)
        self.b3 = NXExpression(number_type, 'b3', 1.1 * self.b_a_tail.value, mm_unit)
        self.b5 = NXExpression(number_type, 'b5', 0.75 * self.b3.value, mm_unit)
        self.z_hole1 = NXExpression(integer_type, 'z_hole1', 2, nd_unit)
        self.z_bolt = NXExpression(integer_type, 'z_bolt', 6, nd_unit)
        self.r2 = NXExpression(number_type, 'r2', 0.8, mm_unit)
        self.r3 = NXExpression(number_type, 'r3', 0.8, mm_unit)
        self.z_tooth_hirth = NXExpression(integer_type, 'z_tooth_hirth', 250, nd_unit)
        self.gamma_hirth = NXExpression(number_type, 'gamma_hirth', 360 / (self.z_tooth_hirth.value * 4), deg_unit)
        self.r_hirth = NXExpression(number_type, 'r_hirth', self.Dr_bolt.value / 2 + self.d_bolt.value / 2 +
                                    self.t1.value + 50, mm_unit)
        self.a_hirth = NXExpression(number_type, 'a_hirth', self.r_hirth.value *
                                    np.sin(np.radians(self.gamma_hirth.value)), mm_unit)
        self.l_hirth = NXExpression(number_type, 'l_hirth', self.r_hirth.value *
                                    np.cos(np.radians(self.gamma_hirth.value)), mm_unit)
        self.alpha_hirth = NXExpression(number_type, 'alpha_hirth', 45, deg_unit)
        self.t_hirth = NXExpression(number_type, 't_hirth', self.a_hirth.value /
                                    np.tan(0.5 * np.radians(self.alpha_hirth.value)), mm_unit)
        self.beta_hirth = NXExpression(number_type, 'beta_hirth',
                                       np.degrees(np.arcsin(self.t_hirth.value / self.l_hirth.value)), deg_unit)

        self.l1_tail_ledge = first_stage_tail.l1_tail_ledge
        self.l2_tail_ledge = first_stage_tail.l2_tail_ledge
        self.d1_tail_ledge = first_stage_tail.d1_tail_ledge
        self.d2_tail1_ledge = first_stage_tail.d2_tail_ledge
        self.angle_tail_ledge = first_stage_tail.angle_tail_ledge
        self.t1_tail_ledge = first_stage_tail.t1_tail_ledge
        self.d_hole3 = NXExpression(number_type, 'd_hole3', round(self.t1_tail_ledge.value / 2), mm_unit)


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
        self.b1 = NXExpression(number_type, 'b1', 0.5 * self.b_a_tail.value, mm_unit)
        self.h_head = NXExpression(number_type, 'h_head', 0.5 * (self.D_out.value - self.D_tail_in.value) + 10, mm_unit)
        self.h3 = NXExpression(number_type, 'h3', 4, mm_unit)
        self.l_l = NXExpression(number_type, 'l_l', self.disks_tail_interval.value, mm_unit)
        self.l_r = NXExpression(number_type, 'l_r', 6, mm_unit)
        self.b4 = NXExpression(number_type, 'b4', 0.8 * self.b_a_tail.value, mm_unit)
        self.z_bolt = FirstStageDisk().z_bolt
        self.r2 = NXExpression(number_type, 'r2', 0.8, mm_unit)
        self.r3 = NXExpression(number_type, 'r3', 0.8, mm_unit)
        self.z_tooth_hirth = FirstStageDisk().z_tooth_hirth
        self.gamma_hirth = FirstStageDisk().gamma_hirth
        self.r_hirth = FirstStageDisk().r_hirth
        self.a_hirth = FirstStageDisk().a_hirth
        self.l_hirth = FirstStageDisk().l_hirth
        self.alpha_hirth = FirstStageDisk().alpha_hirth
        self.t_hirth = FirstStageDisk().t_hirth
        self.beta_hirth = FirstStageDisk().beta_hirth
        self.l1_tail_ledge = second_stage_tail.l1_tail_ledge
        self.l2_tail_ledge = second_stage_tail.l2_tail_ledge
        self.d1_tail_ledge = second_stage_tail.d1_tail_ledge
        self.d2_tail1_ledge = second_stage_tail.d2_tail_ledge
        self.angle_tail_ledge = second_stage_tail.angle_tail_ledge
        self.t1_tail_ledge = second_stage_tail.t1_tail_ledge
        self.d_hole3 = NXExpression(number_type, 'd_hole3', round(self.t1_tail_ledge.value / 2), mm_unit)