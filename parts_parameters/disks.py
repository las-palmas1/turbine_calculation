from parts_parameters.blades import first_stage_tail, second_stage_tail, rk_blades
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
        self.b1 = NXExpression(number_type, 'b1', 0.5 * self.b_a_tail.value, mm_unit)
        self.h_head = NXExpression(number_type, 'h_head', 0.5 * (self.D_out.value - self.D_tail_in.value) + 5, mm_unit)
        self.h3 = NXExpression(number_type, 'h3', 4, mm_unit)
        self.b2 = NXExpression(number_type, 'b2', 1.4 * self.b1.value, mm_unit)
        self.Dr_bolt = NXExpression(number_type, 'Dr_bolt', 0.5 * self.D_out.value, mm_unit)
        self.d_bolt = NXExpression(number_type, 'd_bolt', 8.5, mm_unit)
        self.t1 = NXExpression(number_type, 't1', 2, mm_unit)
        self.l_l = NXExpression(number_type, 'l_l', 20, mm_unit)
        self.l_r = NXExpression(number_type, 'l_r', 7, mm_unit)
        self.b3 = NXExpression(number_type, 'b3', 1.2 * self.b2.value, mm_unit)
        self.b4 = NXExpression(number_type, 'b4', 0.7 * self.b3.value, mm_unit)
        self.b5 = NXExpression(number_type, 'b5', 0.9 * self.b3.value, mm_unit)
        self.z_hole1 = NXExpression(integer_type, 'z_hole1', 2, nd_unit)
        self.z_bolt = NXExpression(integer_type, 'z_bolt', 6, nd_unit)
