from parts_parameters.blades import rk_blades, sa_blades, stage_tails
import parts_parameters.func as func
import os
from parts_parameters.disks import FirstStageDisk, SecondStageDisk, FirstStageDiskComputationModel
from parts_parameters.standart_parts.SingleRowRadialBallBearingGOST8338_75 import BallBearing
from parts_parameters.standart_parts.SlottedRoundNutsGOST11871_88 import SlottedRoundNutGOST
from parts_parameters.standart_parts.MultilegLockWashersGOST11872_89 import MultilegLockWasherGOST
from parts_parameters.standart_parts.HexagonBoltsGOST7798_70 import HexagonBoltGOST
from parts_parameters.standart_parts.HexagonNutsGOST5915_70 import HexagonNutGOST
from parts_parameters.standart_parts.LockWasherGOST6402_70 import LockWasherGOST
from parts_parameters.standart_parts.TabLockWashersGOST13463_77 import TabLockWasherGOST
from parts_parameters.standart_parts.RadialShortCylindricalRollerBearingsGOST8328_75 import RollerBearing

models_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), 'model')

models_name_arr = []

for n, i in enumerate(rk_blades):
    func.create_expressions_file('st%s_rk_blade' % (n + 1), os.path.join(models_dir, 'st%s_rk_blade.prt' % (n + 1)),
                                 func.dict_combination(i, stage_tails[n]), models_name_arr)

for n, i in enumerate(sa_blades):
    func.create_expressions_file('st%s_sa_blade' % (n + 1), os.path.join(models_dir, 'st%s_sa_blade.prt' % (n + 1)), i,
                                 models_name_arr)

func.create_expressions_file('st1_disk', os.path.join(models_dir, 'st1_disk.prt'), FirstStageDisk().__dict__,
                             models_name_arr)

func.create_expressions_file('st2_disk', os.path.join(models_dir, 'st2_disk.prt'), SecondStageDisk().__dict__,
                             models_name_arr)

func.create_expressions_file('st1_disk_model_for_computation',
                             os.path.join(models_dir, 'st1_disk_model_for_computation.prt'),
                             FirstStageDiskComputationModel().__dict__, models_name_arr)

func.create_expressions_file('ball_bearing', os.path.join(models_dir, 'ball_bearing.prt'),
                             BallBearing(d=90, series='light').__dict__, models_name_arr)

slotted_round_nut1 = SlottedRoundNutGOST(d=85, chamfer=2.5, series='narrow')
func.create_expressions_file('slotted_round_nut1', os.path.join(models_dir, 'slotted_round_nut1.prt'),
                             slotted_round_nut1.__dict__, models_name_arr)

multileg_lock_washer1 = MultilegLockWasherGOST(d=slotted_round_nut1.d_in.value, nut_h=slotted_round_nut1.h.value,
                                               nut_D=slotted_round_nut1.D.value)
func.create_expressions_file('multileg_lock_washer1', os.path.join(models_dir, 'multileg_lock_washer1.prt'),
                             multileg_lock_washer1.__dict__, models_name_arr)

bolt_external_stator_body = HexagonBoltGOST(6, 12, 12, 1)
bolt_oil_part4 = HexagonBoltGOST(6, 12, 12, 1)
bolt_oil_part6 = HexagonBoltGOST(6, 12, 12, 1)
bolt_oil_tube = HexagonBoltGOST(6, 13, 13, 1)
bolt_tripot = HexagonBoltGOST(6, 15, 12, 1)
spring_lock_washer_d6 = LockWasherGOST(6)
lock_washer_external_stator_body = TabLockWasherGOST(6, 9, 3)
lock_washer_oil_tube = TabLockWasherGOST(6, 9, 3)
nut_disk_bolt = HexagonNutGOST(8, 6.5)

func.create_expressions_file('bolt_external_stator_body', os.path.join(models_dir, 'bolt_external_stator_body.prt'),
                             bolt_external_stator_body.__dict__, models_name_arr)
func.create_expressions_file('bolt_oil_part4', os.path.join(models_dir, 'bolt_oil_part4.prt'),
                             bolt_oil_part4.__dict__, models_name_arr)
func.create_expressions_file('bolt_oil_part6', os.path.join(models_dir, 'bolt_oil_part6.prt'),
                             bolt_oil_part6.__dict__, models_name_arr)
func.create_expressions_file('bolt_oil_tube', os.path.join(models_dir, 'bolt_oil_tube.prt'),
                             bolt_oil_tube.__dict__, models_name_arr)
func.create_expressions_file('bolt_tripot', os.path.join(models_dir, 'bolt_tripot.prt'),
                             bolt_tripot.__dict__, models_name_arr)
func.create_expressions_file('spring_lock_washer_d6', os.path.join(models_dir, 'spring_lock_washer_d6.prt'),
                             spring_lock_washer_d6.__dict__, models_name_arr)
func.create_expressions_file('lock_washer_external_stator_body',
                             os.path.join(models_dir, 'lock_washer_external_stator_body.prt'),
                             lock_washer_external_stator_body.__dict__, models_name_arr)
func.create_expressions_file('lock_washer_oil_tube', os.path.join(models_dir, 'lock_washer_oil_tube.prt'),
                             lock_washer_oil_tube.__dict__, models_name_arr)
func.create_expressions_file('nut_disk_bolt', os.path.join(models_dir, 'nut_disk_bolt.prt'),
                             nut_disk_bolt.__dict__, models_name_arr)

slotted_round_nut2 = SlottedRoundNutGOST(d=45, chamfer=2, series='narrow')
func.create_expressions_file('slotted_round_nut2', os.path.join(models_dir, 'slotted_round_nut2.prt'),
                             slotted_round_nut2.__dict__, models_name_arr)

multileg_lock_washer2 = MultilegLockWasherGOST(d=slotted_round_nut2.d_in.value, nut_h=slotted_round_nut2.h.value,
                                               nut_D=slotted_round_nut2.D.value)
func.create_expressions_file('multileg_lock_washer2', os.path.join(models_dir, 'multileg_lock_washer2.prt'),
                             multileg_lock_washer2.__dict__, models_name_arr)

roller_bearing = RollerBearing(d=50, series='light')
func.create_expressions_file('roller_bearing', os.path.join(models_dir, 'roller_bearing.prt'),
                             roller_bearing.__dict__, models_name_arr)

lock_washer_disk_bolt = TabLockWasherGOST(8, 6.35, 3)
func.create_expressions_file('lock_washer_disk_bolt', os.path.join(models_dir, 'lock_washer_disk_bolt.prt'),
                             lock_washer_disk_bolt.__dict__, models_name_arr)

lock_washer_tripot = TabLockWasherGOST(6, 11.9, 1.5)
func.create_expressions_file('lock_washer_tripot', os.path.join(models_dir, 'lock_washer_tripot.prt'),
                             lock_washer_tripot.__dict__, models_name_arr)