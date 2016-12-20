from parts_parameters.blades import rk_blades, sa_blades, stage_tails
import parts_parameters.func as func
import os
from parts_parameters.disks import FirstStageDisk, SecondStageDisk, FirstStageDiskComputationModel
from parts_parameters.standart_parts.SingleRowRadialBallBearingGOST8338_75 import BallBearing
from parts_parameters.standart_parts.SlottedRoundNutsGOST11871_88 import SlottedRoundNutGOST
from parts_parameters.standart_parts.MultilegLockWashersGOST11872_89 import MultilegLockWasherGOST

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

multileg_lock_washer = MultilegLockWasherGOST(d=slotted_round_nut1.d_in.value, nut_h=slotted_round_nut1.h.value,
                                              nut_D=slotted_round_nut1.D.value)
func.create_expressions_file('multileg_lock_washer1', os.path.join(models_dir, 'multileg_lock_washer1.prt'),
                             multileg_lock_washer.__dict__, models_name_arr)

