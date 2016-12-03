from parts_parameters.blades import rk_blades, sa_blades, stage_tails
import parts_parameters.func as func
import os
from parts_parameters.disks import FirstStageDisk

models_path = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), 'model')

models_name_arr = []

for n, i in enumerate(rk_blades):
    func.create_expressions_file('st%s_rk_blade' % (n + 1), os.path.join(models_path, 'st%s_rk_blade.prt' % (n + 1)),
                                 func.dict_combination(i, stage_tails[n]), models_name_arr)

for n, i in enumerate(sa_blades):
    func.create_expressions_file('st%s_sa_blade' % (n + 1), os.path.join(models_path, 'st%s_sa_blade.prt' % (n + 1)), i,
                                 models_name_arr)

func.create_expressions_file('st1_disk', os.path.join(models_path, 'st1_disk.prt'), FirstStageDisk().__dict__,
                             models_name_arr)
