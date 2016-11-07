from parts_parameters.rk_blades import rk_blades
import parts_parameters.func as func
import os

models_path = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), 'model')

for n, i in enumerate(rk_blades):
    func.create_expressions_file('st%s_rk_blade' % (n + 1), os.path.join(models_path, 'st%s_rk_blade.prt' % (n + 1)), i)
