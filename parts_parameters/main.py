from parts_parameters.blades import rk_blades, sa_blades
import parts_parameters.func as func
import os

models_path = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), 'model')

for n, i in enumerate(rk_blades):
    func.create_expressions_file('st%s_rk_blade' % (n + 1), os.path.join(models_path, 'st%s_rk_blade.prt' % (n + 1)), i)

for n, i in enumerate(sa_blades):
    func.create_expressions_file('st%s_sa_blade' % (n + 1), os.path.join(models_path, 'st%s_sa_blade.prt' % (n + 1)), i)
