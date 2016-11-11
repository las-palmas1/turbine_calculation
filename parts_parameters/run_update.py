import subprocess as sp
import os

run_managed_path = r'C:\Program Files\Siemens\NX 10.0\UGII\run_managed.exe'

course_work_dir = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
update_exe_path = os.path.join(course_work_dir, 'Model Parameters Update', 'Model Parameters Update', 'bin', 'Debug',
                               'Model Parameters Update.exe')


pr = sp.Popen([run_managed_path, update_exe_path])
pr.wait()
