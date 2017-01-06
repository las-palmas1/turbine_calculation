from report.report_data_update_lib import *
import os

engine = get_TV7_117_object(os.path.join(os.path.dirname(os.path.dirname(__file__)), 'cycle_calculation',
                                         'cycle_calculation_results'))

turbine = get_turbine_object(os.path.join(os.path.dirname(os.path.dirname(__file__)), 'average_streamline_calculation',
                                          'average_streamline_calculation_results'))

par_dict = get_parameter_dict(engine, turbine)

num_ins = NumInserter(r'report_template\first_doc.tex', r'report\report.tex', par_dict)

num_ins.process_document()