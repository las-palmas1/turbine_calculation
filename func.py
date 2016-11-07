import os

number_type = 'Number'
mm_unit = 'MilliMeter'
rad_unit = 'Radian'
m_unit = 'Meter'
deg_unit = 'Degrees'
integer_type = 'Integer'
nd_unit = 'nd'


class NXExpression:
    def __init__(self, exp_type, exp_name, exp_value, exp_unit):
        self.type = exp_type
        self.name = exp_name
        self.value = exp_value
        self.unit = exp_unit


def dict_combination(*args):
    result = dict()
    for i1 in args:
        for i2 in i1:
            result[i2] = i1[i2]
    return result


def number_into_string(exp_value, exp_type):
    if exp_type != 'Point':
        return str(exp_value)
    else:
        return 'Point(' + str(exp_value[0]) + ',' + str(exp_value[1]) + ',' + str(exp_value[2]) + ')'


def add_expressions_into_file(file, expressions_dict):
    for i in expressions_dict:
        file.write(expressions_dict[i].type + ' ' + expressions_dict[i].name + ' ' +
                   number_into_string(expressions_dict[i].value, expressions_dict[i].type) + ' ' +
                   expressions_dict[i].unit + '\n')


def create_expressions_file(exp_filename, part_filename, exp_dict):
    current_dir = os.path.dirname(__file__)
    abs_exp_filename = os.path.join(current_dir, 'output', 'current', exp_filename)
    abs_path_filename = os.path.join(current_dir, 'output', 'model_paths', exp_filename)
    path_file = open(abs_path_filename, 'w')
    path_file.write(part_filename + '\n')
    path_file.close()
    exp_file = open(abs_exp_filename, 'w')
    add_expressions_into_file(exp_file, exp_dict)
    exp_file.close()
    old_versions_dir = os.path.join(os.path.dirname(__file__), 'output', 'old')
    filename_arr = os.listdir(old_versions_dir)
    if filename_arr.count(exp_filename) == 0:
        file = open(os.path.join(old_versions_dir, exp_filename), 'w')
        file.close()
