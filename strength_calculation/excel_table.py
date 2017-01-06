from strength_calculation.post_processing import two_calculations
import xlwt

book = xlwt.Workbook()

# создание листа с результатами первого расчета
sheet1 = book.add_sheet('Sheet 1')
sheet1.write(0, 0, 'i')
sheet1.write(0, 1, 'Радиус, мм')
sheet1.write(0, 2, 'sigma_r, МПа')
sheet1.write(0, 3, 'sigma_t, МПа')

# создание листа с результатами второго расчета
sheet2 = book.add_sheet('Sheet 2')
sheet2.write(0, 0, 'i')
sheet2.write(0, 1, 'Радиус, мм')
sheet2.write(0, 2, 'sigma_r, МПа')
sheet2.write(0, 3, 'sigma_t, МПа')

# создание листа с пересчтитанными результами
sheet3 = book.add_sheet('Sheet 3')
sheet3.write(0, 0, 'i')
sheet3.write(0, 1, 'Радиус, мм')
sheet3.write(0, 2, 'sigma_r, МПа')
sheet3.write(0, 3, 'sigma_t, МПа')

# создание листа с осредненными результами
sheet4 = book.add_sheet('Sheet 4')
sheet4.write(0, 0, 'N')
sheet4.write(0, 1, 'Радиус, мм')
sheet4.write(0, 2, 'sigma_r, МПа')
sheet4.write(0, 3, 'sigma_t, МПа')
sheet4.write(0, 4, 'sigma_eq, МПа')


for i in range(two_calculations.part_number):
    sheet1.write(2 * i + 1, 0, i + 1)
    sheet1.write(2 * i + 2, 0, i + 1)
    sheet1.write(2 * i + 1, 1, round(two_calculations.r_arr[i] * 1e3, 1))
    sheet1.write(2 * i + 2, 1, round(two_calculations.r_arr[i + 1] * 1e3, 1))
    sheet1.write(2 * i + 1, 2, round(two_calculations.first_calc.sigma_r_arr[i][0] / 1e6, 1))
    sheet1.write(2 * i + 2, 2, round(two_calculations.first_calc.sigma_r_arr[i][1] / 1e6, 1))
    sheet1.write(2 * i + 1, 3, round(two_calculations.first_calc.sigma_t_arr[i][0] / 1e6, 1))
    sheet1.write(2 * i + 2, 3, round(two_calculations.first_calc.sigma_t_arr[i][1] / 1e6, 1))

    sheet2.write(2 * i + 1, 0, i + 1)
    sheet2.write(2 * i + 2, 0, i + 1)
    sheet2.write(2 * i + 1, 1, round(two_calculations.r_arr[i] * 1e3, 1))
    sheet2.write(2 * i + 2, 1, round(two_calculations.r_arr[i + 1] * 1e3, 1))
    sheet2.write(2 * i + 1, 2, round(two_calculations.second_calc.sigma_r_arr[i][0] / 1e6, 1))
    sheet2.write(2 * i + 2, 2, round(two_calculations.second_calc.sigma_r_arr[i][1] / 1e6, 1))
    sheet2.write(2 * i + 1, 3, round(two_calculations.second_calc.sigma_t_arr[i][0] / 1e6, 1))
    sheet2.write(2 * i + 2, 3, round(two_calculations.second_calc.sigma_t_arr[i][1] / 1e6, 1))

    sheet3.write(2 * i + 1, 0, i + 1)
    sheet3.write(2 * i + 2, 0, i + 1)
    sheet3.write(2 * i + 1, 1, round(two_calculations.r_arr[i] * 1e3, 1))
    sheet3.write(2 * i + 2, 1, round(two_calculations.r_arr[i + 1] * 1e3, 1))
    sheet3.write(2 * i + 1, 2, round(two_calculations.sigma_r_arr_real[i][0] / 1e6, 1))
    sheet3.write(2 * i + 2, 2, round(two_calculations.sigma_r_arr_real[i][1] / 1e6, 1))
    sheet3.write(2 * i + 1, 3, round(two_calculations.sigma_t_arr_real[i][0] / 1e6, 1))
    sheet3.write(2 * i + 2, 3, round(two_calculations.sigma_t_arr_real[i][1] / 1e6, 1))


for i in range(two_calculations.part_number + 1):
    sheet4.write(i + 1, 0, i + 1)
    sheet4.write(i + 1, 1, round(two_calculations.r_arr[i] * 1e3, 1))
    sheet4.write(i + 1, 2, round(two_calculations.sigma_r_arr_av[i] / 1e6, 1))
    sheet4.write(i + 1, 3, round(two_calculations.sigma_t_arr_av[i] / 1e6, 1))
    sheet4.write(i + 1, 4, round(two_calculations.sigma_eq_arr[i] / 1e6, 1))


book.save(r'output\CalculationResults.xls')