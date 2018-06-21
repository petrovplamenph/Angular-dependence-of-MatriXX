"""
The application is intended for investigation of MatriXX angular dependency
and creation of correction factor set to correct for it.
"""
import statistics
import math
import re
import csv

import matplotlib.pyplot as plt

import tkinter as tk
from tkinter import ttk
from tkinter.filedialog import askopenfilenames, asksaveasfilename


LARGE_FONT = ("Verdana", 12)
NORM_FONT = ("Verdana", 10)
SMALL_FONT = ("Verdana", 8)


def search_cord_start(xy_cord):
    """
    Find at the index of the coordinate corresponding
    to the coordinate -11.810 cm in order to exclude
    the doses at positions at which there are no detectors
    """
    for idx in range(len(xy_cord)):
        if abs(xy_cord[idx] + 11.810) < 0.01:
            return idx


def extract_matrix_values(file_name, Flag=False, deviation=False):
    """
    Parse ASCII file to extract the dose values at  X.Y
    """
    matrix = []

    with open(file_name, encoding='ISO-8859-14') as f:
        lines = f.readlines()
        """
        File contains meta information that we don't need for the calculations.
        Values start after `Y[cm]` & end with closing of the `asciibody` tag.
        """
        for i in range(len(lines)):
            if lines[i].startswith('Y[cm]'):
                values_start_index = i + 1

            if lines[i].startswith('</asciibody>'):
                values_end_index = i
            if lines[i].startswith('Data Factor:'):
                data_factor = re.findall('[0-9]*[.]?[0-9]+', lines[i])
                data_factor = float(data_factor[0])
        # Check if the file corresponds to reference distribution.
        splitted_values = lines[values_start_index - 2].split('\t')
        x_cordinates = [float(value.strip()) for value in splitted_values[1:-1]]

        y_cordinates = []
        # Lines to be traversed
        iterations = values_end_index - values_start_index

        # Iterate `n` times over each line with values we'd want to extract:
        for i in range(iterations):
            """
            Every line of the file looks like this:
            '-11.810   \t    992 \t    988 \t  <... more values ...>  983 \t\n'
            We want list of the values (without the coordinate - the first
            value 11.810). The expected result looks like this:
            [992, 987, 988, <.. more values..>, 983]
            """
            index = values_start_index + i  # index of the line to be processed

            splitted_values = lines[index].split('\t')
            # Remove empty spaces
            result = [data_factor * float(value.strip()) for value in splitted_values[1:-1]]
            # First value is coordinate & last is newline symbol - remove them
            matrix.append(result)
            y_cordinates.append(float(splitted_values[0]))
            # If the the file corresponds to  reference dose distribution
            # exclude the doses outside the boundary of the MatriXX.
        if Flag:
            x_start = search_cord_start(x_cordinates)
            y_start = search_cord_start(y_cordinates)
            # Exclude all doses outside the boundaries of the MatriX detectors.
            matrix = [matrix[i][x_start: x_start + 32] for i in range(y_start, y_start + 32)]
    if deviation:
        return (matrix, x_cordinates[x_start: x_start + 32], y_cordinates[y_start: y_start + 32])
    return matrix


def angle_gen():
    """
    Create a list of all angles the correction factors correspond to.
    """
    angles = []
    for idx in range(len(measured_angles)):
        ang = measured_angles[idx]
        ang = math.radians(ang)
        angles.append(ang)
    angles.append(math.radians(360))
    return (angles)


def madc(matrix, mean_value):
    """
    Calculate the absolute mean deviation of numbers in list.
    """
    abs_deviation = 0
    for idx in range(len(matrix)):
        abs_deviation += (abs(matrix[idx] - mean_value) / len(matrix))
    return abs_deviation


def clean(m1):
    """
    Helper function to convert list of lists to one list
    """
    matrix = []
    for row in m1:
        matrix = matrix + row
    return matrix


def loadSTD(f_names):
    """
    Calculate the relative mean absolute deviations
    and standard deviations of the doses.
    """
    std_doses = []
    std_reslt = []
    mad_result = []
    # Parse the files
    std_ang_data = extract_matrix_values(f_names[0], True, True)
    x_cord = std_ang_data[1]
    y_cord = std_ang_data[2]
    std_doses.append(std_ang_data[0])
    for file_std in f_names[1:]:
        std_ang_data = extract_matrix_values(file_std, True)
        std_doses.append(std_ang_data)
    # x and y contain the x and y coordinates in terms of number of column and number of row
    x = []
    y = []
    for row in range(len(std_doses[0])):
        for col in range(len(std_doses[0][0])):
            temp = []
            for number in range(len(std_doses)):
                temp.append(std_doses[number][row][col])
            mean_value_row_col = statistics.mean(temp)
            if mean_value_row_col == 0:
                mean_dev_row_col = 0
                stdev_data = 0
            else:     
                mean_dev_row_col = madc(temp, mean_value_row_col) * 100 / mean_value_row_col
                stdev_data = statistics.stdev(temp) * 100 / mean_value_row_col
            std_reslt.append(stdev_data)
            mad_result.append(mean_dev_row_col)
            x.append(x_cord[col])
            y.append(y_cord[row])

    mean_std = round(statistics.mean(std_reslt), 4)
    mean_mad = round(statistics.mean(mad_result), 4)
    # Plot the data.
    a = plt.subplot2grid((4, 9), (0, 0), rowspan=4, colspan=4)
    plt.xlabel("columns (X axis)[cm] \n  Average STD:" + str(mean_std) + " %")
    plt.ylabel('rows (Y axis)[cm]')
    a2 = plt.subplot2grid((4, 9), (0, 4), rowspan=4, colspan=4)
    plt.xlabel('\n  Average MAD:' + str(mean_mad) + " %")
    cbar_ax = plt.subplot2grid((4, 9), (0, 8), rowspan=4, colspan=1)
    plt.xlabel("  levels  \n" + "%" + " of average value")
    minial = min(min(std_reslt), min(mad_result)) * 0.92
    maxial = max(max(std_reslt), max(mad_result)) * 1.08
    a_std = a.scatter(x, y, c=std_reslt, vmin=minial, vmax=maxial, cmap='nipy_spectral')
    a_mad = a2.scatter(x, y, c=mad_result, vmin=minial, vmax=maxial, cmap='nipy_spectral')
    plt.colorbar(a_std, cax=cbar_ax)
    a.set_xlim([x_cord[0] - 1, x_cord[-1] + 1])
    a.set_ylim([y_cord[0] - 1, y_cord[-1] + 1])
    a2.set_xlim([x_cord[0] - 1, x_cord[-1] + 1])
    a2.set_ylim([y_cord[0] - 1, y_cord[-1] + 1])
    plt.suptitle("Standard relative deviations                                    Mean relative absolute deviations", fontsize=10)

    plt.show()


def natural_keys(text):
    """
    Keys to sort by, the list of files
    """
    sort_by = re.findall('[0-9]*[.]?[0-9]+', text)
    if sort_by == []:
        popupmsg('The filename should end with angle in degrease \n Example name: abc_80.opg')
    return float(sort_by[-1])


def sort_ang(f_names, Flag):
    """
    Sort the files by the angle they correspond to in ascending order.
    """

    f_names.sort(key=natural_keys)

    ang_list = []
    for item in f_names:
        ang = float(re.findall('[0-9]*[.]?[0-9]+', item)[-1])
        ang_list.append(ang)
    if Flag:
        global cal_filenames
        global ref_angles
        ref_angles = ang_list
        cal_filenames = f_names
    else:
        global meas_filenames
        global measured_angles
        measured_angles = ang_list
        meas_filenames = f_names


def load( obj, controller):
    """
    Read the data from all imported files and store it in global variables
    """
    if len(measured_angles) != len(ref_angles):
        popupmsg('The number of imported files corresponding to reference distribution \n should match the number of files corresponding to reference distribution')
        return
    if measured_angles != ref_angles:
        popupmsg('The filename should end with angle in degrease \n Example name: abc_80.opg.\nThe angles for the reference and measured doses do not match')
        return
    global data_mesh
    global data_cal

    data_mesh = []
    data_cal = []
    for i in range(0, len(measured_angles)):
        file_1 = meas_filenames[i]
        file_2 = cal_filenames[i]
        """

        `matrix_1` represent the measured data
        `matrix_2` - the TPS data
        """
        matrix_1 = extract_matrix_values(file_1)
        matrix_2 = extract_matrix_values(file_2, True)
        if len(matrix_1) != 32 or len(matrix_1[0]) != 32:
            popupmsg("Imported reference distribution instead of measured")
            return
        data_mesh.append(matrix_1)
        data_cal.append(matrix_2)
    controller.show_frame(PageSymmetry)


def ang_sym_chek2(m1, m2):
    """
    Calculate the ration of the doses stored in the 2D array m1
    at coordinates row = i,column =j to the doses stored in the 2D array m2
    at coordinates row = i,column = -j-1 .
    """

    matrix = []
    for row_index in range(len(m2)):

        for value_index in range(len(m2[0])):
            new_value = (m1[row_index][value_index] /
                         m2[row_index][-1 - value_index])
            matrix.append(new_value)
    return matrix


def ang_sym_chek(Flag=False):
    """
    Calculate the ratios of doses at angle tita to doses at angle -tita
    """

    sym_chek_ration = []
    if Flag:
        data_to_chek = data_mesh
        output_str = 'measured doses'
    else:
        data_to_chek = data_cal
        output_str = 'reference doses'

    data_at_zero = ang_sym_chek2(data_to_chek[0], data_to_chek[0])
    sym_chek_ration.append(data_at_zero)
    end = len(data_to_chek)

    for idx in range(1, end):
        data_to_chek_sym = data_to_chek[-idx]
        new_data = ang_sym_chek2(data_to_chek[idx], data_to_chek_sym)
        sym_chek_ration.append(new_data)
    sym_chek_ration.append(data_at_zero)
    mean_up = []
    mean_down = []
    list_mean = []
    angles = angle_gen()
    std_cf_up = []
    std_cf_down = []
    hist_data = []
    stop_iterr = len(sym_chek_ration)
    for i in range(stop_iterr):
        d = sym_chek_ration[i]
        mean_v = statistics.mean(d)
        mean_dev_ang = madc(d, mean_v)
        list_mean.append(mean_v)
        stdev_data = statistics.stdev(d)
        std_cf_up.append(mean_v + stdev_data)
        std_cf_down.append(mean_v - stdev_data)
        mean_up.append(mean_v + mean_dev_ang)
        mean_down.append(mean_v - mean_dev_ang)
        if i != (stop_iterr - 1):
            hist_data = hist_data + d

    lower_limit = min(std_cf_down) - 0.01
    upper_limit = max(std_cf_up) + 0.01
    mean_hist = statistics.mean(hist_data)
    hist_std = round(statistics.stdev(hist_data), 5)
    hist_mad = round(madc(hist_data, mean_hist), 5)
    f = plt.figure()

    a = f.add_subplot(121)
    plt.xlabel(output_str + 'D(theta)/D(-theta) distribution \n,' +
               "STD: " + str(hist_std) + ",  MAD: " + str(hist_mad))
    plt.ylabel('Count')
    a2 = f.add_subplot(122, projection='polar')
    plt.xlabel(output_str + 'D(theta)/D(-theta) angular distribution,')
    a.hist(hist_data, bins=101)

    a2.set_theta_zero_location('N')
    a2.set_theta_direction(-1)

    a2.axis([0, 2 * math.pi, lower_limit, upper_limit])

    a2.plot(angles, mean_up, color='g')
    a2.plot(angles, mean_down, label='MAD', color='g')

    a2.plot(angles, std_cf_up, '--', linewidth=0.85, color='Navy')
    a2.plot(angles, std_cf_down, '--', linewidth=0.85, label='STD', color='Navy')

    a2.plot(angles, list_mean, label='mean', color='Red')
    a2.plot(angles, list_mean, 'ro', markersize=2, color='Black')
    plt.legend(loc=(0, 1))
    plt.suptitle("D(theta)/D(-theta)", fontsize=10)
    plt.show()


def find_angles_index(ang_list):
    """
    Given a list of angles find their indexes in the data set of all angles
    """

    angles_idx = []
    for angle in ang_list:
        idx = measured_angles.index(float(angle))
        angles_idx.append(idx)
    return angles_idx


def calulate_norm(Flag):
    """
    Assign to the measured and reference doses at angle 0
    two different variables.
    Iterate over each value of the measured doses,
    get the value at the same position from the reference doses,
    then perform the nessesery calculations
    to calculate the normalization factor.
    """

    mesh = data_mesh[0]
    cal = data_cal[0]
    matrix = []
    global norm_fact
    if Flag:
        mesh_sym = symetry(mesh)
        mesh_sym_sym = ang_symetry(mesh_sym, mesh_sym)
        cal_sym = symetry(cal)
        cal_sym_sym = ang_symetry(cal_sym, cal_sym)

        for col_index in range(len(mesh_sym_sym)):

            for value_index in range(len(mesh_sym_sym[col_index])):

                # `c` refers to `condition
                c1 = col_index == 0 and value_index == 0
                c2 = col_index == 31 and value_index == 31
                c3 = col_index == 31 and value_index == 0
                c4 = col_index == 0 and value_index == 31
                if (c1 or c2 or c3 or c4):
                    continue

                matrix_mesh_value = mesh_sym_sym[col_index][value_index]
                matrix_cal_value = cal_sym[col_index][value_index]
                new_value = matrix_mesh_value / matrix_cal_value
                matrix.append(new_value)

        norm_fact = statistics.mean(matrix)
        matrix = []
        for dummy_row in range(len(cal_sym_sym)):
            row = []
            for dummy_col in range(len(cal_sym_sym)):
                row.append(norm_fact)
            matrix.append(row)
        norm_fact = matrix

    else:
        for col_index in range(len(mesh)):
            row = []
            for value_index in range(len(mesh)):

                matrix_mesh_value = mesh[col_index][value_index]
                matrix_cal_value = cal[col_index][value_index]
                new_value = matrix_mesh_value / matrix_cal_value
                row.append(new_value)
            matrix.append(row)
        norm_fact = matrix


def symetry(m):
    """
    Symetrisize the values of the matrix on the y axis.
    """
    matrix = []
    for row_index in range(len(m)):
        result_row = []

        for value_index in range(len(m[0])):
            # symerization of matrix1 only by the row_index value (y axis)
            new_value = (float(m[row_index][value_index]) +
                         float(m[-row_index - 1][value_index])) / 2
            result_row.append(new_value)

        matrix.append(result_row)

    return matrix


def ang_symetry(m1, m2):
    """
    Symetrisize the doses at angle theta and - theta.
    """
    matrix = []
    for row_index in range(len(m2)):
        result_row = []

        for value_index in range(len(m2[0])):
            new_value = (m1[row_index][value_index] +
                         m2[row_index][-1 - value_index]) / 2
            result_row.append(new_value)

        matrix.append(result_row)
    return matrix


def calculate_value(v_mesh, v_cal, row, col):

    """
    Calculate each individual normalized correction factor.
    """
    return round(v_cal * norm_fact[row][col] / v_mesh, 5)


def slinding_window(matrix, row_pos, col_pos, rows, cols):
    """
    The function determines the biggest possible symmetrical sliding window
    over which the doses can be averaged and the averages them.
    The values at the MatriXX corners are excluded from the averaging,
    because the Omni-Pro software uses interpretation to determine the doses at
    the corners.
    """
    rows = math.floor(rows / 2)
    cols = math.floor(cols / 2)
    num_col = len(matrix[0]) - 1
    num_row = len(matrix) - 1

    row_top_boundry = rows + row_pos
    if num_row < row_top_boundry:
        rows_up = num_row - row_pos
    else:
        rows_up = rows
    row_bottom_boundry = row_pos - rows

    if 0 > row_bottom_boundry:
        rows_down = row_pos
    else:
        rows_down = rows
    rows = min(rows_down, rows_up)
    col_top_boundry = cols + col_pos
    if num_col < col_top_boundry:
        cols_up = num_col - col_pos
    else:
        cols_up = cols

    col_bottom_boundry = col_pos - cols
    if 0 > col_bottom_boundry:
        cols_down = col_pos
    else:
        cols_down = cols
    cols = min(cols_down, cols_up)
    start_row = row_pos - rows
    end_row = row_pos + rows + 1
    start_col = col_pos - cols
    end_col = col_pos + cols + 1
    summ = 0
    samples = 0
    for row in range(start_row, end_row):
        for col in range(start_col, end_col):
            c1 = (col == 0 and row == 0)
            c2 = (col == 31 and row == 31)
            c3 = (col == 31 and row == 0)
            c4 = (col == 0 and row == 31)
            if (c1 or c2 or c3 or c4):
                continue
            else:
                summ += matrix[row][col]
                samples += 1

    return summ / samples


def calculations(m, how):
    """
    Iterate over all values in the array and apply the sliding window function
    """
    matrix = []
    rows = float(how[0])
    cols = float(how[1])
    for row_index in range(len(m)):
        result_row = []

        for value_index in range(len(m[row_index])):
            # `c` refers to `condition`
            c1 = row_index == 0 and value_index == 0
            c2 = row_index == 31 and value_index == 31
            c3 = row_index == 31 and value_index == 0
            c4 = row_index == 0 and value_index == 31

            if (c1 or c2 or c3 or c4):
                matrix_value = m[row_index][value_index]
            else:
                matrix_value = slinding_window(m, row_index, value_index, rows, cols)
            new_value = matrix_value
            result_row.append(new_value)

        matrix.append(result_row)

    return matrix


def calculations_at_0(m_mesh, m_cal):
    """
    If the user decides that correction at angle zero is not needed
    ,than to all correction factors is assign value of 1.
    Correction factors at the MatriXX corners are calculated, because
    the Omni Pro software uses interpretation to determine the doses at
    the corners.

    """
    matrix = []

    for row_index in range(len(m_cal)):
        result_row = []

        for value_index in range(len(m_cal[row_index])):

            # `c` refers to `condition
            c1 = row_index == 0 and value_index == 0
            c2 = row_index == 31 and value_index == 31
            c3 = row_index == 31 and value_index == 0
            c4 = row_index == 0 and value_index == 31
            # Matrix_2 ia allready symmetrized on both cordinates
            if (c1 or c2 or c3 or c4):
                matrix_mesh_value = m_mesh[row_index][value_index]
                matrix_cal_value = m_cal[row_index][value_index]
                new_value = calculate_value(matrix_mesh_value, matrix_cal_value, row_index, value_index)
            else:
                new_value = 1.
            result_row.append(new_value)

        matrix.append(result_row)

    return matrix


def do_col_calculations(m1):
    """
    Iterate over all values in the 2D array m1,then average all values on
    the same column and assign them to a new array.
    """
    matrix = []

    for col_index in range(len(m1)):
        result_row = []

        for value_index in range(len(m1[col_index])):
            matrix_1_value = 0
            # `c` refers to `condition`
            c1 = col_index == 0 and value_index == 0
            c2 = col_index == 31 and value_index == 31
            c3 = col_index == 31 and value_index == 0
            c4 = col_index == 0 and value_index == 31

            if (c1 or c2 or c3 or c4):
                matrix_1_value = m1[col_index][value_index]
            elif (value_index == 0 or value_index == 31):
                for i in range(1, 31):
                    matrix_1_value += m1[i][value_index] / 30
            else:
                for i in range(0, 32):
                    matrix_1_value += m1[i][value_index] / 32

            new_value = matrix_1_value
            result_row.append(new_value)

        matrix.append(result_row)

    return matrix


def do_row_calculations(m1):
    """
    Iterate over all values in the 2D array m1,then average all values on
    the same row and assign them to a new array.
    """
    matrix = []

    for row_index in range(len(m1)):
        result_row = []

        for value_index in range(len(m1[row_index])):
            matrix_1_value = 0
            if ((row_index == 0 and value_index == 0) or
                    (row_index == 31 and value_index == 31) or
                    (row_index == 31 and value_index == 0) or
                    (row_index == 0 and value_index == 31)):
                matrix_1_value = m1[row_index][value_index]
            elif (row_index == 0 or row_index == 31):
                for i in range(1, 31):
                    matrix_1_value += m1[value_index][i] / 30
            else:
                for i in range(0, 32):
                    matrix_1_value += m1[value_index][i] / 32
            new_value = matrix_1_value
            result_row.append(new_value)

        matrix.append(result_row)

    return matrix


def calculate_cfs(TPS, mesurment):
    """
    Iterate over 2D arrays and calculate all concretion factories
    at the angle the data corresponds to.
    """
    matrix = []
    for row in range(len(TPS)):
        result_row = []
        for col in range(len(TPS[row])):
            v_mesh = mesurment[row][col]
            v_cal = TPS[row][col]
            new_value = calculate_value(v_mesh, v_cal, row, col)
            result_row.append(new_value)
        matrix.append(result_row)
    return matrix


def calirbrate(obj, controller, Flag=False):
    """
    Process the user input than
    calculate all correction factors
    according to the user chosen settings.
    """
    try:
        norm_fact
    except NameError:
        popupmsg("Normalization factor not defined")
        return
    if Flag:
        LUT2 = []
    else:
        global LUT
        LUT = []

    mesh_calib = obj.mesh_how_to_calibrate.get()

    if 'ROWS' in mesh_calib:
        mesh_calib = 'ROWS'
        if 'COLUMNS' in mesh_calib:
            mesh_calib = 'both'
    elif 'COLUMNS' in mesh_calib:
        mesh_calib = 'COLUMNS'
    else:
        mesh_calib = re.findall('[0-9]*[.]?[0-9]+', mesh_calib)

    mesh_diff_calib = obj.mesh_how_to_calibrate_difrently.get()
    if 'ROWS' in mesh_diff_calib:
        mesh_diff_calib = 'ROWS'
        if 'COLUMNS' in mesh_diff_calib:
            mesh_diff_calib = 'both'
    elif 'COLUMNS' in mesh_diff_calib:
        mesh_diff_calib = 'COLUMNS'
    else:
        mesh_diff_calib = re.findall('[0-9]*[.]?[0-9]+', mesh_diff_calib)

    mesh_ang_diff = obj.mesh_ang_to_calibrate_difrently.get()
    mesh_idx_diff = re.findall('[0-9]*[.]?[0-9]+', mesh_ang_diff)
    if mesh_idx_diff != []:
        mesh_idx_diff = find_angles_index(mesh_idx_diff)

    cal_calib = obj.cal_how_to_calibrate.get()
    if 'ROWS' in cal_calib:
        cal_calib = 'ROWS'
        if 'COLUMNS' in cal_calib:
            cal_calib = 'both'
    elif 'COLUMNS' in cal_calib:
        cal_calib = 'COLUMNS'
    else:
        cal_calib = re.findall('[0-9]*[.]?[0-9]+', cal_calib)

    cal_diff_calib = obj.cal_how_to_calibrate_difrently.get()
    if 'ROWS' in cal_diff_calib:
        cal_diff_calib = 'ROWS'
        if 'COLUMNS' in cal_diff_calib:
            cal_diff_calib = 'both'
    elif 'COLUMNS' in cal_diff_calib:
        cal_diff_calib = 'COLUMNS'
    else:
        cal_diff_calib = re.findall('[0-9]*[.]?[0-9]+', cal_diff_calib)

    cal_ang_diff = obj.cal_ang_to_calibrate_difrently.get()
    cal_idx_diff = re.findall('[0-9]*[.]?[0-9]+', cal_ang_diff)
    if len(cal_ang_diff) != 0:
        cal_idx_diff = find_angles_index(cal_idx_diff)

    cond0 = not(obj.corr_at_0.get())

    for idx in range(len(data_mesh)):
        if idx == 0 and cond0:
            factors_0 = calculations_at_0(data_mesh[0], data_cal[0])
            if Flag:
                LUT2.append(factors_0)
            else:
                LUT.append(factors_0)
            continue
        mesurment = data_mesh[idx]
        if idx in mesh_idx_diff:
            if (obj.mesh_y_sym_other.get()):
                mesurment = symetry(mesurment)
            if (obj.mesh_ang_sym_other.get()):
                if idx == 0:
                    mesurment = ang_symetry(mesurment, mesurment)
                else:
                    sym_mesurment = data_mesh[-idx]
                    if (obj.mesh_y_sym_other):
                        sym_mesurment = symetry(sym_mesurment)
                    mesurment = ang_symetry(mesurment, sym_mesurment)
            if mesh_diff_calib == 'both':
                mesurment = do_row_calculations(mesurment)
                mesurment = do_col_calculations(mesurment)
            elif mesh_diff_calib == 'ROWS':
                mesurment = do_row_calculations(mesurment)
            elif mesh_diff_calib == 'COLUMNS':
                mesurment = do_col_calculations(mesurment)
            else:
                mesurment = calculations(mesurment, mesh_diff_calib)

        else:
            if (obj.mesh_y_sym.get()):
                mesurment = symetry(mesurment)
            if (obj.mesh_ang_sym.get()):
                if idx == 0:
                    mesurment = ang_symetry(mesurment, mesurment)
                else:
                    sym_mesurment = data_mesh[-idx]
                    if (obj.mesh_y_sym.get()):
                        sym_mesurment = symetry(sym_mesurment)
                    mesurment = ang_symetry(mesurment, sym_mesurment)
            if mesh_calib == 'both':
                mesurment = do_row_calculations(mesurment)
                mesurment = do_col_calculations(mesurment)
            elif mesh_calib == 'ROWS':
                mesurment = do_row_calculations(mesurment)
            elif mesh_calib == 'COLUMNS':
                mesurment = do_col_calculations(mesurment)
            else:
                mesurment = calculations(mesurment, mesh_calib)

        TPS = data_cal[idx]
        if idx in cal_idx_diff:
            if (obj.cal_y_sym_other.get()):
                TPS = symetry(TPS)
            if (obj.cal_ang_sym_other.get()):
                if idx == 0:
                    TPS = ang_symetry(TPS, TPS)
                else:
                    sym_TPS = data_cal[-idx]
                    if (obj.cal_y_sym_other.get()):
                        sym_TPS = symetry(sym_TPS)
                    TPS = ang_symetry(TPS, sym_TPS)
            if cal_diff_calib == 'both':
                TPS = do_row_calculations(TPS)
                TPS = do_col_calculations(TPS)
            elif cal_diff_calib == 'ROWS':
                TPS = do_row_calculations(TPS)
            elif cal_diff_calib == 'COLUMNS':
                TPS = do_col_calculations(TPS)
            else:
                TPS = calculations(TPS, cal_diff_calib)

        else:
            if (obj.cal_y_sym.get()):
                TPS = symetry(TPS)
            if (obj.cal_ang_sym.get()):
                if idx == 0:
                    TPS = ang_symetry(TPS, TPS)
                else:
                    sym_TPS = data_cal[-idx]
                    if (obj.cal_y_sym.get()):
                        sym_TPS = symetry(sym_TPS)
                    TPS = ang_symetry(TPS, sym_TPS)
            if cal_calib == 'both':
                TPS = do_row_calculations(TPS)
                TPS = do_col_calculations(TPS)
            elif cal_calib == 'ROWS':
                TPS = do_row_calculations(TPS)
            elif cal_calib == 'COLUMNS':
                TPS = do_col_calculations(TPS)
            else:
                TPS = calculations(TPS, cal_calib)
        cfs = calculate_cfs(TPS, mesurment)
        if Flag:
            LUT2.append(cfs)
        else:
            LUT.append(cfs)
        if idx == 0:
            factors_0 = cfs
    if Flag:
        LUT2.append(factors_0)
    else:
        LUT.append(factors_0)
    if not Flag:
        controller.show_frame(PageCalibrationGraphs)
    else:
        calculate_LUTs_ratio(LUT2)
    return


def calculate_LUTs_ratio(LUT2):
    """
    Iterate over all values in two 3D(angle,Y,X) arrays
    then calculate the ratio values with the same coordinates
    """
    LUT_ration = []
    row_itter = len(LUT[0])
    col_itter = len(LUT[0][0])
    for angle in range(len(LUT)):
        angle_ration = []
        for row in range(row_itter):
            for col in range(col_itter):
                new_value = LUT[angle][row][col] / LUT2[angle][row][col]
                angle_ration.append(new_value)
        LUT_ration.append(angle_ration)
    correction_fact_ration_distribution(LUT_ration)


def first_file_row(energy, beam_quality, LUT_name):
    """
    Add header information to the csv file
    """
    name = asksaveasfilename(filetypes=[('CSV files', '*.csv')], title='Export correction factors', defaultextension='.csv')

    with open(name, "w", newline='') as csv_file:
        fieldnames = ['Linac: '+ LUT_name + ' model, ' + str(energy) +', '+ str(beam_quality)  + ', 32, 32']
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames, delimiter=';')
        writer.writeheader()
        csv_file.close()
    return name


def transform_matrix(matrix, angle, name):
    """
    Fill the correction factor a given angle in the csv file
    """
    with open(name, "a", newline='') as csv_file:
        csv_app = csv.writer(csv_file)
        csv_app.writerow([angle])
        csv_app.writerows(matrix)
        csv_file.close()


def Export_cfs(obj):
    """
    Transform the matrix in csv file and fills
    the angles the correction factors correspond to.
    """
    quality = obj.quality.get()
    energy = obj.energy.get()
    LUT_name = str(obj.LUT_name.get())
    energy = re.findall('[0-9]*[.]?[0-9]+', energy)
    quality = re.findall('[0-9]*[.]?[0-9]+', quality)
    if (len(energy) != 1) or (len(quality) != 1) or (len(LUT_name) == 0):
        popupmsg('Fill the fileds: beam quality, beam energy and LUT name')
        return 
    f_name = first_file_row(energy[0], quality[0], LUT_name)
    for idx in range(len(ref_angles)):
        matrix = LUT[idx]
        angle = ref_angles[idx]
        transform_matrix(matrix, angle, f_name)
    transform_matrix(LUT[0], 360.0, f_name)


def correction_fact_distribution(obj):
    """
    Calculate the correction factor mean value,
    standard deviation and absolute standard deviation
    as a function of the gantry angle and plot then.
    Create a histogram plot of the correction factor distribution.
    """
    mean_up = []
    mean_down = []
    list_mean = []
    angles = angle_gen()
    std_cf_up = []
    std_cf_down = []
    hist_data = []
    stop_iterr = len(LUT)
    for i in range(stop_iterr):
        d = clean(LUT[i])
        mean_v = statistics.mean(d)
        mean_dev_ang = madc(d, mean_v)
        list_mean.append(mean_v)
        stdev_data = statistics.stdev(d)
        std_cf_up.append(mean_v + stdev_data)
        std_cf_down.append(mean_v - stdev_data)
        mean_up.append(mean_v + mean_dev_ang)
        mean_down.append(mean_v - mean_dev_ang)
        if i != (stop_iterr - 1):
            hist_data = hist_data + d

    lower_limit = min(std_cf_down) - 0.01
    upper_limit = max(std_cf_up) + 0.01

    f = plt.figure()

    a = f.add_subplot(121)
    plt.xlabel('Correction factors distribution')
    plt.ylabel('Count')

    a2 = f.add_subplot(122, projection='polar')
    plt.xlabel('Correction factors angle distribution')
    a.hist(hist_data, bins=101)
    a2.set_theta_zero_location('N')
    a2.set_theta_direction(-1)

    a2.axis([0, 2 * math.pi, lower_limit, upper_limit])

    a2.plot(angles, mean_up, color='g')
    a2.plot(angles, mean_down, label='MAD(CF)', color='g')

    a2.plot(angles, std_cf_up, '--', linewidth=0.85, color='Navy')
    a2.plot(angles, std_cf_down, '--', linewidth=0.85, label='STD', color='Navy')

    a2.plot(angles, list_mean, label='mean(CF)', color='Red')
    a2.plot(angles, list_mean, 'ro', markersize=2, color='Black')
    plt.legend(loc=(0, 1))
    plt.show()


def correction_fact_ration_distribution(LUT_ration):
    """
    Calculate the correction factor ratio,mean value,
    standard deviation and absolute standard deviation
    as a function of the gantry angle and plot then.
    Create a histogram plot of the ratios distribution and calculate .
    """
    mean_up = []
    mean_down = []
    list_mean = []
    angles = angle_gen()
    std_cf_up = []
    std_cf_down = []
    hist_data = []
    stop_iterr = len(LUT_ration)
    for i in range(stop_iterr):
        d = LUT_ration[i]
        mean_v = statistics.mean(d)
        mean_dev_ang = madc(d, mean_v)
        list_mean.append(mean_v)
        stdev_data = statistics.stdev(d)
        std_cf_up.append(mean_v + stdev_data)
        std_cf_down.append(mean_v - stdev_data)
        mean_up.append(mean_v + mean_dev_ang)
        mean_down.append(mean_v - mean_dev_ang)
        if i != (stop_iterr - 1):
            hist_data = hist_data + d

    lower_limit = min(std_cf_down) - 0.01
    upper_limit = max(std_cf_up) + 0.01
    mean_hist = statistics.mean(hist_data)
    hist_std = round(statistics.stdev(hist_data), 5)
    hist_mad = round(madc(hist_data, mean_hist), 5)
    f = plt.figure()

    a = f.add_subplot(121)
    plt.xlabel('Correction factors ration distribution \n' + 'STD: ' + str(hist_std) + ",   MAD:" + str(hist_mad))
    plt.ylabel('Count')
    a2 = f.add_subplot(122, projection='polar')
    plt.xlabel('Correction factors ratio angular distribution')
    a.hist(hist_data, bins=101)

    a2.set_theta_zero_location('N')
    a2.set_theta_direction(-1)

    a2.axis([0, 2 * math.pi, lower_limit, upper_limit])

    a2.plot(angles, mean_up, color='g')
    a2.plot(angles, mean_down, label='MAD', color='g')

    a2.plot(angles, std_cf_up, '--', linewidth=0.85, color='Navy')
    a2.plot(angles, std_cf_down, '--', linewidth=0.85, label='STD', color='Navy')

    a2.plot(angles, list_mean, label='mean(ration)', color='Red')
    a2.plot(angles, list_mean, 'ro', markersize=2, color='Black')
    plt.legend(loc=(0, 1))
    plt.show()


#################BEGIN THE GUI##########################################


def tutorial():
    """
    Application help menu
    """
    app_intro = ("""                                            With this application, it is possible to:

            - evaluate the consistency of the measured/reference doses

            - evaluate the change in the dose distribution at different gantry angles

            - investigate the MatriXX angular dependency

            - compare different calibrations

            - create correction factor sets in the format required by Omni Pro I'mRT software
         """)
    file_specification = ("""                                                                          Files specifications: 

            The measured data at each angle should be imported into OmniPro-I’mRT software and exported as ASCII file. The reference dose distributions 
        should be interpolated to 7.619 mm grid (If OmniPro-I'mRT is used for the interpolation the ages are aligned to the 2 central lines) and 
        exported as ASCII files. In case exact interpolation to 7.619 mm grid is not possible, small deviations are allowed, in this case 7.62 mm grid 
        is recommended.   
            Each file name should contain information about the angle it corresponds to and no numbers are allowed after the angle. For example, the  
        ASCII file corresponding to measurements at angle 185˚ can be named “meas_185.opg". The data separator must be TAB. 
        The minimum information a file must contain is shown below: data factor and X and Y coordinates of the detectors.

        Separator:          [TAB]
        Data Factor:        ...

        <asciibody>

        X  ...
        Y
        .  ...
        .  ...
        .  ...

                         """)

    evaluete_dev = ("""                          Evaluating  the consistency the measured/reference doses  
                  or evaluating  the change in the dose distribution at different gantry angles.

    Click on Stat Page "Evaluate STD and MAD".Then import the files of interest.
    After importing the files, relative standard(STD) and absolute deviations(MAD)
    will be automatically calculated for each chamber.
    The average value of the deviations is also calculated.
        """)

    file_import_cfs = ("""                                                         Importing files for calibration

    Select all files (press ctrl then chose the files) then click OPEN. 
    The reference and measured dose distributions should correspond to the same set of angles
        """)
    Symmetry_tut = ("""                                                 Investigate doses difference for gantry angles 0˚-180˚ and 180˚ - 360˚

    Click "Proceed correction factor calculation". Then import the files you want to investigate. 
    The reference and measured dose distributions must correspond to the same set of angles. 
    Clicking the buttons "Evaluate measured dose angle symmetry" or "Evaluate measured dose angle symmetry" 
    will display the ratios of dose at position row = i, column = j and angle = theta and their standard and absolute mean deviations 
    to dose at position row = i ,column = number of columns + 1 - j, angle = -theta. The mean value of the rations along 
    with their standard deviation(STD) and absolute mean deviation(MAD) as function of the gantry angle will also be displayed.
        """)
    calibration_tut_intro = ("""                                                                         Calibration menu 

    To enter the calibration menu reference and measured dose distributions have to be imported. 
    The reference and measured dose distributions must correspond to the same set of angles. 
    Once the files are imported it is not necessary to import them again. 
    Measured and reference doses are processed separately.  
    First moving average parametres have to be chosen. 
    All other option by default are set inactive and are not mandatory.
    The top buttons for symmetrization regard to the angles not listed in the field "Enter angles to be calibrated differently". 
    The bottom buttons regard to the angles listed in the field "Enter angles to be calibrated differently".
        """)

    cfs_calculation_tut = ("""                                                              Calibration menu  
                                                 How calibration factors are calculated   

    Each factor is calculated from the ratio of: 
    reference dose at row = i, column = j, angle = theta to measured dose at  
    position row = i, column = j, angle = theta. Then the ratio is multiplied by  
    normalization factor.
    Before the calculation, measured and the reference doses can be processed differently.
        """)

    norm_fact_tut = ("""                                                                                    Calibration menu  
                                                                                Normalization factor 

    At the top left corner, there two are possible ways to calculate normalization factor:
    1. Correction factor individual for each detector. 
        Individual correction factor is calculated for every detector by the ration of 
    reference dose at position row = i, column = j, angle = 0˚to measured dose at  
    position row = i, column = j, angle = 0˚. 

    2.The same for all detectors. 
        The doses at angle are 0˚ are separately processed and then their ration calculates. 
    The calculation includes 4 steps: 
    a)Symmetrization on the X axis: 
    The doses at position row = i, column = j and doses dose at position row = i ,column = total number of columns+ 1 - j are averaged. 
    b)Symmetrization on the Y axis: 
    The doses at position row = i, column = j and doses dose at position row = number of row + 1 - i ,column = j are averaged. 
    c)The ratio of reference dose at position row = i, column = j to measured dose at position row = i, column = j is calculated. 
    d)All ratios are averaged in order to obtain a signal correction factor for all detectors. 
    *The doses at the corners of the MatriXX (there are no detectors and the dose is interpolated via the MatriXX software) are excluded 
    from the normalization factor calculatio, for them normalization factor is calculated as (1. Correction factor 
    individual for each detector). 
        """)

    norm_fact_err = ("""                             Error Normalization factor not defined

       Prior to calibration normalization factor have to be chosen. 
       In the top left corner (near the exit button) chose calibration factor.
        """)
    moving_avg_tut = ("""                                                                          Calibration menu
                                                                        Chose moving average 

    The top field "chose moving average" are mandatory. 
    The bottom fields are mandatory only if in the filed "Enter angles to be calibrated differently" 
    angles are enterd. 
    The sliding window helps to smooth out the doses by filtering out the “noise” from random fluctuations. 
    Sliding window over which the doses would be averaged can be chosen. Once the window is  
    beyond the borders of the MatriXX, the window will get smaller symmetrically on both sides. 
    For example, sliding window "row height: 1 ,column weight: 3", means that doses at positions row: i, column: j 
    is calculated by averaging the doses at positions: (i, j-1)+(i, j)+(i,j+1). 
    If even number is "N" inputted as a parameter in the sliding window the calculation will proceed as  
    N-1 is the input (because the number of neighbor detectors on both sides will be uneven if N is even). 
    If "row height: 1 ,column weight: 1" is chosen, the window will include only the detector itself. 
    If "ROWS" is entered in the filed then the doses corresponding to a detector on a given row.
    will be the average dose of all doses on the given row. 
    If "COLUMNS" is entered in the filed then the doses will be averaged over each column. 
    If both "ROWS" and "COLUMNS" are entered doses would be averaged over all detectors. 
    "ROWS" with sliding window for columns or "COLUMNS" with sliding window for rows is unavailabe and such input would result in 
    executing only in avereging the dose over row/column.
    *The doses at the corners of the MatriXX (there are no detectors and the dose is interpolated via the MatriXX software)
     are exclude from the window and are processed separately.
     *Enter only integer numbers in this field
        """)

    y_axis_sym_tut = ("""                                                        Calibration menu
                                                 Symmetrize on the Y axis

    The dose at position row = i,column = j ,angle = theta would be calculated as the average of 
    the doses at position row = i, column = j ,angle = theta, and doses  
    at position row = total number of rows + 1 - i ,column = j, angle=theta.""") 

  
    ang_sym_tut =("""                                                               Calibration menu
                                                          Angle symmetrization 

    The dose at position row = i, column = j ,angle = theta would be calculated as the average of 
    the doses at position row = i, column = j ,angle = theta, and doses  
    at position row = i ,column = total number of columns + 1- j, angle=-theta. """)

  
    diff_calib_tut = ("""                                                       Calibration menu
                                                 Filed:Enter angles to be calibrated differently  

    This field is optional. 
    Separate the angles in the input filed with commas (Example input:90,91,92). 
    Determine how to process the doses at the listed angles by the buttons and fields  
    under filed "Enter angles to be calibrated differently".""")

  

    extport_tut = ("""                                                       Export correction factor set 

     To export the factors as a csv file in the required by Omni-Pro I'mRT format filling the,
    fileds: beam quality, beam energy and LUT name is mandatory, then click export.
    *use comma as decimal separator 
     """)

  

    dist_graph_tut = ("""                                                                          Show graph (correction factors distribution)

    After calculating the correction factors click "Show graph" to see graphical representation of their distribution. 
    Histogram of all correction factors will be displayed, along their mean value, standard and mean absolute deviations as function of the gantry angles. 
    *correction factors at angle zero are also included in the histogram.""")

    compare_tut = ("""                                                                  Compare to other Calibration 

    To compare correction factors calculated using different parameters click the button "Compare to other Calibration". 
    First the calibration parameters have to be defined. Next click "compare". 
    The comparison is done by calculating the ratio: 
    correction factor from set1 at position row = i, column = j, angle = theta 
    correction factor from set2 at position row = i, column = j, angle = theta. 
    Histogram of all ratios will be displayed 
    (under the histogram standard deviation (STD)and mean absolute deviation (MAD) of all ratios will be displayed). 
    , along their mean value, standard and mean absolute deviations as function of the gantry angles. 
    *ratios at angle zero are also calculated.""")
    txt_list = [app_intro, file_specification, evaluete_dev, file_import_cfs,
                Symmetry_tut, calibration_tut_intro, cfs_calculation_tut, 
                norm_fact_tut, moving_avg_tut, y_axis_sym_tut, ang_sym_tut,
                 diff_calib_tut, extport_tut, dist_graph_tut, compare_tut]
    title_list =['App introduction', 'Files specifications', 'Dose deviations', 'Importing files for calibration',
                 'Evaluate dose angle symmetry', 'Calibration menu', 'How correction factors are calculated ', 
                 'Normalization factor', 'Moving average', 'Symmetrize on the Y axis', 'Angle symmetrization',
                 'Angles to be calibrated differently', 'Export file', 'Correction factors distribution', 'Compare to other Calibration']
    
    last_page = len(title_list) - 1

    def leavemini(what, index):
        """
        Destroy opend page and open the next page.
        """
        what.destroy()
        StepByStepTut(index + 1)

    def StepByStepTut(index):
        """
        Foce a "step by step" turturil behavior
        """

        tut = tk.Tk()
        tut.wm_title('Tutorial')
        label = ttk.Label(tut, text=txt_list[index], font=NORM_FONT)
        label.pack(side="top", fill="x", pady=10)
        if index == last_page:
            B = ttk.Button(tut, text="Done!",
                           command=tut.destroy)
        else:
            B = ttk.Button(tut, text="Next!",
                           command=lambda: leavemini(tut, index))
        B.pack(side="bottom")
        tut.minsize(700, 200)
        tut.eval('tk::PlaceWindow %s center' % tut.winfo_pathname(tut.winfo_id()))

        tut.mainloop()

    def More_tut(what):
        """
        Help menu
        """
        what.destroy()
        tut = tk.Tk()
        tut.wm_title('Help')
        for idx in range(last_page + 1):
            B = ttk.Button(tut, text=title_list[idx], width=33,
                           command=lambda idx=idx: popupmsg(txt_list[idx]))
            B.pack()
        tut.mainloop()

    tut = tk.Tk()
    tut.wm_title("Tutorial")
    label = ttk.Label(tut, text="What do you need help with?", font=NORM_FONT)
    label.pack(side="top", fill="x", pady=10)
    B1 = ttk.Button(tut, text="Overview of the application", width=30,
                    command=lambda: leavemini(tut, -1))
    B1.pack()


    B2 = ttk.Button(tut, text="Normalization factor not defined", width=30, command=lambda: popupmsg(norm_fact_err))
    B2.pack()

    B3 = ttk.Button(tut, text="More", width=30, command=lambda:More_tut(tut))
    B3.pack()
    tut.mainloop()


def popupmsg(msg):
    """
    Pop-up message
    """
    popup = tk.Tk()
    popup.wm_title("!")
    label = ttk.Label(popup, text=msg, font=NORM_FONT)
    label.pack(side="top", fill="x", pady=10)
    B1 = ttk.Button(popup, text="Okay", command=lambda: popup.destroy())
    B1.pack()
    popup.eval('tk::PlaceWindow %s center' % popup.winfo_pathname(popup.winfo_id()))
    popup.mainloop()


class CorrGui(tk.Tk):
    """
    GUI for the application
    """
    def __init__(self, *args, **kwargs):
        tk.Tk.__init__(self, *args, **kwargs)
        img = tk.PhotoImage(file='icon.gif')
        self.tk.call('wm', 'iconphoto', self._w, img)
        tk.Tk.wm_title(self, 'Correction of MatriXX angular dependency ')
        container = tk.Frame(self)
        container.pack(side="top", fill="both", expand=True)

        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)
        menubar = tk.Menu(container)

        menubar.add_command(label="Exit", command=lambda: app.destroy())
        choice = tk.Menu(menubar, tearoff=1)
        choice.add_command(label="Individual for each detector",
                           command=lambda: calulate_norm(False))
        choice.add_command(label="The same for all detectors",
                           command=lambda: calulate_norm(True))

        menubar.add_command(label="Tutorial", command=tutorial)
        tk.Tk.config(self, menu=menubar)
        menubar.add_cascade(label="Normalization factor", menu=choice)
        self.frames = {}
        for F in (StartPage, HomePage, PageDeviations,
                  PageImportData, PageSymmetry, PageCompareCalibrations,
                  PageCalibrateCustom, PageCalibrationGraphs):
            frame = F(container, self)
            self.frames[F] = frame
            frame.grid(row=0, column=0, sticky="nsew")
        self.show_frame(StartPage)

    def show_frame(self, cont):
        """
        Show the nessesery page
        """
        frame = self.frames[cont]
        frame.tkraise()


class StartPage(tk.Frame):
    """
    Start page
    """

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)

        label = tk.Label(self, text=("""App for angular correction of MatriXX use at your own risk
         no promise of waranty."""), font=NORM_FONT)
        label.place(x=200, y=5)
        button1 = ttk.Button(self, text="Agree",
                             command=lambda: controller.show_frame(HomePage))
        button1.place(x=330, y=50)
        button2 = ttk.Button(self, text="Disagree",
                             command=lambda: quit())
        button2.place(x=430, y=50)


        label2 = tk.Label(self, text=("""
        The program is intended for investigation of MatriXX angular dependency and creation of correction factor set to correct for it. 
    The measured data at each angle should  imported into OmniPro-I’mRT software and exported as ASCII file.
    The reference dose distributions should be interpolated to 7.619 mm grid (If OmniPro-I'mRT is used for the interpolation the ages are 
    aligned to the 2 central lines) and exported as ASCII files.In case exact interpolation to 7.619 mm grid is not possible, small
    deviations are allowed, in this case 7.62 mm grid is recommended.  
    Each file name should contain information about the angle it corresponds to and no numbers are allowed after the angle. 
    For example the ASCII file corresponding to measurements at angle 185˚ can be named “meas_185.opg”."""), font=SMALL_FONT)
        label2.place(x = 3, y = 90)


class HomePage(tk.Frame):
    """
    Home page
    You can go to :PageDeviations, or PageImportData
    """

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = tk.Label(self, text="Home", font=LARGE_FONT)
        label.pack(pady=10, padx=10)
        frame2 = tk.Frame(self)
        frame2.pack(side="top")

        button1 = ttk.Button(frame2, text="Evaluate STD and MAD",
                             command=lambda: controller.show_frame(PageDeviations), width=26)
        button1.pack(side="left", pady=10, padx=10)
        button2 = ttk.Button(frame2, text="Correction factors calculation ",
                             command=lambda: controller.show_frame(PageImportData), width=26)
        button2.pack(side="left", pady=10, padx=10)

        label2 = tk.Label(self, text="To investigate the mean and absolute deviation of doses/detectors responses \n at each detectors position click: Evaluate STD and MAD", font=NORM_FONT)
        label2.pack(side="top", pady=10, padx=10)


class PageDeviations(tk.Frame):
    """
    Page one
    You can go back to HomePage
    Calculate the deviations of the doses at each position
    """
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = ttk.Label(self, text="Evaluate STD and MAD", font=LARGE_FONT)
        label.pack(pady=10, padx=10)

        mlabel = tk.Label(self, text=("""
        Evaluate the consistency of measured/reference doses or evaluate the change in the dose distribution at different gantry angles.
        After importing the files of interest relative standard(STD) and absolute deviations(MAD) will be calculated for each chamber.
        The average value of the deviations is also calculated."""),font=SMALL_FONT)
        mlabel.pack()
        frame2 = tk.Frame(self)
        frame2.pack(side="top")
        button_sel_mesh = ttk.Button(frame2, text='Import files',
                                     command=lambda: self.SelectFilesSTD(), width=16)
        button_sel_mesh.pack(side="left", padx=10, pady=10)
        button1 = ttk.Button(frame2, text="Go to Home page",
                             command=lambda: controller.show_frame(HomePage), width=16)
        button1.pack(side="left", pady=10, padx=10)

    def SelectFilesSTD(self):
        """
        Open the files to be processed for:
        dose relative mean and absolute dose deviations
        at each detector position
        """
        files = askopenfilenames(filetypes=(('ASCII', '*.opg'),
                                            ('All files', '*.*')),
                                 title='Import files'
                                 )
        InputFileList = list(self.tk.splitlist(files))
        loadSTD(InputFileList)


class PageImportData(tk.Frame):
    """
    PageImportData
    You can go to Page: HomePage or PageSymmetry
    Impot the files to be investigated
    """
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = ttk.Label(self, text="Import Data", font=LARGE_FONT)
        label.pack(pady=10, padx=10)
        button_back = ttk.Button(self, text='Back',
                                 command=lambda: controller.show_frame(HomePage))
        button_back.pack(side="top", padx=10, pady=10)
        frame2 = tk.Frame(self)
        frame2.pack(side="top")
        button_select_mesh = ttk.Button(frame2, text='Import measurements',
                                        command=lambda: self.selectFiles(), width=20)
        button_select_mesh.pack(side="left", padx=10, pady=10)
        button_select_cal = ttk.Button(frame2, text='Import reference',
                                       command=lambda: self.selectFiles(True), width=20)
        button_select_cal.pack(side="left", padx=10, pady=10)
        button1 = ttk.Button(self, text="Load data",
                             command=lambda: load(self, controller))
        button1.pack(pady=10)
        label2 = ttk.Label(self, text="Import all measured distributions together(not one by one),by pressing ctrl and selecting the files of interest!", font=NORM_FONT)
        label2.pack(pady=10, padx=10)
        label3 = ttk.Label(self, text="Import all reference distributions together(not one by one),by pressing ctrl and selecting the files of interest!", font=NORM_FONT)
        label3.pack(pady=10, padx=10)
        label3 = ttk.Label(self, text="The reference and measured dose distributions should correspond to the the same set of angles!", font=NORM_FONT)
        label3.pack(pady=10, padx=10)

    def selectFiles(self, Flag=False):
        """
        Open the files needed for correction factor calculation
        """
        if Flag:
            txt = "Select reference doses"
        else:
            txt = "Select measured doses"

        files = askopenfilenames(filetypes=(('ASCII', '*.opg'),
                                            ('All files', '*.*')),
                                 title=txt
                                 )
        InputFileList = list(self.tk.splitlist(files))
        InputFileList = sort_ang(InputFileList, Flag)


class PageSymmetry(tk.Frame):
    """
    PageSymmetry
    you can go to :HomePage or PageCalibrateCustom.
    You can evaluate doses angle symmetry or continue to correction factors calculation
    """
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)

        label = ttk.Label(self, text="Page Symmetry", font=LARGE_FONT)
        label.pack(pady=10)

        button1 = ttk.Button(self, text="Back to Home",
                             command=lambda: controller.show_frame(HomePage))
        button1.pack(side="top", pady=10)
        frame2 = tk.Frame(self)
        frame2.pack(side="top")
        button2 = ttk.Button(frame2, text="Evaluate reference dose angle symmetry",
                             command=lambda: ang_sym_chek())
        button2.pack(side="left", pady=10, padx=10)
        button3 = ttk.Button(frame2, text="Evaluate measured dose angle symmetry",
                             command=lambda: ang_sym_chek(True))
        button3.pack(side="left", pady=10, padx=10)
        button4 = ttk.Button(self, text="Calibrate",
                             command=lambda: controller.show_frame(PageCalibrateCustom))
        button4.pack(side="top", pady=10)

        label2 = ttk.Label(self, text=("""
            You can choose to proceed to calibration or first to see whether doses differ from each other for gantry angles:
        0˚-180˚ and 180˚ - 360˚.
            Clicking the buttons  "Evaluate measured dose angle symmetry" or "Evaluate measured dose angle symmetry"
        will display the ratios of dose at position row = i,column = j and angle = theta 
        to dose dose at position row = i, column = number of columns + 1 - j, angle = -theta"""), font=NORM_FONT)
        label2.pack(pady=10)


class PageCalibrateCustom(tk.Frame):
    """
    PageCalibrateCustom
    You can go to :HomePage or PageCalibrationGraphs.
    Chose options for the calibration
    """
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)

        self.corr_at_0 = tk.BooleanVar()

        self.mesh_ang_sym = tk.BooleanVar()
        self.mesh_y_sym = tk.BooleanVar()
        self.mesh_ang_sym_other = tk.BooleanVar()
        self.mesh_y_sym_other = tk.BooleanVar()

        self.cal_ang_sym = tk.BooleanVar()
        self.cal_y_sym = tk.BooleanVar()
        self.cal_ang_sym_other = tk.BooleanVar()
        self.cal_y_sym_other = tk.BooleanVar()

        self.corr_at_0.set(False)

        self.mesh_ang_sym.set(False)
        self.mesh_y_sym.set(False)
        self.mesh_ang_sym_other.set(False)
        self.mesh_y_sym_other.set(False)

        self.cal_ang_sym.set(False)
        self.cal_y_sym.set(False)
        self.cal_ang_sym_other.set(False)
        self.cal_y_sym_other.set(False)

        label = ttk.Label(self, text="Calibration", font=LARGE_FONT)
        label.pack(pady=10, padx=10)
        frame1 = tk.Frame(self)
        frame1.pack()
        button1 = ttk.Button(frame1, text="Back",
                             command=lambda: controller.show_frame(HomePage))
        button1.pack(side="left", padx=10)
        button2 = tk.Radiobutton(frame1, text="Correct at angle zero",
                                  variable=self.corr_at_0, value=True, indicatoron=0)
        button2.pack(side="left", padx=10)
        button3 = tk.Radiobutton(frame1, text="Don't correct at angle zero",
                                  variable=self.corr_at_0, value=False, indicatoron=0)
        button3.pack(side="left", padx=10)
        button4 = ttk.Button(frame1, text="Calibrate",
                             command=lambda: calirbrate(self, controller))
        button4.pack(side="left", padx=10)

        frame2 = tk.Frame(self)
        frame2.config(highlightbackground="white", highlightthickness=1)
        frame2.pack(side="left")

        left_frame_label = ttk.Label(frame2, text="Measured doses", font=LARGE_FONT)
        left_frame_label.pack(pady=5)
        left_button1 = ttk.Radiobutton(frame2, text="Symmetrize on Y axis",
                                       variable=self.mesh_y_sym_other, value=True)
        left_button1.pack(side="bottom", pady=5, anchor="w")
        left_button2 = ttk.Radiobutton(frame2, text="Don't symmetrize on Y axis",
                                       variable=self.mesh_y_sym_other, value=False)
        left_button2.pack(side="bottom", pady=5, anchor="w")
        left_button3 = ttk.Radiobutton(frame2, text="Angle symmetrize",
                                       variable=self.mesh_ang_sym_other, value=True)
        left_button3.pack(side="bottom", pady=5, anchor="w")
        left_button4 = ttk.Radiobutton(frame2, text="Don't angle symmetrize",
                                       variable=self.mesh_ang_sym_other, value=False)
        left_button4.pack(side="bottom", pady=5, anchor="w")
        self.mesh_how_to_calibrate_difrently = tk.StringVar()
        left_1_Entry = ttk.Entry(frame2, textvariable=self.mesh_how_to_calibrate_difrently, width=37)
        left_1_Entry.insert(string="row height:     ,column weight: ", index=0)
        left_1_Entry.pack(side="bottom")
        left_1_Entry.focus_set()
        left_label_1 = ttk.Label(frame2, text="Chose moving average (enter " + "ROWS" + " to average doses of row \n         or " + "COLUMNS" + " to average over column)", font=SMALL_FONT)
        left_label_1.pack(side="bottom", pady=5)
        self.mesh_ang_to_calibrate_difrently = tk.StringVar()
        left_2_Entry = ttk.Entry(frame2, textvariable=self.mesh_ang_to_calibrate_difrently, width=37)
        left_2_Entry.insert(string="Angles: ", index=0)
        left_2_Entry.pack(side="bottom")
        left_2_Entry.focus_set()
        left_label_2 = ttk.Label(frame2, text="Optional(Enter angles to be calibrated differently)", font=SMALL_FONT)
        left_label_2.pack(side="bottom", pady=5)
        left_button11 = ttk.Radiobutton(frame2, text="Symmetrize on Y axis",
                                        variable=self.mesh_y_sym, value=True)
        left_button11.pack(side="bottom", pady=5, anchor="w")
        left_button21 = ttk.Radiobutton(frame2, text="Don't symmetrize on Y axis",
                                        variable=self.mesh_y_sym, value=False)
        left_button21.pack(side="bottom", pady=5, anchor="w")
        left_button31 = ttk.Radiobutton(frame2, text="Angle symmetrize",
                                        variable=self.mesh_ang_sym, value=True)
        left_button31.pack(side="bottom", pady=5, anchor="w")
        left_button41 = ttk.Radiobutton(frame2, text="Don't angle symmetrize",
                                        variable=self.mesh_ang_sym, value=False)
        left_button41.pack(side="bottom", pady=5, anchor="w")
        self.mesh_how_to_calibrate = tk.StringVar()
        left_11_Entry = ttk.Entry(frame2, textvariable=self.mesh_how_to_calibrate, width=37)
        left_11_Entry.insert(string="row height:     ,column weight: ", index=0)
        left_11_Entry.pack(side="bottom")
        left_11_Entry.focus_set()

        left_label = ttk.Label(frame2, text="Chose moving average (enter " + "ROWS" + " to average doses of row \n        or " + "COLUMNS" + " to average over column)", font=SMALL_FONT)
        left_label.pack(side="bottom", pady=5, padx=10)

        frame3 = tk.Frame(self)
        frame3.config(highlightbackground="white", highlightthickness=1)
        frame3.pack(side="right")
        right_frame_label = ttk.Label(frame3, text="Reference doses", font=LARGE_FONT)
        right_frame_label.pack(pady=5)
        right_button1 = ttk.Radiobutton(frame3, text="Symmetrize on Y axis",
                                        variable=self.cal_y_sym_other, value=True)
        right_button1.pack(side="bottom", pady=5, anchor="w")
        right_button2 = ttk.Radiobutton(frame3, text="Don't symmetryzie on Y axis",
                                        variable=self.cal_y_sym_other, value=False)
        right_button2.pack(side="bottom", pady=5, anchor="w")
        right_button3 = ttk.Radiobutton(frame3, text="Angle symmetrize",
                                        variable=self.cal_ang_sym_other, value=True)
        right_button3.pack(side="bottom", pady=5, anchor="w")
        right_button4 = ttk.Radiobutton(frame3, text="Don't angle symmetrize",
                                        variable=self.cal_ang_sym_other, value=False)
        right_button4.pack(side="bottom", pady=5, anchor="w")
        self.cal_how_to_calibrate_difrently = tk.StringVar()
        right_1_Entry = ttk.Entry(frame3, textvariable=self.cal_how_to_calibrate_difrently, width=37)
        right_1_Entry.insert(string="row height:     ,column weight: ", index=0)
        right_1_Entry.pack(side="bottom")
        right_1_Entry.focus_set()

        right_label_1 = ttk.Label(frame3, text="Chose moving average (enter " + "ROWS" + " to average doses of row \n        or " + "COLUMNS" + " to average over column)", font=SMALL_FONT)
        right_label_1.pack(side="bottom", pady=5, anchor="w")
        self.cal_ang_to_calibrate_difrently = tk.StringVar()
        right_2_Entry = ttk.Entry(frame3, textvariable=self.cal_ang_to_calibrate_difrently, width=37)
        right_2_Entry.insert(string="Angles:", index=0)
        right_2_Entry.pack(side="bottom")
        right_2_Entry.focus_set()
        right_label_2 = ttk.Label(frame3, text="Optional(Enter angles to be calibrated differently)", font=SMALL_FONT)
        right_label_2.pack(side="bottom", pady=5)
        right_button11 = ttk.Radiobutton(frame3, text="Symmetrize on Y axis",
                                         variable=self.cal_y_sym, value=True)
        right_button11.pack(side="bottom", pady=5, anchor="w")
        right_button21 = ttk.Radiobutton(frame3, text="Don't symmetrize on Y axis",
                                         variable=self.cal_y_sym, value=False)
        right_button21.pack(side="bottom", pady=5, anchor="w")
        right_button31 = ttk.Radiobutton(frame3, text="Angle symmetrize",
                                         variable=self.cal_ang_sym, value=True)
        right_button31.pack(side="bottom", pady=5, anchor="w")
        right_button41 = ttk.Radiobutton(frame3, text="Don't angle symmetrize",
                                         variable=self.cal_ang_sym, value=False)
        right_button41.pack(side="bottom", pady=5, anchor="w")
        self.cal_how_to_calibrate = tk.StringVar()
        right_11_Entry = ttk.Entry(frame3, textvariable=self.cal_how_to_calibrate, width=37)
        right_11_Entry.insert(string="row height:     ,column weight: ", index=0)
        right_11_Entry.pack(side="bottom")
        right_11_Entry.focus_set()

        left_label = ttk.Label(frame3, text="Chose moving average (enter " + "ROWS" + " to average doses of row \n          or " + "COLUMNS" + " to average over column)", font=SMALL_FONT)
        left_label.pack(side="bottom", pady=5, padx=10)


class PageCompareCalibrations(tk.Frame):
    """
    PageCompareCalibrations
    You can go to Pages:PageCalibrationGraphs
    Chose options for another correction factor sets
    and compere with the set calculated previously
    """
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)

        self.corr_at_0 = tk.BooleanVar()

        self.mesh_ang_sym = tk.BooleanVar()
        self.mesh_y_sym = tk.BooleanVar()
        self.mesh_ang_sym_other = tk.BooleanVar()
        self.mesh_y_sym_other = tk.BooleanVar()

        self.cal_ang_sym = tk.BooleanVar()
        self.cal_y_sym = tk.BooleanVar()
        self.cal_ang_sym_other = tk.BooleanVar()
        self.cal_y_sym_other = tk.BooleanVar()

        self.corr_at_0.set(False)

        self.mesh_ang_sym.set(False)
        self.mesh_y_sym.set(False)
        self.mesh_ang_sym_other.set(False)
        self.mesh_y_sym_other.set(False)

        self.cal_ang_sym.set(False)
        self.cal_y_sym.set(False)
        self.cal_ang_sym_other.set(False)
        self.cal_y_sym_other.set(False)

        label = ttk.Label(self, text="Compare Calibrations", font=LARGE_FONT)
        label.pack(pady=10, padx=10)
        frame1 = tk.Frame(self)
        frame1.pack()
        button1 = ttk.Button(frame1, text="Back",
                             command=lambda: controller.show_frame(PageCalibrationGraphs))
        button1.pack(side="left", padx=10)
        button2 = tk.Radiobutton(frame1, text="Correct at angle zero",
                                  variable=self.corr_at_0, value=True, indicatoron=0)
        button2.pack(side="left", padx=10)
        button3 = tk.Radiobutton(frame1, text="Don't correct at angle zero",
                                  variable=self.corr_at_0, value=False, indicatoron=0)
        button3.pack(side="left", padx=10)
        button4 = ttk.Button(frame1, text="Compare",
                             command=lambda: calirbrate(self, controller, True))
        button4.pack(side="left", padx=10)


        frame2 = tk.Frame(self)
        frame2.config(highlightbackground="white", highlightthickness=1)
        frame2.pack(side="left")

        left_frame_label = ttk.Label(frame2, text="Measured doses", font=LARGE_FONT)
        left_frame_label.pack(pady=5)
        left_button1 = ttk.Radiobutton(frame2, text="Symmetrize on Y axis",
                                       variable=self.mesh_y_sym_other, value=True)
        left_button1.pack(side="bottom", pady=5, anchor="w")
        left_button2 = ttk.Radiobutton(frame2, text="Don't symmetrize on Y axis",
                                       variable=self.mesh_y_sym_other, value=False)
        left_button2.pack(side="bottom", pady=5, anchor="w")
        left_button3 = ttk.Radiobutton(frame2, text="Angle symmetrize",
                                       variable=self.mesh_ang_sym_other, value=True)
        left_button3.pack(side="bottom", pady=5, anchor="w")
        left_button4 = ttk.Radiobutton(frame2, text="Don't angle symmetrize",
                                       variable=self.mesh_ang_sym_other, value=False)
        left_button4.pack(side="bottom", pady=5, anchor="w")
        self.mesh_how_to_calibrate_difrently = tk.StringVar()
        left_1_Entry = ttk.Entry(frame2, textvariable=self.mesh_how_to_calibrate_difrently, width=37)
        left_1_Entry.insert(string="row height:     ,column weight: ", index=0)
        left_1_Entry.pack(side="bottom")
        left_1_Entry.focus_set()
        left_label_1 = ttk.Label(frame2, text="Chose moving average (enter " + "ROWS" + " to average doses of row \n         or " + "COLUMNS" + " to average over column)", font=SMALL_FONT)
        left_label_1.pack(side="bottom", pady=5)
        self.mesh_ang_to_calibrate_difrently = tk.StringVar()
        left_2_Entry = ttk.Entry(frame2, textvariable=self.mesh_ang_to_calibrate_difrently, width=37)
        left_2_Entry.insert(string="Angles: ", index=0)
        left_2_Entry.pack(side="bottom")
        left_2_Entry.focus_set()
        left_label_2 = ttk.Label(frame2, text="Optional(Enter angles to be calibrated differently)", font=SMALL_FONT)
        left_label_2.pack(side="bottom", pady=5)
        left_button11 = ttk.Radiobutton(frame2, text="Symmetrize on Y axis",
                                        variable=self.mesh_y_sym, value=True)
        left_button11.pack(side="bottom", pady=5, anchor="w")
        left_button21 = ttk.Radiobutton(frame2, text="Don't symmetrize on Y axis",
                                        variable=self.mesh_y_sym, value=False)
        left_button21.pack(side="bottom", pady=5, anchor="w")
        left_button31 = ttk.Radiobutton(frame2, text="Angle symmetrize",
                                        variable=self.mesh_ang_sym, value=True)
        left_button31.pack(side="bottom", pady=5, anchor="w")
        left_button41 = ttk.Radiobutton(frame2, text="Don't angle symmetrize",
                                        variable=self.mesh_ang_sym, value=False)
        left_button41.pack(side="bottom", pady=5, anchor="w")
        self.mesh_how_to_calibrate = tk.StringVar()
        left_11_Entry = ttk.Entry(frame2, textvariable=self.mesh_how_to_calibrate, width=37)
        left_11_Entry.insert(string="row height:     ,column weight: ", index=0)
        left_11_Entry.pack(side="bottom")
        left_11_Entry.focus_set()

        left_label = ttk.Label(frame2, text="Chose moving average (enter " + "ROWS" + " to average doses of row \n        or " + "COLUMNS" + " to average over column)", font=SMALL_FONT)
        left_label.pack(side="bottom", pady=5, padx=10)

        frame3 = tk.Frame(self)
        frame3.config(highlightbackground="white", highlightthickness=1)
        frame3.pack(side="right")
        right_frame_label = ttk.Label(frame3, text="Reference doses", font=LARGE_FONT)
        right_frame_label.pack(pady=5)
        right_button1 = ttk.Radiobutton(frame3, text="Symmetrize on Y axis",
                                        variable=self.cal_y_sym_other, value=True)
        right_button1.pack(side="bottom", pady=5, anchor="w")
        right_button2 = ttk.Radiobutton(frame3, text="Don't symmetryize on Y axis",
                                        variable=self.cal_y_sym_other, value=False)
        right_button2.pack(side="bottom", pady=5, anchor="w")
        right_button3 = ttk.Radiobutton(frame3, text="Angle symmetrize",
                                        variable=self.cal_ang_sym_other, value=True)
        right_button3.pack(side="bottom", pady=5, anchor="w")
        right_button4 = ttk.Radiobutton(frame3, text="Don't angle symmetrize",
                                        variable=self.cal_ang_sym_other, value=False)
        right_button4.pack(side="bottom", pady=5, anchor="w")
        self.cal_how_to_calibrate_difrently = tk.StringVar()
        right_1_Entry = ttk.Entry(frame3, textvariable=self.cal_how_to_calibrate_difrently, width=37)
        right_1_Entry.insert(string="row height:     ,column weight: ", index=0)
        right_1_Entry.pack(side="bottom")
        right_1_Entry.focus_set()

        right_label_1 = ttk.Label(frame3, text="Chose moving average (enter " + "ROWS" + " to average doses of row \n        or " + "COLUMNS" + " to average over column)", font=SMALL_FONT)
        right_label_1.pack(side="bottom", pady=5, anchor="w")
        self.cal_ang_to_calibrate_difrently = tk.StringVar()
        right_2_Entry = ttk.Entry(frame3, textvariable=self.cal_ang_to_calibrate_difrently, width=37)
        right_2_Entry.insert(string="Angles:", index=0)
        right_2_Entry.pack(side="bottom")
        right_2_Entry.focus_set()
        right_label_2 = ttk.Label(frame3, text="Optional(Enter angles to be calibrated differently)", font=SMALL_FONT)
        right_label_2.pack(side="bottom", pady=5)
        right_button11 = ttk.Radiobutton(frame3, text="Symmetrize on Y axis",
                                         variable=self.cal_y_sym, value=True)
        right_button11.pack(side="bottom", pady=5, anchor="w")
        right_button21 = ttk.Radiobutton(frame3, text="Don't symmetrize on Y axis",
                                         variable=self.cal_y_sym, value=False)
        right_button21.pack(side="bottom", pady=5, anchor="w")
        right_button31 = ttk.Radiobutton(frame3, text="Angle symmetrize",
                                         variable=self.cal_ang_sym, value=True)
        right_button31.pack(side="bottom", pady=5, anchor="w")
        right_button41 = ttk.Radiobutton(frame3, text="Don't angle symmetrize",
                                         variable=self.cal_ang_sym, value=False)
        right_button41.pack(side="bottom", pady=5, anchor="w")
        self.cal_how_to_calibrate = tk.StringVar()
        right_11_Entry = ttk.Entry(frame3, textvariable=self.cal_how_to_calibrate, width=37)
        right_11_Entry.insert(string="row height:     ,column weight: ", index=0)
        right_11_Entry.pack(side="bottom")
        right_11_Entry.focus_set()

        left_label = ttk.Label(frame3, text="Chose moving average (enter " + "ROWS" + " to average doses of row \n          or " + "COLUMNS" + " to average over column)", font=SMALL_FONT)
        left_label.pack(side="bottom", pady=5, padx=10)


class PageCalibrationGraphs(tk.Frame):
    """
    You can go to: PageCalibrateCustom or PageCompareCalibrations
    Export the correction factors as a csv file or diplay some 
    statistics reqrding the correction factors.
    """
    def __init__(self, parent, controller):
        
        tk.Frame.__init__(self, parent)

        label = ttk.Label(self, text="Calculated correction factors", font=LARGE_FONT)
        label.pack(padx=10)
        frame1 = tk.Frame(self)
        frame1.pack()
        button1 = ttk.Button(frame1, text="Back",
                             command=lambda: controller.show_frame(PageCalibrateCustom))
        button1.pack(side="left", pady=10, padx=10)
        button2 = ttk.Button(frame1, text="Show graph",
                             command=lambda: correction_fact_distribution(self))
        button2.pack(side="left", pady=10, padx=10)
        button3 = ttk.Button(frame1, text="Compare to other calibration",
                             command=lambda: controller.show_frame(PageCompareCalibrations))
        button3.pack(side="left", pady=10, padx=10)

        frame2 = tk.Frame(self)
        frame2.pack()

        self.energy = tk.StringVar()
        Entry_energy = ttk.Entry(frame2, textvariable=self.energy, width=20)
        Entry_energy.insert(string="Beam energy:", index=0)
        Entry_energy.pack(side="left", pady=10, padx=10)
        Entry_energy.focus_set()
        
        self.quality = tk.StringVar()
        Entry_quality = ttk.Entry(frame2, textvariable=self.quality, width=20)
        Entry_quality.insert(string="Beam quality:", index=0)
        Entry_quality.pack(side="left", pady=10, padx=10)
        Entry_quality.focus_set()
        
        label = ttk.Label(frame2, text="Enter LUT name:")
        label.pack(pady=10, padx=10)
        
        self.LUT_name = tk.StringVar()
        Entry_name = ttk.Entry(frame2, textvariable=self.LUT_name, width=20)
        Entry_name.pack(side="left", pady=0, padx=10)

        button4 = ttk.Button(self, text="Export LUT",
                             command=lambda: Export_cfs(self))
        button4.pack(pady=10, padx=10)

if __name__ == "__main__":
    app = CorrGui()
    app.geometry("800x600")
    app.eval('tk::PlaceWindow %s center' % app.winfo_pathname(app.winfo_id()))
    app.mainloop()
