from __future__ import division
__author__ = 'Jacob'
import math
import numpy as np
import matplotlib.pyplot as pyplot


def singal_to_noise(flux, s, exposure_time, number_of_images=1, gain=10):
    return flux * exposure_time * math.sqrt(gain / (((flux / number_of_images) + s) * exposure_time)) #The SNR formula


def mag_error_calc(snr_value): # Egn iii
    return 1.0875 / snr_value


def mag_calc(flux, zeropoint): #Egn i
    return -2.5 * math.log10(flux) + zeropoint


def abs_mag(distnace, apparent_mag): #Uses equation v, sovled for Absolute Mag
    return (-5) * math.log10(distnace) + 5 + apparent_mag


def zero_point(apparent_mag, flux): #Egn i solved for zeropoint
    return apparent_mag + 2.5 * math.log10(flux)


def find_standard_flux(m1, m2, f2): #Gives the F1 from Eqn vi
    return (10 ** (-2.5 * (m1 - m2))) / f2


def star_distance(flux, luminosity): #Egn iv solved for distance
    return math.sqrt(luminosity / flux)


def star_luminosity(flux, distance): #Egn iv solved for luminosity
    return flux * distance * distance


def mag_difference(apparent_mag_2, flux_1, flux_2): #Egn vi solved for m1
    return (-2.5) * math.log10(flux_1 / flux_2) + apparent_mag_2


def new_mags(flux, zero_point): #Eqn i
    return (-2.5) * math.log10(flux) + zero_point


two_star_data_file = 'two_star_data'
v_mag_data = 'v_mag_data'
b_mag_data = 'b_mag_data'
b_standard_data = 'standard_other_b.txt'
v_standard_data = 'standard_other_v.txt'
standard_data = "standard.txt"
standard_other_mag = 'standard_other_mag.txt'

with open(two_star_data_file, 'r') as f:
    for line in f:
        lines = line.strip() #removes the return at end
        col = line.split() # splits the line into separate columns
        print singal_to_noise(float(col[4]), math.pi * (float(col[2]) ** 2) * float(col[5]), int(col[6]), 4)


of_v = open('v_mag_error.txt', 'w')
with open(v_mag_data, 'r') as f_v:
    for line in f_v:
        l = line.strip() #removes the return at end
        col = line.split() # splits the line into separate columns
        SNR = singal_to_noise(float(col[5]), float(col[3]) * float(col[6]), 1.0/60.0)
        of_v.write(str(SNR) + '\t' + str(mag_error_calc(SNR)) + '\n')
of_v.close()

of_b = open('b_mag_error.txt', 'w')
with open(b_mag_data, 'r') as f_b:
    for line in f_b:
        li = line.strip() #removes the return at end
        col = line.split() # splits the line into separate columns
        SNR = singal_to_noise(float(col[5]), float(col[3]) * float(col[6]), 1.0/60.0)
        of_b.write(str(SNR) + '\t' + str(mag_error_calc(SNR)) + '\n')
of_b.close()

zero_point_file = open('zero_point.txt', 'w')
with open(standard_other_mag, 'r') as standard_b:
    standard_b.readline()
    for line in standard_b:
        lin = line.strip() #removes the return at end
        col = line.split() #splits the line into separate columns
        b_flux = find_standard_flux(float(col[7]), float(col[1]), float(col[2]))
        v_flux = find_standard_flux(float(col[6]), float(col[3]), float(col[4]))
        zero_point_file.write(col[0] + " " + str(zero_point(b_flux, float(col[2]))) + '\t' + str(zero_point(v_flux, float(col[4]))) +  '\n')
zero_point_file.close()

b_zero_point_avg = 0
v_zero_point_avg = 0
b_zero_points = []
v_zero_points = []
file_one = open('zero_point.txt', 'r')

with file_one as zp:
    for lines in zp:
        line = lines.strip()
        col = lines.split()
        print col[1]
        print col[2]
        b_zero_points.append(float(col[1]))
        v_zero_points.append(float(col[2]))

file_one.close()

for index in v_zero_points:
    v_zero_point_avg += index

for index in b_zero_points:
    b_zero_point_avg += index

b_zero_point_avg /= len(b_zero_points)
v_zero_point_avg /= len(v_zero_points)

std_dev_b = 0
std_dev_v = 0

#Finds Standard Deviation for B
diffsquared_b = 0
sum_diffsquared_b = 0
for val in b_zero_points:
    diffsquared_b = (val-b_zero_point_avg)**2
    sum_diffsquared_b = diffsquared_b + sum_diffsquared_b
    std_dev_b = (((sum_diffsquared_b)/(len(b_zero_points) - 1))**(1/2))


#Finds Standard Deviation for V
diffsquared_v = 0
sum_diffsquared_v = 0
for val in v_zero_points:
    diffsquared_v = (val-v_zero_point_avg)**2
    sum_diffsquared_v = diffsquared_v + sum_diffsquared_v
    std_dev_v = (((sum_diffsquared_v)/(len(v_zero_points) - 1))**(1/2))

corrected_data_file = open('corrected_data_b.txt', 'w')
with open('b_mag_data') as bdf:
    for lines in bdf:
        line = lines.strip() #removes the return at end
        col = lines.split() #splits the line into separate columns
        corrected_data_file.write(str(new_mags(float(col[5]), b_zero_point_avg)) + '\n')
corrected_data_file.close()

corrected_data_file_v = open('corrected_data_v.txt', 'w')
with open('v_mag_data') as vdf:
    for lines in vdf:
        line = lines.strip() #removes the return at end
        col = lines.split() #splits the line into separate columns
        corrected_data_file_v.write(str(new_mags(float(col[5]), v_zero_point_avg)) + '\n')
corrected_data_file_v.close()

ordered_v = []

with open('corrected_data_v.txt', 'r') as data:
    for lines in data:
        line = lines.strip() #removes the return at end
        col = lines.split() #splits the line into separate columns
        ordered_v.append(float(col[0]))

ordered_v = sorted(ordered_v, reverse=True) #proer ordering for V

pyplot.errorbar(ordered_v, ordered_v, xerr=std_dev_b, yerr=std_dev_v)
pyplot.show()

print b_zero_point_avg
print v_zero_point_avg
print std_dev_v
print std_dev_b