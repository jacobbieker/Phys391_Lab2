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
    return apparent_mag - (-2.5 * math.log10(flux))


def find_standard_flux(m1, m2, f2): #Gives the F1 from Eqn vi
    return (10 ** ((-2.5) * (m1 - m2))) / f2


def star_distance(flux, luminosity): #Egn iv solved for distance
    return math.sqrt(luminosity / flux)


def star_luminosity(flux, distance): #Egn iv solved for luminosity
    return flux * distance * distance


def mag_difference(apparent_mag_2, flux_1, flux_2): #Egn vi solved for m1
    return (-2.5) * math.log10(flux_1 / flux_2) + apparent_mag_2


def new_mags(flux, zero_point): #Eqn i
    return (-2.5) * math.log10(flux) + zero_point


def bv(Fb,Fv): #Calculate the BV value
    return Fb / Fv


def quad_error(b_error, v_error):
    return math.sqrt((b_error ** 2) + (v_error ** 2))

def quad_array_error(error_array):
    sum_squared = 0
    for value in error_array:
        sum_squared += (value ** 2)
    return math.sqrt(sum_squared)

def quad_quad_error(b_error, v_error, b_zp_error, v_zp_error):
    return math.sqrt((b_error ** 2) + (v_error ** 2) + (v_zp_error ** 2) + (b_zp_error ** 2))

two_star_data_file = 'two_star_data'
v_mag_data = 'v_mag_data'
b_mag_data = 'b_mag_data'
b_standard_data = 'standard_other_b.txt'
v_standard_data = 'standard_other_v.txt'
standard_data = "standard.txt"
standard_other_mag = 'standard_other_mag.txt'

snr_array = []
with open(two_star_data_file, 'r') as f:
    for line in f:
        lines = line.strip() #removes the return at end
        col = line.split() # splits the line into separate columns
        snr_array.append(singal_to_noise(float(col[4]), math.pi * (float(col[2]) ** 2) * float(col[5]), int(col[6]), 4))

for index, value in enumerate(snr_array):
    if index > 0:
        print (str(index) + " Factor: " + str(value) + " Other: " + str(snr_array[index - 1]) + " This Big: " + str(value/snr_array[index - 1]))

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
        zero_point_file.write(col[0] + " " + str(zero_point(float(col[6]), float(col[2]))) + '\t' + str(zero_point(float(col[7]), float(col[4]))) +  '\n')
zero_point_file.close()

'''
Error calculations in V zero point and B zero point
'''
b_zero_point_error = []
v_zero_point_error = []
with open('b_mag_error.txt', 'r') as error_b:
    for line in error_b:
        line = line.strip()
        col = line.split()
        b_zero_point_error.append(float(col[1]))

b_zero_point_error_avg = 0
for index in b_zero_point_error:
    b_zero_point_error_avg += index

b_zero_point_error_avg /= len(b_zero_point_error)

std_dev_b_zp = 10
#Finds Standard Deviation for B_zp_error
diffsquared_b_zp = 0
sum_diffsquared_b_zp = 0
for val in b_zero_point_error:
    diffsquared_b_zp = (val-b_zero_point_error_avg)**2
    sum_diffsquared_b_zp = diffsquared_b_zp + sum_diffsquared_b_zp
    std_dev_b_zp = (((sum_diffsquared_b_zp)/(len(b_zero_point_error) - 1))**(1/2))

zp_error_b = std_dev_b_zp

with open('v_mag_error.txt', 'r') as error_v:
    for line in error_v:
        line = line.strip()
        col = line.split()
        v_zero_point_error.append(float(col[1]))

v_zero_point_error_avg = 0
for index in v_zero_point_error:
    v_zero_point_error_avg += index

v_zero_point_error_avg /= len(v_zero_point_error)
std_dev_v_zp = 10
#Finds Standard Deviation for B_zp_error
diffsquared_v_zp = 0
sum_diffsquared_v_zp = 0
for val in b_zero_point_error:
    diffsquared_v_zp = (val-v_zero_point_error_avg)**2
    sum_diffsquared_v_zp = diffsquared_v_zp + sum_diffsquared_v_zp
    std_dev_v_zp = (((sum_diffsquared_v_zp)/(len(v_zero_point_error) - 1))**(1/2))

zp_error_v = std_dev_v_zp

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

ordered_b_flux = []
ordered_v_flux = []

print("B-ZP-AVG: " + str(b_zero_point_avg) + " +- " + str(std_dev_b) + " ZP Error: " + str(zp_error_b) + " Total Error: " + str(std_dev_b+std_dev_b_zp) + "\n" + " V-ZP: " + str(v_zero_point_avg) + " +- " + str(std_dev_v) + " ZP Error: " + str(std_dev_v_zp) +" Total Error: " + str(std_dev_v+zp_error_v))
corrected_data_file = open('corrected_data_b.txt', 'w')
with open('b_mag_data') as bdf:
    for lines in bdf:
        line = lines.strip() #removes the return at end
        col = lines.split() #splits the line into separate columns
        ordered_b_flux.append(float(col[5]))
        corrected_data_file.write(str(new_mags(float(col[5]), b_zero_point_avg)) + " +- " + str(zp_error_b) + '\n')
corrected_data_file.close()

corrected_data_file_v = open('corrected_data_v.txt', 'w')
with open('v_mag_data') as vdf:
    for lines in vdf:
        line = lines.strip() #removes the return at end
        col = lines.split() #splits the line into separate columns
        ordered_v_flux.append(float(col[5]))
        corrected_data_file_v.write(str(new_mags(float(col[5]), v_zero_point_avg)) + " +- " + str(zp_error_v) + '\n')
corrected_data_file_v.close()

ordered_v = []
ordered_b = []
with open('corrected_data_v.txt', 'r') as data:
    for lines in data:
        line = lines.strip() #removes the return at end
        col = lines.split() #splits the line into separate columns
        ordered_v.append(float(col[0]))

with open('corrected_data_b.txt', 'r') as data:
    for lines in data:
        line = lines.strip() #removes the return at end
        col = lines.split() #splits the line into separate columns
        ordered_b.append(float(col[0]))

bv_array = []
#Now find B-V before reording V
for index, item in enumerate(ordered_v):
    print (str(index) + " BF " + str(ordered_b_flux[index]) + " VF " + str(ordered_v_flux[index]) + " BV " + str(bv(ordered_b_flux[index], ordered_v_flux[index])))
    bv_array.append(ordered_b[index] - ordered_v[index])

bv_error = []
v_mag_error_array = []
b_mag_error_array = []
#Finds the error in the proper ordering of V and B mags
with open('v_mag_error.txt', 'r') as v_mag_error:
    for lines in v_mag_error:
        line = lines.strip() #removes the return at end
        col = lines.split() #splits the line into separate columns
        v_mag_error_array.append(float(col[1]))

for index, value in enumerate(v_mag_error_array):
    v_mag_error_array[index] = value #+ v_zero_point_error[index]

with open('b_mag_error.txt', 'r') as b_mag_error:
    for lines in b_mag_error:
        line = lines.strip() #removes the return at end
        col = lines.split() #splits the line into separate columns
        b_mag_error_array.append(float(col[1]))

for index, value in enumerate(b_mag_error_array):
    b_mag_error_array[index] = value #+ b_zero_point_error[index]

#actually finding the error
bv_error_file = open('bv_error.txt', 'w')
for index, value in enumerate(v_mag_error_array):
    bv_error.append(quad_quad_error(b_mag_error_array[index], v_mag_error_array[index], std_dev_b_zp, std_dev_v_zp))

for value in bv_error:
    bv_error_file.write(str(value) + '\n')

bv_error_file.close()

pyplot.errorbar(bv_array, ordered_v, xerr=bv_error, yerr=v_mag_error_array, linestyle='None')
pyplot.xlim(-0.5, 2)
pyplot.gca().invert_yaxis()
pyplot.xticks(np.arange(-0.5, 2, 0.1), rotation=270)
pyplot.yticks(np.arange(6, 14, 0.5))
pyplot.xlabel("B-V")
pyplot.ylabel("V")
pyplot.title("M39")
pyplot.show()


def distance_to_cluster(apparent_mag_star, absolute_mag_Sun): #Egn v solved for d
    return 10 ** ((apparent_mag_star - absolute_mag_Sun + 5) / 5)

medium_distance = distance_to_cluster(11.1, 4.8)
closest_distance = distance_to_cluster(11.1 - 0.103783671155, 4.8)
farthest_distance = distance_to_cluster(11.1 + 0.103783671155, 4.8)

print("Med: " + str(medium_distance) + " Small: " + str(closest_distance) + " Large: " + str(farthest_distance))


def percent_error(experimental, theoretical): #Calculates the percentage error in a measurement
    return ((abs(experimental - theoretical)) / theoretical) * 100

#Found distance for M39 is 245.275876 parsecs
medium_distance_error = percent_error(medium_distance, 245.275876)
closest_distance_error = percent_error(closest_distance, 245.275876)
farthest_distance_error = percent_error(farthest_distance, 245.275876)

print("Med: " + str(medium_distance_error) + " Small: " + str(closest_distance_error) + " Large: " + str(farthest_distance_error))

luminosity_ratio_12 = 5910. * (245.275876 ** 2) # Eqn iv solved for L for the brightest man sequence star

print luminosity_ratio_12
masses = math.log(luminosity_ratio_12, 3.5) #Calculates the Mass from the relationship that L = M ^3.5
print ("Solar Masses: " + str(masses))
age = luminosity_ratio_12/masses #Calculates the age based off of the luminosity and the mass
print ("Age: " + str(age))