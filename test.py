from __future__ import division
__author__ = 'Jacob'
import math
luminosity_ratio_12 = 61949 * (245.275876 ** 2) # Eqn iv solved for L for the brightest man sequence star

print luminosity_ratio_12
masses = math.log(luminosity_ratio_12, 3.5) #Calculates the Mass from the relationship that L = M ^3.5
print ("Solar Masses: " + str(masses))
age = luminosity_ratio_12/masses #Calculates the age based off of the luminosity and the mass
print ("Age: " + str(age))
