# Project-3---Simulating-a-Radio-Interferometer

""" Radio Astronomy Project """
""" Project 3: Simulating a Radio Interferometer """

# A project that creates an iPython notebook that simulates the response of an interferometer to an astrophysical signal and uses this data to construct an image of the sky. 



""" Import modules """ 

import numpy as np
from scipy.constants import c, pi
import math
import matplotlib.pyplot as plt
import astropy
import astropy; from astropy import coordinates, wcs 
import astropy; from astropy.coordinates import SkyCoord



""" Plot of the VLA in D-Configuration """

# For this project we are using the EVLA interferometer. The EVLA is made up of 27 radio telescopes each with diameter 25m. For this similation we use the EVLA in its D shape configuaration arranged in a Y-shape.
# To plot the EVLA in its D-configuration I have taken the x,y,z positions outlined below from data from the NRAO Very Large Array. [1]


# VLA x,y,z coordinates given in nanoseconds
           #x      #y        #z
west = [[76.69, 11.67],   #-108.36       #W1
        [49.29, -123.87], #-67.42        #W2
        [96.46, -248.46], #-136.94       #W3
        [156.49, -407.06], #-225.51      #W4
        [228.83, -597.84], #-331.98      #W5
        [311.96, -817.22], #-454.39      #W6
        [405.70, -1064.49], #-592.36     #W7
        [509.53, -1338.54], #-745.23     #W8
        [623.12, -1638.19]] #-912.51     #W9

east = [[151.26, 23.33], #-218.44        #E1
        [37.71, 135.65], #-50.59         #E2
        [73.37, 271.95], #-103.23        #E3
        [118.76, 445.77], #-170.46       #E4
        [173.02, 653.27], #-250.51       #E5
        [235.66, 893.16], #-343.18       #E6
        [305.29, 1163.76], #-448.46      #E7
        [381.68, 1463.33], #-565.35      #E8
        [465.79, 1790.89]] #-692.95      #E9

north = [[2.24, 0.05],      #1.71        #N1
        [-100.24, -15.93], #152.45       #N2
        [-174.91, -27.56], #262.39       #N3
        [-249.59, -39.15], #372.31       #N4
        [-361.68, -56.66], #537.09       #N5
        [-495.22, -77.43], #733.79       #N6
        [-645.82, -100.90], #955.52      #N7
        [-812.58, -126.88], #1200.98     #N8
        [-995.39, -155.53]] #1469.71     #N9

# Converting the VLA x,y,z coordinates to number of wavelengths

W = 5*np.array(west) # convert nanoseconds to seconds (x10e-9) and multiply by frequency 5GHz (prefixes cancel) for conversion
E = 5*np.array(east)
N = 5*np.array(north)

# Slicing the arrays into columns 
 
W_x = W[:,0]
W_y = W[:,1]
#W_z = W[:,2]

E_x = E[:,0]
E_y = E[:,1]
#E_z = E[:,2]

N_x = N[:,0]
N_y = N[:,1]
#N_z = N[:,2]

# Plotting the EVLA x,y,z coordinates to illustrate it in its D-configuration 
 
plt.figure(figsize=(5,5)) 
plt.title('EVLA in its D-configuration') # title 
plt.ylabel('South-North') # axes labelling 
plt.xlabel('West-East')
plt.xticks([]) # removing the scales on the axes 
plt.yticks([]) 
plt.plot(W_y, -W_x, 'kx', E_y, -E_x, 'kx', N_y, -N_x, 'kx')
plt.show()



""" Calculation of the Field of View """

# We can calculate the field of view for one telescope to give the total field of view for the interferometer 
# The equation used for this is θ = 1.02λ/D

D = 25 # diameter of each dish is 25m
freq = 5*1e9 # operating frequency is 5GHz
wave = c/freq # recieving wavelength calculated from given frequency
print('wavelength =', round(wave,3),'m') # print statement check rounded to 3sf 

FoV = ((1.02*wave)*(180/pi))/D # equation calculation and conversion into degrees
print('Field of view= ', round(FoV,3), 'degrees') # print statement check rounded to 3sf 



""" Calculation of offset between sources in degrees """

# To create the map of the field of view we need to calculate the offset between the two sources. We are given the coordinates of the two sources in right ascension RA and declination DEC measured in hours minutes seconds (hms) and degrees arcminutes arcseconds respectively. 

# Right ascension and declination are astronomical coordinates which specify the direction of a point on the celestial sphere in the equatorial coordinate system. [2]

# Looking at the coordinates given we can see that source 2 is offset from source 1 (centre of field of view) by RA 10 seconds and DEC 3 arcminutes, shown below. 

# Source 1:
# Position: J 05 00 00 +45 00 00
# Flux density: 3.6 Jy

# Source 1 should be at the centre of the field-of-view

# Source 2:
# Position: J 05 00 10 +45 03 00
# Flux density: 5.8 Jy
 
# Offset: Right Ascention (RA)  =  10 seconds (05 00 00 -> 05 00 10) 
#            Declination (Dec)  =  3 arcminutes (+45 00 00 -> +45 03 00) 

# Now we can convert the calculated offset into degrees:

# delta RA
RA = 10 * pi / 43200 * (180/pi) # 10 second offset converted to degrees 
print('RA = ', round(RA,2)) # print statement check 

# delta D
Dec = 3 * pi / (60*180) * (180/pi) # 3 arcminute offset converted to degrees 
print('Dec =', round(Dec,2)) # print statement check 

# Flux values of sources
F1 = 3.6 # source 1 flux
F2 = 5.8 # source 2 flux 

# We can write the sources as an n by 3 array for the n sources to be in the sky with the first index being delta RA, second being delta DEC and third being flux

sources = [[0, 0, F1],     # source 1
           [RA, Dec, F2]]  # source 2 

source = np.array(sources) 



""" Primary beam calculation and plot """

# In an interferometer, going further away from the centre of the field of view corresponds to a delay from that which obtains at the centre. This leads to a primary beam correction where the brightness distribution on the sky needs to be modified by the beam shape (primary beam) of the individual antennas. 

# The primary beam descibes the sensitivity of the telescope as a function of direction. As we are dealing with circular dishes then the best way to describe the function is as a sinc function hence a gaussian distribution. [3]

# The primary beam FWHM is the “field-of-view” of a single-pointing interferometric image. 

# FWHM = 2σ √(2 ln⁡2 ) where σ = stand_dev.

# The primary beam is associated with the collecting area, which is just the Fourier Transform of the aperture.

# To calculate the dirty image via FFT, we require operation of order (N^2 log2N) 

FWHM = 256 
stand_dev = FWHM/(2*np.sqrt(2*np.log(2))) # width of gaussian calculation 
print('standard deviation = ', stand_dev)                                   
x = y = np.arange(-128,128,1) # centre the gaussian (start,stop,incr)
xx,yy = np.meshgrid(x,y) # evaluate functions on a grid

z = np.hypot(xx,yy) # image returned in which each pixel contains its own distance from the point specified with np.arange
gauss_image = np.exp(-z*z/(stand_dev**2)) # generates a gaussian of width w

plt.title('Primary Beam') # title
plt.xlabel('Right Ascension (J5000)') # axes labelling 
plt.ylabel('Declination (J5000)')
plt.xticks([]) # removing the scales on the axes 
plt.yticks([])
plt.imshow(np.abs(gauss_image), origin = 'lower', cmap='gray') # plots the gaussian which we can approx as the primary beam  
plt.show()



""" Calculation of sky xy coordinates """

# A 2D array of pixels is created and the 2 sources are plotted on the grid 
# These positions are calculated below and scaled using the field of view calculated above 

# u,v are E-W, N-S spatial frequencies [wavelengths]

xysky = np.zeros((256,256))   # set size of the image to be 256x256 pixels

x = source[:,0]     # slice the array into columns
y = source[:,1]     
z = source[:,2]

for i in range(len(source)):
   
    # Scaling
    x_pixel = int((x[i]*256/FoV)+128) # multiply by 256/FoV to set the scale
    y_pixel = int((y[i]*256/FoV)+128)    # add 128 to centre it
    xysky[y_pixel, x_pixel] = z[i]
   
    print(z[i], x_pixel, y_pixel) # print flux and coordinates of sources then inputted them into the array below
    
    xysky[128,128]=3.6 # source one set flux = 3.6 Jy
    xysky[219,204]=5.8 # source two set flux = 5.8 Jy



""" Plot Field of View """

# Multiplying the primary beam by the 2D array of sources gives a map of the sky in the interferometers field of view

sky = xysky*gauss_image  # multiply the 2D array by the primary beam 

plt.title('Field of View') # title 
plt.xlabel('Right Ascension (J5000)') # axes labelling 
plt.ylabel('Declination (J5000)')
plt.imshow(sky, origin='lower', cmap='gray')
plt.show()

# print positions and flux densities of the sources in field-of-view
print('Source 1 Flux Density =', round(sky[128,128], 2), 'Jy') 
print('Source 2 Flux Density =', round(sky[219,204], 2), 'Jy')



""" Map of Visibilities  """

# The convolution theorem states that the Fourier Transform of a product of two functions is the convolution of the Fourier Transforms of the individual functions. Fourier theory states that any signal or image can be expressed as the sum of sine curves, so by doing a fourier transform on the signal we can put it into its sinusoidal components which still contains all of the information of the original signal. [5]

# So we must convolve the fourier transform of the sampling function (synthesised beam) with the fourier transform of the sky brightness image to give the fourier transform of the observed sky which we can then take the fourier transform of to get the sky as observed by the interferometer.

# We obtain by the synthesised beam by filling in our sampled u,v points with real part = 1 and imaginary part = 0. 

# Then taking a fourier transform of the xy sky with the point sources will give a map of the visibilities. These visibilities will be imaginary and real parts and we can plot both. 

fftsky = np.fft.fft2(sky) # take the fourier transform of the sky

plt.title('Fourier Transform of Field of View (Real)') # title
plt.imshow(np.abs(np.fft.fftshift(fftsky.real)), origin = 'lower', cmap='gray') # plot the real part
plt.xlabel('v [wavelengths]') # axes labelling 
plt.ylabel('u [wavelengths]')
plt.show()

plt.title('Fourier Transform of Field of View (Imaginary)') # title
plt.xlabel('v [wavelengths]') # axes labelling 
plt.ylabel('u [wavelengths])')
plt.imshow(np.abs(np.fft.fftshift(fftsky.imag)), origin = 'lower', cmap='gray') # plot the imaginary part
plt.show()



""" Baseline calculation  """

# The baseline of a telescope is the difference beween two telescope antenna coordinates. For example for the baseline between two antennas we have a vector: B = (Bx, By, Bz) = (x2−x1, y2− y1, z2− z1). [4]

# This vector difference in positions can point in any direction in space. The part of the baseline we want to use in calculating u,v is the phase center direction. This is the component perpendicular to the direction.

# The phase center direction as a unit vector = (HA, DEC), where HA is the hour angle and DEC is the declination. 

# Hour angle and declination are coordinates used in the equatorial coordinate system to give the direction of a point on the celestial sphere. The hour angle of a point is the angle between two planes: one containing Earth's axis and the zenith. Declination is measured north or south of the celestial equator, along the hour circle passing through the point in question. [2]

# The angle may be converted between degrees and time where 24h = 360° and 60 arcminutes = 1 degree. 

# The spatial frequencies u,v are the distances expressed in wavelength units, so we can get the u,v coordinates from the matrix coordinate transformation dotted with the baseline vector. 

# This will give : 

# u = cosH*delta(x) - sinHsinL*delta(y)
# v = sinHsinD*delta(x) + (cosHsinDsinL + cosDcosL)*delta(y)

# Note that u,v depend on the hour angle, so as the Earth rotates and the source appears to move across the sky, the array samples different u,v at different times.  

# An array of n telescopes will have n(n-1)/2 baselines to calculate the baselines between every combination of any two antenna (labelled antA and antB) 

# To do this, iterate a for loop over all antenna in the array

# The reverse baselines also exists as the complex conjugate of each baseline which will be acounted for by plotting (-u, -v) 

all_antenna = west + east + north      # create array with all antenna positions
baselines = []      # create empty array to put in all the baselines
for antA in range(len(all_antenna)):
    for antB in range(len(all_antenna)):
        if antA > antB:
            baselines.append(np.subtract(all_antenna[antA], all_antenna[antB]))   # finds delta x,y (AB) 
            baselines.append(np.subtract(all_antenna[antB], all_antenna[antA]))   # finds delta x,y (BA) 

baselines = np.array(baselines)     # create array with all the calculated baselines in 

print("Number of baselines is ", int(len(baselines)/2)) # check how many baselines have been calculated



""" Matrix calculation for each of the Hour Angles """

# A synthesis imaging radio instrument consists of a number of radio elements which represent measurement points in u,v space.  Now we can convert the array of dishes on the ground to a set of points in u,v space. [6]

# As the Earth rotates, 15 degrees of the sky will be swept by the EVLA. This will be between -12.5 < theta < 12.5 

# Data is taken every 30 seconds which gives 120 hour angles at an increments of 0.125 degrees

H = []

start = -2.5*(pi/180) # convert to radians 
stop  = 12.5*(pi/180)
incre = 0.125*(pi/180)
 
H = np.arange(start, stop, incre) # says to start at -2.5 deg, finish at 12.5 deg in increments of 0.125 deg 

H = np.array(H)

Lat  = (34.07*pi/180) # latitude in radians

# Declination of the source is D = 45 degrees
Dec  = pi/4 # declination in radians

# The first step is to determine a consistent coordinate system.  Antennas are typically measured in units such as meters along the ground.  So we can convert from (x, y) to (u, v) is done through a rotation matrix.

M  = [] # create an empty array for the values of Matrix M

for i in range(len(H)):
    M.append([[np.cos(H[i]), -np.sin(H[i])*np.sin(Lat)],
              [np.sin(H[i])*np.sin(Dec), np.cos(H[i])*np.sin(Dec)*np.sin(Lat)+ np.cos(Dec)*np.cos(Lat)]])
             
M=np.array(M)

UV = []     # create an empty array for the values of u and v 

for i in range(len(M)):
    for j in range(len(baselines)):
       UV.append(np.dot(M[i], baselines[j])) # dot product of the matrix and delta x, delta y
UV = np.array(UV)   # create an array with all the u and v values 


""" Creation of uv Coverage """

# To create a diagram illustrating snapshot u,v coverage we plot the 2D array
# An array of zeros is created and each u,v coordinate is scaled to fit this array
# Each coordinate is represented by a 1 and each empty space is represented by a 0

u = UV[:,0] # slice the array 
v = UV[:,1] 

uv_coverage = np.zeros((256,256))   # set size of the image to be 256x256 pixels

for i in range(len(UV)):
    u_pixel = int(u[i]*(256/8000)+128) # image scaling 
    v_pixel = int(v[i]*(256/8000)+128) 
    uv_coverage[v_pixel, u_pixel] += 1

plt.title('Map of UV Coverage') # title 
plt.xlabel('u [wavelengths]') # axes labelling 
plt.ylabel('v [wavelengths]')
plt.imshow(uv_coverage, origin = 'lower', cmap = 'gray', extent=[-0.05, 0.05,-0.05, 0.05])
plt.show()



""" Fourier Transform of the u,v coverage """

# Calulate the dirty beam by taking the fourier transform of the visibilities

synthesisedbeam = np.fft.fft2(uv_coverage) # fourier transform 

plt.title('Synthesised Beam') # title
plt.xlabel('Righ Ascension Offset (J5000)') # axes labelling 
plt.ylabel('Declination Offset (J5000)')
plt.imshow(np.abs(np.fft.fftshift(synthesisedbeam)), origin='lower', cmap= 'gray', extent=[-0.5, 0.5,-0.5, 0.5])
plt.show()



""" Fourier Transform of Observed Sky """

# To produce the fourier transform of the observed sky we must perform the calculation explained above (uv sampled = fftsky * uv coverage) which multiplys the fourier transform of the sky by the u,v coverage

fftobssky = (fftsky*uv_coverage) # fourier transform 

plt.title('Sampled u,v Plane') # title
plt.xlabel('v [wavelengths]') # axes labelling 
plt.ylabel('u [waveengths]')
plt.imshow(np.abs(fftobssky), origin='lower', cmap='gray', extent=[-0.05, 0.05,-0.05, 0.05])
plt.show()



""" Dirty Map """

# Take the fourier transform of the sampled u,v plane to produce dirty map
# This gives the image observed by the interferometer

obssky = np.fft.fft2(fftobssky) # fourier transform 

plt.xlabel('Right Ascension Offset (J5000)') # axes labelling 
plt.ylabel('Declination Offset (J5000)')
plt.title('"Dirty Map"') # title
plt.imshow(np.abs(obssky), cmap= 'gray', extent=[-0.1, 0.1,-0.1, 0.1] )

# The next step after this would be to CLEAN the image which attempts to deconvolve the dirty beam from the dirty map. [



""" References """ 
# [1] An Introduction to the NRAO Very Large Array, https://science.nrao.edu/facilities/vla/docs/manuals/oss2016A/ant_positions.pdf
# [2] U.S. Naval Observatory, Nautical Almanac Office, P. Kenneth Seidelmann (ed.). Explanatory Supplement to the Astronomical Almanac. University Science Books, Mill Valley, CA. p. 724. ISBN 0-935702-68-7., (1992)
# [3] Jackson N., Principles of Interferometry, University of Manchester, Jodrell Bank Observatory, (2016) 
# [4] https://www.cv.nrao.edu/course/astr534/Interferometers1.html
# [5] Thompson A.R., Moran J.M., Swenson G.W., ISBN 0-471-25492-4, Wiley, (2001) 
# [6] Taylor G.B., Carilli C., Perley R.A. (eds), Astronomical Society of the Pacific Conf. Ser, vol. 180, (1998)
# [7] Garrett M., Basics of radio interferometry & Aperture Synthesis,  https://www.strw.leidenuniv.nl/radioastronomy/lib/exe/fetch.php?media=radio_astronomy_lec_4_ma_garrett.pdf
