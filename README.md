# A1-Image-Processing-Algorithm
Year 3 Lab, Cycle 3, A1 Image-Processing Galaxy Detection
Date: 06/01/2021
Authors: Sulaiman Rasool, Georgios Alevras
----------------------------------------------------------

Python Version: 3.8.2
Dependencies Used: 
    => Math, OS (Included as standard with Python 3.8.2)
    => Astropy: 4.2
    => Numpy: 1.19.1
    => Matplotlib: 3.3.1
==========================================


Instructions to run the A1 Astronomical image processing / galaxy detection code, and obtain results required by the lab-script.


-----------------------------------------------------------------------------------------------------------------------------------
The .zip file initially contains the following items:

1) 'galaxy_detection.py' - python file:
	This contains all the code we developed for the lab, inlculuding all image processing and galaxy detection. 
	All the lab results are obtained using this piece of code.

2) 'mosaic.fits':
	This is the .fits file that has the original astronomical data before any processing is occured.

3) 'mosaic_clean_image.fits':
	This is an altered verion of the .fits file mentioned above, with all the bleeding, bright star and anomalous 
	data removed. This is the working image used by the galaxy detection algorithm.

4) 'Numpy Files' folder:
	Here the code will save a numpy file (.npy) which contains all the catalogued data of the galaxies that have been detected
	This file can then be moved to the working directory disscused below so that analysis and plots can be made on the data. 
-----------------------------------------------------------------------------------------------------------------------------------


******* How to run the data **********

In order for the code to run, .fits and .npy files need to be in the same working directory as the python file.
Initially before the detection algorithm runs, two plots are printed, press 'X' to close these in order for the algorithm to commence.

The script serves 2 purposes: 
	a) it either runs the detection algorithm to find a given number of galaxies and saves this catalogue as a .npy
		file, 
	b) it loads the catalogued data produced previously and plots the number count.


1) To run the detection algorithm:
	==> Keep lines: 437 - 456 commented, and specify at line: 395 the galaxy_counter to your desired value (how many galaxies to find)
2) To observe the data:
	==> Move .npy file produced from the 'Numpy Files' folder into the working directory (where the python file is)
	==> Comment out line 434 and uncomment lines: 437 - 456


N.B.: When running the script, 5 files will be produced:
	1) catalogued image with galaxies (catalogued_image_galaxies.fits)
	2) catalogued image with masking for galaxies (catalogued_image_masking.fits)
	3) numpy file (catalogue_data.npy inside numpy folder)
	4) background_histogram.png (histogram of the noise)
	5) mosaic_clean_image_heatmap.png (heatmap)
To run the code again; ensure these 5 files are deleted, so that new can be produced. 
