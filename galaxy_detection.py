import math
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import os


"""
    Year 3 Lab, Cycle 3, A1 Image-Processing Galaxy Detection
    Date: 06/01/2021
    Authors: Sulaiman Rasool, Georgios Alevras
    
    Python Version: 3.8.2
    Dependencies Used: 
        => Math, OS (Included as standard with Python 3.8.2)
        => Astropy: 4.2
        => Numpy: 1.19.1
        => Matplotlib: 3.3.1
"""


def circle(x, y, h, k, r):
    """
        Returns boolean: whether point in square is within a circle of radius = square length
    """
    if (x-h)**2 + (y-k)**2 <= r**2:  # Uses circle equation to find a point
        return True
    else:
        return False


def get_error(counts, pixels, st_error):
    """
        Returns the magnitude error for each galaxy (following normal error propagation) and the fact that counts
        follow a Poisson distribution
    """
    flux_error = np.sqrt(counts + (pixels*st_error)**2)
    mag_error = 2.5*np.log10(flux_error)
    return mag_error


class Galaxy:
    """
        An object (OOP) representing a galaxy.
    """
    def __init__(self, ID, x, y, pixel_count, counts, flux, magnitude, background):
        self.id = ID
        self.x = x
        self.y = y
        self.pixels = pixel_count
        self.counts = counts
        self.magnitude_error = get_error(self.counts, self.pixels, st_error=4.663e-3)
        self.flux = flux
        self.magnitude = magnitude
        self.background = background

    def __repr__(self):
        return 'Galaxy(ID = ' + str(self.id) + ', X-position = ' + str(self.x) + ', Y-position = ' + str(self.y) + \
               ', Pixels = ' + str(self.pixels) + ', Counts = ' + str(self.counts) + ', Count Error = ' + \
               str(self.magnitude_error) + ', Flux = ' + str(self.flux) + ', Magnitude = ' + str(self.magnitude) + \
               ', Background = ' + str(self.background) + ')'


def get_data(FILE_NAME, save=False):
    """
    :param FILE_NAME: filename of image
    :param save: decide whether to save new clean image
    :return: header, data, zero-point
    """
    image = fits.open(FILE_NAME)
    hd = image[0].header
    zp = hd['MAGZPT']
    dt = image[0].data
    dt = dt[::-1]
    dt = fix_image(dt)
    dt = dt[::-1]
    image[0].data = dt
    if save:
        image.writeto('mosaic_clean_image.fits')
    image.close()
    dt = dt[::-1]
    return hd, dt, zp


def fix_image(dt):
    """
        This function removes the original bright objects - more are removed at a later stage.
    """

    """Removal of edge slices of picture, all with pixel count: 3421"""
    dt = dt[99:4496]  # Removing top and bottom slices of data with pixel = 3421
    dt = np.array(dt)  # Removing left and right slices of data with pixel = 3421
    cols_left = [i for i in range(0, 150)]
    cols_right = [i for i in range(2470-150, 2570-150)]
    dt = np.delete(dt, cols_left, axis=1)
    dt = np.delete(dt, cols_right, axis=1)

    """Removing bright objects (total of 6)"""
    square = dt[965:1607, 975:1617]
    m, n = -1, -1
    for i in range(965, 1607 - 1):
        m += 1
        for j in range(975, 1617 - 1):
            n += 1
            b = circle(i, j, 1286, 1296, 320)
            if b:
                square[m][n] = np.random.normal(1000, 15.78, 1)[0]
            else:
                pass
        n = 0
    dt[965:1607, 975:1617] = square

    square_2 = dt[1108:1274, 542:708]
    m, n = -1, -1
    for i in range(542, 708 - 1):
        m += 1
        for j in range(1108, 1274 - 1):
            n += 1
            b = circle(i, j, 625, 1181, 55)
            if b:
                square_2[m][n] = np.random.normal(1000, 15.78, 1)[0]
            else:
                pass
        n = 0
    dt[1108:1274, 542:708] = square_2

    square_3 = dt[1659:1815, 745:901]
    m, n = -1, -1
    for i in range(1659, 1815 - 1):
        m += 1
        for j in range(745, 901 - 1):
            n += 1
            b = circle(i, j, 1737, 823, 50)
            if b:
                square_3[m][n] = np.random.normal(1000, 15.78, 1)[0]
            else:
                pass
        n = 0
    dt[1659:1815, 745:901] = square_3

    square_4 = dt[2150:2400, 680:816]
    m, n = -1, -1
    for i in range(2150, 2400 - 1):
        m += 1
        for j in range(680, 816 - 1):
            n += 1
            b = circle(i, j, 2233, 755, 47)
            if b:
                square_4[m][n] = np.random.normal(1000, 15.78, 1)[0]
            else:
                pass
        n = 0
    dt[2150:2400, 680:816] = square_4

    square_5 = dt[684:820, 1915:2051]
    m, n = -1, -1
    for i in range(684, 820 - 1):
        m += 1
        for j in range(1915, 2051 - 1):
            n += 1
            b = circle(i, j, 752, 1983, 50)
            if b:
                square_5[m][n] = np.random.normal(1000, 15.78, 1)[0]
            else:
                pass
        n = 0
    dt[684:820, 1915:2051] = square_5

    square_6 = dt[1177:1297, 2028:2148]
    m, n = -1, -1
    for i in range(1177, 1297 - 1):
        m += 1
        for j in range(2028, 2148 - 1):
            n += 1
            b = circle(i, j, 1237, 2088, 60)
            if b:
                square_6[m][n] = np.random.normal(1000, 15.78, 1)[0]
            else:
                pass
        n = 0
    dt[1177:1297, 2028:2148] = square_6

    """Removing bleeding of bright objects and other artefacts (total of 10)"""
    bleeding_1 = dt[0:4397, 1270:1298]
    for i in range(len(bleeding_1)):
        for j in range(len(bleeding_1[0])):
            bleeding_1[i][j] = np.random.normal(1000, 15.78, 1)[0]
    dt[0:4397, 1270:1298] = bleeding_1

    bleeding_2 = dt[1092:1311, 620:630]
    for i in range(len(bleeding_2)):
        for j in range(len(bleeding_2[0])):
            bleeding_2[i][j] = np.random.normal(1000, 15.78, 1)[0]
    dt[1092:1311, 620:630] = bleeding_2

    bleeding_3 = dt[1675:1808, 818:827]
    for i in range(len(bleeding_3)):
        for j in range(len(bleeding_3[0])):
            bleeding_3[i][j] = np.random.normal(1000, 15.78, 1)[0]
    dt[1675:1808, 818:827] = bleeding_3

    bleeding_4 = dt[2153:2290, 750:758]
    for i in range(len(bleeding_4)):
        for j in range(len(bleeding_4[0])):
            bleeding_4[i][j] = np.random.normal(1000, 15.78, 1)[0]
    dt[2153:2290, 750:758] = bleeding_4

    bleeding_5 = dt[4029:4087, 952:1502]
    for i in range(len(bleeding_5)):
        for j in range(len(bleeding_5[0])):
            bleeding_5[i][j] = np.random.normal(1000, 15.78, 1)[0]
    dt[4029:4087, 952:1502] = bleeding_5

    bleeding_6 = dt[4112:4200, 863:1557]
    for i in range(len(bleeding_6)):
        for j in range(len(bleeding_6[0])):
            bleeding_6[i][j] = np.random.normal(1000, 15.78, 1)[0]
    dt[4112:4200, 863:1557] = bleeding_6

    bleeding_7 = dt[4240:4294, 1240:1325]
    for i in range(len(bleeding_7)):
        for j in range(len(bleeding_7[0])):
            bleeding_7[i][j] = np.random.normal(1000, 15.78, 1)[0]
    dt[4240:4294, 1240:1325] = bleeding_7

    bleeding_8 = dt[4358:4396, 1142:1388]
    for i in range(len(bleeding_8)):
        for j in range(len(bleeding_8[0])):
            bleeding_8[i][j] = np.random.normal(1000, 15.78, 1)[0]
    dt[4358:4396, 1142:1388] = bleeding_8

    bleeding_9 = dt[0:87, 2010:2320]
    for i in range(len(bleeding_9)):
        for j in range(len(bleeding_9[0])):
            bleeding_9[i][j] = np.random.normal(1000, 15.78, 1)[0]
    dt[0:87, 2010:2320] = bleeding_9

    bleeding_10 = dt[0:311, 2219:2320]
    for i in range(len(bleeding_10)):
        for j in range(len(bleeding_10[0])):
            bleeding_10[i][j] = np.random.normal(1000, 15.78, 1)[0]
    dt[0:311, 2219:2320] = bleeding_10

    return dt


def get_background_hist(dt):
    """
    :param dt: the data of the image
    :return: background count, standard deviation (noise)
    """

    plt.figure(1)  # Background histogram
    data_1d = np.ndarray.flatten(dt)
    data_1d_c = data_1d[data_1d < 3500]
    data_1d_c = data_1d_c[data_1d_c > 3350]
    mean = round(float(np.mean(data_1d_c)), 2)
    stdv = round(float(np.std(data_1d_c)), 2)
    y, x, _ = plt.hist(x=data_1d_c, bins=149, color='m', alpha=0.7, rwidth=0.9, label='Number Count \nμ = ' + str(mean)
                                                                                      + '\nσ = ' + str(stdv))
    x_avg = x[list(y).index(max(y))]
    plt.legend()
    plt.xlabel('Number Count', fontname='Times New Roman', fontsize=12)
    plt.ylabel('Frequency', fontname='Times New Roman', fontsize=12)
    plt.ticklabel_format(scilimits=(0, 5))
    plt.savefig('background_histogram.png')

    plt.figure(2)  # Python heatmap of the clean image
    plt.imshow(dt, cmap='hot', interpolation='nearest')
    plt.xlabel('Width Pixel Number', fontname='Times New Roman', fontsize=11)
    plt.ylabel('Height Pixel Number', fontname='Times New Roman', fontsize=11)
    plt.title('FITS Image Heatmap', fontname='Times New Roman', fontsize=14)
    plt.savefig('mosaic_clean_image_heatmap.png')

    plt.show()

    return x_avg, stdv


def fill_catalogue(dt, cutoff, background, zero_p, store_data, x_h, y_h, save=False):
    """
        Algorithm detects galaxies in clean image, produces a catalogue with galaxies and their attributes.
    """
    catal = []
    dt_w = dt.copy()  # Working data, which will include the masking
    dt_g = dt_w.copy()  # Galaxy data, only includes galaxies

    brightest = np.amax(dt_w)  # Get the next brightest point (ie the next galaxy)
    loop_counter = 0
    galaxy_counter = 0

    #  Keeps iterating (looking for the next galaxy, while the next brightest pixel is above the cutoff, at 1.8 std)
    while brightest > cutoff:
        r_x, r_y = 1, 1  # Initial square of 2x2 pixels around brightest spot
        pos_0, pos_1 = np.where(dt_w == brightest)  # Co-ordinates of brightest spot
        x = pos_0[0]
        y = pos_1[0]
        if r_x + 1 < x < x_h - r_x - 1 and r_y + 1 < y < y_h - r_y - 1:
            r_x += 1
            r_y += 1
            # Boundary conditions: checks if boundary of square has reached background levels
            b_top, b_bottom, b_left, b_right = dt_w[x - r_x][y], dt_w[x + r_x][y], dt_w[x][y - r_y], dt_w[x][y + r_y]
            max_p = np.max([b_top, b_bottom, b_right, b_left])
            avg = np.average([b_top, b_bottom, b_right, b_left])
            if avg > cutoff-35:  # Will expand square if boundary values are above this limit
                while max_p > cutoff:  # Keeps expanding the square (galaxy boundary) until it reaches the background
                    r_x += 1
                    r_y += 1
                    if r_x < x < x_h - r_x and r_y < y < y_h - r_y:
                        b_top, b_bottom, b_left, b_right = dt_w[x - r_x][y], dt_w[x + r_x][y], dt_w[x][y - r_y], \
                                                           dt_w[x][y + r_y]
                        max_p = np.max([b_top, b_bottom, b_right, b_left])
                    else:  # Break the loop (galaxy boundaries found) once background levels reached
                        r_x -= 1
                        r_y -= 1
                        break
                square = dt_w[x - r_x:x + r_x, y - r_y:y + r_y]  # Part of image that is the galaxy
                m, n = -1, -1
                pixels = 0
                counts = 0
                #  Loop creates a circle inscribing square boundaries
                for i in range(x - r_x, x + r_x):
                    m += 1
                    for j in range(y - r_y, y + r_y):
                        n += 1
                        if circle(i, j, x, y, r_x):
                            if square[m][n] > 3448:
                                counts += square[m][n]  # Measures count values of galaxy
                                pixels += 1  # Measures the total number of pixels in galaxy
                            square[m][n] = 2000  # Masks area of galaxy with a low pixel value (below cut-off)
                    n = -1
                # Get the local background for each galaxy
                stfu_b = 2*r_x  # square with boundary increased
                if r_x + stfu_b < x < x_h - r_x - stfu_b and r_y + stfu_b < y < y_h - r_y - stfu_b:
                    square_b = dt_w[x - r_x-stfu_b:x + r_x+stfu_b, y - r_y-stfu_b:y + r_y+stfu_b]
                    counts_b = square_b[square_b > 2000]
                    counts_b = counts_b[counts_b < cutoff - 100]
                    bg = np.average(counts_b)  # Obtains local background
                    if not math.isnan(bg):
                        background = bg  # Assign background as local background
                    else:
                        pass
                if counts > 0 and pixels > 0:
                    flux = counts - pixels * background  # Gets the flux with local background adjusted
                    magnitude = zero_p - 2.5 * np.log10(flux)  # Determines magnitude given image's zero-point
                    # Appends galaxy with specific attributes to catalogue
                    catal.append(Galaxy(galaxy_counter, x, y, pixels, counts, flux, magnitude, background))
                    galaxy_counter += 1
                loop_counter += 1
                dt_g[x - r_x:x + r_x, y - r_y:y + r_y] = square
                stfu = 3  # shift the f-value up
                if r_x + stfu < x < x_h - r_x - stfu and r_y + stfu < y < y_h - r_y - stfu:
                    dt_w[x - r_x - stfu:x + r_x + stfu, y - r_y - stfu:y + r_y + stfu] = 2000
                brightest = np.amax(dt_w)  # Gets the next brightest spot of image
            else:  # When square boundaries were below cutoff value, make a small galaxy in-place
                # Trace bright adjacent points to locate dim galaxies
                square = dt_w[x - 3:x + 3, y - 3:y + 3]
                m, n = -1, -1
                pixels = 0
                counts = 0
                for i in range(x - 2, x + 2):
                    m += 1
                    for j in range(y - 2, y + 2):
                        n += 1
                        if circle(i, j, x, y, 2):
                            if square[m][n] > 3448:
                                counts += square[m][n]
                                pixels += 1
                            square[m][n] = 2000
                    n = -1
                if pixels == 0:
                    pixels = 1
                    counts = 3450
                stfu = 2  # shifts the f-value up
                if r_x + stfu < x < x_h - r_x - stfu and r_y + stfu < y < y_h - r_y - stfu:
                    dt_w[x - r_x - stfu:x + r_x + stfu, y - r_y - stfu:y + r_y + stfu] = 2000
                if counts/pixels > cutoff - 25:
                    flux = counts - pixels * background  # Gets the flux with local background adjusted
                    magnitude = zero_p - 2.5 * np.log10(flux)  # Determines magnitude given image's zero-point
                    # Appends galaxy with specific attributes to catalogue
                    catal.append(Galaxy(galaxy_counter, x, y, pixels, counts, flux, magnitude, background))
                    galaxy_counter += 1
                    loop_counter += 1
                    dt_w[x - 3:x + 3, y - 3:y + 3] = square
                    brightest = np.amax(dt_w)
                else:
                    dt_w[x:x + 1, y:y + 1] = 0
                    brightest = np.amax(dt_w)
                    loop_counter += 1
        else:
            dt_w[x:x + 1, y:y + 1] = 0
            brightest = np.amax(dt_w)
            loop_counter += 1
        print('Galaxies added: ', galaxy_counter)
        if galaxy_counter == 10000:  # Break algorithm at specified total number of galaxies (if wanted)
            break

    if save:
        # Saves the catalogue in a Numpy file, to not require re-running the script each time
        array_to_save = np.array(catal)
        file = store_data
        np.save(os.path.join('Numpy Files', file), array_to_save)

    image = fits.open('mosaic.fits')
    dt_w = dt_w[::-1]
    image[0].data = dt_w
    image.writeto('catalogued_image_masking.fits')  # Save working image with masking
    image.close()
    image = fits.open('mosaic.fits')
    dt_g = dt_g[::-1]
    image[0].data = dt_g
    image.writeto('catalogued_image_galaxies.fits')  # Save working image with galaxies
    image.close()
    plt.imshow(dt_w, cmap='Blues', interpolation='nearest')
    plt.xlabel('Width Pixel Number', fontname='Times New Roman', fontsize=11)
    plt.ylabel('Height Pixel Number', fontname='Times New Roman', fontsize=11)
    plt.title('FITS Image Heatmap', fontname='Times New Roman', fontsize=14)
    plt.show()
    return catal


if __name__ == '__main__':
    header, data, zero_point = get_data('mosaic.fits', save=False)  # Obtain header, data and zero-point form image

    avg_background, std = get_background_hist(data)  # Get background statistics: mean and standard deviation
    limit = avg_background + 1.8*std  # Define the cut-off for a source as 1.8 standard deviations from the background

    file_to_save = 'catalogue_data.npy'  # Numpy file to save catalogue
    """ The following block of code is used to test code and validate it on small sub-section of the image """
    # x_lower, x_higher, y_lower, y_higher = 1650, 1950, 400, 650
    # data = data[y_lower:y_higher, x_lower:x_higher]

    """ Running the code on the entire image once testing and validation are complete """
    catalogue = fill_catalogue(data, limit, avg_background, zero_point, file_to_save, 4397, 2320, save=True)

    """ Plotting of results (number count plot) """
    # catalogue = np.load(file_to_save, allow_pickle=True)
    # galaxy_mags = np.array([g.magnitude for g in catalogue][50:])  # Filters first 50 bright objects
    # magnitudes = np.linspace(min(galaxy_mags), (max(galaxy_mags)), 30)
    # mag_errors = [g.magnitude_error for g in catalogue][50:]
    # number = [np.log10(np.count_nonzero(galaxy_mags < m)) for m in magnitudes]  # Obtain number count (N(<m))
    # count_errors = [np.log10(np.sqrt(n)) for n in number]  # Count errors using Poisson Distribution
    # plt.plot(magnitudes, number, 'x', label='Data Obtained')
    # plt.errorbar(magnitudes, number, yerr=count_errors, linestyle='')
    # theoretical_x = [11, 12, 13, 14, 15]
    # theoretical_y = [0.6*x - 5.85 for x in theoretical_x]
    # plt.plot(theoretical_x, theoretical_y, label='Theoretical Data (gradient = 0.6)')
    # h_assymptote = [1.01*number[-1], number[-1]*1.01]
    # plt.plot([10, 21], h_assymptote, '--', label='Total number of galaxies: ' + str(round(10**number[-1], 0)))
    # plt.xlabel('Magnitude [m]', fontname='Times New Roman', fontsize=12)
    # plt.ylabel('Log$_{10}$ of Number Count [log$_{10}$(N(m))]', fontname='Times New Roman', fontsize=12)
    # plt.legend()
    # plt.xlim(10, 21)
    # plt.ticklabel_format(scilimits=(0, 5))
    # plt.savefig('count_plot.png')
    # plt.show()
