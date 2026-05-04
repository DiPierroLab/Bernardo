import math
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import poisson
from sklearn.metrics import r2_score
from scipy.special import factorial


def bp_notation(num):
    units = ["bp", "Kbp", "Mbp"]
    unit = units[int(np.log10(num) // 3)]
    val = np.power(10, np.log10(num) % 3 - np.log10(num)) * num
    return val, unit


def exp(x, a, b):
    return a * np.exp(-b * x)


def site_distribution(ax, instances, bincount=200):
    """Creates a histogram restriction site frequency along a chromosome"""
    counts, bins = np.histogram(instances, bins=bincount)
    ax.plot(bins[1:], counts)
    ax.set_xlabel("Position (bp)")
    ax.set_ylabel("Site frequency (binned)")
    ax.set_title("Site positions")


def avg_distance_between_sites(ax, instances, bincount=200):
    """Creates a plot showing the average distance between sites at locations along a chromosome"""
    counts, bins = np.histogram(instances, bins=bincount)
    # Find bin size
    binsize = bins[1] - bins[0]
    # Remove bins with zero instances
    bins = bins[1:][counts != 0]
    counts = counts[counts != 0]
    ax.plot(bins, binsize / counts)
    ax.set_xlabel("Position (bp)")
    ax.set_ylabel("Average site distance (binned)")
    ax.set_title("Site distances")


def distances_scatter(ax, distances, xscale="linear", yscale="linear", bincount=300):
    """Creates a histogram restriction site frequency along a chromosome"""
    counts, bins = np.histogram(distances, bins=bincount)
    print(bins)
    ax.scatter(bins[:-1], counts, s=5, color="lightcoral")
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    ax.set_xlabel("Distance between sites (bp)")
    ax.set_ylabel("Frequency (binned)")
    ax.set_title("Distance distribution")


# Find an autocorrelation function that doesn't take forever but that says what formula it uses.
def autocorrelation(ax, distances):
    """Plots autocorrelation """
    acorr = np.correlate(distances, distances, 'same')
    ax.plot(range(len(acorr)), acorr)
    ax.set_title("ACF")
    ax.set_xlabel("Lag")
    ax.set_ylabel("ACF Coefficient")


def fit_distances(ax, func, distances, bincount=300):
    """Plots of a fit of the distances distribution and returns parameters"""
    counts, bins = np.histogram(distances, bins=bincount)
    x = bins[:-1]
    popt, pcov = curve_fit(func, x, counts, maxfev=6000000)
    y = func(x, *popt)
    ax.plot(x, y, label=func.__name__ + " fit", color="navy", linestyle="dashed")
    ax.set_ylim([1, max(counts)])

    r_squared = r2_score(counts, func(x, *popt))
    return popt, pcov, r_squared


def cdf_distances_distribution(ax, distances, bincount=300):
    counts, bins = np.histogram(distances, bins=bincount)
    # To fit to a cdf, find the integral of the distances distribution
    # Use that as inputs to fit A(1-e^(lx))
    integrated = [None] * len(counts)
    integrated[0] = counts[0]
    for k in range(1, len(counts)):
        integrated[k] = integrated[k-1] + counts[k]

    ax.scatter(bins[:-1], integrated, s=5, color="lightcoral")

    def cdf(x, A, l):
        return A/l * (1 - np.exp(-l * x))

    popt, pcov = curve_fit(cdf, bins[:-1], integrated, maxfev=600000)
    x = bins[:-1]
    y = cdf(x, *popt)
    ax.plot(x, y, color="blue", linestyle="dashed", label="CDF Fit")
    ax.set_ylabel("Sum of instances within distance")
    ax.set_xlabel("Distance between sites (bp)")
    ax.set_title("CDF fit")
    return popt, pcov


def shifted_returns(ax, distances, lag=1):
    """Creates a shifted return scatter"""
    if lag == 0:
        ax.scatter(distances, distances, s=5, color="lightcoral")
    else:
        ax.scatter(distances[:-lag], distances[lag:], s=5, color="lightcoral")
    ax.set_xlabel("Shifted up " + str(lag))
    ax.set_ylabel("Shifted back " + str(lag))


def shifted_return_density(ax, distances, lag=1, bincount=300):
    """Plots the density of points in a shifted return scatter as a function of distance from the origin"""
    # Get the distance from the origin of each point in the shift. Bin these distances and plot them.
    if lag == 0:
        shifted_distances = np.sqrt(np.square(distances) + np.square(distances))
    else:
        shifted_distances = np.sqrt(np.square(distances[lag:]) + np.square(distances[:-lag]))

    counts, bins = np.histogram(shifted_distances, bincount)
    ax.plot(bins[1:], counts)
    ax.set_xlabel("Shift dist from origin")
    ax.set_ylabel("Density")


def poisson(x, avg):
    return avg ** x * np.exp(-avg) / factorial(x)


def poisson_distribution(ax, instances, regioncount=300):
    """Plots the poisson distribution of sites per region of a set length"""
    region_counts, region_bins = np.histogram(instances, bins=regioncount)
    print(region_counts)

    counts, bins = np.histogram(region_counts, bins=range(min(region_counts), max(region_counts)))
    # print(counts)
    # Normalize the data
    counts = np.divide(counts, sum(counts))
    ax.plot(bins[:-1], counts)

    # Plot a poisson fit
    popt, pcov = curve_fit(poisson, bins[:-1], counts)
    ax.plot(bins[:-1], poisson(bins[:-1], *popt), linestyle="dashed", color="navy", label="Poisson Fit")

    val, units = bp_notation(round(region_bins[1])-round(region_bins[0]))
    ax.set_ylabel("Frequency")
    ax.set_xlabel("Sites per " + str(val) + units)
    ax.set_xlim([0, max(region_counts)])
    return popt, pcov, r2_score(counts, poisson(bins[:-1], *popt))


def ideal_poisson(ax, restriction_length, region_length, xdata=np.array(range(100))):
    """Plots the ideal poisson distribution of sites per region of a set length"""
    # There is a 4^-n probability that a given sequence of length n will be a match for another sequence of length n,
    # and there are 2 * regionlength possible sequences of length n (read in both directions).
    avg_sites_per_segment = (4 ** -restriction_length) * 2 * region_length
    ax.plot(xdata, poisson(xdata, avg_sites_per_segment), linestyle="dotted", color="lightblue", label="Ideal Poisson")