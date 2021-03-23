import numpy as np

from scipy.stats import norm

def gaussian_kernel(kernlen, nsig):
    """Returns a 2D Gaussian kernel."""

    x = np.linspace(-nsig, nsig, kernlen+1)
    kern1d = np.diff(norm.cdf(x))
    kern2d = np.outer(kern1d, kern1d)

    return kern2d / kern2d.sum()
