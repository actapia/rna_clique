import functools

import numpy as np
import scipy.stats, scipy.linalg
import matplotlib as mpl

from matplotlib import pyplot as plt
from numpy.typing import ArrayLike

from collections import namedtuple
from numbers import Real

# @dataclass
# class Ellipsoid:
#     center: np.ndarray # p vector
#     axes: np.ndarray   # p x p

# See Johnson and Whichern's Applied Multivariate Statistical Analysis, Chapters
# 1--5 for insights into the math used here.

Ellipsoid = namedtuple("Ellipsoid", ["center", "axes"])

def get_stat_distance_ellipsoid(data: ArrayLike, dist: Real) -> Ellipsoid:
    """Get the ellipsoid representing points at a given statistical distance."""
    eigs = np.linalg.eig(np.cov(data.T))
    return Ellipsoid(
        np.mean(data, axis=0),
        dist * np.sqrt(eigs.eigenvalues) * eigs.eigenvectors
    )
    
def get_multivariate_normal_density_ellipsoid(
        data: ArrayLike,
        density: Real = 0.95
) -> Ellipsoid:
    """Get the ellipsoid containing given probability density at the data mean.

    This assumes that the data are multivariate normal.

    Parameters:
        data:            A n x p array, the data matrix.
        density (float): Probability density to be contained in the ellipsoid.

    Returns:
        An ellipsoid containg the given probability density at the data mean.
    """
    return get_stat_distance_ellipsoid(
        data,
        np.sqrt(scipy.stats.chi2.ppf(density, data.shape[1]))
    )

def get_multivariate_normal_confidence_ellipsoid(
        data: ArrayLike,
        conf: Real = 0.95,
        use_chi2: bool = False
) -> Ellipsoid:
    """Get a confidence ellipsoid for the population mean from the data.

    This assumes that the population distribution is multivariate normal.

    For large sample sizes, the chi-squared distribution can be used instead of
    the F-distribution. In general (and especially for small sample sizes), the
    F-distribution should be used.

    Parameters:
        data:            A n x p array, the data matrix.
        conf (float):    Confidence level for the ellipsoid.
        use_chi2 (bool): Whether to use chi-squared instead of F distribution.

    Returns:
        A confidence ellipsoid for the population mean with given confidence.
    """
    if use_chi2:
        dist2 = scipy.stats.chi2.ppf(conf, data.shape[1])
    else:
        dist2 = (
            (
                data.shape[1] * (data.shape[0] - 1)
            )/(data.shape[0] - data.shape[1])
        ) * scipy.stats.f.ppf(
            conf,
            data.shape[1],
            data.shape[0] - data.shape[1]
        )
    dist2 /= data.shape[0]
    return get_stat_distance_ellipsoid(data, np.sqrt(dist2))

get_confidence_ellipsoid = get_multivariate_normal_confidence_ellipsoid

def conf_ellipsoid(conf):
    return functools.partial(get_confidence_ellipsoid, conf=conf)

def draw_ellipse(
        ellipse: Ellipsoid,
        axes: bool= False,
        **kwargs
):
    """Draw an ellipse (2D Ellipsoid).

    Parameters:
        ellipse:       The ellipse to draw.
        axes (bool):   Whether to draw the axes of the ellipse.
    """
    if axes:
        plt.quiver(
            *np.broadcast_to(
                ellipse.center,
                (ellipse.center.shape[0], ellipse.center.shape[0])
            ).T,
            ellipse.axes[0],
            ellipse.axes[1],
            angles='xy',
            scale_units='xy',
            scale=1
        )
    ellipse_patch = mpl.patches.Ellipse(
        tuple(ellipse.center),
        width=2*np.linalg.norm(ellipse.axes[:,0]),
        height=2*np.linalg.norm(ellipse.axes[:,1]),
        angle=np.rad2deg(np.arctan2(*reversed(ellipse.axes[:,0]))),
        **kwargs
    )
    plt.gca().add_artist(ellipse_patch)
