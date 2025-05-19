import functools

import numpy as np
import scipy.stats, scipy.linalg
import matplotlib as mpl

from matplotlib import pyplot as plt
from numpy.typing import ArrayLike

from collections import namedtuple
from typing import Optional
from numbers import Real

# @dataclass
# class Ellipsoid:
#     center: np.ndarray # p vector
#     axes: np.ndarray   # p x p

# See Johnson and Whichern's Applied Multivariate Statistical Analysis, Chapters
# 1--5 for insights into the math used here.

Ellipsoid = namedtuple("Ellipsoid", ["center", "axes"])

@functools.cache
def uv_sphere(seg: int = 100) -> np.ndarray:
    """Returns a grid of points for a unit UV sphere.

    Neighboring points in the array are connected.

    Parameters:
        seg (int): Number of points in each UV dimension.

    Returns:
        A grid of seg**2 points forming a unit UV sphere.
    """
    # Based on https://matplotlib.org/stable/gallery/mplot3d/surface3d_2.html
    u = np.linspace(0, 2 * np.pi, seg)
    v = np.linspace(0, np.pi, seg)
    return np.array(
        [
            np.outer(np.cos(u), np.sin(v)),
            np.outer(np.sin(u), np.sin(v)),
            np.outer(np.ones(np.size(u)), np.cos(v)),
        ]
    )

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
        ax: Optional[mpl.axes.Axes] = None,
        **kwargs
):
    """Draw an ellipse (2D Ellipsoid).

    Parameters:
        ellipse:       The ellipse to draw.
        axes (bool):   Whether to draw the axes of the ellipse.
        ax:            Matplotlib Axes on which to draw the ellipse.
    """
    if ax is None:
        ax = plt.gca()
    if axes:
        ax.quiver(
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
    ax.add_artist(ellipse_patch)

def draw_3d_ellipsoid(
        ellipsoid: Ellipsoid,
        axes: bool = False,
        ax: Optional[mpl.axes.Axes] = None,
        seg: int = 100,
        **kwargs
):
    """Draw a 3D ellipsoid.

    Parameters:
        ellipsoid:   The ellipsoid to draw.
        axes (bool): Whether to draw the axes of the ellipsoid.
        ax:          Matplotlib Axes on which to draw the ellipsoid.
    """
    if ax is None:
        ax = plt.gca()
    sph = uv_sphere(seg)
    if axes:
        ax.quiver(
            *np.broadcast_to(
                ellipsoid.center,
                (ellipsoid.center.shape[0], ellipsoid.center.shape[0])
            ).T,
            ellipsoid.axes[0],
            ellipsoid.axes[1],
            ellipsoid.axes[2],
            color="black"
        )
    ax.plot_surface(
        *(
            np.matmul(
                ellipsoid.axes,
                sph.transpose(1, 0, 2)
            ).transpose(1, 0, 2) + ellipsoid.center.to_numpy().reshape(
                (3, 1, 1)
            )
        ),
        **kwargs
    )

def draw_ellipsoid(ellipsoid: Ellipsoid, *args, **kwargs):
    """Draw an ellipsoid (2D or 3D)."""    
    match ellipsoid.axes.shape[0]:
        case 2:
            draw_ellipse(ellipsoid, *args, **kwargs)
        case 3:
            draw_3d_ellipsoid(ellipsoid, *args, **kwargs)
        case _:
            raise ValueError("Cannot draw ellipsoid in > 3 dimensions.")


        
