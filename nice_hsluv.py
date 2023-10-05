import hsluv as hsluv_
from collections import namedtuple
from nice_colorsys import spaces, rgb
hsluv = namedtuple("hsluv", ["hue", "saturation", "lightness"])
spaces["hsluv"] = hsluv
hsluv.to_rgb = lambda x: rgb(*hsluv_.hsluv_to_rgb(x))
spaces["rgb"].to_hsluv = lambda x: hsluv(*hsluv_.rgb_to_hsluv(x))
for s in (set(spaces) - {"rgb", "hsluv"}):
    spaces[s].to_hsluv = lambda x: x.to_rgb().to_hsluv()
for s in (set(spaces) - {"rgb", "hsluv"}):
    setattr(
        hsluv,
        f"to_{s}",
        (lambda _s: lambda x: getattr(x.to_rgb(), f"to_{_s}")())(s)
    )
__all__ = ["hsluv"]