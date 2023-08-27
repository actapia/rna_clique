import colorsys
from itertools import permutations
from collections import namedtuple

hls = namedtuple("hls", ["hue", "lightness", "saturation"])
hsv = namedtuple("hsv", ["hue", "saturation", "value"])
rgb = namedtuple("rgb", ["red", "green", "blue"])
yiq = namedtuple("yiq", ["y", "i", "q"])
chrominance = namedtuple("chrominance", ["i", "q"])
yiq.luma = property(lambda x: x.y)
yiq.chrominance = property(lambda x: chrominance(x.i, x.q))
spaces = {"hls": hls, "hsv": hsv, "rgb": rgb, "yiq": yiq}
__all__ = list(spaces)
for s in spaces:
    spaces[s].to_rgb = (lambda _s: lambda x: rgb(*(getattr(colorsys, f"{_s}_to_rgb")(*x))))(s)
non_rgb = set(spaces) - {"rgb"}
for s in non_rgb:
    setattr(rgb, f"to_{s}", (lambda _s: lambda x: spaces[_s](*(getattr(colorsys, f"rgb_to_{_s}")(*x))))(s))
for s1, s2 in permutations(set(spaces) - {"rgb"}, 2):
    setattr(
        spaces[s1],
        f"to_{s2}",
        (lambda _s2: lambda x: spaces[_s2](*(getattr(x.to_rgb(), f"to_{_s2}")())))(s2)
    )
