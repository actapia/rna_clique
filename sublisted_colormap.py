import matplotlib as mpl

class SublistedColormap(mpl.colors.ListedColormap):
    """ListedColormap that supports slicing."""
    def __getitem__(self, idx):
        try:
            return SublistedColormap(
                self.colors[slice(idx.start, idx.stop, idx.step)]
            )
        except AttributeError:
            return self.colors[idx]

    @classmethod
    def from_listedcolormap(cls, colormap: mpl.colors.ListedColormap):
        """Create a SublistedColormap from a ListedColormap"""
        return SublistedColormap(colormap.colors)

def convert_to_slcm(cm: mpl.colors.Colormap) -> mpl.colors.Colormap:
    """Converts ListedColormaps to SublistedColormaps and leaves others alone.

    Parameters:
        cm: A Colormap to possibly convert.

    Returns:
        If cm is a ListedColormap, a corresponding SublistedColormap, else cm.
    """
    if isinstance(cm, mpl.colors.ListedColormap):
        return SublistedColormap.from_listedcolormap(cm)
    else:
        return cm

colormaps = {k : convert_to_slcm(v) for (k, v) in mpl.colormaps.items()}
