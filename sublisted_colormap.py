import matplotlib as mpl

class SublistedColormap(mpl.colors.ListedColormap):
    def __getitem__(self, idx):
        try:
            return SublistedColormap(
                self.colors[slice(idx.start, idx.stop, idx.step)]
            )
        except AttributeError:
            return self.colors[idx]

    @classmethod
    def from_listedcolormap(cls, colormap):
        return SublistedColormap(colormap.colors)

def convert_to_slcm(cm):
    if isinstance(cm, mpl.colors.ListedColormap):
        return SublistedColormap.from_listedcolormap(cm)
    else:
        return cm

colormaps = {k : convert_to_slcm(v) for (k, v) in mpl.colormaps.items()}
