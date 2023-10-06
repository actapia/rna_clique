import numpy as np

def tril_jagged(mat):
    return [r[:(i+1)].tolist() for (i, r) in enumerate(np.tril(mat))]