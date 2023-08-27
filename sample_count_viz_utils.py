import pandas as pd

def read_sample_count_test_results(path):
    return pd.read_csv(
        path,
        sep=" ",
        names=[
            "samples",
            "components",
            "big",
            "ideal",
            "n"
        ]
    )