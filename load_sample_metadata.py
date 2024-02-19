from pathlib import Path
import pandas as pd

def get_sample_metadata(metadata_path: Path) -> pd.DataFrame:
    sample_info = pd.read_csv(
        metadata_path,
        sep=" ",
        names=[
            "url",
            "name",
            "genotype",
            "endophyte"
        ]
    )
    sample_info["endophyte"] = sample_info["endophyte"] == "infected"
    return sample_info.sort_values(
        ["name"]
    ).reset_index(drop=True).sort_values(["genotype", "endophyte"])
