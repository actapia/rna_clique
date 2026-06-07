import re
import string

lower_set = set(string.ascii_lowercase) | {"_"}
paren_re = re.compile(r"\(.*\)")
underscore_re = re.compile(r"_+")

def column_rename(s: str, remove_paren: bool = True) -> str:
    """Rename a Pandas column to make it easier to access as an attribute.

    Parameters:
        s (str):             The original column name.
        remove_paren (bool): Remove text within parentheses.

    Returns:
        New string based on the original using letters and underscores only.
    """
    if remove_paren:
        s = paren_re.sub("",s)
    return underscore_re.sub(
        "_",
        "".join(
            c for c in s.lower().strip().replace(" ", "_") if c in lower_set
        )
    )
