import re
import string

lower_set = set(string.ascii_lowercase) | {"_"}
paren_re = re.compile(r"\(.*\)")
underscore_re = re.compile(r"_+")

def column_rename(s, remove_paren=True):
    if remove_paren:
        s = paren_re.sub("",s)
    return underscore_re.sub(
        "_",
        "".join(
            c for c in s.lower().strip().replace(" ", "_") if c in lower_set
        )
    )
