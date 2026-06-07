import typing

def unoptional(t):
    fields = typing.get_args(t)
    type_ = typing.get_origin(t)
    if type_ is typing.Union:
        non_none = [f for f in fields if f is not type(None)]            
        if len(fields) == 2 and len(non_none) == 1:
            # Handle "Optional" types.
            non_none = non_none[0]
            return unoptional(non_none)
    return t

def get_type_name(t):
    if isinstance(t, type):
        type_name = t.__qualname__
    else:
        type_name = str(t)
    if t.__module__ != "builtins":
        type_name = "{}.{}".format(t.__module__, type_name)
    return type_name

def column_to_text(c):
    sp = c.split("_")
    sp[0] = sp[0][0].upper() + sp[0][1:]
    return " ".join(sp)
