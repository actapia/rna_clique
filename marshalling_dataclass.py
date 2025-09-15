import inspect
import typing
from collections.abc import Collection, Mapping
from dataclasses import (
    Field,
    _recursive_repr,
    dataclass,
    MISSING,
    KW_ONLY,
    _get_field
)

# Based on dataclasses.field function.
def custom_field(
        cls,
        *,        
        default=MISSING,
        default_factory=MISSING,
        init=True,
        repr=True,
        hash=None,
        compare=True,
        metadata=None,
        kw_only=MISSING,
        **kwargs,
):
    if default is not MISSING and default_factory is not MISSING:
        raise ValueError('cannot specify both default and default_factory')
    return cls(default, default_factory, init, repr, hash, compare,
               metadata, kw_only,  **kwargs)

def typing_to_cast(t):
    fields = typing.get_args(t)
    type_ = typing.get_origin(t)
    if type_ is not None:
        if isinstance(type_, type):
            if issubclass(type_, tuple):            
                return lambda x: type_(
                    typing_to_cast(f)(v) for (f, v) in zip(fields, x)
                )
            if issubclass(type_, Collection):
                if issubclass(type_, Mapping):
                    if len(fields) == 2:
                        return lambda x: type_(
                            tuple(
                                typing_to_cast(f)(v)
                                for (f, v) in zip(fields, y)
                            ) for y in x.items()
                        )
                elif len(fields) == 1:
                    return lambda x: type_(map(typing_to_cast(fields[0]), x))
                raise ValueError(f"Cannot automatically (un)marshal {t!r}.")
        elif type_ is typing.Union:
            non_none = [f for f in fields if f is not type(None)]            
            if len(fields) == 2 and len(non_none) == 1:
                # Handle "Optional" types.
                non_none = non_none[0]
                return lambda x: typing_to_cast(
                    non_none
                )(x) if x is not None else None
    return t
        
class MarshallingField(Field):
    __slots__ = Field.__slots__ + ("_marshal", "_unmarshal")

    def __init__(self, *args, marshal=False, unmarshal=True, **kwargs):
        super().__init__(*args, **kwargs)
        self._unmarshal = unmarshal
        self._marshal = marshal

    def _unmarshal_helper(self):
        if self._unmarshal:
            if self._unmarshal is True:
                return self.type
            return self._unmarshal
        return lambda x: x

    @property
    def unmarshal(self):
        return typing_to_cast(self._unmarshal_helper())

    def _marshal_helper(self):
        if self._marshal:
            return self._marshal
        return lambda x: x    

    @property
    def marshal(self):
        return typing_to_cast(self._marshal_helper())
        
    def _repr_field(self, f):
        att = getattr(self, f)
        if f.startswith("_"):
            return str(att)
        else:
            return repr(att)

    def _repr_fields(self):
        return {f: self._repr_field(f) for f in self.__slots__}

    @_recursive_repr
    def _repr_fields1(self):
        return self._repr_fields()

    def __repr__(self):
        return "{}({})".format(
            type(self).__name__,
            ",".join(
                map(
                    "=".join,
                    self._repr_fields1().items()
                )
            )
        )

def marshalling_field(marshal=False, **kwargs):
    return custom_field(MarshallingField, marshal=marshal, **kwargs)

def make_optional(f):
    def inner(x):
        return f(x) if x is not None else None
    return inner

def marshalling_dataclass(
        optional=False,
        marshal_type_optional=True,
        unmarshal_type_optional=True,
        hide_none=False,
        unmarshal_init=False
):
    def inner(c):
        annotations = inspect.get_annotations(c)
        for name, type_ in annotations.items():
            if not type_ is KW_ONLY:
                f = _get_field(c, name, type_, MISSING)
                if optional:
                    fields = typing.get_args(type_)
                    type_ = typing.get_origin(type_)
                    if type_ is typing.Union and any(
                            t is type(None) for t in fields
                    ):
                        if f.default is MISSING:
                            f.default = None
                        try:
                            if marshal_type_optional and f._marshal:
                                f._marshal = typing.Optional[f._marshal]
                            if unmarshal_type_optional and callable(
                                    f._unmarshal
                            ):
                                f._unmarshal = make_optional(f._unmarshal)
                        except AttributeError:
                            pass
                    
                setattr(c, name, f)
        class wrapped(dataclass(c)):
            def __init__(self, *args, **kwargs):
                super().__init__(*args, **kwargs)
                if unmarshal_init:
                    self._unmarshal()

            def _unmarshal(self):
                for name, field in self.__dataclass_fields__.items():
                    try:                    
                        setattr(
                            self,
                            name,
                            field.unmarshal(getattr(self, name))
                        )
                    except AttributeError:
                        pass

            @classmethod
            def from_marshalled_representation(cls, d):
                res = cls(**d)
                if not unmarshal_init:
                    res._unmarshal()
                return res

            def marshal(self, hide_none=hide_none):
                res = {}
                for name, field in self.__dataclass_fields__.items():
                    try:
                        value = field.marshal(getattr(self, name))
                    except AttributeError:
                        value = getattr(self, name)
                    if value is not None or not hide_none:
                        res[name] = value
                return res

        for attr in ['__doc__', '__name__', '__qualname__', '__module__']:
            setattr(wrapped, attr, getattr(c, attr))
        return wrapped
    return inner
