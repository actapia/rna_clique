import abc
import inspect
import typing

from collections.abc import Collection, Mapping
from dataclasses import (
    Field,
    _recursive_repr,
    dataclass,
    MISSING,
    KW_ONLY,
    _get_field,
    _MISSING_TYPE
)
from typing import Callable, Any, Optional, TypeVar

from .identity import id_

S = TypeVar("S")
T = TypeVar("T")
FieldT = TypeVar("FieldT", bound=Field)

# Based on dataclasses.field function.
def custom_field(
        cls: FieldT,
        *,        
        default=MISSING,
        default_factory: Callable[[], Any] | _MISSING_TYPE = MISSING,
        init: bool = True,
        repr: bool = True,
        hash: Optional[bool] = None,
        compare: bool = True,
        metadata: Optional[Mapping] =None,
        kw_only: bool | _MISSING_TYPE = MISSING,
        **kwargs,
) -> FieldT:
    """Create a field belonging to a custom type.

    This function is provided for the sake of conveniently specifying defaults
    for a field, much like the dataclasses.field function on which this function
    was modeled.

    By default, a created field will have no default value (specified with
    default or default_factory), will appear in the dataclass __init__ function
    and __repr__, will be used for comparisons and in __hash__, has no metadata,
    and is not a keyword-only argument in __init__.

    Parameters:
        cls (type):      Field subclass for which to make an instance.
        default:         Default value to be used.
        default_factory: Nullary function for setting a default value.
        init (bool):     Should appear in the dataclass's __init__?
        repr (bool):     Should appear in the dataclass's __repr__?
        hash (bool):     Should be used in the dataclass's __hash__?
        compare (bool):  Should be used for comparisons on the dataclass?
        metadata (dict): Arbitrary data to store with the field.
        kw_only (bool):  Should be keyword-only in dataclass's __init__?

    Returns:
        A field of the specified type, constructed with given kwargs.
    """
    if default is not MISSING and default_factory is not MISSING:
        raise ValueError('cannot specify both default and default_factory')
    return cls(default, default_factory, init, repr, hash, compare,
               metadata, kw_only,  **kwargs)

def typing_to_cast(t) -> Callable:
    """Convert a type annotation into a function to convert values to that type.

    This function attempts to convert a type annotation into a function that can
    be used to cast values to that type. The function depends on the destination
    type being constructable with just one value.

    Currently, this function handles annotations for typical unary Collection
    generics, e.g., list[int], set[float], but not dict[int, float]. This
    includes nested generics like list[set[float]]. Of course, non-collections
    are handled as well, provided it is possible to construct them with just one
    argument.

    The function also handles Optional annotations, handling the case where the
    argument is None automatically.

    Parameters:
        t: The type annotation for which to get a cast function.

    Returns:
        A function that convert values to match the provided type annotation.
    """    
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
    """Dataclass Field with functions for marshalling and unmarshalling values.

    The marshal and unmarshal functions can be accessed via the properties
    marshal and unmarshal, respectively.
    """
    __slots__ = Field.__slots__ + ("_marshal", "_unmarshal")

    def __init__(
            self,
            *args,
            marshal = None,
            unmarshal: bool | Callable = True,
            **kwargs
    ):
        """Construct a MarshallingField with specified marshallability.

        The marshal parameter can be a unary function, None, or a type
        annotation. In the first case, the provided function will be used to
        marshal values for the dataclass's marshalled representation. In the
        second case, no marshalling will be performed; the value will be kept as
        is in the marshalled representation. In the last case, the field will
        attempt to create a marshalling function from the annotation using
        type_to_cast.

        The unmarshal parameter can either be a boolean value, a unary function,
        or type annotation. A value of False indicates that no unmarshalling
        should be performed; values should be kept as they are. A value of True
        indicates that a casting function based on the MarshallingField's type
        attribute should be used (see typing_to_cast). If unmarshal is a
        function, that function will be used to marshal values for the marshalled
        representation. Finally, if unmarshal is a type annotation, the field
        will attempt to create an unmarshalling function for that annotation
        using type_to_cast.

        Parameters:
            marshal:   Function for marshalling values.
            unmarshal: Whether to unmarshal values, or a function to do so.
        """
        super().__init__(*args, **kwargs)
        self._unmarshal = unmarshal
        self._marshal = marshal

    def _unmarshal_helper(self):
        if self._unmarshal:
            if self._unmarshal is True:
                return self.type
            return self._unmarshal
        return id_

    @property
    def unmarshal(self):
        return typing_to_cast(self._unmarshal_helper())

    def _marshal_helper(self):
        if self._marshal is not None:
            return self._marshal
        return id_

    @property
    def marshal(self):
        return typing_to_cast(self._marshal_helper())

    # These function are based on Field's implementation.
        
    def _repr_field(self, f) -> str:
        """Get the representation of an attribute to use in the field's repr."""
        att = getattr(self, f)
        if f.startswith("_"):
            return str(att)
        else:
            return repr(att)

    def _repr_fields(self) -> dict[str, str]:
        """Get representations of attributes to use in the field's repr."""
        return {f: self._repr_field(f) for f in self.__slots__}

    @_recursive_repr
    def _repr_fields1(self) -> dict[str, str]:
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

def marshalling_field(
        marshal: Optional[Callable] = None,
        **kwargs
) -> MarshallingField:
    """Create a MarshallingField using sensible default values.
    
    By default, a created field will not perform marshalling, will have no
    default value (specified with default or default_factory), will appear in
    the dataclass __init__ function and __repr__, will be used for comparisons
    and in __hash__, has no metadata, and is not a keyword-only argument in
    __init__.

    Parameters:
        marshal: Function to use for marshalling values of the field.

    Returns:
        A marshalling field constructed with the given parameters.
    """
    return custom_field(MarshallingField, marshal=marshal, **kwargs)

def make_optional(f: Callable[[S], T]) -> Callable[[Optional[S]], Optional[T]]:
    """Wrap a unary function to allow it to pass None.

    The returned function always returns None if the argument is None.

    Parameters:
        f: A unary function.

    Returns:
        A wrapped version of f that return None when the argument is None.
    """
    def inner(x):
        return f(x) if x is not None else None
    return inner

class MarshallingDataclassBase(abc.ABC):
    """Abstract base class for marshalling dataclasses, for type annotations.

    Ordinary dataclasses do not derive from a specific class, which makes them
    difficult to include in annotations. This class addresses that problem for
    marshalling dataclasses.
    """
    pass

def marshalling_dataclass(
        optional: bool = False,
        marshal_type_optional: bool = True,
        unmarshal_type_optional: bool = True,
        hide_none: bool = False,
        unmarshal_init: bool = False
) -> Callable[[type], type[MarshallingDataclassBase]]:
    """Make a function to make a marshalling dataclass from a class.

    This function "wraps" a class to make a marshalling dataclass. A marshalling
    dataclass is similar to an ordinary dataclass, and, in fact, every
    marshalling dataclass is a subclass of a dataclass created from the wrapped
    class. A marshalling dataclass differs from an ordinary dataclass in that
    the former can be marshalled or unmarshalled. That is, an instances of a
    marshalled dataclass can be converted to or from a dict containing a
    system-neutral representation of the data, using only types suitable for
    serialization as JSON, YAML, etc.

    To perform marshalling and unmarshalling of a marshalling dataclass, each
    field needs to be associated with a function for marshalling and
    unmarshalling values of that field. For "primitive" types like int, this can
    just be the identity function. More complex types will require specialized
    marshalling and unmarshalling functions.

    Many types in Python's standard library can be constructed from a single
    argument of a simple type. For example, fractions.Fraction and pathlib.Path
    can both be constructed from strings. To simplify usage of marshalling
    dataclasses, unmarshalling for these cases can be handled automatically
    based on type annotations. Since many types support conversion to multiple
    different simple types, marshalling cannot be handled automatically this
    way.

    If a field cannot be marshalled or unmarshalled automatically, a
    marshalling_field will need to be constructed explicitly and provided with
    marshal and unmarshal attributes.

    When manually specifying the marshal or unmarshal attributes of a field, a
    type annotation can be provided to try to generate a marshalling or
    unmarshalling function automatically using type_to_cast. See the
    documentation for MarshallingField and type_to_cast for more information.

    If a type annotation will not work, marshal and unmarshal functions can
    always be provided explicitly. In both cases, the expected function must be
    unary.

    This function accepts various parameters that affect the creation and
    behavior of the marshalling dataclass.

    If the optional parameter is True, then fields with Optional types will
    automatically have None as their default values. Otherwise, the marshalling
    dataclass will use the same behavior as an ordinary dataclass, giving no
    default values to fields unless they are explicitly specified.

    If marshal_type_optional is also True, then for all MarshallingFields
    annotated with Optional types, the _marshal attribute will be made Optional
    automatically. Likewise, if unmarshal_type_optional and optional are both
    true, then for all MarshallingFields that are annotated with Optional types
    and speciy a custom unmarshalling function, the function will automatically
    be wrapped to pass through None values.

    Marshalling of marshalling dataclass is performed using the class's marshal
    method. This method returns a dict, the marshalled representation of the
    instance's data. Optionally, the marshalled representation can exclude None
    values. the hidee_none parameter to this function specifies whether to do
    this by default. Regardless of the setting, the default behavior can also be
    overridden by passing a value for hide_none to the marshal method itself.

    By default, no unmarshalling is performed when an instance of a marshalling
    datacalss is constructed. This enables easy copying of marshalling
    dataclasses by passing the return value of the asdict method as keyword
    arguments to __init__. To perform unmarshalling while constructing an
    instance, the from_marshalled_representation classmethod should normally be
    used instead. In some cases, it might make more sense to have unmarshalling
    performed at __init__. To enable this, pass True for the unmarshal_init
    parameter.

    Parameters:        
        optional (bool):                Use None as default default for Optional
                                        fields.
        marshal_type_optional (bool):   Add None marshalling passthrough for
                                        Optional fields, if optional is also
                                        True.
        unmarshal_type_optional (bool): Add None unmarshalling passthrough for
                                        Optional fields, if optional is also
                                        True.
        hide_none (bool):               Exclude fields with None value from
                                        marshalled representation.
        unmarshal_init (bool):          Unmarshal values in __init__.

    Returns:
        Function mapping classes to marshalling dataclasses with given options.
    """
    def inner(c: type) -> type[MarshallingDataclassBase]:
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
        class wrapped(dataclass(c), MarshallingDataclassBase):
            def __init__(self, *args, **kwargs):
                super().__init__(*args, **kwargs)
                if unmarshal_init:
                    self._unmarshal()

            def _unmarshal(self):
                """Unmarshal attributes of the class, updating them."""
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
            def from_marshalled_representation(cls, d: dict):
                """Make an instance of the dataclass from a marshalled dict.

                This function constructs an instance of the marshalling
                dataclass using values from a marshalled representation stored
                in a dict. The values are then unmarshalled immediately after
                construction using the unmarshalling functions specified
                (implicitly or explicitly) in the definition of the dataclass.

                Parameters:
                    d (dict): dict containing marshalled representation of data.

                Returns:
                    Instance of this dataclass with same data, unmarshalled.
                """
                res = cls(**d)
                if not unmarshal_init:
                    res._unmarshal()
                return res

            def marshal(self, hide_none: bool = hide_none) -> dict:
                """Make a marshalled representation of this dataclass instance.

                This function uses the marshal functions associated with each
                ofthe dataclass's fields to create a "marshalled representation"
                of the dataclass instance.

                Parameters:
                    hide_node (bool): Exclude fields with None values.

                Returns:
                    A dict, a marshalled representation of this instance.
                """
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
