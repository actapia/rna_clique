import operator
import functools

from collections import Counter
from collections.abc import MutableMapping, Mapping, Set

# def as_frozen_multiset(c):
#     return frozenset(Counter(c).items())

class FrozenMultiset(Set):
    def __init__(self, existing=None):
        self._set = frozenset()
        try:
            self._set = frozenset(existing._set)
        except AttributeError:
            self._set = frozenset(Counter(existing).items())

    def __contains__(self, k):
        return k in self._set

    def __iter__(self):
        return iter(self._set)

    def __len__(self):
        return len(self._set)

    def __repr__(self):
        return "FrozenMultiset({})".format(repr(set(self._set)))

    def __hash__(self):
        return hash(self._set)

    @classmethod
    def from_counts(cls, counts):
        res = cls()
        try:
            res._set = frozenset(counts.items())
        except AttributeError:
            res._set = frozenset(counts)
        return res

class MultisetKeyDict:
    def __init__(self, mapping = None):
        self._dict = {}
        if mapping is not None:
            try:
                self._dict = {k: v for (k, v) in mapping.multiset_iter()}
            except AttributeError:
                try:
                    self._dict = {FrozenMultiset(k): v for (k, v) in mapping.items()}
                except AttributeError:
                    for k, v in mapping:
                        self._dict[FrozenMultiset(k)] = v

    def __getitem__(self, k):
        return self._dict[FrozenMultiset(k)]

    def __setitem__(self, k, v):
        self._dict[FrozenMultiset(k)] = v

    def __delitem__(self, k):
        del self._dict[FrozenMultiset(k)]

    def set_iter(self):
        for k, v in self._dict.items():
            yield frozenset(x[0] for x in k), v

    def multiset_iter(self):
        yield from self._dict.items()

    def __iter__(self):
        yield from self.set_iter()
 
    def items(self):
        return self.multiset_iter()

    def __len__(self):
        return len(self._dict)

    def __contains__(self, k):
        return FrozenMultiset(k) in self._dict

    def set_keys(self):
        for k, _ in self.set_iter():
            yield k

    def multiset_keys(self):
        yield from self._dict.keys()

    def keys(self):
        return self.multiset_keys()

    def values(self):
        yield from self._dict.values()

    def _dict_op(self, op, other):
        new = MultisetKeyDict()
        new._dict = op(self._dict, other._dict)
        return new

    def __or__(self, other):
        return self._dict_op(operator.or_, other)

    def key_elements(self):
        return frozenset.union(*self.set_keys())
