"""For n = 10, 11, 12 cards with three properties, check that every
configuration with the maximum number of sets contains a magic square.

For each n, the script enumerates all n-card configurations, counts their
sets, checks whether they contain any of the 39 magic squares (affine planes
of F_3^3), and records the maximum number of sets overall and among the
configurations with no magic square.
"""
from itertools import combinations, product

points = list(product(range(3), repeat=3))


def third_card(a, b):
    """The card that forms a set with a and b."""
    return tuple((-x - y) % 3 for x, y in zip(a, b))


def mask(cards):
    """A collection of cards as a bitmask over the 27 points."""
    return sum(1 << points.index(c) for c in cards)


# All 117 sets and all 39 magic squares, as bitmasks.
sets = {mask([a, b, third_card(a, b)]) for a, b in combinations(points, 2)}
planes = set()
for a, b, c in combinations(points, 3):
    if c != third_card(a, b):
        planes.add(mask([tuple((a[k] + i * (b[k] - a[k]) + j * (c[k] - a[k])) % 3
                               for k in range(3))
                         for i in range(3) for j in range(3)]))
assert len(sets) == 117 and len(planes) == 39

if __name__ == "__main__":
    for n in (10, 11, 12):
        max_sets = 0
        max_sets_no_magic_square = 0
        for config in combinations(range(27), n):
            m = sum(1 << i for i in config)
            num_sets = sum(1 for s in sets if m & s == s)
            max_sets = max(max_sets, num_sets)
            if not any(m & p == p for p in planes):
                max_sets_no_magic_square = max(max_sets_no_magic_square, num_sets)
        print("n = %d: max sets %d, max sets with no magic square %d"
              % (n, max_sets, max_sets_no_magic_square))
