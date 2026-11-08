"""Check that every maximizing configuration of n cards with 3 properties,
for n = 10, 11, 12, contains a complete magic square.

Cards with 3 properties are the 27 points of F_3^3, and three cards form a set
exactly when they sum to 0 (a line of F_3^3). A magic square is one of the 39
affine planes of F_3^3. For each n, the script

  1. enumerates the n-point configurations of F_3^3,
  2. counts their sets,
  3. checks whether they contain any of the 39 affine planes, and
  4. records the maximum number of sets overall, M(n,3), and the maximum
     number of sets among configurations that contain no magic square.

If the second maximum is smaller than M(n,3), then every maximizing
configuration contains a magic square.

Fixed-set symmetry reduction (default): the affine maps x -> Ax + b with A
invertible permute the 27 cards, send sets to sets and planes to planes, and
can move any set onto any other. So every configuration with at least one set
has the same number of sets, and contains a magic square exactly when, some
configuration containing the fixed set {(0,0,0), (0,0,1), (0,0,2)}. It
therefore suffices to enumerate the C(24, n-3) configurations containing this
set; no maximizing configuration is lost. Use --full to enumerate all C(27, n)
configurations instead (slower, same maxima).

Usage:
    python set_search_3d_magic_square.py            # n = 10, 11, 12
    python set_search_3d_magic_square.py --full 10  # no symmetry reduction
"""
import argparse
import time
from itertools import combinations, product
from math import comb

POINTS = list(product(range(3), repeat=3))
INDEX = {p: i for i, p in enumerate(POINTS)}
FIXED_SET = [(0, 0, 0), (0, 0, 1), (0, 0, 2)]


def third_card(a, b):
    """The unique card that forms a set with a and b, namely -a-b."""
    return tuple((-x - y) % 3 for x, y in zip(a, b))


def to_mask(cards):
    """Bitmask with bit i set for each card POINTS[i]."""
    m = 0
    for c in cards:
        m |= 1 << INDEX[c]
    return m


def all_sets():
    """The 117 sets (lines) of F_3^3, as bitmasks."""
    lines = {to_mask([a, b, third_card(a, b)]) for a, b in combinations(POINTS, 2)}
    assert len(lines) == 117
    return sorted(lines)


def all_planes():
    """The 39 affine planes (magic squares) of F_3^3, as bitmasks."""
    planes = set()
    for a, b, c in combinations(POINTS, 3):
        if third_card(a, b) == c:
            continue  # collinear: a set, not a plane
        u = [(b[k] - a[k]) % 3 for k in range(3)]
        v = [(c[k] - a[k]) % 3 for k in range(3)]
        plane = [tuple((a[k] + i * u[k] + j * v[k]) % 3 for k in range(3))
                 for i in range(3) for j in range(3)]
        planes.add(to_mask(plane))
    assert len(planes) == 39
    return sorted(planes)


SETS = all_sets()
PLANES = all_planes()


def search(n, full=False):
    """Return (M(n,3), max sets without a magic square, number of maximizing
    configurations, how many of them contain a magic square, configurations
    checked). With full=False, only configurations containing FIXED_SET
    are enumerated."""
    if full:
        base, free, k = 0, [1 << i for i in range(27)], n
    else:
        base = to_mask(FIXED_SET)
        free = [1 << i for i in range(27) if not base >> i & 1]
        k = n - 3
    best = best_no_plane = 0
    num_best = num_best_with_plane = 0
    checked = 0
    for chosen in combinations(free, k):
        m = base | sum(chosen)
        checked += 1
        count = sum(1 for s in SETS if m & s == s)
        if count < best and count <= best_no_plane:
            continue  # cannot change either maximum
        has_plane = any(m & p == p for p in PLANES)
        if count > best:
            best, num_best, num_best_with_plane = count, 0, 0
        if count == best:
            num_best += 1
            num_best_with_plane += has_plane
        if not has_plane and count > best_no_plane:
            best_no_plane = count
    return best, best_no_plane, num_best, num_best_with_plane, checked


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("n", type=int, nargs="*", default=[10, 11, 12])
    parser.add_argument("--full", action="store_true",
                        help="enumerate all C(27,n) configurations")
    args = parser.parse_args()

    print("$n$ & $M(n,3)$ & Max sets, no magic square & Configurations checked \\\\ \\hline")
    for n in args.n:
        start = time.time()
        best, best_no_plane, num_best, num_with_plane, checked = search(n, args.full)
        expected = comb(27, n) if args.full else comb(24, n - 3)
        assert checked == expected
        print("$%d$ & $%d$ & $%d$ & $%s$ \\\\" % (n, best, best_no_plane, f"{checked:,}".replace(",", "{,}")))
        print("%% n=%d: %d maximizing configurations checked, %d contain a magic square; %.1f s"
              % (n, num_best, num_with_plane, time.time() - start))
        if best_no_plane < best:
            print("%% n=%d: every maximizing configuration contains a magic square" % n)
        else:
            print("%% n=%d: some maximizing configuration has no magic square" % n)
