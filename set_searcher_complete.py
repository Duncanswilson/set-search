from tqdm import tqdm
from itertools import combinations 
import logging

logging.basicConfig(level=logging.DEBUG, format='%(message)s')

###############################################################################
# 1. Basic GF(3) arithmetic helpers
###############################################################################

def add_mod3(a, b, c):
    """Component-wise addition mod 3 of two 4D tuples."""
    return tuple((a[i] + b[i]+c[i]) % 3 for i in range(4))

def is_set(a,b,c):
    """Returns true if three elements form a set."""
    return add_mod3(a,b,c) == (0, 0, 0, 0)
        
def all_points_ag43():
    """Return a list of all 81 points (4D tuples) in AG(4,3)."""
    pts = []
    for x1 in range(3):
        for x2 in range(3):
            for x3 in range(3):
                for x4 in range(3):
                    pts.append((x1,x2,x3,x4))
    return pts


if __name__ == "__main__":
    # https://docs.python.org/3/library/itertools.html#itertools.combinations
    all_points = all_points_ag43()
    max_sets = 0
    max_config = [] 
    for comb in tqdm(combinations(all_points, 6)):
        num_sets = 0
        for possible_set in combinations(comb, 3):
            if is_set(possible_set[0], possible_set[1], possible_set[2]):
                num_sets +=1
        if num_sets > max_sets:
            max_sets = num_sets
            max_config = comb
            print("New maximum number of sets found", max_sets)
            print(max_config)
            print("")
    print("Maximum number of sets found", max_sets)
    print("Best configuration found", max_config)
