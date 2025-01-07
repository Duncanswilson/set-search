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
    for comb in tqdm(combinations(all_points, 8)):
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
        if num_sets == 8:
            num_same_elements = 0
            for i in range(4):
                if(comb[0][i] == comb[1][i] == comb[2][i] == comb[3][i] == comb[4][i] == comb[5][i] == comb[6][i] == comb[7][i]):
                    num_same_elements+=1
            if(num_same_elements != 2):
                print("Non-subset of F_3^2 found")
                print(comb)
                break 
        if num_sets > 8:
            print("More than 8 sets found")
            print(comb)
            break 
