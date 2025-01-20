from itertools import combinations 
import time

def add_mod3(a, b, c):
    """Component-wise addition mod 3 of two 3D tuples."""
    return tuple((a[i] + b[i]+c[i]) % 3 for i in range(3))

def is_set(a,b,c):
    """Returns true if three elements form a set."""
    return add_mod3(a,b,c) == (0, 0, 0)

if __name__ == "__main__":
    all_points = [(x,y,z) for x in range(3) for y in range(3) for z in range(3)]
    max_sets = 0
    for n in range(3, 28):
        begin_time = time.time()
        for comb in combinations(all_points, n):
            
            num_sets = 0
            for possible_set in combinations(comb, 3):
                if is_set(possible_set[0], possible_set[1], possible_set[2]):
                    num_sets +=1
            if num_sets > max_sets:
                max_sets = num_sets
        end_time = time.time()
       
        print("Max sets for %d cards: %d"%(n, max_sets))