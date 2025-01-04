from tqdm import tqdm
import itertools
import logging

logging.basicConfig(level=logging.DEBUG, format='%(message)s')

###############################################################################
# 1. Basic GF(3) arithmetic helpers
###############################################################################

def add_mod3(a, b):
    """Component-wise addition mod 3 of two 4D tuples."""
    return tuple((a[i] + b[i]) % 3 for i in range(4))

def sub_mod3(a, b):
    """Component-wise subtraction mod 3 of two 4D tuples."""
    return tuple((a[i] - b[i]) % 3 for i in range(4))

def scalar_mul_mod3(v, alpha):
    """Multiply a 4D tuple v by scalar alpha in GF(3)."""
    return tuple((alpha * v[i]) % 3 for i in range(4))

def mat_vec_mul_mod3(M, v):
    """
    Multiply a 4x4 matrix M (list of 4 lists) by a 4D vector v, all mod 3.
    Returns a 4D tuple.
    """
    # M is a 4x4, M[r] is row r
    # v is a 4D tuple
    res = [0,0,0,0]
    for r in range(4):
        val = 0
        for c in range(4):
            val += M[r][c]*v[c]
        res[r] = val % 3
    return tuple(res)

def identity_matrix_4():
    """Return the 4x4 identity matrix over GF(3)."""
    return [[1 if r==c else 0 for c in range(4)] for r in range(4)]

def copy_matrix(M):
    """Return a (deep) copy of a 4x4 matrix M."""
    return [row[:] for row in M]

###############################################################################
# 2. Points and lines in AG(4,3)
###############################################################################

def all_points_ag43():
    """Return a list of all 81 points (4D tuples) in AG(4,3)."""
    pts = []
    for x1 in range(3):
        for x2 in range(3):
            for x3 in range(3):
                for x4 in range(3):
                    pts.append((x1,x2,x3,x4))
    return pts

def canonical_directions_gf3_4():
    """
    Return a list of 40 canonical direction vectors in GF(3)^4,
    one per equivalence class d ~ 2d.
    """
    all_nonzero = []
    for x1 in range(3):
        for x2 in range(3):
            for x3 in range(3):
                for x4 in range(3):
                    if (x1,x2,x3,x4) != (0,0,0,0):
                        all_nonzero.append((x1,x2,x3,x4))
    used = set()
    directions = []
    for v in all_nonzero:
        v2 = scalar_mul_mod3(v, 2)  # 2 ~ -1 mod 3
        canon = min(v, v2)          # pick lexicographically smaller as canonical
        if canon not in used:
            directions.append(canon)
            used.add(canon)
            used.add(v2)
    return directions

def all_lines_ag43():
    """
    Return a set of all 1080 lines in AG(4,3), each line as a 3-tuple of points
    (sorted in lex order).
    """
    pts = all_points_ag43()
    dirs = canonical_directions_gf3_4()
    line_set = set()
    
    for p in pts:
        for d in dirs:
            p1 = add_mod3(p, d)
            p2 = add_mod3(p, scalar_mul_mod3(d, 2))
            triple = sorted([p, p1, p2])
            line_set.add(tuple(triple))
    
    return line_set
def has_geometric_potential(partial_subset, new_point):
    """Check if adding new_point maintains geometric properties we want"""
    if len(partial_subset) < 2:
        return True
        
    # Check if new_point forms at least one line with existing points
    forms_line = False
    for p1, p2 in itertools.combinations(partial_subset, 2):
        if points_form_line(p1, p2, new_point):
            forms_line = True
            break
    
    return forms_line

def points_form_line(p1, p2, p3):
    """Check if three points are collinear in AG(4,3)"""
    # Vector from p1 to p2
    v12 = sub_mod3(p2, p1)
    # Vector from p1 to p3
    v13 = sub_mod3(p3, p1)
    # Check if vectors are parallel (one is scalar multiple of other)
    return are_parallel_mod3(v12, v13)
###############################################################################
# 3. Building a canonical label for a subset of points
#
#    The idea: two subsets that differ by an affine transform in AG(4,3)
#    should map to the *same* canonical label (a sorted tuple of points).
#
#    We'll do a two-step approach:
#       (a) Translate so that the "smallest" point in the subset goes to 0.
#       (b) If there is a second distinct point, attempt to map that point to (1,0,0,0).
#           (We build a linear transform that sends p1 -> e1).
#
#    The result is then sorted lexicographically, so we can store it in a set/dict.
###############################################################################
def get_orbit_representatives():
    """Get representatives of orbits under the action of the affine group"""
    reps = set()
    for pts in itertools.combinations(all_points_ag43(), 3):
        canon = canonical_label(pts)
        if canon not in reps:
            reps.add(canon)
    return reps

def find_basis_for_vector(v):
    """
    Given a nonzero vector v in GF(3)^4, produce a 4x4 invertible matrix M
    whose first column is v.
    """
    all_pts = all_points_ag43()
    basis = [v]
    for candidate in all_pts:
        if candidate == (0,0,0,0):
            continue
        if is_linear_combination_mod3(candidate, basis):
            continue
        basis.append(candidate)
        if len(basis) == 4:
            break
    
    M = [[0]*4 for _ in range(4)]
    for c in range(4):
        for r in range(4):
            M[r][c] = basis[c][r]
    return M

def is_linear_combination_mod3(v, vec_list):
    """
    Check if 'v' in GF(3)^4 is a linear combination of vectors in vec_list
    (also in GF(3)^4).
    
    We'll do a small system solve, or a trivial brute force if the list is short.
    """
    # For a small list, we can brute force:
    # But let's do a quick approach using a rank-based method:
    
    # We'll put vectors [vec_list, v] into a matrix and check rank.
    mat = [list(x) for x in vec_list]  # each x is 4D tuple
    # Convert v to a row as well
    original_len = len(mat)
    mat.append(list(v))
    
    # Compute rank mod 3
    r = rank_mod3(mat)
    # If rank didn't increase after adding 'v', then v was a combination.
    return (r == original_len)

def rank_mod3(matrix_rows):
    """
    Compute the rank of a matrix (given as a list of row lists) mod 3
    using a standard Gaussian elimination approach, but mod 3.
    """
    mat = [row[:] for row in matrix_rows]  # copy
    R = len(mat)
    C = len(mat[0]) if R>0 else 0
    row_idx = 0
    
    for col in range(C):
        # Find pivot row
        pivot = -1
        for r in range(row_idx, R):
            if mat[r][col] % 3 != 0:
                pivot = r
                break
        if pivot < 0:
            continue  # no pivot in this column
        # swap pivot row to row_idx
        if pivot != row_idx:
            mat[row_idx], mat[pivot] = mat[pivot], mat[row_idx]
        
        # pivot value
        pv = mat[row_idx][col] % 3
        inv_pv = inv_mod3(pv)  # multiply by this to normalize pivot to 1
        
        # normalize pivot row
        for cc in range(col, C):
            mat[row_idx][cc] = (mat[row_idx][cc]*inv_pv) % 3
        
        # eliminate below
        for rr in range(row_idx+1, R):
            if mat[rr][col] != 0:
                factor = mat[rr][col] % 3
                for cc in range(col, C):
                    mat[rr][cc] = (mat[rr][cc] - factor*mat[row_idx][cc])%3
        
        row_idx += 1
        if row_idx == R:
            break
    
    return row_idx  # number of pivots found

def inv_mod3(x):
    """Return inverse of x mod 3, assuming x in {1,2}."""
    # 1 -> 1, 2 -> 2 (because 2*2 = 4 = 1 mod 3)
    return 1 if x==1 else 2

def invert_4x4_matrix_mod3(M):
    """
    Invert a 4x4 matrix mod 3. Return the inverse matrix (4x4).
    We assume M is invertible. If not, you get an error or nonsense.
    """
    # We'll augment M with identity and do row ops
    mat = [row[:] + row_id[:] for row, row_id in zip(M, identity_matrix_4())]
    # mat is now 4x8
    # Do Gaussian elimination
    R = 4; C = 8
    # forward
    pivot_row = 0
    for col in range(4):
        # find pivot in or below pivot_row
        pivot = -1
        for r in range(pivot_row, 4):
            if mat[r][col] % 3 != 0:
                pivot = r
                break
        if pivot<0:
            raise ValueError("Matrix not invertible mod 3 (no pivot found).")
        if pivot != pivot_row:
            mat[pivot_row], mat[pivot] = mat[pivot], mat[pivot_row]
        pv = mat[pivot_row][col] % 3
        inv_pv = inv_mod3(pv)
        # scale pivot row
        for cc in range(col, C):
            mat[pivot_row][cc] = (mat[pivot_row][cc]*inv_pv)%3
        # eliminate below
        for rr in range(pivot_row+1, 4):
            if mat[rr][col] != 0:
                factor = mat[rr][col]
                for cc in range(col, C):
                    mat[rr][cc] = (mat[rr][cc] - factor*mat[pivot_row][cc])%3
        pivot_row += 1
    
    # back
    for col in reversed(range(4)):
        # pivot row is col in this scheme
        pr = col
        # pivot is mat[pr][col] = 1
        for rr in range(pr):
            factor = mat[rr][col]
            for cc in range(col, C):
                mat[rr][cc] = (mat[rr][cc] - factor*mat[pr][cc])%3
    
    # Now the right half of mat is M^-1
    invM = []
    for r in range(4):
        invM.append(mat[r][4:8])
    return invM

def canonical_label(subset):
    """
    Return a tuple of points that is the canonical representative of 'subset'
    under affine transformations of AG(4,3).
    """
    if not subset:
        return ()
    
    # 1) find p0 = min(subset)
    p0 = min(subset)
    # 2) translate
    translated = [sub_mod3(p, p0) for p in subset]
    
    # 3) find p1
    translated.sort()
    if translated[0] != (0,0,0,0):
        raise RuntimeError("Logic error: after translation, first point must be 0?")
    if len(translated) == 1:
        return tuple(translated)
    
    p1 = None
    for pt in translated:
        if pt != (0,0,0,0):
            p1 = pt
            break
    
    if p1 is None:
        return tuple(translated)
    
    try:
        M_cols = find_basis_for_vector(p1)
        B = invert_4x4_matrix_mod3(M_cols)
    except:
        return tuple(sorted(translated))
    
    transformed = [mat_vec_mul_mod3(B, pt) for pt in translated]
    transformed.sort()
    return tuple(transformed)

###############################################################################
# 4. Example: enumerating subsets up to canonical form
#
#    We'll do it in a naive way for small N.  For each combination of size 3N,
#    we compute its canonical form.  If we've never seen that form, we count
#    how many lines are contained, and track the maximum.
#
#    For big N, you'd do a more careful *backtracking* with incremental
#    canonical forms/pruning, not a raw combination loop.
###############################################################################

def lines_in_subset_count(subset, all_lines):
    """Count how many lines (each a 3-tuple) are fully contained in 'subset'."""
    subs = set(subset)
    count = 0
    for line in all_lines:
        # line is a 3-tuple
        if line[0] in subs and line[1] in subs and line[2] in subs:
            count += 1
    return count

def are_parallel_mod3(v1, v2):
    """Check if two vectors are parallel in GF(3)^4 (one is scalar multiple of other)"""
    # Find first nonzero component in v1
    for i, x in enumerate(v1):
        if x != 0:
            # If v1[i] = a and v2[i] = b, then check if v1 = (b/a)v2
            scalar = 1 if x == v2[i] else (2 if x == (2 * v2[i]) % 3 else 0)
            if scalar == 0:
                return False
            # Check if this scalar works for all components
            return all((x * scalar) % 3 == y for x, y in zip(v1, v2))
    return False

def points_form_line(p1, p2, p3):
    """Check if three points are collinear in AG(4,3)"""
    # Vector from p1 to p2
    v12 = sub_mod3(p2, p1)
    # Vector from p1 to p3
    v13 = sub_mod3(p3, p1)
    # Check if vectors are parallel (one is scalar multiple of other)
    return are_parallel_mod3(v12, v13)

def max_lines_in_3N_subset_backtrack(N):
    """
    Use backtracking with incremental canonical forms to find the max # of lines
    in any subset of size 3N in AG(4,3).
    """
    points = all_points_ag43()
    all_lines = list(all_lines_ag43())
    seen_canonical = set()
    best_result = {'count': 0, 'subset': None}
    
    # Initialize progress tracking
    progress = {'combinations_checked': 0}
    progress_bar = tqdm(desc=f"Searching N={N}")

    def has_geometric_potential(partial_subset, new_point):
        """
        Modified to be less restrictive - a point has potential if either:
        1. It forms a line with existing points
        2. The subset size is small enough that we need more points anyway
        """
        if len(partial_subset) < 4:  # Always accept early points
            return True
            
        # Count how many lines this point could potentially participate in
        potential_lines = 0
        for p1, p2 in itertools.combinations(partial_subset, 2):
            if points_form_line(p1, p2, new_point):
                potential_lines += 1
                if potential_lines >= 2:  # If point could form multiple lines, it's promising
                    return True
        
        return potential_lines > 0  # Accept if it forms at least one line

    def count_partial_lines(subset):
        return sum(1 for line in all_lines if
                  sum(1 for p in line if p in subset) == 2)

    def estimate_max_possible_lines(partial_subset, remaining_points, points_needed):
        """
        More accurate estimation of maximum possible lines we could form
        by adding points_needed more points to partial_subset.
        """
        current_size = len(partial_subset)
        current_lines = lines_in_subset_count(partial_subset, all_lines)
        
        # Count partial lines (2 points from subset)
        partial_lines = count_partial_lines(partial_subset)
        
        # Estimate additional lines possible with new points
        # This is a very rough upper bound based on maximum possible lines
        # through each new point we could add
        max_new_lines_per_point = min(
            current_size * (current_size - 1) // 2,  # Max lines through existing points
            points_needed * (points_needed - 1) // 2  # Max lines between new points
        )
        
        return current_lines + partial_lines + max_new_lines_per_point

    def can_extend_to_solution(partial_subset, remaining_points, target_lines):
        """
        Modified to use better estimation of potential lines
        """
        current_size = len(partial_subset)
        points_needed = 3*N - current_size
        if points_needed < 0:
            return False
                
        # If we already have enough lines, great!
        current_lines = lines_in_subset_count(partial_subset, all_lines)
        if current_lines >= target_lines:
            return True
        
        # Otherwise, estimate maximum possible lines we could achieve
        max_possible = estimate_max_possible_lines(
            partial_subset, remaining_points, points_needed
        )
        
        return max_possible >= target_lines

    def backtrack(partial_subset, start_idx):
        if len(partial_subset) == 3*N:
            progress['combinations_checked'] += 1
            if progress['combinations_checked'] % 1000 == 0:
                progress_bar.update(1000)
                progress_bar.set_postfix({'found': best_result['count']})
            
            canform = canonical_label(partial_subset)
            if canform not in seen_canonical:
                seen_canonical.add(canform)
                line_count = lines_in_subset_count(partial_subset, all_lines)
                if line_count > best_result['count']:
                    best_result['count'] = line_count
                    best_result['subset'] = partial_subset[:]
                    print(f"\nNew best found! Lines: {line_count}")
                    print(f"Subset: {partial_subset}")
            return

        # Handle first point
        if len(partial_subset) == 0:
            partial_subset.append((0,0,0,0))
            backtrack(partial_subset, 1)
            partial_subset.pop()
            return

        # Handle second point (use canonical representatives)
        if len(partial_subset) == 1:
            canonical_seconds = [
                (1,0,0,0), (1,1,0,0), (1,1,1,0), (1,1,1,1),
                (2,0,0,0), (2,1,0,0), (2,1,1,0)
            ]
            for p in canonical_seconds:
                partial_subset.append(p)
                backtrack(partial_subset, 2)
                partial_subset.pop()
            return

        if not can_extend_to_solution(partial_subset, points[start_idx:], best_result['count']):
            return
            
        for i in range(start_idx, len(points)):
            new_point = points[i]
            if not has_geometric_potential(partial_subset, new_point):
                continue
            partial_subset.append(new_point)
            backtrack(partial_subset, i + 1)
            partial_subset.pop()

    backtrack([], 0)
    progress_bar.close()
    return best_result['count']

###############################################################################
# 5. Demonstration
###############################################################################

if __name__ == "__main__":
    for testN in [1,2,3,4]:
        print(f"\n=== Searching for N={testN} (subset size={3*testN}) ===")
        maxc = max_lines_in_3N_subset_backtrack(testN)
        print(f"Max # of lines in a {3*testN}-point subset = {maxc}")