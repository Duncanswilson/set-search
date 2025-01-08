from itertools import combinations, product
from typing import List, Tuple, Set

def generate_deck(num_properties: int) -> List[Tuple[int, ...]]:
    """Generate a deck of SET cards with given number of properties."""
    return [tuple(x) for x in product(range(3), repeat=num_properties)]

def forms_set(cards: List[Tuple[int, ...]]) -> bool:
    """Check if three cards form a SET.
    
    For each property, all cards must have either:
    - All same values, or
    - All different values
    """
    if len(cards) != 3:
        return False
    
    # Check each property position
    for values in zip(*cards):
        unique_vals = len(set(values))
        if unique_vals == 2:  # If 2 same and 1 different, not a SET
            return False
        # Must be either all same (1) or all different (3)
        if unique_vals not in (1, 3):
            return False
    return True

def find_sets(cards: List[Tuple[int, ...]], current: List[Tuple[int, ...]], added_card: Tuple[int, ...]) -> Set[Tuple[Tuple[int, ...], ...]]:
    """Find new SETs formed by adding a card to existing selection."""
    sets = set()
    # Look at pairs from current cards
    for card1 in current:
        for card2 in current:
            if card1 < card2:  # Avoid duplicates
                if forms_set([card1, card2, added_card]):
                    sets.add(tuple(sorted([card1, card2, added_card])))
    return sets

def get_cube_index(card: Tuple[int, ...]) -> int:
    """Get cube index (0-8) following paper's layout.
    Uses first two properties to determine cube placement."""
    return card[0] * 3 + card[1]  # No modulo needed since values already 0-2

def consecutive_maximization(deck: List[Tuple[int, ...]], n: int) -> List[Tuple[int, ...]]:
    """Implement Consecutive Maximization Method."""
    if n < 3:
        return deck[:n]
    
    selected = []
    
    # Start with first card from first cube
    selected.append(deck[0])
    first_cube = get_cube_index(deck[0])
    
    # Find second card from different cube
    for card in deck[1:]:
        if get_cube_index(card) != first_cube:
            selected.append(card)
            second_cube = get_cube_index(card)
            break
    
    # Find third card that forms a SET with first two cards
    for card in deck:
        if card not in selected and forms_set([selected[0], selected[1], card]):
            selected.append(card)
            break
    
    # For remaining selections
    remaining = [card for card in deck if card not in selected]
    used_cubes = {get_cube_index(card) for card in selected}
    
    while len(selected) < n:
        # Try to select from new cube when possible
        best_card = None
        
        # First try to find a card from a new cube that forms most SETs
        max_sets = -1
        for card in remaining:
            cube = get_cube_index(card)
            if cube not in used_cubes:
                num_sets = len(find_sets(deck, selected, card))
                if num_sets > max_sets:
                    max_sets = num_sets
                    best_card = card
        
        # If no card from new cube, find card that forms most SETs
        if not best_card:
            for card in remaining:
                num_sets = len(find_sets(deck, selected, card))
                if num_sets > max_sets:
                    max_sets = num_sets
                    best_card = card
        
        selected.append(best_card)
        used_cubes.add(get_cube_index(best_card))
        remaining.remove(best_card)
    
    return selected

def count_internal_sets(cards: List[Tuple[int, ...]]) -> int:
    """Count number of internal SETs."""
    return len([combo for combo in combinations(cards, 3) if forms_set(list(combo))])

# Example usage
if __name__ == "__main__":
    deck = generate_deck(4)
    for n in [9, 12, 15, 18, 21, 24, 27]:
        selected = consecutive_maximization(deck, n)
        sets = count_internal_sets(selected)
        print(f"N={n}: Found {sets} internal SETs")