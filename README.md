# set-search
This repo is the codebase for the paper:  **Investigation into Maximal Number of Sets for 12 Cards**

It offers multiple implementations for exploring the space of possible `n` card arrangements to find the maximal number of sets for the card game Set. 

## `set_searcher_complete.py` 
Iterates over all possible configurations of boards for `n` cards with four properties (i.e. the normal Set deck) looking for the maximal number of sets. The script prints the board and the max sets found upon completion. 

## `set_search_3d_simple.py`
Iterates over all possible configurations of boards for `3 <= n < 27` cards with three properties looking for the maximal number of sets. The script prints the board and the max sets found upon completion.

## `set_search_3d_table.py`
Code to reproduce Table 3 in the paper.

## `set_search_3d_magic_square.py`
Computer check for Appendix B of the paper: for `n = 10, 11, 12`, every maximizing configuration of `n` cards with three properties contains a complete magic square (one of the 39 affine planes of F_3^3). For each `n`, the script enumerates the `n`-card configurations, counts their sets, checks whether they contain any of the 39 planes, and records the maximum number of sets overall and among configurations with no magic square.

By default it uses a fixed-set symmetry reduction: affine maps of F_3^3 send sets to sets and planes to planes and can move any set onto any other, so it suffices to enumerate the `C(24, n-3)` configurations containing the set `{(0,0,0), (0,0,1), (0,0,2)}`. No maximizing configuration is lost. This takes about 10 seconds. Pass `--full` to enumerate all `C(27, n)` configurations instead (about two minutes for all three; same maxima). Requires only Python 3.

## `cmm.py` 
Code for running the Consecutive Maximization Algorithm from [Vinci09](https://www.setgame.com/sites/default/files/teacherscorner/SETPROOF.pdf).



