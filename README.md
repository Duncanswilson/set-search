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
Checks the computer-assisted step in Appendix B of the paper: for `n = 10, 11, 12` cards with three properties, every configuration with the maximum number of sets contains a magic square. For each `n`, it enumerates all `n`-card configurations, counts their sets, checks whether they contain any of the 39 magic squares (affine planes), and records the maximum number of sets overall and among the configurations with no magic square.

## `cmm.py` 
Code for running the Consecutive Maximization Algorithm from [Vinci09](https://www.setgame.com/sites/default/files/teacherscorner/SETPROOF.pdf).



