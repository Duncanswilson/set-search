# set-search
This repo is the codebase for the paper:  **Investigation into Maximal Number of Sets for 12 Cards**

It offers multiple implementations for exploring the space of possible `n` card arrangements to find the maximal number of sets for the card game Set. 

## `set_searcher_complete.py` 
Iterates over all possible configurations of boards for `n` cards with four properties (i.e. the normal Set deck) looking for the maximal number of sets. The script prints the board and the max sets found upon completion. 

## `set_search_3d_simple.py`
Iterates over all possible configurations of boards for `3 <= n < 27` cards with three properties looking for the maximal number of sets. The script prints the board and the max sets found upon completion.

## `set_search_3d_table.py`
Code to reproduce Table 3 in the paper.

## `cmm.py` 
Code for running the Consecutive Maximization Algorithm from [Vinci09](https://www.setgame.com/sites/default/files/teacherscorner/SETPROOF.pdf).



