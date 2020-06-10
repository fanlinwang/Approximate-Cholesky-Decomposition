## Baseline
#### V=1000, E=5000
get column(L137-141): 11.76%
random number init: 5.88%
get column (merge/sort)(L153-162): 35.28%
compute(L167-185): 23.52%
add new edge (L187-207): 11.76%
#### V=10000, E=160000
get column(L137-141): 8.17%
random number init: 4.9%
get column (merge/sort)(L153-162): 22.21%
compute(L167-185): 30.72%
add new edge (L187-207): 32.61%

## VerMg
#### V=1000, E=5000
merge and sort(L400-411): 5%
sort(L392+L415): 15%
random number init: 5%
get column (L420-445): 15%, search and insert (L440-444): 10.52%
compute(L450-468): 25%
add new edge(L470-484): 25%
#### V=10000, E=160000
merge and sort (L400-411): 22.24%
sort(L392+L415): 1%
random number init: 4.81%
get column (L420-445): 26.71%, search and insert (L440-444): 9.26%
compute(L450-468): 24.66%
add new edge(L470-484): 17.52%

## VerMgSearch
#### V=1000, E=5000
bitset_search (L606-608): 9.09%
#### V=10000, E=160000
bitset_search (L606-608): 7.96%

## VerMgSearchSIMD
#### V=1000, E=5000
merge and sort(L1215-1227): 21.04%
sort(L1207+L1231): 0%
compute(L1270-1287): 21.05% (binsearch 0%)
SIMD(L1301-1306): 10.52%
SIMD(L1310-1342): 21.04%
add new edge(L1345-1353): 21.04%

#### V=10000, E=160000
merge and sort(L1215-1227): 28.84%
sort(L1207+L1231): 0.63%
compute(L1270-1287): 30.55% (binsearch 1.68%)
SIMD(L1301-1306): 3.36%
SIMD(L1310-1342): 19.33%
add new edge(L1345-1353): 13.68%


