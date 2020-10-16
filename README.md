# crystal
Proof of concept crystal algebra calculations



usage: sage crystal.sage [-h] -m mutation_method -t type -r rank [-k root] [-i inversions] [-n N] [word [word ...]]

positional arguments:

  word                  start calculations from word (needs to be reduced word corresponding to longest element of Weyl group)

required arguments:

  -m mutation_method, --method mutation_method
                        mutation method, one of: 20200916, 20201012, old
  -t type, --type type  type of algebra
  -r rank, --rank rank  rank of algebra
  -k root, --root root  root of algebra

optional arguments:

  -i inversions, --inversions inversions
                        process only decompositions with set number of inversions
  -n N, --num N         process only N random decompositions
