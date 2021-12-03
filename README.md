# crystal
Proof of concept crystal algebra calculations

- crystal.sage - calculation of crystal half-potential 
- crystal-compare.sage - comparison of half-potential calculation methods (can be used to batch compare 
all words for given algebra type and rank)
- crystal_[abcd]n.sage - calculation of crystal half-potential using determinant method (slow algorithm)
- crystal_[abcd]n_comp.sage - comparison of combinatorial method half-potential calculation with determinant method
- crystal_newton_poly.sage - calculation of Newton polytope for half-potential
- crystal_GHKK_support.sage - calculation of X-cluster Gross-Hacking-Keel-Kontsevich potential support

### USAGE
    sage crystal.sage [-h] -m mutation_method -t type -r rank [-k root] [-i inversions] [-n N] [word [word ...]]

    positional arguments:

	word                  start calculations from word (needs to be reduced word corresponding to longest element of Weyl group)

    required arguments:

	-m mutation_method, --method mutation_method 
		mutation method, one of: 
			new - mutation algorithm for A_n, B_n, C_n, D_n
			old - simple mutation algorithm for A_n
                        
	-t type, --type type  type of algebra
  
	-r rank, --rank rank  rank of algebra
  
	-k root, --root root  root of algebra

    optional arguments:

	-i inversions, --inversions inversions  
                        process only decompositions with set number of inversions
                        
	-n N, --num N         process only N random decompositions
  

