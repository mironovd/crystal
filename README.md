# crystal
Proof of concept crystal algebra calculations

- crystal.sage - calculation of crystal half-potential 
- crystal-compare.sage - comparison of half-potential calculation methods (can be used to batch compare 
all words for given algebra type and rank)
- crystal_GHKK_support.sage - calculation of X-cluster Gross-Hacking-Keel-Kontsevich potential support and checking if Newton polytopes are void

### USAGE
    usage: sage crystal_GHKK_support.sage [-h] -m mutation_method -t type -r rank [-k root] [-i inversions] [-n N] [-p] [-P] [-v] [word [word ...]]

    Crystal calculations

    positional arguments:
        word                  start calculations from word 
                                (needs to be reduced word corresponding to longest element of Weyl group)

    optional arguments:
        -h, --help            show this help message and exit

    required arguments:
        -m mutation_method, --method mutation_method
                mutation method, one of:  new, old
                        new - mutation algorithm for A_n, B_n, C_n, D_n
                        old - simple mutation algorithm for A_n

        -t type, --type type  type of algebra
        -r rank, --rank rank  rank of algebra
        -k root, --root root  root of algebra

    optional arguments:
        -n N, --num N         process only N random decompositions
        -p, --polymake        use Polymake to check Newton polytopes
        -P, --no-polymake     do not use Polymake to check Newton polytopes
        -v, --verbose         use Polymake to check Newton polytopes


