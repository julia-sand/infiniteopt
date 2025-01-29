This library runs a direct numerical integration model using InfiniteOpt: an optimisation library for julia

There are two main cases covered: 
1. KL: minimising the Kullback-Leibler divergence
2. EP: minimising the mean entropy production

To run a model and save a csv of the results use 

>julia --project=. -e 'using Pkg; Pkg.instantiate()'
>
>julia --project=. infiniteopt/kl/model_setup.jl 

use --h to check the passable



