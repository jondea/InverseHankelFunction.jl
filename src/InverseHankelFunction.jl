module InverseHankelFunction

import SpecialFunctions: hankelh1, gamma
import LambertW: lambertw

export invhankelh1, PassingThrough

export hankelh1n,
       HankelH1N

export diffhankelh1,
       diff2hankelh1,
       diff3hankelh1

export invhankelh1n,
       invhankelh1n_sortedvec,
       diffinvhankelh1n_sortedvec,
       invhankelh1n_adaptive_solve

export hankel_arg_asymptotic_scale,
       small_arg_hankelh1,
       inv_small_arg_hankelh1,
       large_arg_hankelh1,
       inv_large_arg_hankelh1

include("hankelh1n.jl")

include("diffhankel.jl")

include("invhankelh1n.jl")

include("hankel_asymptotics.jl")

include("invhankel.jl")

end # module
