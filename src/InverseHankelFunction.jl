module InverseHankelFunction

import SpecialFunctions: hankelh1, gamma
import LambertW: lambertw

export invhankelh1, PassingThrough

export normalisedhankelh1,
       invnormalisedhankelh1,
       invnormalisedhankelh1_sortedvec

export hankel_arg_asymptotic_scale,
       small_arg_hankelh1,
       inv_small_arg_hankelh1,
       large_arg_hankelh1,
       inv_large_arg_hankelh1

include("diffhankel.jl")

include("invnormalisedhankel.jl")

include("hankel_asymptotics.jl")

include("invhankel.jl")

end # module
