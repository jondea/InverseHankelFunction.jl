
"""
    diffhankelh1(ν, z [, h[, hm1]])

Returns the first derivative of hankelh1 with respect to `z`.
Can optionally pass in a precomputed `H_ν(z)` and `H_{ν-1}(z)`.
"""
diffhankelh1(ν, z, h=hankelh1(ν, z), hm1=hankelh1(ν-1, z)) = hm1 - ν/z*h


"""
    diff2hankelh1(ν, z [, h[, hm1]])

Returns the first and second derivatives of hankelh1 with respect to `z` as tuple
Can optionally pass in a precomputed `H_ν(z)` and `H_{ν-1}(z)`.
"""
function diff2hankelh1(ν, z, h=hankelh1(ν, z), hm1=hankelh1(ν-1, z))
    hm2 = hankelh1(ν-2, z)
    hp = hm1 - ν/z*h
    hm1p = hm2 - (ν-1)/z*hm1
    hpp = hm1p - ν/z*hp + ν/z^2*h
    return (hp, hpp)
end

"""
    diff3hankelh1(ν, z [, h[, hm1]])

Returns the first, second and third derivatives of hankelh1 with respect to `z` as tuple
Can optionally pass in a precomputed `H_ν(z)` and `H_{ν-1}(z)`.
"""
function diff3hankelh1(ν, z, h=hankelh1(ν, z), hm1=hankelh1(ν-1, z))
    hm2 = hankelh1(ν-2, z)
    hm3 = hankelh1(ν-3, z)
    hp = hm1 - ν/z*h
    hm1p = hm2 - (ν-1)/z*hm1
    hm2p = hm3 - (ν-2)/z*hm2
    hpp = hm1p - ν/z*hp + ν/z^2*h
    hm1pp = hm2p - (ν-1)/z*hm1p + (ν-1)/z^2*hm1
    hppp = hm1pp - ν/z*hpp + 2ν/z^2*hp - 2ν/z^3*h
    return (hp, hpp, hppp)
end
