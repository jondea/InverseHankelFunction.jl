
"""
    diffhankelh1(ν, z)

Returns the first derivative of hankelh1 with respect to `z`
"""
diffhankelh1(ν, z) = hankelh1(ν-1, z) - ν/z*hankelh1(ν, z)


"""
    diff2hankelh1(ν, z)

Returns the first and second derivatives of hankelh1 with respect to `z` as tuple
"""
function diff2hankelh1(ν, z)
    h = hankelh1(ν, z)
    hm1 = hankelh1(ν-1, z)
    hm2 = hankelh1(ν-2, z)
    hp = hm1-ν/z*h
    hm1p = hm2-(ν-1)/z*hm1
    hpp = hm1p -ν/z*hp + ν/z^2*h
    return (hp, hpp)
end

"""
    diff3hankelh1(ν, z)

Returns the first, second and third derivatives of hankelh1 with respect to `z` as tuple
WARNING: BROKEN
"""
function diff3hankelh1(ν, z)
    h = hankelh1(ν, z)
    hm1 = hankelh1(ν-1, z)
    hm2 = hankelh1(ν-2, z)
    hm3 = hankelh1(ν-3, z)
    hp = hm1-ν/z*h
    hm1p = hm2-(ν-1)/z*hm1
    hpp = hm1p -ν/z*hp + ν/z^2*h
    hppp = 1.0
    return (hp, hpp, hppp)
end
