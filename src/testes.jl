# version for a sorted vector of hns, reuses the z from the previous solve for the next one
function testes(hns::AbstractVector{<:Real})
    indices_gt1 = findall(x -> x> 1, hns)
    indices_le1 = findall(x -> x<=1, hns)

    sortperm_gt1 = sortperm(hns[indices_gt1])
    sortperm_le1 = sortperm(hns[indices_le1], rev=true)

    zs_gt1 = (hns[indices_gt1][sortperm_gt1])
    zs_le1 = (hns[indices_le1][sortperm_le1])

    zs = similar(hns, complex(eltype(hns)))

    zs[indices_gt1][sortperm_gt1] .= zs_gt1
    zs[indices_le1][sortperm_le1] .= zs_le1

    return zs
end

testes([3.4, 2.1, 0.9, 1.1])
