# The decomposition tree of the cheapest-repair family, and the product cover
# of what does not decompose (src/Diagnostics.jl). Pure `F` combinatorics, no
# solver: build `fmin` as a Vector{Vector{Int}} of equal-size fact sets,
# decompose it, and check the tree shape; then cover a family that does not
# factor and check that every layer is exact and the layers are all of it.

using Resolver.Diagnostics: decompose, tree_shape, flatten_menus, product_menus,
    subsets, DecompNode, product_cover, cover_selections

# each ground `used` is the sorted union of the family's facts
ground(F) = sort!(unique!(reduce(vcat, F; init = Int[])))
# all k-subsets of 1:n, the uniform (threshold) family
ksubsets(n, k) = subsets(collect(1:n), k)

@testset "decompose: validated fixtures" begin
    shape(F) = tree_shape(decompose(F, ground(F)))

    # two disjoint edges — an OR of two AND-of-menus
    @test shape([[1, 2], [3, 4]]) ==
        "OR(AND(menu[1], menu[2]), AND(menu[3], menu[4]))"
    # K_{2,2} — a product of two choose-one menus
    @test shape([[1, 3], [1, 4], [2, 3], [2, 4]]) ==
        "AND(OR(menu[1], menu[2]), OR(menu[3], menu[4]))"
    # K_3 — every 2-subset of three, the tame threshold
    @test shape([[1, 2], [1, 3], [2, 3]]) == "T2of[1, 2, 3]"
    # K_4 — every 2-subset of four
    @test shape(ksubsets(4, 2)) == "T2of[1, 2, 3, 4]"
    # a threshold factor beside a menu factor
    @test shape([sort!(vcat(t, [m])) for t in [[1, 2], [1, 3], [2, 3]]
                 for m in [4, 5]]) ==
        "AND(T2of[1, 2, 3], OR(menu[4], menu[5]))"
    # C_5 — an odd cycle, genuinely incompressible
    @test shape([[1, 2], [2, 3], [3, 4], [4, 5], [1, 5]]) ==
        "PRIME[1, 2, 3, 4, 5]"
    # three singletons — one choose-one menu (an OR of singleton leaves)
    @test shape([[1], [2], [3]]) == "OR(menu[1], menu[2], menu[3])"

    # k > 2 thresholds — detection is |F| == binomial(N, k), general in k
    @test shape(ksubsets(4, 3)) == "T3of[1, 2, 3, 4]"
    @test shape(ksubsets(5, 3)) == "T3of[1, 2, 3, 4, 5]"
    # a k=3 threshold factor beside a menu factor
    @test shape([sort!(vcat(t, [m])) for t in ksubsets(4, 3) for m in [5, 6]]) ==
        "AND(T3of[1, 2, 3, 4], OR(menu[5], menu[6]))"
end

@testset "decompose: flatten agrees with product_menus" begin
    # where the family is a product of choose-one menus, the tree flattens to
    # exactly the factors `analyse` uses today; where it is not, product_menus
    # returns nothing and flatten declines too (threshold/prime/OR-of-edges)
    for F in ([[1, 3], [1, 4], [2, 3], [2, 4]],   # K_{2,2}, a product
              [[1], [2], [3]],                     # one menu
              [[1, 3, 5], [1, 3, 6], [2, 3, 5], [2, 3, 6]])  # three-menu product
        used = sort!(unique!(reduce(vcat, F; init = Int[])))
        factors = product_menus(F, used)
        @test factors !== nothing
        @test flatten_menus(decompose(F, used)) == factors
    end
    # non-products: product_menus declines and so does flatten
    for F in ([[1, 2], [1, 3], [2, 3]],            # K_3, a threshold
              [[1, 2], [2, 3], [3, 4], [4, 5], [1, 5]],  # C_5, prime
              [[1, 2], [3, 4]])                    # OR of two edges
        used = sort!(unique!(reduce(vcat, F; init = Int[])))
        @test product_menus(F, used) === nothing
        @test flatten_menus(decompose(F, used)) === nothing
    end
end

@testset "cover: exact products, layer by layer" begin
    cover(F) = product_cover([sort(m) for m in F], ground(F))
    # every selection of every layer, as the members it claims to be
    reached(layers) = Set(cover_selections(layers))

    # a product is one layer and is itself: nothing to cover, nothing left
    K22 = [[1, 3], [1, 4], [2, 3], [2, 4]]
    @test cover(K22) == [[[[1], [2]], [[3], [4]]]]

    # the path :1–:3–:2–:4, prime: the 2×1 rectangle {1,2}×{3} is the largest
    # there is, and the member it misses is the layer after it — as two menus
    # of one entry, since nothing couples them: a coupling is made only where
    # it reaches a member the split cannot
    P4 = [[1, 3], [2, 3], [2, 4]]
    @test cover(P4) == [[[[1], [2]], [[3]]], [[[2]], [[4]]]]

    # a coupled pair — the family holds three of the four combinations of
    # {1,2}×{3,4} — is one menu of compound options beside a menu of singles,
    # which is what a factor the family couples has to print as
    coupled = [[1, 3, 5], [1, 4, 5], [2, 3, 5],
               [1, 3, 6], [1, 4, 6], [2, 3, 6]]
    @test cover(coupled) == [[[[1, 3], [1, 4], [2, 3]], [[5], [6]]]]

    # every layer is exact and the layers are complete, on each of these and
    # on the k=3 threshold besides: what a cover offers is the family and
    # nothing else
    for F in (K22, P4, coupled, ksubsets(4, 3), ksubsets(5, 3),
              [[1, 2], [2, 3], [3, 4], [4, 5], [1, 5]])
        members = Set(Vector{Int}[sort(m) for m in F])
        layers = cover(F)
        sels = cover_selections(layers)
        @test allunique(sels)          # no repair said twice
        @test Set(sels) == members     # exact, and all of it
    end
end
