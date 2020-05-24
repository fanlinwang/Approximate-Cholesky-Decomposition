using Laplacians
using SuiteSparse
using LinearAlgebra
using SparseArrays
using Statistics
import Laplacians.LLMatOrd
import Laplacians.approxChol
import Laplacians.LDLsolver
import Laplacians.approxchol_lapGiven
import Laplacians.LLord
import Laplacians.LDLinv
import Laplacians.LLcol

function generate(n::Integer, m::Integer, uniform::Bool)
    s = Set(2:n)
    t = Set(1)
    edges = Set()
    cur_vertex = 1;
    while  !(isempty(s))
        num = rand(1:n)
        if (num in s)
            pop!(s, num)
            push!(t, num)
            push!(edges, tuple(cur_vertex, num))
        end
        cur_vertex = num
    end
    
    for i in 1:m-n+1
        u = rand(1:n)
        v = rand(1:n)
        while tuple(u,v) in edges || tuple(v,u) in edges || u==v
            u = rand(1:n)
            v = rand(1:n)
        end
        push!(edges, tuple(u,v))
    end
    
    u=[];v=[];w=[]
    for edge in edges
        if uniform
            push!(w, 1.0, 1.0)
        else
            weight = rand(Float64)
            push!(w, weight, weight)
        end
        push!(u, edge[1], edge[2])
        push!(v, edge[2], edge[1])
    end
#     println(u,v,w)
#     println(typeof(w))
    mat = sparse(u,v,w,n,n)
    return convert(SparseMatrixCSC{Float64, Int64}, mat)
end

a = generate(parse(Int, ARGS[1]), parse(Int, ARGS[2]), true)
llmat = LLMatOrd(a)
@time approxChol(llmat)
# @benchmark approxChol(llmat)

