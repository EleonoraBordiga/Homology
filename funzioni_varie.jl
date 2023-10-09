# Funzioni definite in precedenza che verranno richiamate


# Facce di un simplesso 

delete(s::AbstractSet, x) = delete!(copy(s), x)

push(s::AbstractSet, x...) = push!(copy(s), x...)
push(s::AbstractSet) = s

function faces(s::AbstractSet)
    (delete(s,v) for v ∈ s)
end

function allsubfaces(simplexes,acc)
    if isempty(simplexes)
        for t in acc
          if length(t)==0
             delete!(acc,t)
          end
        end
        return acc
    else
        s,ss... = simplexes
        f = faces(s)
        return allsubfaces((f...,ss...),push(acc,f...))
    end
end

allfaces(s) = allsubfaces((s,),Set((s,)))


function SimplicialComplex(simplexes::Tuple)
    S = Set()
    [push!(S,s) for t in simplexes for s in allfaces(t)]
    return S
end

# Complesso simpliciale del toro 

T = SimplicialComplex(( 
                        Set([:a,:b,:f]),
                        Set([:b,:c,:f]),
                        Set([:a,:c,:g]),
                        Set([:a,:d,:f]),
                        Set([:c,:f,:g]),
                        Set([:a,:d,:g]),
                        Set([:d,:e,:f]),
                        Set([:f,:g,:h]),
                        Set([:d,:g,:j]),
                        Set([:e,:f,:i]),
                        Set([:f,:h,:i]),
                        Set([:h,:i,:j]),
                        Set([:g,:h,:j]),
                        Set([:d,:e,:j]),
                        Set([:a,:e,:i]),
                        Set([:b,:i,:j]),
                        Set([:a,:e,:j]),
                        Set([:a,:b,:i]),
                        Set([:b,:c,:j]),
                        Set([:a,:c,:j]) 
))

# Complesso simpliciale del nastro di Moebius

M = SimplicialComplex((
                       Set([:a,:d,:e]),
                       Set([:a,:b,:e]),
                       Set([:b,:c,:e]),
                       Set([:c,:e,:f]),
                       Set([:c,:d,:f]),
                       Set([:a,:d,:f])
))

# Complesso simpliciale della bottiglia di Klein

K = SimplicialComplex(( 
                        Set([:a,:b,:f]),
                        Set([:b,:c,:f]),
                        Set([:a,:c,:g]),
                        Set([:a,:d,:f]),
                        Set([:c,:f,:g]),
                        Set([:a,:e,:g]),
                        Set([:d,:e,:f]),
                        Set([:f,:g,:h]),
                        Set([:e,:g,:j]),
                        Set([:e,:f,:i]),
                        Set([:f,:h,:i]),
                        Set([:h,:i,:j]),
                        Set([:g,:h,:j]),
                        Set([:d,:e,:j]),
                        Set([:a,:e,:i]),
                        Set([:b,:i,:j]),
                        Set([:a,:d,:j]),
                        Set([:a,:b,:i]),
                        Set([:b,:c,:j]),
                        Set([:a,:c,:j]) 
))

# Complesso simpliciale del piano proiettivo

P = SimplicialComplex(( 
                        Set([:a,:b,:g]),
                        Set([:b,:c,:g]),
                        Set([:c,:d,:h]),
                        Set([:a,:f,:g]),
                        Set([:c,:g,:h]),
                        Set([:d,:e,:h]),
                        Set([:e,:f,:g]),
                        Set([:g,:h,:i]),
                        Set([:e,:h,:k]),
                        Set([:e,:g,:j]),
                        Set([:g,:i,:j]),
                        Set([:i,:j,:k]),
                        Set([:e,:f,:k]),
                        Set([:e,:d,:j]),
                        Set([:c,:j,:k]),
                        Set([:a,:f,:k]),
                        Set([:c,:d,:j]),
                        Set([:b,:c,:k]),
                        Set([:a,:b,:k]),
                        Set([:h,:i,:k]) 
))


# Bordo di una catena

function addchains(c...)
    d = mergewith(+,c...)
    [pop!(d,k) for (k,v) in pairs(d) if iszero(v)]
    return d
end

delete(s::Vector,i) = deleteat!(copy(s),i)
 
function boundary(c::Dict)
  d = Dict()
  for s in keys(c)               
    v = get(c,s,0)                
    for i in 1:length(s)
        ŝ = delete(s,i)           
        if iseven(i)              
            d = addchains(d, Dict((ŝ) => -v)) 
        else 
            d = addchains(d,Dict((ŝ) => v))
        end
    end
  end
  return d
end


