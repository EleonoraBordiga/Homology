
# Funzione che, dato un complesso simpliciale K e un intero p, crea una base (ordinata) dello 
# spazio vettoriale delle p-catene C_p(K) (ovvero un vettore contenente tutti i p-simplessi ordinati di K).
# Scriveremo la matrice bordo rispetto a queste basi di C_p(K) e di C_(p-1)(K).  

function basis(K::AbstractSet,p::Int)
    return [sort!(collect(s)) for s in K if length(s) == p+1]
end


# Funzione che, dato un complesso simpliciale K e un intero p, conta il numero di p-simplessi di K,
# ovvero la dimensione dello spazio vettoriale delle p-catene C_p(K).

dimPchains(K::AbstractSet,p::Int) = length(basis(K,p))


# Funzione che, dato un complesso simpliciale K e un suo p-simplesso ordinato s, crea un vettore v le cui
# componenti sono i coefficienti, rispetto alla base delle p-1 catene p_basis(K,p-1), del bordo della catena elementare
# associata ad s.
# Useremo questa funzione per costruire le colonne della matrice dp.

function linearComb(K::AbstractSet,s::Vector)
    p = length(s)-1
    v = Vector(undef, dimPchains(K,p-1))
    for i in 1:dimPchains(K,p-1)
        v[i] = get(boundary(Dict(s=>1)),basis(K,p-1)[i],0)
    end
    return v
end

# Funzione che, dato un complesso simpliciale K, calcola la sua dimensione, ovvero la massima dimensione 
# dei simplessi che lo compongono.

function dimSimplicialComplex(K::AbstractSet)
    return maximum([length(s)-1 for s in K])
end

# Funzione che, dato un complesso simpliciale K e un intero p, calcola la matrice
# associata all'operatore (lineare) bordo d_p:C_p(K) -> C_(p-1)(K).
# Nello spazio vettoriale di partenza consideriamo la base p_basis(K,p), mentre in quello di arrivo 
# p_basis(K,p-1).
# Inizialmente sono trattati, separatamente, i casi particolari di p per i quali l'operatore bordo 
# è nullo, ovvero p <= 0 o p >= k+1, dove indichiamo con k la dimensione del complesso simpliciale. 
# Nel caso in cui p < 0 o p > k+1 non viene fornita alcuna matrice.
# Tuttavia, per poter calcolare l'omologia in grado 0 e in grado k, abbiamo bisogno di una matrice
# bordo nulla delle dimensioni "giuste": d_0 dovrà essere moltiplicata a destra per d_1, quindi
# il numero di colonne di d_0 dovrà essere pari al numero di righe di d_1 (mentre il numero di righe 
# di d_0 è irrilevante, poniamo 1 in modo arbitrario).
# Analogamente d_(k+1) dovrà essere moltiplicata a sinistra per d_k, quindi il numero di righe di d_(k+1)
# dovrà essere pari al numero di colonne di d_k (per il numero di colonne vale la considerazione del caso
# precedente).

function dp(K::AbstractSet,p::Int)

  dim_p = dimPchains(K,p)
  dim_pminus1 = dimPchains(K,p-1)

  if p > dimSimplicialComplex(K)+1 || p < 0
        return println("Operatore nullo")

  elseif p == dimSimplicialComplex(K)+1
        return dp = zeros(Rational,dim_pminus1,1)

  elseif p == 0
        return dp = zeros(Rational,1,dimPchains(K,0))
  end
 
  dp = Matrix{Rational}(undef,dim_pminus1,dim_p) 
 
  for i in 1:dim_p
        dp[:,i] = linearComb(K,basis(K,p)[i])
  end
 
  return dp
end


# Esempio: matrici bordo del toro

dp(T,0)
dp(T,1)
dp(T,2)