# Operazioni elementari su righe/colonne di una matrice

function rowSwap(A,i,j)
    A[i,:], A[j,:] = A[j,:], A[i,:]
end

function colSwap(A,i,j)
    A[:,i], A[:,j] = A[:,j], A[:,i]
end

function scaleCol(A,i,c)
   A[:,i] = c.*A[:,i] 
end

function scaleRow(A,i,c)
   A[i,:] = c.*A[i,:]
end
 
function colCombine(A,i,j,c)
   A[:,i] += c.* A[:,j]
end

function rowCombine(A,i,j,c)
   A[i,:] += c.*A[j,:]
end
# Funzione che riduce a scala una matrice A per colonne e corrispondentemente B per righe.
# Ad ogni operazione elementare su una colonna di A, viene effettuata la medesima operazione sulla riga di B 
# con lo stesso indice.

function simultaneousReduce(A::Matrix,B::Matrix)
   if size(A,2) != size(B,1)
      error("Le matrici hanno le dimensioni sbagliate.")
   end
 
   numRows, numCols = size(A) 
 
   i,j = 1,1
   while true
      if i > numRows || j > numCols
         break
      end
 
      if A[i,j] == 0
         nonzeroCol = j
         while nonzeroCol <= numCols && A[i,nonzeroCol] == 0
            nonzeroCol += 1
         end
 
         if nonzeroCol == numCols+1
             i += 1
            continue
         end
         colSwap(A, j, nonzeroCol)
         rowSwap(B, j, nonzeroCol)
      end
 
      pivot = A[i,j]
      scaleCol(A, j, 1//1 / pivot)
      scaleRow(B, j, 1//1 / pivot)
 
      for otherCol in range(1, numCols)
         if otherCol == j
            continue
         end

         if A[i, otherCol] != 0
            c = -A[i, otherCol]
            colCombine(A, otherCol, j, c)
            rowCombine(B, j, otherCol, -c)
         end
      end
      i += 1; j+= 1
   end

   return A,B
end

# Funzione che riduce a scala per righe una matrice B (la chiamiamo dopo simultaneousReduce(A,B) perchè quest'ultima
# funzione termina appena la matrice A è ridotta in forma a scala, mentre B potrebbe non esserlo ancora).

function rowReducing(B)

   numRows, numCols = size(B)

    i, j = 1, 1
    while true
        if i > numRows || j > numCols
            break
        end

        if B[i, j] == 0
            nonzeroRow = i
            while nonzeroRow <= numRows && B[nonzeroRow, j] == 0
                nonzeroRow += 1
            end

            if nonzeroRow == numRows+1
                j += 1
                continue
            end

            rowSwap(B, i, nonzeroRow)
         end

        pivot = B[i, j]
        scaleRow(B, i, 1//1 / pivot)

        for otherRow in range(1,numRows)
            if otherRow == i
                continue
            end
 
            if B[otherRow, j] != 0
                c = -B[otherRow, j]
                rowCombine(B, otherRow, i, c)
            end
         end

        i += 1; j += 1
      end

    return B
end

# Esempio

K = SimplicialComplex((Set([0,1,2]),Set([1,2,3])))
B = dp(K,2)

A = dp(K,1)

# Funzione che conta il numero di colonne non nulle di una matrice A.

function numPivotCols(A)
   z = zeros(size(A,1))
   return size(A,2) - count([all(A[:, j] == z) for j in range(1,size(A,2))])
end


# Funzione che conta il numero di righe non nulle di una matrice A.

function numPivotRows(A)
   z = zeros(size(A,2))
   return size(A,1) - count([all(A[i, :] == z) for i in range(1,size(A,1))])
end

# Funzione che, dato un complesso simpliciale K e un intero p, calcola il p-esimo numero di Betti di K.

function bettiNumber(K::AbstractSet,p::Int)

   if p < 0 || p > dimSimplicialComplex(K)
      return 0
   end

   A, B = copy(dp(K,p)), copy(dp(K,p+1))

   simultaneousReduce(A,B)
   rowReducing(B)
 
   kernelDim = dimPchains(K,p) - numPivotCols(A)
   imageDim = numPivotRows(B)
 
   return kernelDim - imageDim
end

# Esempi:

# Toro

for p in 0:2
println("numero di Betti del toro #$p:\n", bettiNumber(T,p), "\n")
end

# Nastro di Moebius

for p in 0:2
   println("numero di Betti del nastro di Moebius #$p:\n", bettiNumber(M,p), "\n")
end

# Bottiglia di Klein

for p in 0:2
   println("numero di Betti della bottiglia di Klein #$p:\n", bettiNumber(K,p), "\n")
end

# Piano proiettivo

for p in 0:2
   println("numero di Betti del piano proiettivo #$p:\n", bettiNumber(P,p), "\n")
end