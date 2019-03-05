module FarooqHamidInverse

#=
This is a matrix invese algorithm based on Ahmad Farooq and Khan Hamid of King Khalid University 
2010 paper titled An efficient and simple algorithm for Matrix Inversion. 
The algorithm by itself has a time complexity of O(n^3). However it is extremely easy to 
implement and is faster than Gauss-Jordan elimination in matrices of size less than 500x500. 

Other advantage of this algorithm is that it is capable of finding pseudo-inverse of non square
matrix. 
=#

function fhInv(A)
    B = A
    det = 1.0
    siz = size(B)[1]
    for p in 1:siz
        pivot = B[p,p]
        det = det * pivot
        for i in 1:siz
            B[i, p] =  -B[i, p]/pivot
        end
        for i in 1:siz
            if i!=p
                for j in 1:siz
                    if j!=p
                        B[i, j] = B[i, j] + B[p, j]*B[i,p]
                    end
                end
            end
        end
        for i in 1:siz
            B[p, i] = B[p, i]/pivot
        end
        B[p, p] = 1/pivot
    end
    return B
end

function fhInv!(B)
    det = 1.0
    siz = size(B)[1]
    for p in 1:siz
        pivot = B[p,p]
        det = det * pivot
        for i in 1:siz
            B[i, p] =  -B[i, p]/pivot
        end
        for i in 1:siz
            if i!=p
                for j in 1:siz
                    if j!=p
                        B[i, j] = B[i, j] + B[p, j]*B[i,p]
                    end
                end
            end
        end
        for i in 1:siz
            B[p, i] = B[p, i]/pivot
        end
        B[p, p] = 1/pivot
    end
    return B
end