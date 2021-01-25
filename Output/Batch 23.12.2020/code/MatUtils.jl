using LinearAlgebra # For the matrix manipulations and preprocessing.
using Random # for the generation of a random cp Matrix
using Kronecker

############## read and write

"""Saves an input matrix to a specified path as a .txt-file """
function saveMat2txt(M, saveDir, saveName)
    open(saveDir*saveName*".txt", "w") do f
        n = size(M)[1]
        for i = 1:n
            rowStr = string(M[i,:])
            rowStr = replace(rowStr, "," => "")
            rowStr = replace(rowStr, "[" => "")
            rowStr = replace(rowStr, "]" => "")

            write(f,   rowStr*"  \n")
        end
    end
end

"""Loads from specified path of a .txt-file a square matrix."""
function loadMatfromtxt(loadPath)
    local n, Mflat
    open(loadPath, "r") do f
        count = 1
        Mflat = []
        for ln in eachline(f)
            splitln = split(ln)
            for Mijstr in splitln
                append!(Mflat, parse(Float64, Mijstr))
            end
            n = length(splitln)
        end
    end
    M = reshape(Mflat,(n,n))
    return M
end

############## Testing properties

function testNN(M)
    isNN = minimum(M) > 0
    return isNN
end


function testPSD(M)
    eigV = eigvals(M)
    isPSD1 = sum(imag(eigV)) < 1.0e-10
    isPSD2 = minimum(real(eigvals(M))) > - 1.0e-10
    return isPSD1 && isPSD1
end

function testDD(M)
    isDD = all(sum(M-Diagonal(diag(M)), dims = 2) .<= diag(M))
    return isDD
end

############## Generating

"""Vector input, output a Circulant matrix of said vector"""
function circmatrix(v)
    hcat([circshift(v, k) for k = -length(v)+1:0]...)
end

"""Normalization function"""
makediagone(A) = diagm(1 ./ sqrt.(diag(A))) * A * diagm( 1 ./ sqrt.(diag(A)))


""" Randomly generates a cp-matrix with size: n and cp-rank at most r.""" # Randomly generated matrices:
function genCPmatrix(n,r)
    A = zeros(Float64, (n,n))
    factor_dict = Dict()
    a_ℓ = Dict()
    for ℓ in 1:r
        a_ℓ[ℓ] = rand(Float64, (n, 1))
        #factor_dict[ℓ]
        A = A + a_ℓ[ℓ]* transpose(a_ℓ[ℓ])
    end
    A = makediagone(A)
    return A, a_ℓ
end


""" Generates a nonnegatvie symmetric matrix""" # every symetric nonnegative diagonally dominant matrix is cp:
function genSNNmatrix(n)
    B = rand(Float64, (n,n)) # nonnegative
    A = B + transpose(B) # symmetrize
    A = makediagone(A)
    return A
end

"""gives a diagonal offset to a matrix"""
function DominateDiagonal(A,λ)
    n = size(A)[1]
    A = A + λ *I(n)
    A = makediagone(A)
    return A
end

"""Keep offseting the diagonal until the matrix is PSD"""
function genDDSNNmatrix(n)
    A = genSNNmatrix(n)
    while ~testDD(A)
        A = A + 0.1*I(n)
    end
    A = makediagone(A)
    return A
end

################### Misc.:

function n2div4(M)
    return floor(size(M)[1]^2 / 4)
end

function extName(srttxt)
    return string(srttxt[1:(end-4)])
end





################### TO DO:


# Known negative examples:
""" Generates a tridiagonal matrix"""
function DNNclass(n)
    A = Tridiagonal(dl, d, du)
    A[(1,n)] = 1
    A[(n,1)] = 1
    return A
end

# CPSD examples :
function CPSDclass(n)
    dl = ones(1,n-1);
    du = dl;
    d = float32(rand(0:10, (1,n)));

     A = Tridiagonal(dl, d, du)
    A[(1,n)] = 1
    A[(n,1)] = 1
    return A
end

"""Generate bipartite matrix"""
function BipartiteGraphclass(p,q,a,b)
    Atop = hcat((a + q)*I(p), ones(p,q) )
    Abot = hcat(ones(q,p), (b + p)*I(q) )
    A = vcat(Atop, Abot)
    A = makediagone(A)
    return A
end

""" Generates cp-matrix via Kronecker products of cp-matrices with size: n and cp-rank at most r."""
function genVariantCPmatrix(n,r)
    A = zeros(Float64, (n^3,n^3))
    factor_dict = Dict()
    a_ℓ = Dict()
    b_ℓ = Dict()
    c_ℓ = Dict()
    d_ℓ = Dict()
    for ℓ in 1:r
        a_ℓ[ℓ] = rand(Float64, (n, 1))
        b_ℓ[ℓ] = rand(Float64, (n, 1))
        c_ℓ[ℓ] = rand(Float64, (n, 1))

        d_ℓ[ℓ] =  a_ℓ[ℓ], b_ℓ[ℓ], c_ℓ[ℓ]
        #factor_dict[ℓ]
        A = A + kron(kron(a_ℓ[ℓ]* transpose(a_ℓ[ℓ]), b_ℓ[ℓ]* transpose(b_ℓ[ℓ])), c_ℓ[ℓ]* transpose(c_ℓ[ℓ]))
    end
    A = makediagone(A)
    return A, d_ℓ
end
