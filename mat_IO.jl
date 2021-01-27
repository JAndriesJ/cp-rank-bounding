## file management: read and write
"""Check if dir exitst if not make one."""
function check_save_dir(save_dir, save_sub_dir_name)
    saveSubDir = save_dir*save_sub_dir_name*"\\"
    if ~isdir(saveSubDir)
        mkdir(saveSubDir)
    end
    return saveSubDir
end

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
    M = real(reshape(Mflat,(n,n)))
    return M
end


"""cuts every thing in the string after the '.' """
function cut_ext(srttxt)
    ext_ind = findlast(isequal('.'), srttxt)
    return string(srttxt[1:(ext_ind-1)])
end


##

## Testing properties of the matrix
"""Test if matrix is nonnegative"""
function testNN(A)
    isNN = minimum(A) > -1.0e-10
    return isNN
end

"""Test if matrix M is PSD"""
function testPSD(A)
    eigV = eigvals(A)
    isPSD1 = sum(imag(eigV)) < 1.0e-10
    isPSD2 = minimum(real(eigV)) > - 1.0e-10
    return isPSD1 && isPSD1
end
