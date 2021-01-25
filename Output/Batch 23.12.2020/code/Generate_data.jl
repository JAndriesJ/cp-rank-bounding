# This script is a colllection of completely positive matrices used for the testing of cp-rank computation.
using Dates # used in the output file naming
using Random # for the generation of a random cp Matrix
Random.seed!(1234)
include("MatUtils.jl")

# from Completely positive matrices Berman Shaked-Monderer  World scientific
M6 = makediagone([  8  12 16 4 6   8;
                    12 20 28 6 10 14;
                    16 28 40 8 14 20;
                    4   6  8 2  3  4;
                    6  10 14 3  5  7;
                    8  14 20 4  7 10])

# Matrices from:  I.M. Bomze, W. Schachinger, R. Ullrich. From seven to eleven: Completely positive matrices with high cp-rank. Linear Algebra and its Applications 459 (2014), 208 – 221.
M7       = circmatrix([531.0, 81, 150, 150, 81, 531, 926]/926)

M7tilde  = circmatrix([108, 27, 4, 4, 27, 108, 163.0]/163)

M8tilde  = makediagone([541.0 880 363 24 55 11 24 0;
                       880 2007 1496 363 48 22 22 24;
                       363 1496 2223 1452 363 24 22 11;
                       24 363 1452 2325 1584 363 48 55;
                       55 48 363 1584 2325 1452 363 24;
                       11 22 24 363 1452 2223 1496 363;
                       24 22 22 48 363 1496 2007 880;
                       0 24 11 55 24 363 880 541])

M9tilde  = makediagone([2548 1628 363 60 55 55 60 363 1628;
                       1628 2548 1628 363 60 55 55 60 363;
                       363 1628 2483 1562 363 42 22 55 60;
                       60 363 1562 2476 1628 363 42 55 55;
                       55 60 363 1628 2548 1628 363 60 55;
                       55 55 42 363 1628 2476 1562 363 60;
                       60 55 22 42 363 1562 2483 1628 363;
                       363 60 55 55 60 363 1628 2548 1628;
                       1628 363 60 55 55 60 363 1628 2548])

M11tilde = makediagone([781 0 72 36 228 320 240 228 36 96 0;
                             0 845 0 96 36 228 320 320 228 36 96;
                            72 0 827 0 72 36 198 320 320 198 36;
                            36 96 0 845 0 96 36 228 320 320 228;
                            228 36 72 0 781 0 96 36 228 240 320;
                            320 228 36 96 0 845 0 96 36 228 320;
                            240 320 198 36 96 0 745 0 96 36 228;
                            228 320 320 228 36 96 0 845 0 96 36;
                            36 228 320 320 228 36 96 0 845 0 96;
                            96 36 198 320 240 228 36 96 0 745 0;
                            0 96 36 228 320 320 228 36 96 0 845])


saveDir = "C:\\Users\\andries\\.julia\\dev\\CP-Rank-Bounding\\Data\\"

if false
    cp_matrices = [M6, M7, M7tilde, M8tilde, M9tilde, M11tilde]
    cp_matrices_names = ["M6", "M7", "M7tilde", "M8tilde", "M9tilde", "M11tilde"]
    saveSubDir = saveDir*"CPmats\\"
    for ind in 1:6
        saveMat2txt(cp_matrices[ind], saveDir, cp_matrices_names[ind])
    end
end

if false
    saveSubDir = saveDir*"CPrandMat\\"
    local count
    count =  0
    for  n in 6:9
        for ind in 1:10
            count = count + 1
            r = rand(4:30)
            M, a_ℓ = genCPmatrix(n,r)
            Mtext  = "($count)-$n Rcp_$r"
            saveMat2txt(M, saveSubDir, Mtext)
        end
    end
end

if true
    local count
    count = 0
    saveSubDir = saveDir*"DDSNN_Mat\\"
    for  n in 6:9
        for ind in 1:10
           count = count + 1
            M = genDDSNNmatrix(n)
            Mtext  = "($count)-$n DDSNN"
            saveMat2txt(M, saveSubDir, Mtext)
        end
    end
end

if true
    saveSubDir = saveDir*"SNN_Mat\\"
    local count
    count =  0
    for  n in 6:9
        for ind in 1:10
            count = count + 1
            M = genSNNmatrix(n)
            Mtext  = "($count)-$n SNN"
            saveMat2txt(M, saveSubDir, Mtext)
        end
    end
end
