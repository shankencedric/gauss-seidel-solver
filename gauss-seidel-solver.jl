using LinearAlgebra

function getParameters()
    println("Square matrix size n:")
    n = parse(Int, readline())

    println("Matrix elements (elements separated by ' ', enter for next row):")
    
    matrixA = Matrix{Float64}(undef, n, n)
    for row_idx in 1:n
        row_elems_str = split(strip(readline()), " ")
        if (length(row_elems_str) != n) 
            println("Row is not of length n. Try again.\n")
            return getParameters()
        end
        row_elems = [parse(Float64, elem) for elem in row_elems_str]
        matrixA[row_idx, : ] = row_elems
    end

    println("Output vector (n x 1) elements, separated by ' ':")
    b_elems_str = split(strip(readline()), " ")
    if (length(b_elems_str) != n) 
        println("Vector is not of length n. Try again.\n")
        return getParameters()
    end
    vectorb = [parse(Float64, elem) for elem in b_elems_str]

    println("Vector name of desired ('x', 'w', etc.):")
    xName = readline()

    println("Initial guess for $xName:")
    x0_str = split(strip(readline()), " ")
    if (length(x0_str) != n) 
        println("Vector is not of length n. Try again.\n")
        return getParameters()
    end
    x0 = [parse(Float64, elem) for elem in x0_str]

    println("Max iterations (default is k = 1000):")
    k_input = readline()
    k = length(inp) > 0 ? parse(Int, k_input) : 1000

    println("Max tolerance (default is '-3' for 10^-3):")
    tol_input = readline()
    tol = 10^(length(tol_input) > 0 ? parse(Int, tol_input) : -3)

    return [n, matrixA, xName, vectorb, x0, k, tol]
end

function gaussSeidel(n::Int, matrixA::Matrix, xName::String, vectorb::Vector, x0::Vector, k::Int, tol::Float64)
    x_k = deepcopy(x0)

    for k_idx in 1:k
        x_kminus1 = deepcopy(x_k) 
        
        for i in 1:n
            sum1 = 0
            for j in 1:i-1
                sum1 += matrixA[i, j] * x_k[j]
            end
            sum2 = 0
            for j in i+1:n
                sum2 += matrixA[i, j] * x_k[j]
            end
            x_k[i] = (vectorb[i] - sum1 - sum2) / matrixA[i, i]
        end
        
        # Check for convergence
        if (norm(x_k - x_kminus1, Inf) / norm(x_k, Inf)) < tol
            println("Converged in $k_idx iterations.\n")
            for i in 1:n
                print("$xName")
                print("$i = ")
                println(x_k[i])
            end
            println("$xName = $x_k")
            return x_k
        end
    end
    
    println("Maximum iterations reached.\n")
    for i in 1:n
        print("$xName")
        print("$i = ")
        println(x_k[i])
    end
    println("$xName = $x_k")
    return x_k
end

# Parameters (can call getParameters instead but nah)
n = 3
A = [1 1e-4 1e-16; 11e-4 -5 11e-4; 1e-16 1e-4 1]
xName = "w"
b = [1.0; 0.0; 1.0]
x0 = [1.0; 1.0; 1.0]
k = 1000
tol = 10^-1000

# Gauss Seidel execution
gaussSeidel(n, A, xName, b, x0, k, tol)
