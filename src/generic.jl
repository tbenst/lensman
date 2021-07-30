function vecvec2mat(vecvec)
    dim1 = length(vecvec)
    dim2 = length(vecvec[1])
    my_array = zeros(Int64, dim1, dim2)
    for i in 1:dim1
        for j in 1:dim2
            my_array[i,j] = vecvec[i][j]
        end
    end
    return my_array
end

"apply each function to args"
funprod(functions...) = (args...)->map(x->x(args...), functions)

zeroToOne(x) = (x.-minimum(x))/(minimum(x) + maximum(x))