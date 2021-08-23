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

"""
Dictionary puns for variable assignment.
    
mydict = Dict(:x=>1,:y=>2,:z=>3)
@pun (x, y, z) = mydict
@pun x = mydict

expands to
(x, y, z) = (mydict[:x], mydict[:y], mydict[:z])
x = mydict[:x]
"""
macro pun(ex)
    @assert ex.head == :(=) "no `=` found in expression."
    vars = ex.args[1] # :((x, y, z))
    d = ex.args[2] # "mydict"
    if typeof(vars)==Symbol
        var_sym = QuoteNode(vars)
        esc(:($vars = $d[$var_sym]))
    else # Expr of Array
        vars_sym = map(Symbol ∘ string, vars.args) # [:x, :y, :z]
        esc(:($vars = [$d[k] for k in $vars_sym]))
    end
end

"""
Variable puns for dictionary assignment.
    
mydict = Dict()
w = 0; x = 1; y = 2; z = 3
@assign mydict = (x, y, z)
@assign mydict = w

expands to
(mydict[:x], mydict[:y], mydict[:z]) = (x, y, z)
mydict[:w] = w
"""
macro assign(ex)
    @assert ex.head == :(=) "no `=` found in expression."
    d = ex.args[1] # "mydict"
    vars = ex.args[2] # :((x, y, z))
    if typeof(vars)==Symbol
        esc(Meta.parse("$d[:$vars] = $vars"))
    else # Expr of Array
        vars_sym = map(Symbol ∘ string, vars.args) # [:x, :y, :z]
        left = "("
        for k in vars_sym
            left *= "$d[:$k], "
        end
        left *= ")"
        left = Meta.parse(left)
        esc(:($left = $vars))
    end
end

"No operation."
function noop(args...; kwargs...)
    nothing
end