using Observables
##
a = Observable([2])
y = map(a) do val
    val.+1
end

on(y) do val
    println("y changed to $val")
end

on(y) do val
    println("nuh-uh")
end

a[] .+= 1