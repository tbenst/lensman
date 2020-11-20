using Transducers
using Base: iterate

struct IntStream
    length::Int
    sleepTime::Real
end

function Base.iterate(iter::IntStream)
    sleep(iter.sleepTime)
    return (1,2)
end

function length(is::IntStream)
    return is.length
end

function Base.iterate(iter::IntStream, state)
    if state > iter.length
        return nothing
    else
        sleep(iter.sleepTime)
        return (state,state+1)
    end
end

x = IntStream(10,0.1)

function printPass(x)
    println(x)
    x
end

# is there a no-op alternative to collect that still causes this to run?
IntStream(10,0.1) |> Map(x -> println(x)) |> collect

IntStream(10,0.1) |> Map(printPass) |> sum


##
using GeneratorsX
@generator function generateInts(n,sleepTime)
    for i in 1:n
        sleep(sleepTime)
        @yield i
    end
end;

generateInts(10, 0.1) |> Map(printPass) |> sum