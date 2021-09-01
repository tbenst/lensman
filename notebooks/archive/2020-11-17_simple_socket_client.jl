##
using Sockets
``
HOSTNAME = "171.65.101.131" # B113
PORT = 1050
pv = connect(HOSTNAME, PORT)
# pv = connect(PORT) # localhost
println("client: connected")
@async while isopen(pv)
    # write(stdout, readline(pv, keep=true))
    # @show readline(pv, keep=true)
    # @show dims = read(pv, Int64)
    bSize = read(pv, sizeof(Int16)*4)
    # reshape fails unless an Int64
    C, Z, H, W = Int64.(reinterpret(Int16,bSize))
    @show C, Z, H, W
    bArray = read(pv, sizeof(Int16)*C*H*W)
    @show size(bArray)
    @show typeof(C)
    @show typeof(1)
    # array = reshape(reinterpret(Int16, bArray),1, 1, 16,16)
    array = permutedims(reshape(reinterpret(Int16, bArray),C,Z,H,W),[1,2,4,3])
    # @show array = reshape(reinterpret(Int16, bArray),C,Z,H,W)
    @show size(array)
    @show array[1,1,:,:]
    # @show readline(pv, keep=true)
    # write(stdout, readline(pv, keep=true))
    break
end

nothing