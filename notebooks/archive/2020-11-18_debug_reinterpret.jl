## pure julia TCP working!
using Sockets
``

PORT = 1050

@async begin
    server = listen(PORT)
    while true
        sock = accept(server)
        println("server: got connection")

        # send array shape (C, Z, H, W)
        C, Z, H, W = 1, 1, 16, 16
        arr = Array{Int16}(undef,4)
        arr[:] .= [C, Z, H, W]
        write(sock, arr)
        println("server: sent sizes")

        # send mock array
        arr = Array{Int16}(undef, H, W)
        arr[:] .= collect(1:H*W)
        write(sock, arr)
        println("server: sent array")
    end
end

sleep(0.5)
##

client = connect(PORT)
# client = connect(HOSTNAME, PORT)
println("client: connected")
bSize = read(client, 8)
C, Z, H, W  = reinterpret(Int16,bSize)
@show "client", C, Z, H, W
bArray = read(client, 2* C * H * W ) # read one frame

println("client: received a frame!")
@show "client", typeof(bArray)
@show "client", size(bArray)
array = reshape(reinterpret(Int16, bArray),1, 1, 16,16)
# array = copy(array)
@show "client", typeof(array)
@show "client", size(array)
permArray = permutedims(array, [1,2,4,3]) # may copy
@show "client", typeof(permArray)
@show "client", permutedims(array, [1,2,4,3])[1,1,:,:]

nothing