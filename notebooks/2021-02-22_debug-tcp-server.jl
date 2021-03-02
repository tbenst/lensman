using Sockets
@async begin
    server = listen(2001)
    while true
        sock = accept(server)
        write(sock, collect(1:100))
        write(sock, collect(200:500))
    end
end
##

HOSTNAME = "localhost"
PORT = 2001
pv = connect(HOSTNAME, PORT)
println("connected")
data = []

while isopen(pv)
    println("open")
    push!(data, read(pv,10))
    println("pushed")
end
