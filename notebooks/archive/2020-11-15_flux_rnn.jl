using Flux

m = Chain(LSTM(3,2))
data = rand(3,5,7)
inputs = [data[:,:,t] for t in 1:size(data,3)]
output = m.(inputs)
size(output), size(output[1])