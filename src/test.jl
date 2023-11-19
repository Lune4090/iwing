include("test1.jl")
include("test2.jl")

using .b, .c

tmp = sc("hello")
c.fc()

fbd(tmp)