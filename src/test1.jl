module b
include("test2.jl")
using .c
export fb
function fb(var::sc)
    print(var.x)
end
end

using .c
function fbd(var::sc)
    print(var.x)
end