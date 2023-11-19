module c

export fc, sc

struct sc
    x::String
end

function fc()
    println("This is module c's func")
end

end