ab = 1

for i in 1:1
    nanka = Symbol("a" * "b")
    println(@eval Symbol("a" * "b"))
    println(@eval $(Symbol("a" * "b")))
    println(eval(Symbol("a" * "b")))
end

for i in 1:3
    nanka = Symbol("a" * "b")
    println(@eval nanka)
end

dump(Symbol("a" * "b"))
dump($(Symbol("a" * "b")))
