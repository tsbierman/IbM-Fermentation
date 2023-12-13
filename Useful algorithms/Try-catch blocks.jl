try
    throw(ErrorException("Not going well"))
catch e
    if isa(e, ErrorException)
        println(typeof(getfield(e, :msg)))
    else
        println("What")
    end
end
