function read_user_vars(
    freestream::NamedTuple{Float64},
    dx::Union{Missing,Float64},
    dt::Union{Missing,Float64},
    Re::Float64,
    T::Union{Missing, Float64})

    # dx 
    if ismissing(dx) == true
        dx = 2.0/Re # default to a grid Re of 2
    end

    # Freestream
    xkey = :Ux in keys(freestream)
    ykey = :Uy in keys(freestream)
    anglekey = :inclination in keys(freestream)

    if xkey == false
        Ux = t -> t^0.0 # other way to say Ux = 1.0. Here t is a dummy variable, like x in f(x)
    else
        Ux = freestream.Ux
    end

    if ykey == false
        Uy = t -> 0.0 # other way to say Uy = 0.0

        #=
        ty -> 0.0*ty^0.0 ==> Does not work when trying to extract Uy from freestream.
        I believe that this is due the fact that implicit functions cannot be zero out (t -> t*0.0 ==> WRONG).
        Hence just assign Uy = t -> 0.0, if ykey == false.
        =#
    else
        Uy = freestream.Uy
    end

    if anglekey == false
        inclination = t -> 0.0 # other way to say inclination = 0.0
    else
        inclination = freestream.inclination
    end

    tnames = (:Ux, :Uy, :inclination)
    tvals = (Ux, Uy, inclination)

    freestream = (;zip(tnames, tvals)...)


    # TODO: FINISH FROM HERE
    # Approximate dt using dx if not provided by user
    if ismissing(dt) == true
        if ismissing(T) == false
            tvect = range(0.0,T, length = 5000)
        else
            tvect = range(0.0,20.0, length = 5000)
        end
        Umax = sqrt(maximum(freestream.Ux.(tvect))^2.0 .+
            maximum(freestream.Uy.(tvect))^2.0 )

        dt = 0.1 * dx/(5.0*Umax) #satisfy CFL 0.2 constraint, with
                                 #safety factor on max velocity

        if ismissing(T) == true
            T = 20.0*dt
        end
    end

    return dx, dt, T, freestream
end
