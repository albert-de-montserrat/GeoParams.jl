export computeViscosity_TauII,
    computeViscosity_EpsII, computeViscosity_TauII!, computeViscosity_EpsII!

"""
    compute viscosity given strain rate 2nd invariant

    τ = 2ηε -> η = τ/2/ε
"""
function computeViscosity_EpsII(εII, v, args)
    τII = computeCreepLaw_TauII(εII, v, args) # gives 
    η = 0.5 * τII / εII
    return η
end

"""
    compute viscosity given strain rate 2nd invariant

    τ = 2ηε -> η = 2ε/τ
"""
@inline function computeViscosity_TauII(τII, v, args)
    εII = computeCreepLaw_EpsII(τII, v, args)
    η = 2 * εII / τII
    return η
end

@inline function computeViscosity_TauII(τII::T, v::Tuple, args) where {T}
    return computeViscosity(computeViscosity_TauII, τII, v, args, Val(length(v)))
end
@inline function computeViscosity_EpsII(εII::T, v::Tuple, args) where {T}
    return computeViscosity(computeViscosity_EpsII, εII, v, args, Val(length(v)))
end

@inline @generated function computeViscosity(
    fn::F, CII::T, v::Tuple, args::NamedTuple, ::Val{N}
) where {F,T,N}
    quote
        η = 0.0
        Base.Cartesian.@nexprs $N i -> η += 1 / fn(CII, v[i], args)
        return 1 / η
    end
end

@inline function computeViscosity_TauII!(τII::T, v::Tuple, args) where {T}
    Threads.@Threads for I in eachindex(τII)
        computeViscosity(
            computeViscosity_TauII,
            τII,
            v,
            (; zip(keys(args), getindex.(values(args), I))...),
            Val(length(v)),
        )
    end
end

@inline function computeViscosity_EpsII!(τII::T, v::Tuple, args) where {T}
    Threads.@Threads for I in eachindex(τII)
        computeViscosity(
            computeViscosity_EpsII,
            τII,
            v,
            (; zip(keys(args), getindex.(values(args), I))...),
            Val(length(v)),
        )
    end
end

@inline function computeViscosity_EpsII!(εII::T, v::Tuple, args) where {T}
    return computeViscosity(computeViscosity_EpsII, εII, v, args, Val(length(v)))
end

# OPTION 1
strainCircuit(TauII, v, args) = strainCircuit(TauII, v, args, Val(length(v)))

@inline @generated function strainCircuit(TauII, v, args, ::Val{N}; n=1) where {N}
    quote
        c = 0.0
        Base.Cartesian.@nexprs $N i ->
            c += if v[i] isa Tuple
                1 / strainCircuit(TauII, v[i], args, Val(length(v[i])); n=-1)
            else
                computeCreepLaw_EpsII(TauII, v[i], args)^n
            end
        return c
    end
end

function viscosityCircuit_TauII(τII::T, v, args) where {T}
    εII = strainCircuit(τII, v, args)
    return η = 2 * εII / τII
end

# v = (DiffusionCreep(), DislocationCreep())
# v = (DiffusionCreep(), (DiffusionCreep(), DislocationCreep()))

# EpsII = 1e-10
# TauII = 20e6
# P, T, f, d = 20e6, 300.0, 1.0, 1.0
# args = (P=P, T=T, f=f, d=d)

# @benchmark viscosityCircuit_TauII($TauII, $v, $args)
# @benchmark strainCircuit($TauII, $v, $args)

# strainCircuit(TauII, v, args)

# computeCreepLaw_TauII(EpsII, DiffusionCreep(), args)
# @btime computeCreepLaw_TauII($EpsII, $(DiffusionCreep()), $args)

# computeCreepLaw_TauII(EpsII, DislocationCreep(), args)
# computeCreepLaw_EpsII(TauII, DiffusionCreep(), args)
# computeCreepLaw_EpsII(TauII, DislocationCreep(), args)

# v = (DiffusionCreep(), DislocationCreep())

# computeViscosity_EpsII2(EpsII, v, args)
# computeViscosity_EpsII(EpsII, v, args)

# @benchmark computeViscosity_EpsII2($EpsII, $v, $args)
# @benchmark computeViscosity_EpsII($EpsII, $v, $args)
