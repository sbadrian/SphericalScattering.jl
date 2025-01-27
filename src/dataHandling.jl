

abstract type Field end

struct FarField <: Field
    locations#::Vector
end

struct ElectricField <: Field
    locations#::Vector
end

struct MagneticField <: Field
    locations#::Vector
end

struct ScalarPotential <: Field
    locations
end

abstract type Excitation end

#abstract type Parameter end

struct Parameter
    nmax::Int
    relativeAccuracy::AbstractFloat
end

Parameter() = Parameter(-1, 1e-12)
