"""
    scatteredfield(sphere::Sphere, excitation::PlaneWave, quantity::Field; parameter::Parameter=Parameter())

Compute the electric field scattered by a sphere, for an incident uniform field.
"""
function scatteredfield(sphere::Sphere, excitation::UniformField, quantity::Field; parameter::Parameter=Parameter())

    sphere.embedding == excitation.embedding || error("Excitation and sphere are not in the same medium.") # verify excitation and sphere are in the same medium

    F = zeros(fieldType(quantity), size(quantity.locations))

    # --- compute field in Cartesian representation
    for (ind, point) in enumerate(quantity.locations)
        F[ind] = scatteredfield(sphere, excitation, point, quantity; parameter=parameter)
    end

    return F
end



"""
    scatteredfield(sphere::DielectricSphere, excitation::UniformField, point, quantity::ElectricField; parameter::Parameter=Parameter())

Compute the electric field scattered by a Dielectric sphere, for an incident uniform field with polarization in x-direction.

The point and the returned field are in Cartesian coordinates.
"""
function scatteredfield(
    sphere::DielectricSphere, excitation::UniformField, point, quantity::ElectricField; parameter::Parameter=Parameter()
)

    ε0 = sphere.embedding.ε
    ε1 = sphere.filling.ε
    E0 = field(excitation, point, quantity)

    R = sphere.radius
    r = norm(point)

    if r <= R
        return (3 * ε0 / (ε1 + 2 * ε0) - 1) * E0 # Sihvola&Lindell 1988, (8)
    end

    return E0 * (-(ε1 - ε0) / (ε1 + 2 * ε0) * R^3 / r^3) + 3 * (ε1 - ε0) / (ε1 + 2 * ε0) * R^3 * point / r^5 * dot(E0, point) #Sihvola&Lindell 1988, (9)

end

"""
    scatteredfield(sphere::DielectricSphere, excitation::UniformField, point, quantity::ScalarPotential; parameter::Parameter=Parameter())

Compute the scalar potential scattered by a Dielectric sphere, for an incident uniform field with polarization in x-direction.

The point and the returned field are in Cartesian coordinates.
"""
function scatteredfield(
    sphere::DielectricSphere, excitation::UniformField, point, quantity::ScalarPotential; parameter::Parameter=Parameter()
)

    ε0 = sphere.embedding.ε
    ε1 = sphere.filling.ε
    Φ0 = field(excitation, point, quantity)

    R = sphere.radius
    r = norm(point)

    if r <= R
        return (3 * ε0 / (ε1 + 2 * ε0) - 1) * Φ0 # Sihvola&Lindell 1988, (8)
    end

    return -Φ0 * ((ε1 - ε0) / (ε1 + 2 * ε0) * R^3 / r^3) #Sihvola&Lindell 1988, (9)

end



"""
    scatteredfield(sphere::DielectricSphereThinLayerPotentialJump, excitation::UniformField, point, quantity::ElectricField; parameter::Parameter=Parameter())

Compute the electric field scattered by a dielectric sphere with a thin coating,
where the displacement field in the coating is only in radial direction.
We assume an an incident uniform field with polarization in z-direction.

The point and the returned field are in Cartesian coordinates.
"""
function scatteredfield(
    sphere::DielectricSphereThinLayerPotentialJump,
    excitation::UniformField,
    point,
    quantity::ElectricField;
    parameter::Parameter=Parameter(),
)

    ẑ = SVector(0.0, 0.0, 1.0)
    # Currently, we only support a constant field in z-direction
    @assert excitation.direction == ẑ

    point_sph = cart2sph(point)
    θ = point_sph[2]

    cosθ = dot(ẑ, point) / norm(point)

    @assert cosθ ≈ cos(point_sph[2]) atol = 1e-6

    r = norm(point)

    E0 = excitation.amplitude

    A, K = scatterCoeff(sphere, excitation)

    if r > sphere.radius
        E = SVector((-2 * A / r^3) * cosθ, (+A / r^3) * (-sin(θ)), 0.0)
    else
        E = SVector((-K * cosθ, K * sin(θ), 0.0))
    end

    return E
end

function scatteredfield(
    sphere::DielectricSphereThinLayerPotentialJump,
    excitation::UniformField,
    point,
    quantity::ScalarPotentialJump;
    parameter::Parameter=Parameter(),
)

    ẑ = SVector(0.0, 0.0, 1.0)
    # Currently, we only support a constant field in z-direction
    @assert excitation.direction == ẑ

    point_sph = cart2sph(point)
    @show point_sph
    θ = point_sph[2]

    cosθ = dot(ẑ, point) / norm(point)
    @show cosθ
    @show cos(point_sph[2])
    @assert cosθ ≈ cos(point_sph[2]) atol = 1e-6

    ~, K = scatterCoeff(sphere, excitation)

    return sphere.thickness * (sphere.filling.ε / sphere.thinlayer.ε) * K * cosθ

end

"""
    scatteredfield(sphere::DielectricSphereThinLayerPotentialJump, excitation::UniformField, point, quantity::ScalarPotential; parameter::Parameter=Parameter())

Compute the scalar potential scattered by a dielectric sphere with a thin coating,
where the displacement field in the coating is only in radial direction.
We assume an an incident uniform field with polarization in z-direction.

The point and the returned field are in Cartesian coordinates.
"""
function scatteredfield(
    sphere::DielectricSphereThinLayerPotentialJump,
    excitation::UniformField,
    point,
    quantity::ScalarPotential;
    parameter::Parameter=Parameter(),
)

    ẑ = SVector(0.0, 0.0, 1.0)
    # Currently, we only support a constant field in z-direction
    @assert excitation.direction == ẑ

    cosθ = dot(ẑ, point) / norm(point)
    r = norm(point)

    R = sphere.radius

    A, K = scatterCoeff(sphere, excitation)

    if r > R
        return (+A / r^2) * cosθ # Jones 1995, (C.1a)
    else
        return (-K * r * cosθ)
    end
end

function scatterCoeff(sp::DielectricSphereThinLayerPotentialJump, ex::UniformField)
    R = sp.radius
    Δ = sp.thickness
    εₘ = sp.thinlayer.ε
    εₑ = sp.embedding.ε
    εᵢ = sp.filling.ε

    # We divide the second equation by
    # εₑ to improve the conditioning
    b = [R; 1.0] .* ex.amplitude

    A = [
        R^(-2)   R-Δ * (εᵢ / εₘ)
        -2R^(-3)     εᵢ/εₑ
    ]

    Ainv = [
        εᵢ/εₑ     -(R - Δ * (εᵢ / εₘ))
        2R^(-3)     R^(-2)
    ] ./ (εᵢ / εₑ * R^(-2) + (R - Δ * (εᵢ / εₘ)) * 2R^(-3))

    if norm(A * Ainv - I) > 2e-6
        println("Matrix inversion is unstable: ", norm(A * Ainv - I))
        println("Condition number is: ", cond(A))
        println("Alternative inversion leads to: ", norm(A * pinv(A) - I))
    end
    x = Ainv * b

    if norm(A * x - b) / norm(b) > eps(eltype(R)) * 10
        print("No stable solution possible: ", norm(A * x - b) / norm(b))
    end

    return x[1], x[2]
end

"""
    scatteredfield(sphere::PECSphere, excitation::UniformField, point, quantity::ElectricField; parameter::Parameter=Parameter())

Compute the electric field scattered by a PEC sphere, for an incident uniform field with polarization in x-direction.
The point and returned field are in Cartesian coordinates.
"""
function scatteredfield(sphere::PECSphere, excitation::UniformField, point, quantity::ElectricField; parameter::Parameter=Parameter())

    E0 = field(excitation, point, quantity)

    R = sphere.radius
    r = norm(point)

    T = eltype(point)
    if r <= R
        return -E0
    end

    return -E0 * R^3 / r^3 + 3 * R^3 / r^5 * point * dot(E0, point) # Griffits, Example 3.8
end

"""
    scatteredfield(sphere::PECSphere, excitation::UniformField, point, quantity::ScalarPotential; parameter::Parameter=Parameter())

Compute the scalar potential scattered by a PEC sphere, for an incident uniform field with polarization in x-direction.
The point and returned field are in Cartesian coordinates.
"""
function scatteredfield(
    sphere::PECSphere, excitation::UniformField, point, quantity::ScalarPotential; parameter::Parameter=Parameter()
)

    Φ0 = field(excitation, point, quantity)

    R = sphere.radius
    r = norm(point)

    T = eltype(point)
    if r <= R
        return -Φ0
    end

    return Φ0 * (-R^3 / r^3) # Griffits, Example 3.8
end

"""
    scatteredfield(sphere::LayeredSphere, excitation::UniformField, point, quantity::ElectricField; parameter::Parameter=Parameter())

Compute the electric field scattered by a layered dielectric sphere, for an incident uniform field with polarization in x-direction.
The point and returned field are in Cartesian coordinates.
"""
function scatteredfield(
    sphere::LayeredSphere{LN,LR,LC},
    excitation::UniformField{FC,FT,FR},
    point,
    quantity::ElectricField;
    parameter::Parameter=Parameter(),
) where {LN,LR,LC,FC,FT,FR}
    # using `Sihvola and Lindell, 1988, Transmission line analogy for calculating the effective permittivity of mixtures with spherical multilayer scatterers`

    E0 = excitation.amplitude

    r = norm(point)
    a = sphere.radii
    dir = excitation.direction
    n = length(a)
    perms = getfield.(vcat(sphere.embedding, sphere.filling), 1)

    T = promote_type(LR, LC, FC, FT, FR)

    Bk = zeros(SMatrix{2,2,T}, n)
    B = Matrix(I, 2, 2)

    for k in range(n; stop=1, step=-1)
        B11 = perms[k + 1] + 2 * perms[k]
        B12 = 2 * (perms[k + 1] - perms[k]) * a[k]^-3
        B21 = (perms[k + 1] - perms[k]) * a[k]^3
        B22 = 2 * perms[k + 1] + perms[k]

        Bk[k] = 1 / (3 * perms[k]) * ([B11 B12; B21 B22])
        B = Bk[k] * B
    end

    Ck = zeros(T, n + 1)
    Dk = zeros(T, n + 1)

    Ck[1] = 1
    Dk[n + 1] = 0
    Dk[1] = B[2, 1] / B[1, 1]
    Ck[n + 1] = 1 / B[1, 1]

    for k in range(n; stop=2, step=-1)
        Bk[k]
        [Ck[k + 1], Dk[k + 1]]
        Ck[k], Dk[k] = Bk[k] * [Ck[k + 1], Dk[k + 1]]
    end

    # find out in what layer `point` is
    pos = 0

    for k in 1:n
        if r < a[k]
            pos = k
        end
    end

    C = Ck[pos + 1]
    D = Dk[pos + 1]

    return C * E0 * dir - D * E0 * dir / r^3 + 3 * D * E0 * point * dot(dir, point) / r^5 - E0 * dir
end

"""
    scatteredfield(sphere::LayeredSphere, excitation::UniformField, point, quantity::ScalarPotential; parameter::Parameter=Parameter())

Compute the scalar potential scattered by a layered dielectric sphere, for an incident uniform field with polarization in x-direction.
The point and returned field are in Cartesian coordinates.
"""
function scatteredfield(
    sphere::LayeredSphere{LN,LR,LC},
    excitation::UniformField{FC,FT,FR},
    point,
    quantity::ScalarPotential;
    parameter::Parameter=Parameter(),
) where {LN,LR,LC,FC,FT,FR}
    # using `Sihvola and Lindell, 1988, Transmission line analogy for calculating the effective permittivity of mixtures with spherical multilayer scatterers`

    Φ0 = field(excitation, point, quantity)

    r = norm(point)
    a = sphere.radii
    dir = excitation.direction
    n = length(a)
    perms = getfield.(vcat(sphere.embedding, sphere.filling), 1)

    T = promote_type(LR, LC, FC, FT, FR)

    Bk = zeros(SMatrix{2,2,T}, n)
    B = Matrix(I, 2, 2)

    for k in range(n; stop=1, step=-1)
        B11 = perms[k + 1] + 2 * perms[k]
        B12 = 2 * (perms[k + 1] - perms[k]) * a[k]^-3
        B21 = (perms[k + 1] - perms[k]) * a[k]^3
        B22 = 2 * perms[k + 1] + perms[k]

        Bk[k] = 1 / (3 * perms[k]) * ([B11 B12; B21 B22])
        B = Bk[k] * B
    end

    Ck = zeros(T, n + 1)
    Dk = zeros(T, n + 1)

    Ck[1] = 1
    Dk[n + 1] = 0
    Dk[1] = B[2, 1] / B[1, 1]
    Ck[n + 1] = 1 / B[1, 1]

    for k in range(n; stop=2, step=-1)
        #Bk[k]
        #[Ck[k+1],Dk[k+1]]
        Ck[k], Dk[k] = Bk[k] * [Ck[k + 1], Dk[k + 1]]
    end

    # find out in what layer `point` is
    pos = 0

    for k in 1:n
        if r < a[k]
            pos = k
        end
    end

    C = Ck[pos + 1]
    D = Dk[pos + 1]

    return C * Φ0 - D * Φ0 / r^3 - Φ0
end

"""
    scatteredfield(sphere::LayeredSpherePEC, excitation::UniformField{FC,FT,FR}, point, quantity::ScalarPotential; parameter::Parameter=Parameter())

Compute the scalar potential scattered by a layered dielectric sphere with PEC core, for an incident uniform field with polarization in x-direction.
The point and returned field are in Cartesian coordinates.
"""
function scatteredfield(
    sphere::LayeredSpherePEC{LN,LD,LR,LC},
    excitation::UniformField{FC,FT,FR},
    point,
    quantity::ScalarPotential;
    parameter::Parameter=Parameter(),
) where {LN,LD,LR,LC,FC,FT,FR}
    # using `Sihvola and Lindell, 1988, Transmission line analogy for calculating the effective permittivity of mixtures with spherical multilayer scatterers`

    Φ0 = field(excitation, point, quantity)

    r = norm(point)
    a = sphere.radii
    dir = excitation.direction
    n = length(a) - 1
    perms = getfield.(vcat(sphere.embedding, sphere.filling), 1)

    T = promote_type(LR, LC, FC, FT, FR)

    Bk = zeros(SMatrix{2,2,T}, n)
    B = Matrix(I, 2, 2)

    for k in range(n; stop=1, step=-1)
        B11 = perms[k + 1] + 2 * perms[k]
        B12 = 2 * (perms[k + 1] - perms[k]) * a[k]^-3
        B21 = (perms[k + 1] - perms[k]) * a[k]^3
        B22 = 2 * perms[k + 1] + perms[k]

        Bk[k] = 1 / (3 * perms[k]) * ([B11 B12; B21 B22])
        B = Bk[k] * B
    end
    Ck = zeros(T, n + 1)
    Dk = zeros(T, n + 1)

    Ck[1] = 1
    Dk[1] = (B[2, 1] + B[2, 2] * a[end]^3) / (B[1, 1] + B[1, 2] * a[end]^3)
    Ck[n + 1] = 1 / (B[1, 1] + B[1, 2] * a[end]^3)
    Dk[n + 1] = a[end]^3 * Ck[n + 1]

    for k in range(n; stop=2, step=-1)
        #Bk[k]
        #[Ck[k+1],Dk[k+1]]
        Ck[k], Dk[k] = Bk[k] * [Ck[k + 1], Dk[k + 1]]
    end
    # find out in what layer `point` is
    pos = 0

    for k in 1:(n + 1)
        if r < a[k]
            pos = k
        end
    end

    if pos == n + 1
        C = 0
        D = 0
        R = 0
    else
        C = Ck[pos + 1]
        D = Dk[pos + 1]
        R = a[pos + 1]
    end

    return C * Φ0 - D * Φ0 / r^3 - Φ0

end

"""
    scatteredfield(sphere::LayeredSpherePEC, excitation::UniformField, point, quantity::ElectricField; parameter::Parameter=Parameter())

Compute the electric field scattered by a layered dielectric sphere with PEC core, for an incident uniform field with polarization in x-direction.
The point and returned field are in Cartesian coordinates.
"""
function scatteredfield(
    sphere::LayeredSpherePEC{LN,LD,LR,LC},
    excitation::UniformField{FC,FT,FR},
    point,
    quantity::ElectricField;
    parameter::Parameter=Parameter(),
) where {LN,LD,LR,LC,FC,FT,FR}
    # using `Sihvola and Lindell, 1988, Transmission line analogy for calculating the effective permittivity of mixtures with spherical multilayer scatterers`

    E0 = excitation.amplitude
    r = norm(point)
    a = sphere.radii
    dir = excitation.direction
    n = length(a) - 1
    perms = getfield.(vcat(sphere.embedding, sphere.filling), 1)

    T = promote_type(LR, LC, FC, FT, FR)

    Bk = zeros(SMatrix{2,2,T}, n)
    B = Matrix(I, 2, 2)

    for k in range(n; stop=1, step=-1)
        B11 = perms[k + 1] + 2 * perms[k]
        B12 = 2 * (perms[k + 1] - perms[k]) * a[k]^-3
        B21 = (perms[k + 1] - perms[k]) * a[k]^3
        B22 = 2 * perms[k + 1] + perms[k]

        Bk[k] = 1 / (3 * perms[k]) * ([B11 B12; B21 B22])
        B = Bk[k] * B
    end

    Ck = zeros(T, n + 1)
    Dk = zeros(T, n + 1)

    Ck[1] = 1
    Dk[1] = (B[2, 1] + B[2, 2] * a[end]^3) / (B[1, 1] + B[1, 2] * a[end]^3)
    Ck[n + 1] = 1 / (B[1, 1] + B[1, 2] * a[end]^3)
    Dk[n + 1] = a[end]^3 * Ck[n + 1]

    for k in range(n; stop=2, step=-1)
        Ck[k], Dk[k] = Bk[k] * [Ck[k + 1], Dk[k + 1]]
    end
    # find out in what layer `point` is
    pos = 0

    for k in 1:(n + 1)
        if r < a[k]
            pos = k
        end
    end

    if pos == n + 1
        C = 0
        D = 0
        R = 0
    else
        C = Ck[pos + 1]
        D = Dk[pos + 1]
        R = a[pos + 1]
    end

    return C * E0 * dir - D * E0 * dir / r^3 + 3 * D * E0 * point * dot(dir, point) / r^5 - E0 * dir

end

fieldType(F::ElectricField) = SVector{3,Complex{eltype(F.locations[1])}}
fieldType(F::ScalarPotential) = Complex{eltype(F.locations[1])}
fieldType(F::ScalarPotentialJump) = Complex{eltype(F.locations[1])}
