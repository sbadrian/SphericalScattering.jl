
"""
    scatteredfield(sphere::Sphere, excitation::PlaneWave, quantity::Field; parameter::Parameter=Parameter())
    
Compute the electric field scattered by a PEC sphere, for an incident plane wave.
"""
function scatteredfield(sphere::Sphere, excitation::PlaneWave, quantity::Field; parameter::Parameter=Parameter())

    sphere.embedding == excitation.embedding || error("Excitation and sphere are not in the same medium.") # verify excitation and sphere are in the same medium

    T = typeof(excitation.frequency)

    F = zeros(SVector{3,Complex{T}}, size(quantity.locations))

    # --- rotate coordinates
    # rotate!(points, -excitation.rotation)

    # --- compute field in Cartesian representation
    for (ind, point) in enumerate(quantity.locations)
        F[ind] = scatteredfield(sphere, excitation, point, quantity; parameter=parameter)
    end

    # --- rotate resulting field
    # rotate!(F, excitation.rotation)

    return F
end



"""
    scatteredfield(sphere::Sphere, excitation::PlaneWave, point, quantity::ElectricField; parameter::Parameter=Parameter())

Compute the electric field scattered by a PEC or dielectric sphere, for an incident plane wave
travelling in +z-direction with E-field polarization in x-direction.

The point and the returned field are in Cartesian coordinates.
"""
function scatteredfield(sphere::Sphere, excitation::PlaneWave, point, quantity::Field; parameter::Parameter=Parameter())

    point_sph = cart2sph(point) # [r ϑ φ]
    r = point_sph[1]

    T = typeof(excitation.frequency)

    eps = parameter.relativeAccuracy

    if !(quantity isa FarField)
        sphere isa PECSphere && r < sphere.radius && return SVector{3,Complex{T}}(0.0, 0.0, 0.0) # inside the sphere the field is 0
    end

    if quantity isa MagneticField
        A₀ = 1 / impedance(sphere, r) * excitation.amplitude
    else
        A₀ = excitation.amplitude
    end

    Fr = Complex{T}(0.0)
    Fϑ = Complex{T}(0.0)
    Fϕ = Complex{T}(0.0)

    δF = T(Inf)
    n = 0

    sinϑ = abs(sin(point_sph[2]))  # note: theta only defined from from 0 to pi
    cosϑ = cos(point_sph[2])       # ok for theta > pi
    sinϕ = sin(point_sph[3])
    cosϕ = cos(point_sph[3])

    # first two values of the Associated Legendre Polynomial
    plm = Vector{T}()
    push!(plm, -sinϑ)
    push!(plm, -T(3.0) * sinϑ * cosϑ)

    try
        while δF > eps || n < 10
            n += 1

            coeffs = scatterCoeff(sphere, excitation, n)        


            if quantity isa FarField
                k = wavenumber(excitation)
                Nn_ϑ, Nn_ϕ, Mn_ϑ, Mn_ϕ = expansion(sphere, excitation, plm, cosϑ, sinϑ, n)
                aₙ = coeffs[1]
                bₙ = coeffs[2]
            else
                k = wavenumber(sphere, excitation, r)

                kr = k*r
                s = sqrt(π / 2 / kr)
    
                if r >= sphere.radius
                    Nn_r, Nn_ϑ, Nn_ϕ, Mn_ϑ, Mn_ϕ = expansion(sphere, excitation, plm, kr, s, cosϑ, sinϑ, n)
                    aₙ = coeffs[1]
                    bₙ = coeffs[2]
                else
                    Nn_r, Nn_ϑ, Nn_ϕ, Mn_ϑ, Mn_ϕ = expansion_dielectric_inner(sphere, excitation, plm, kr, s, cosϑ, sinϑ, n)
                    aₙ = coeffs[3]
                    bₙ = coeffs[4]
                end
            end

            if quantity isa ElectricField
                ΔFr = +(cosϕ / (im * kr^2)) * aₙ * Nn_r
                ΔFϑ = -(cosϕ / kr) * (aₙ * Nn_ϑ + bₙ * Mn_ϑ)
                ΔFϕ = +(sinϕ / kr) * (aₙ * Nn_ϕ + bₙ * Mn_ϕ)
            elseif quantity isa MagneticField
                ΔFr = +(sinϕ / (im * kr^2)) * bₙ * Nn_r
                ΔFϑ = -(sinϕ / kr) * (aₙ * Mn_ϑ + bₙ * Nn_ϑ)
                ΔFϕ = -(cosϕ / kr) * (aₙ * Mn_ϕ + bₙ * Nn_ϕ)
            else # FarField
                # See Jin, (7.4.44)-(7.4.45)
                # We deviate from Jin by replacing exp(-im*kr) / kr with 1/k
                ΔFr = T(0)
                ΔFϑ = -im * cosϕ * 1 / k * im^n * (aₙ * Nn_ϑ + bₙ * Mn_ϑ)
                ΔFϕ = +im * sinϕ * 1 / k * im^n * (aₙ * Nn_ϕ + bₙ * Mn_ϕ)
            end

            Fr += ΔFr
            Fϑ += ΔFϑ
            Fϕ += ΔFϕ

            δF = (abs(ΔFr) + abs(ΔFϑ) + abs(ΔFϕ)) / (abs(Fr) + abs(Fϑ) + abs(Fϕ)) # relative change

            n > 1 && push!(plm, (T(2.0) * n + 1) * cosϑ * plm[n] / n - (n + 1) * plm[n - 1] / n) # recurrence relationship for next associated Legendre polynomials
        end
    catch

    end

    return convertSpherical2Cartesian(A₀ .* SVector(Fr, Fϑ, Fϕ), point_sph)
end



"""
    scatterCoeff(sphere::PECSphere, excitation::PlaneWave, n::Int, ka)

Compute scattering coefficients for a plane wave travelling in +z-direction with polarization in x-direction.
"""
function scatterCoeff(sphere::PECSphere, excitation::PlaneWave, n::Int)

    T = typeof(excitation.frequency)

    k = wavenumber(excitation)

    ka = k*sphere.radius

    s = sqrt(π / 2 / ka)

    Ĵ  = ka * s * besselj(n + T(0.5), ka)   # Riccati-Bessel function
    Ĥ  = ka * s * hankelh2(n + T(0.5), ka)  # Riccati-Hankel function
    Ĵ2 = ka * s * besselj(n - T(0.5), ka)   # for derivate needed
    Ĥ2 = ka * s * hankelh2(n - T(0.5), ka)  # for derivate needed

    # Use recurrence relationship
    dĴ = (Ĵ2 - n / ka * Ĵ)    # derivative Riccati-Bessel function
    dĤ = (Ĥ2 - n / ka * Ĥ)    # derivative Riccati-Hankel function

    aₙ = -im^(-T(n)) * (dĴ / dĤ) * (2 * n + 1) / (n * (n + 1))  # Jin (7.4.41)
    bₙ = -im^(-T(n)) * (Ĵ / Ĥ) * (2 * n + 1) / (n * (n + 1))    # Jin (7.4.42)

    return aₙ, bₙ
end



function scatterCoeff(sphere::DielectricSphere, excitation::PlaneWave, n::Int)

    f = excitation.frequency
    T = typeof(f)

    ε2 = sphere.embedding.ε
    μ2 = sphere.embedding.μ

    ε1 = sphere.filling.ε
    μ1 = sphere.filling.μ

    c2 = 1/sqrt(ε2*μ2)
    c1 = 1/sqrt(ε1*μ1)

    k2 = 2π * f / c2
    k1 = 2π * f / c1

    εᵣ = ε1 / ε2
    μᵣ = μ1 / μ2

    k₂a = k2*sphere.radius
    s₂ = sqrt(π / 2 / k₂a)

    Ĵ₂  = k₂a * s₂ * besselj(n + T(0.5), k₂a)   # Riccati-Bessel function
    Ĥ₂  = k₂a * s₂ * hankelh2(n + T(0.5), k₂a)  # Riccati-Hankel function
    Ĵ₂2 = k₂a * s₂ * besselj(n - T(0.5), k₂a)   # for derivate needed
    Ĥ₂2 = k₂a * s₂ * hankelh2(n - T(0.5), k₂a)  # for derivate needed

    # Use recurrence relationship
    dĴ₂ = (Ĵ₂2 - n / k₂a * Ĵ₂)    # derivative Riccati-Bessel function
    dĤ₂ = (Ĥ₂2 - n / k₂a * Ĥ₂)    # derivative Riccati-Hankel function

    k₁a = k1*sphere.radius
    s₁ = sqrt(π / 2 / k₁a)

    Ĵ₁  = k₁a * s₁ * besselj(n + T(0.5), k₁a)   # Riccati-Bessel function
    #Ĥ₁  = k₁a * s₁ * hankelh2(n + T(0.5), k₁a)  # Riccati-Hankel function
    Ĵ₁2 = k₁a * s₁ * besselj(n - T(0.5), k₁a)   # for derivate needed
    #Ĥ₁2 = k₁a * s₁ * hankelh2(n - T(0.5), k₁a)  # for derivate needed

    # Use recurrence relationship
    dĴ₁ = (Ĵ₁2 - n / k₁a * Ĵ₁)    # derivative Riccati-Bessel function
    #dĤ₁ = (Ĥ₁2 - n / k₁a * Ĥ₁)    # derivative Riccati-Hankel function

    pF = im^(-T(n)) * (2 * n + 1) / (n * (n + 1))  

    aₙ = pF * (√εᵣ*dĴ₂*Ĵ₁ - √μᵣ*Ĵ₂*dĴ₁) / (√μᵣ*Ĥ₂*dĴ₁ - √εᵣ*dĤ₂*Ĵ₁)     # Jin (7.4.65)
    bₙ = pF * (√μᵣ*dĴ₂*Ĵ₁ - √εᵣ*Ĵ₂*dĴ₁) / (√εᵣ*Ĥ₂*dĴ₁ - √μᵣ*dĤ₂*Ĵ₁)     # Jin (7.4.66)
    cₙ = pF * (im*√εᵣ*μᵣ) / (√μᵣ*Ĥ₂*dĴ₁ - √εᵣ*dĤ₂*Ĵ₁)      # Jin (7.4.66)
    dₙ = pF * (im*√εᵣ*μᵣ) / (√εᵣ*Ĥ₂*dĴ₁ - √μᵣ*dĤ₂*Ĵ₁)      # Jin (7.4.66)

    return aₙ, bₙ, cₙ, dₙ
end



"""
    expansion(sphere::Sphere, excitation::PlaneWave, plm, kr, s, cosϑ, p, n::Int) 

Compute functional dependencies of the Mie series for a plane wave
travelling in +z-direction with polarization in x-direction.
"""
function expansion(sphere::Sphere, excitation::PlaneWave, plm, kr, s, cosϑ, sinϑ, n::Int)

    T = typeof(excitation.frequency)

    Ĥ  = kr * s * hankelh2(n + T(0.5), kr)     # Riccati-Hankel function
    Ĥ2 = kr * s * hankelh2(n - T(0.5), kr)

    dĤ = Ĥ2 - n / kr * Ĥ          # derivative of Riccati-Hankel function

    p = plm[n]

    if abs(cosϑ) < 0.999999
        if n == 1 # derivative of associated Legendre Polynomial
            dp = cosϑ * plm[1] / sqrt(T(1.0) - cosϑ * cosϑ)
        else
            dp = (n * cosϑ * plm[n] - (n + 1) * plm[n - 1]) / sqrt(T(1.0) - cosϑ * cosϑ)
        end

        Mn_ϑ = Ĥ * p / sinϑ
        Mn_ϕ = Ĥ * dp

        Nn_ϑ = dĤ * dp
        Nn_ϕ = dĤ * p / sinϑ

    elseif cosϑ > 0.999999
        aux = (n + T(1.0)) * n / T(2.0)

        Mn_ϑ = -Ĥ * aux
        Mn_ϕ = Mn_ϑ

        Nn_ϑ = -dĤ * aux
        Nn_ϕ = Nn_ϑ

    elseif cosϑ < -0.999999
        aux = (n + T(1.0)) * n / T(2.0) * T(-1.0)^n

        Mn_ϑ = Ĥ * aux
        Mn_ϕ = -Mn_ϑ

        Nn_ϑ = -dĤ * aux
        Nn_ϕ = -Nn_ϑ
    end

    Nn_r = n * (n + 1) * Ĥ * p
    Nn_ϑ = im * Nn_ϑ
    Nn_ϕ = im * Nn_ϕ

    return Nn_r, Nn_ϑ, Nn_ϕ, Mn_ϑ, Mn_ϕ
end



"""
    expansion(sphere::Sphere, excitation::PlaneWave, plm, cosϑ, sinϑ, n::Int)

Compute far-field functional dependencies of the Mie series for a plane wave travelling in -z direction with polarization in x-direction.
"""
function expansion(sphere::Sphere, excitation::PlaneWave, plm, cosϑ, sinϑ, n::Int)

    T = typeof(excitation.frequency)

    p = plm[n]

    if abs(cosϑ) < 0.999999
        if n == 1 # derivative of associated Legendre Polynomial
            dp = cosϑ * plm[1] / sqrt(T(1.0) - cosϑ * cosϑ)
        else
            dp = (n * cosϑ * plm[n] - (n + 1) * plm[n - 1]) / sqrt(T(1.0) - cosϑ * cosϑ)
        end

        Mn_ϑ = p / sinϑ
        Mn_ϕ = dp

        Nn_ϑ = dp
        Nn_ϕ = p / sinϑ

    elseif cosϑ > 0.999999
        aux = (n + T(1.0)) * n / T(2.0)

        Mn_ϑ = -aux
        Mn_ϕ = Mn_ϑ

        Nn_ϑ = -aux
        Nn_ϕ = Nn_ϑ

    elseif cosϑ < -0.999999
        aux = (n + T(1.0)) * n / T(2.0) * T(-1.0)^n

        Mn_ϑ = aux
        Mn_ϕ = -Mn_ϑ

        Nn_ϑ = -aux
        Nn_ϕ = -Nn_ϑ
    end

    return Nn_ϑ, Nn_ϕ, Mn_ϑ, Mn_ϕ
end



"""
    expansion_dielectric_inner(sphere::Sphere, excitation::PlaneWave, plm, kr, s, cosϑ, p, n::Int) 

Compute functional dependencies of the Mie series for a plane wave
travelling in +z-direction with polarization in x-direction.
"""
function expansion_dielectric_inner(sphere::Sphere, excitation::PlaneWave, plm, kr, s, cosϑ, sinϑ, n::Int)

    T = typeof(excitation.frequency)

    Ĵ  = kr * s * besselj(n + T(0.5), kr)     # Riccati-Hankel function
    Ĵ2 = kr * s * besselj(n - T(0.5), kr)

    dĴ = Ĵ2 - n / kr * Ĵ          # derivative of Riccati-Hankel function

    p = plm[n]

    if abs(cosϑ) < 0.999999
        if n == 1 # derivative of associated Legendre Polynomial
            dp = cosϑ * plm[1] / sqrt(T(1.0) - cosϑ * cosϑ)
        else
            dp = (n * cosϑ * plm[n] - (n + 1) * plm[n - 1]) / sqrt(T(1.0) - cosϑ * cosϑ)
        end

        Mn_ϑ = Ĵ * p / sinϑ
        Mn_ϕ = Ĵ * dp

        Nn_ϑ = dĴ * dp
        Nn_ϕ = dĴ * p / sinϑ

    elseif cosϑ > 0.999999
        aux = (n + T(1.0)) * n / T(2.0)

        Mn_ϑ = -Ĵ * aux
        Mn_ϕ = Mn_ϑ

        Nn_ϑ = -dĴ * aux
        Nn_ϕ = Nn_ϑ

    elseif cosϑ < -0.999999
        aux = (n + T(1.0)) * n / T(2.0) * T(-1.0)^n

        Mn_ϑ = Ĵ * aux
        Mn_ϕ = -Mn_ϑ

        Nn_ϑ = -dĴ * aux
        Nn_ϕ = -Nn_ϑ
    end

    Nn_r = n * (n + 1) * Ĵ * p
    Nn_ϑ = im * Nn_ϑ
    Nn_ϕ = im * Nn_ϕ

    return Nn_r, Nn_ϑ, Nn_ϕ, Mn_ϑ, Mn_ϕ
end