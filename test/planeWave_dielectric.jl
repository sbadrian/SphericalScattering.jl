f = 1e8
κ = 2π * f / c   # Wavenumber

# Embedding
μ2 = 𝜇
ε2 = 𝜀

# Filling
ε1 = 𝜀*(8)
μ1 = 𝜇

c2 = 1/sqrt(ε2*μ2)
c1 = 1/sqrt(ε1*μ1)

k2 = 2π * f / c2
k1 = 2π * f / c1

ω = 2*π*f

η2 = sqrt(μ2/ε2)
η1 = sqrt(μ1/ε1)

𝓣k2 = Maxwell3D.singlelayer(wavenumber=k2, alpha=-im*μ2*ω, beta=1/(-im*ε2*ω))
𝓣k1 = Maxwell3D.singlelayer(wavenumber=k1, alpha=-im*μ1*ω, beta=1/(-im*ε1*ω))

𝓚k2 = Maxwell3D.doublelayer(wavenumber=k2)
𝓚k1 = Maxwell3D.doublelayer(wavenumber=k1)

𝐸 = Maxwell3D.planewave(; direction=ẑ, polarization=x̂, wavenumber=k2)

𝒆 = (n × 𝐸) × n
H = (-1/(im*μ2*ω))*curl(𝐸)
𝒉 = (n × H) × n
nx𝒉 = (n × H)

Tk2 = Matrix(assemble(𝓣k2, RT, RT))
Tk1 = Matrix(assemble(𝓣k1, RT, RT))

Kk2_rt = Matrix(assemble(𝓚k2, RT, RT))
Kk1_rt = Matrix(assemble(𝓚k1, RT, RT))

e = Vector(assemble(𝒆, RT))
h = Vector(assemble(𝒉, RT))

Z_PMCHWT = 
[
-(Kk2_rt + Kk1_rt)  (Tk2 + Tk1) ./ η2;
((1/η2)^2 .* Tk2 + (1/η1)^2 .* Tk1) .* η2  (+Kk2_rt + Kk1_rt)
]

eh = [-e; -h .* η2]

mj_PMCHWT = Z_PMCHWT\eh

m = mj_PMCHWT[1:numfunctions(RT)]
j = mj_PMCHWT[(1+numfunctions(RT)):end] ./ η2

function efield(𝓣, j, X_j, 𝓚, m, X_m, pts)
    return potential(MWSingleLayerField3D(𝓣), pts, j, X_j) .+ 
    potential(BEAST.MWDoubleLayerField3D(𝓚), pts, m, X_m)
end

function hfield(𝓣, m, X_m, 𝓚, j, X_j , pts)
    return potential(MWSingleLayerField3D(𝓣), pts, m, X_m) .+ 
    potential(BEAST.MWDoubleLayerField3D(𝓚), pts, j, X_j)
end

function efarfield(𝓣, j, X_j, 𝓚, m, X_m, pts)
    return potential(MWFarField3D(𝓣), pts, j, X_j) .+
    potential(BEAST.MWDoubleLayerFarField3D(𝓚), pts, m, X_m)
end

points_cartNF_inside, points_sphNF = getDefaultPoints(0.5)

EF₂MoM = efield(𝓣k2, j, RT, 𝓚k2, -m, RT, points_cartNF)
EF₁MoM = efield(𝓣k1, -j, RT, 𝓚k1, +m, RT, points_cartNF_inside)

ex = planeWave(; wavenumber=k2, embedding=Medium(ε2, μ2), frequency=f)
sp = DielectricSphere(; radius=spRadius, filling=Medium(ε1, μ1), embedding=Medium(ε2, μ2))

EF₂ = scatteredfield(sp, ex, ElectricField(points_cartNF))
EF₁ = scatteredfield(sp, ex, ElectricField(points_cartNF_inside))

diff_EF₂ = norm.(EF₂ - EF₂MoM) ./ maximum(norm.(EF₂))  # worst case error
diff_EF₁ = norm.(EF₁ - EF₁MoM) ./ maximum(norm.(EF₁))  # worst case error

@test maximum(20 * log10.(abs.(diff_EF₂))) < -27 # dB 
@test maximum(20 * log10.(abs.(diff_EF₁))) < -27 # dB 

##
HF₂MoM = hfield(𝓣k2, (1/η2)^2 .* m, RT, 𝓚k2, -j, RT, points_cartNF)
HF₁MoM = hfield(𝓣k1, -(1/η1)^2 .* m, RT, 𝓚k1, +j, RT, points_cartNF_inside)

ex = planeWave(; wavenumber=k2, embedding=Medium(ε2, μ2), frequency=f)
sp = DielectricSphere(; radius=spRadius, filling=Medium(ε1, μ1), embedding=Medium(ε2, μ2))

HF₂ = scatteredfield(sp, ex, MagneticField(points_cartNF))

diff_HF₂ = norm.(HF₂ - HF₂MoM) ./ maximum(norm.(HF₂))  # worst case error
@test maximum(20 * log10.(abs.(diff_HF₂))) < -27 # dB 

##
function get_spherical_coordinates(fld, pts, ϑ, ϕ)
    retfld = copy(fld)
    for i in eachindex(ϑ)
        for j in eachindex(ϕ)
            retfld[i, j] = SphericalScattering.convertCartesian2Spherical(fld[i,j], pts[i, j])
        end
    end

    return retfld
end

sEF₂ = get_spherical_coordinates(EF₂, points_cartNF, ϑ, ϕ)
sEF₂M = get_spherical_coordinates(EF₂MoM, points_cartNF, ϑ, ϕ)

##
sHF₂ = get_spherical_coordinates(HF₂, points_cartNF, ϑ, ϕ)
sHF₂M = get_spherical_coordinates(HF₂MoM, points_cartNF, ϑ, ϕ)
sHF₁ = get_spherical_coordinates(EF₁, points_cartNF_inside, ϑ, ϕ)
sHF₁M = get_spherical_coordinates(EF₁MoM, points_cartNF_inside, ϑ, ϕ)