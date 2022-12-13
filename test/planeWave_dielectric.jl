f = 1e8
κ = 2π * f / c   # Wavenumber

μ1 = 𝜇
ε1 = 𝜀

ε2 = 𝜀*4
μ2 = 𝜇*3

c1 = 1/sqrt(ε1*μ1)
c2 = 1/sqrt(ε2*μ2)

k1 = 2π * f / c1
k2 = 2π * f / c2

ω = 2*π*f

η1 = sqrt(μ1/ε1)
η2 = sqrt(μ2/ε2)


𝓣k1 = Maxwell3D.singlelayer(wavenumber=k1, alpha=-im*μ1*ω, beta=1/(-im*ε1*ω))
𝓣k2 = Maxwell3D.singlelayer(wavenumber=k2, alpha=-im*μ2*ω, beta=1/(-im*ε2*ω))

𝓚k1 = Maxwell3D.doublelayer(wavenumber=k1)
𝓚k2 = Maxwell3D.doublelayer(wavenumber=k2)

𝐸 = Maxwell3D.planewave(; direction=ẑ, polarization=x̂, wavenumber=k1)

𝒆 = (n × 𝐸) × n
H = (-1/(im*μ1*ω))*curl(𝐸)
𝒉 = (n × H) × n
nx𝒉 = (n × H)

Tk1 = Matrix(assemble(𝓣k1, RT, RT))
Tk2 = Matrix(assemble(𝓣k2, RT, RT))

Kk1_rt = Matrix(assemble(𝓚k1, RT, RT))
Kk2_rt = Matrix(assemble(𝓚k2, RT, RT))

e = Vector(assemble(𝒆, RT))
h = Vector(assemble(𝒉, RT))

Z_PMCHWT = 
[
-(Kk1_rt + Kk2_rt)  (Tk1 + Tk2) ./ η1;
((1/η1)^2 .* Tk1 + (1/η2)^2 .* Tk2) .* η1  (+Kk1_rt + Kk2_rt)
]

eh = [-e; -h .* η1]

mj_PMCHWT = Z_PMCHWT\eh

m = mj_PMCHWT[1:numfunctions(RT)]
j = mj_PMCHWT[(1+numfunctions(RT)):end] ./ η1

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

EF_MoM_1 = efield(𝓣k1, j, RT, 𝓚k1, m, RT, points_cartNF)

ex = planeWave(; wavenumber=k1)
sp = DielectricSphere(; radius=spRadius, filling=Medium(ε2, μ2))
EF = scatteredfield(sp, ex, ElectricField(points_cartNF))