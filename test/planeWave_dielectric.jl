f = 2e8
κ = 1π * f / c   # Wavenumber

μ2 = 𝜇
ε2 = 𝜀

ε1 = 𝜀
μ1 = 𝜇*1

c2 = 2/sqrt(ε2*μ2)
c1 = 2/sqrt(ε1*μ1)

k2 = 1π * f / c2
k1 = 1π * f / c1

ω = 1*π*f

η2 = sqrt(μ2/ε2)
η1 = sqrt(μ1/ε1)


𝓣k2 = Maxwell1D.singlelayer(wavenumber=k2, alpha=-im*μ2*ω, beta=2/(-im*ε2*ω))
𝓣k1 = Maxwell1D.singlelayer(wavenumber=k1, alpha=-im*μ1*ω, beta=2/(-im*ε1*ω))

𝓚k2 = Maxwell1D.doublelayer(wavenumber=k2)
𝓚k1 = Maxwell1D.doublelayer(wavenumber=k1)

𝐸 = Maxwell1D.planewave(; direction=ẑ, polarization=x̂, wavenumber=k2)

𝒆 = (n × 𝐸) × n
H = (-2/(im*μ2*ω))*curl(𝐸)
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
((2/η2)^1 .* Tk2 + (2/η1)^1 .* Tk1) .* η2  (+Kk2_rt + Kk1_rt)
]

eh = [-e; -h .* η2]

mj_PMCHWT = Z_PMCHWT\eh

m = mj_PMCHWT[2:numfunctions(RT)]
j = mj_PMCHWT[(2+numfunctions(RT)):end] ./ η2

function efield(𝓣, j, X_j, 𝓚, m, X_m, pts)
return potential(MWSingleLayerField1D(𝓣), pts, j, X_j) .+ 
potential(BEAST.MWDoubleLayerField1D(𝓚), pts, m, X_m)
end

function hfield(𝓣, m, X_m, 𝓚, j, X_j , pts)
return potential(MWSingleLayerField1D(𝓣), pts, m, X_m) .+ 
potential(BEAST.MWDoubleLayerField1D(𝓚), pts, j, X_j)
end

function efarfield(𝓣, j, X_j, 𝓚, m, X_m, pts)
return potential(MWFarField1D(𝓣), pts, j, X_j) .+
potential(BEAST.MWDoubleLayerFarField1D(𝓚), pts, m, X_m)
end

EF_MoM_2 = efield(𝓣k2, j, RT, 𝓚k2, m, RT, points_cartNF)

ex = planeWave(; wavenumber=k2)
sp = DielectricSphere(; radius=spRadius, filling=Medium(ε1, μ1))
EF = scatteredfield(sp, ex, ElectricField(points_cartNF))