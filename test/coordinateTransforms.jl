
@testset "Coordinate transforms" begin

    # ----- sph2cart
    vec = SVector(1.0, π / 4, π / 4)
    xyz = SphericalScattering.sph2cart(vec)

    @test xyz[1] ≈ 0.5
    @test xyz[2] ≈ 0.5
    @test xyz[3] ≈ 1 / √2

    # ----- convertCartesian2Spherical
    F_cart    = SVector(1.0, 1.0, 1.0)
    point_sph = SVector(1, π / 2, π / 2)
    F_sph     = SphericalScattering.convertCartesian2Spherical(F_cart, point_sph)

    @test F_sph[1] ≈ +1.0
    @test F_sph[2] ≈ -1.0
    @test F_sph[3] ≈ -1.0

    # ----- rotate
    point = SVector(2.0, 2.2, 3.1)
    @test_nowarn SphericalScattering.rotate!(point, SVector(1.0, 2.3))
end
