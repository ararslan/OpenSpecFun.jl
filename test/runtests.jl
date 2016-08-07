using OpenSpecFun
using Base.Test

const SF = OpenSpecFun

# error functions
@testset "error functions" begin
    @test SF.erf(Float16(1)) ≈ 0.84270079294971486934
    @test SF.erf(1) ≈ 0.84270079294971486934
    @test SF.erfc(1) ≈ 0.15729920705028513066
    @test SF.erfc(Float16(1)) ≈ 0.15729920705028513066
    @test SF.erfcx(1) ≈ 0.42758357615580700442
    @test SF.erfcx(Float32(1)) ≈ 0.42758357615580700442
    @test SF.erfcx(Complex64(1)) ≈ 0.42758357615580700442
    @test SF.erfi(1) ≈ 1.6504257587975428760
    @test SF.erfinv(0.84270079294971486934) ≈ 1
    @test SF.erfcinv(0.15729920705028513066) ≈ 1
    @test SF.dawson(1) ≈ 0.53807950691276841914

    @test SF.erf(1+2im) ≈ -0.53664356577856503399-5.0491437034470346695im
    @test SF.erfc(1+2im) ≈ 1.5366435657785650340+5.0491437034470346695im
    @test SF.erfcx(1+2im) ≈ 0.14023958136627794370-0.22221344017989910261im
    @test SF.erfi(1+2im) ≈ -0.011259006028815025076+1.0036063427256517509im
    @test SF.dawson(1+2im) ≈ -13.388927316482919244-11.828715103889593303im

    for elty in [Float32,Float64]
        for x in logspace(-200, -0.01)
            @test_approx_eq_eps SF.erf(SF.erfinv(x)) x 1e-12*x
            @test_approx_eq_eps SF.erf(SF.erfinv(-x)) -x 1e-12*x
            @test_approx_eq_eps SF.erfc(SF.erfcinv(2*x)) 2*x 1e-12*x
            if x > 1e-20
                xf = Float32(x)
                @test_approx_eq_eps SF.erf(SF.erfinv(xf)) xf 1e-5*xf
                @test_approx_eq_eps SF.erf(SF.erfinv(-xf)) -xf 1e-5*xf
                @test_approx_eq_eps SF.erfc(SF.erfcinv(2xf)) 2xf 1e-5*xf
            end
        end
        @test SF.erfinv(one(elty)) == Inf
        @test SF.erfinv(-one(elty)) == -Inf
        @test_throws DomainError SF.erfinv(convert(elty,2.0))

        @test SF.erfcinv(zero(elty)) == Inf
        @test_throws DomainError SF.erfcinv(-one(elty))
    end

    @test SF.erfinv(one(Int)) == SF.erfinv(1.0)
    @test SF.erfcinv(one(Int)) == SF.erfcinv(1.0)
end

# airy
@testset "airy" begin
    @test SF.airy(1.8) ≈ SF.airyai(1.8)
    @test SF.airyprime(1.8) ≈ -0.0685247801186109345638
    @test SF.airyaiprime(1.8) ≈ SF.airyprime(1.8)
    @test SF.airybi(1.8) ≈ 2.595869356743906290060
    @test SF.airybiprime(1.8) ≈ 2.98554005084659907283
    @test_throws SF.AmosException SF.airy(200im)
    @test_throws SF.AmosException SF.airybi(200)
    @test_throws ArgumentError SF.airy(5,one(Complex128))
    z = 1.8 + 1.0im
    for elty in [Complex64,Complex128]
        @test SF.airy(convert(elty,1.8)) ≈ 0.0470362168668458052247
        z = convert(elty,z)
        @test SF.airyx(z) ≈ SF.airyx(0,z)
        @test SF.airyx(0, z) ≈ SF.airy(0, z) * exp(2/3 * z * sqrt(z))
        @test SF.airyx(1, z) ≈ SF.airy(1, z) * exp(2/3 * z * sqrt(z))
        @test SF.airyx(2, z) ≈ SF.airy(2, z) * exp(-abs(real(2/3 * z * sqrt(z))))
        @test SF.airyx(3, z) ≈ SF.airy(3, z) * exp(-abs(real(2/3 * z * sqrt(z))))
        @test_throws ArgumentError SF.airyx(5,z)
    end
    @test_throws MethodError SF.airy(complex(big(1.0)))
end

# bessely0, bessely1, besselj0, besselj1
@testset "bessel" begin
    @test SF.besselj0(Float32(2.0)) ≈ SF.besselj0(Float64(2.0))
    @test SF.besselj1(Float32(2.0)) ≈ SF.besselj1(Float64(2.0))
    @test SF.bessely0(Float32(2.0)) ≈ SF.bessely0(Float64(2.0))
    @test SF.bessely1(Float32(2.0)) ≈ SF.bessely1(Float64(2.0))
    @test SF.besselj0(2) ≈ SF.besselj0(2.0)
    @test SF.besselj1(2) ≈ SF.besselj1(2.0)
    @test SF.bessely0(2) ≈ SF.bessely0(2.0)
    @test SF.bessely1(2) ≈ SF.bessely1(2.0)
    @test SF.besselj0(2.0 + im) ≈ SF.besselj(0, 2.0 + im)
    @test SF.besselj1(2.0 + im) ≈ SF.besselj(1, 2.0 + im)
    @test SF.bessely0(2.0 + im) ≈ SF.bessely(0, 2.0 + im)
    @test SF.bessely1(2.0 + im) ≈ SF.bessely(1, 2.0 + im)

    @test_throws MethodError SF.besselj(1.2,big(1.0))
    @test_throws MethodError SF.besselj(1,complex(big(1.0)))
    @test_throws MethodError SF.besseljx(1,big(1.0))
    @test_throws MethodError SF.besseljx(1,complex(big(1.0)))

    # besselh
    true_h133 = 0.30906272225525164362 - 0.53854161610503161800im
    @test SF.besselh(3,1,3) ≈ true_h133
    @test SF.besselh(-3,1,3) ≈ -true_h133
    @test SF.besselh(3,2,3) ≈ conj(true_h133)
    @test SF.besselh(-3,2,3) ≈ -conj(true_h133)
    @test_throws SF.AmosException SF.besselh(1,0)

    @test_throws MethodError SF.besselh(1,big(1.0))
    @test_throws MethodError SF.besselh(1,complex(big(1.0)))
    @test_throws MethodError SF.besselhx(1,big(1.0))
    @test_throws MethodError SF.besselhx(1,complex(big(1.0)))

    # besseli
    true_i33 = 0.95975362949600785698
    @test SF.besseli(3,3) ≈ true_i33
    @test SF.besseli(-3,3) ≈ true_i33
    @test SF.besseli(3,-3) ≈ -true_i33
    @test SF.besseli(-3,-3) ≈ -true_i33
    @test SF.besseli(Float32(-3),Complex64(-3,0)) ≈ -true_i33
    @test_throws SF.AmosException SF.besseli(1,1000)
    @test_throws DomainError SF.besseli(0.4,-1.0)

    @test_throws MethodError SF.besseli(1,big(1.0))
    @test_throws MethodError SF.besseli(1,complex(big(1.0)))
    @test_throws MethodError SF.besselix(1,big(1.0))
    @test_throws MethodError SF.besselix(1,complex(big(1.0)))


    # besselj
    @test SF.besselj(0,0) == 1
    for i = 1:5
        @test SF.besselj(i,0) == 0
        @test SF.besselj(-i,0) == 0
        @test SF.besselj(-i,Float32(0)) == 0
        @test SF.besselj(-i,Float32(0)) == 0
    end

    j33 = SF.besselj(3,3.)
    @test SF.besselj(3,3) == j33
    @test SF.besselj(-3,-3) == j33
    @test SF.besselj(-3,3) == -j33
    @test SF.besselj(3,-3) == -j33
    @test SF.besselj(3,3f0) ≈ j33
    @test SF.besselj(3,complex(3.)) ≈ j33
    @test SF.besselj(3,complex(3f0)) ≈ j33
    @test SF.besselj(3,complex(3)) ≈ j33

    j43 = SF.besselj(4,3.)
    @test SF.besselj(4,3) == j43
    @test SF.besselj(-4,-3) == j43
    @test SF.besselj(-4,3) == j43
    @test SF.besselj(4,-3) == j43
    @test SF.besselj(4,3f0) ≈ j43
    @test SF.besselj(4,complex(3.)) ≈ j43
    @test SF.besselj(4,complex(3f0)) ≈ j43
    @test SF.besselj(4,complex(3)) ≈ j43

    @test j33 ≈ 0.30906272225525164362
    @test j43 ≈ 0.13203418392461221033
    @test_throws DomainError    SF.besselj(0.1, -0.4)
    @test SF.besselj(0.1, complex(-0.4)) ≈ 0.820421842809028916 + 0.266571215948350899im
    @test SF.besselj(3.2, 1.3+0.6im) ≈ 0.01135309305831220201 + 0.03927719044393515275im
    @test SF.besselj(1, 3im) ≈ 3.953370217402609396im
    @test SF.besselj(1.0,3im) ≈ SF.besselj(1,3im)
    @test_throws SF.AmosException SF.besselj(20,1000im)
    @test_throws MethodError SF.besselj(big(1.0),3im)

    # besselk
    true_k33 = 0.12217037575718356792
    @test SF.besselk(3,3) ≈ true_k33
    @test SF.besselk(-3,3) ≈ true_k33
    true_k3m3 = -0.1221703757571835679 - 3.0151549516807985776im
    @test_throws DomainError SF.besselk(3,-3)
    @test SF.besselk(3,complex(-3)) ≈ true_k3m3
    @test SF.besselk(-3,complex(-3)) ≈ true_k3m3
    @test_throws SF.AmosException SF.besselk(200,0.01)
    # base julia issue #6564
    @test SF.besselk(1.0,0.0) == Inf

    @test_throws MethodError SF.besselk(1,big(1.0))
    @test_throws MethodError SF.besselk(1,complex(big(1.0)))
    @test_throws MethodError SF.besselkx(1,big(1.0))
    @test_throws MethodError SF.besselkx(1,complex(big(1.0)))


    # bessely
    y33 = SF.bessely(3,3.)
    @test SF.bessely(3,3) == y33
    @test SF.bessely(3.,3.) == y33
    @test SF.bessely(3,Float32(3.)) ≈ y33
    @test SF.bessely(-3,3) ≈ -y33
    @test y33 ≈ -0.53854161610503161800
    @test_throws DomainError SF.bessely(3,-3)
    @test SF.bessely(3,complex(-3)) ≈ 0.53854161610503161800 - 0.61812544451050328724im
    @test_throws SF.AmosException SF.bessely(200.5,0.1)
    @test_throws DomainError SF.bessely(0.4,-1.0)
    @test_throws DomainError SF.bessely(0.4,Float32(-1.0))
    @test_throws DomainError SF.bessely(1,Float32(-1.0))
    @test_throws DomainError SF.bessely(Cint(3),Float32(-3.))
    @test_throws DomainError SF.bessely(Cint(3),Float64(-3.))

    @test_throws MethodError SF.bessely(1.2,big(1.0))
    @test_throws MethodError SF.bessely(1,complex(big(1.0)))
    @test_throws MethodError SF.besselyx(1,big(1.0))
    @test_throws MethodError SF.besselyx(1,complex(big(1.0)))


    # besselhx
    for elty in [Complex64,Complex128]
        z = convert(elty, 1.0 + 1.9im)
        @test SF.besselhx(1.0, 1, z) ≈ convert(elty,-0.5949634147786144 - 0.18451272807835967im)
        @test SF.besselhx(Float32(1.0), 1, z) ≈ convert(elty,-0.5949634147786144 - 0.18451272807835967im)
    end

    @test_throws MethodError SF.besselh(1,1,big(1.0))
    @test_throws MethodError SF.besselh(1,1,complex(big(1.0)))
    @test_throws MethodError SF.besselhx(1,1,big(1.0))
    @test_throws MethodError SF.besselhx(1,1,complex(big(1.0)))

    # base julia issue #6653
    for f in (SF.besselj,SF.bessely,SF.besseli,SF.besselk,SF.hankelh1,SF.hankelh2)
        @test f(0,1) ≈ f(0,Complex128(1))
        @test f(0,1) ≈ f(0,Complex64(1))
    end

    # scaled bessel[ijky] and hankelh[12]
    for x in (1.0, 0.0, -1.0), y in (1.0, 0.0, -1.0), nu in (1.0, 0.0, -1.0)
        z = Complex128(x + y * im)
        z == zero(z) || @test SF.hankelh1x(nu, z) ≈ SF.hankelh1(nu, z) * exp(-z * im)
        z == zero(z) || @test SF.hankelh2x(nu, z) ≈ SF.hankelh2(nu, z) * exp(z * im)
        (nu < 0 && z == zero(z)) || @test SF.besselix(nu, z) ≈ SF.besseli(nu, z) * exp(-abs(real(z)))
        (nu < 0 && z == zero(z)) || @test SF.besseljx(nu, z) ≈ SF.besselj(nu, z) * exp(-abs(imag(z)))
        z == zero(z) || @test SF.besselkx(nu, z) ≈ SF.besselk(nu, z) * exp(z)
        z == zero(z) || @test SF.besselyx(nu, z) ≈ SF.bessely(nu, z) * exp(-abs(imag(z)))
    end
    @test_throws SF.AmosException SF.hankelh1x(1, 0)
    @test_throws SF.AmosException SF.hankelh2x(1, 0)
    @test_throws SF.AmosException SF.besselix(-1, 0)
    @test_throws SF.AmosException SF.besseljx(-1, 0)
    @test SF.besselkx(1, 0) == Inf
    @test_throws SF.AmosException SF.besselyx(1, 0)
    @test_throws DomainError SF.besselix(0.4,-1.0)
    @test_throws DomainError SF.besseljx(0.4, -1.0)
    @test_throws DomainError SF.besselkx(0.4,-1.0)
    @test_throws DomainError SF.besselyx(0.4,-1.0)
end
