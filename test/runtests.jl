using OpenSpecFun
using Base.Test

# error functions
@test erf(Float16(1)) ≈ 0.84270079294971486934
@test erf(1) ≈ 0.84270079294971486934
@test erfc(1) ≈ 0.15729920705028513066
@test erfc(Float16(1)) ≈ 0.15729920705028513066
@test erfcx(1) ≈ 0.42758357615580700442
@test erfcx(Float32(1)) ≈ 0.42758357615580700442
@test erfcx(Complex64(1)) ≈ 0.42758357615580700442
@test erfi(1) ≈ 1.6504257587975428760
@test erfinv(0.84270079294971486934) ≈ 1
@test erfcinv(0.15729920705028513066) ≈ 1
@test dawson(1) ≈ 0.53807950691276841914

@test erf(1+2im) ≈ -0.53664356577856503399-5.0491437034470346695im
@test erfc(1+2im) ≈ 1.5366435657785650340+5.0491437034470346695im
@test erfcx(1+2im) ≈ 0.14023958136627794370-0.22221344017989910261im
@test erfi(1+2im) ≈ -0.011259006028815025076+1.0036063427256517509im
@test dawson(1+2im) ≈ -13.388927316482919244-11.828715103889593303im

for elty in [Float32,Float64]
    for x in logspace(-200, -0.01)
        @test_approx_eq_eps erf(erfinv(x)) x 1e-12*x
        @test_approx_eq_eps erf(erfinv(-x)) -x 1e-12*x
        @test_approx_eq_eps erfc(erfcinv(2*x)) 2*x 1e-12*x
        if x > 1e-20
            xf = Float32(x)
            @test_approx_eq_eps erf(erfinv(xf)) xf 1e-5*xf
            @test_approx_eq_eps erf(erfinv(-xf)) -xf 1e-5*xf
            @test_approx_eq_eps erfc(erfcinv(2xf)) 2xf 1e-5*xf
        end
    end
    @test erfinv(one(elty)) == Inf
    @test erfinv(-one(elty)) == -Inf
    @test_throws DomainError erfinv(convert(elty,2.0))

    @test erfcinv(zero(elty)) == Inf
    @test_throws DomainError erfcinv(-one(elty))
end

@test erfinv(one(Int)) == erfinv(1.0)
@test erfcinv(one(Int)) == erfcinv(1.0)

# airy
@test airy(1.8) ≈ airyai(1.8)
@test airyprime(1.8) ≈ -0.0685247801186109345638
@test airyaiprime(1.8) ≈ airyprime(1.8)
@test airybi(1.8) ≈ 2.595869356743906290060
@test airybiprime(1.8) ≈ 2.98554005084659907283
@test_throws OpenSpecFun.AmosException airy(200im)
@test_throws OpenSpecFun.AmosException airybi(200)
@test_throws ArgumentError airy(5,one(Complex128))
z = 1.8 + 1.0im
for elty in [Complex64,Complex128]
    @test airy(convert(elty,1.8)) ≈ 0.0470362168668458052247
    z = convert(elty,z)
    @test airyx(z) ≈ airyx(0,z)
    @test airyx(0, z) ≈ airy(0, z) * exp(2/3 * z * sqrt(z))
    @test airyx(1, z) ≈ airy(1, z) * exp(2/3 * z * sqrt(z))
    @test airyx(2, z) ≈ airy(2, z) * exp(-abs(real(2/3 * z * sqrt(z))))
    @test airyx(3, z) ≈ airy(3, z) * exp(-abs(real(2/3 * z * sqrt(z))))
    @test_throws ArgumentError airyx(5,z)
end
@test_throws MethodError airy(complex(big(1.0)))

# bessely0, bessely1, besselj0, besselj1
@test besselj0(Float32(2.0)) ≈ besselj0(Float64(2.0))
@test besselj1(Float32(2.0)) ≈ besselj1(Float64(2.0))
@test bessely0(Float32(2.0)) ≈ bessely0(Float64(2.0))
@test bessely1(Float32(2.0)) ≈ bessely1(Float64(2.0))
@test besselj0(2) ≈ besselj0(2.0)
@test besselj1(2) ≈ besselj1(2.0)
@test bessely0(2) ≈ bessely0(2.0)
@test bessely1(2) ≈ bessely1(2.0)
@test besselj0(2.0 + im) ≈ besselj(0, 2.0 + im)
@test besselj1(2.0 + im) ≈ besselj(1, 2.0 + im)
@test bessely0(2.0 + im) ≈ bessely(0, 2.0 + im)
@test bessely1(2.0 + im) ≈ bessely(1, 2.0 + im)

@test_throws MethodError besselj(1.2,big(1.0))
@test_throws MethodError besselj(1,complex(big(1.0)))
@test_throws MethodError besseljx(1,big(1.0))
@test_throws MethodError besseljx(1,complex(big(1.0)))

# besselh
true_h133 = 0.30906272225525164362 - 0.53854161610503161800im
@test besselh(3,1,3) ≈ true_h133
@test besselh(-3,1,3) ≈ -true_h133
@test besselh(3,2,3) ≈ conj(true_h133)
@test besselh(-3,2,3) ≈ -conj(true_h133)
@test_throws OpenSpecFun.AmosException besselh(1,0)

@test_throws MethodError besselh(1,big(1.0))
@test_throws MethodError besselh(1,complex(big(1.0)))
@test_throws MethodError besselhx(1,big(1.0))
@test_throws MethodError besselhx(1,complex(big(1.0)))

# besseli
true_i33 = 0.95975362949600785698
@test besseli(3,3) ≈ true_i33
@test besseli(-3,3) ≈ true_i33
@test besseli(3,-3) ≈ -true_i33
@test besseli(-3,-3) ≈ -true_i33
@test besseli(Float32(-3),Complex64(-3,0)) ≈ -true_i33
@test_throws OpenSpecFun.AmosException besseli(1,1000)
@test_throws DomainError besseli(0.4,-1.0)

@test_throws MethodError besseli(1,big(1.0))
@test_throws MethodError besseli(1,complex(big(1.0)))
@test_throws MethodError besselix(1,big(1.0))
@test_throws MethodError besselix(1,complex(big(1.0)))


# besselj
@test besselj(0,0) == 1
for i = 1:5
    @test besselj(i,0) == 0
    @test besselj(-i,0) == 0
    @test besselj(-i,Float32(0)) == 0
    @test besselj(-i,Float32(0)) == 0
end

j33 = besselj(3,3.)
@test besselj(3,3) == j33
@test besselj(-3,-3) == j33
@test besselj(-3,3) == -j33
@test besselj(3,-3) == -j33
@test besselj(3,3f0) ≈ j33
@test besselj(3,complex(3.)) ≈ j33
@test besselj(3,complex(3f0)) ≈ j33
@test besselj(3,complex(3)) ≈ j33

j43 = besselj(4,3.)
@test besselj(4,3) == j43
@test besselj(-4,-3) == j43
@test besselj(-4,3) == j43
@test besselj(4,-3) == j43
@test besselj(4,3f0) ≈ j43
@test besselj(4,complex(3.)) ≈ j43
@test besselj(4,complex(3f0)) ≈ j43
@test besselj(4,complex(3)) ≈ j43

@test j33 ≈ 0.30906272225525164362
@test j43 ≈ 0.13203418392461221033
@test_throws DomainError    besselj(0.1, -0.4)
@test besselj(0.1, complex(-0.4)) ≈ 0.820421842809028916 + 0.266571215948350899im
@test besselj(3.2, 1.3+0.6im) ≈ 0.01135309305831220201 + 0.03927719044393515275im
@test besselj(1, 3im) ≈ 3.953370217402609396im
@test besselj(1.0,3im) ≈ besselj(1,3im)
@test_throws OpenSpecFun.AmosException besselj(20,1000im)
@test_throws MethodError besselj(big(1.0),3im)

# besselk
true_k33 = 0.12217037575718356792
@test besselk(3,3) ≈ true_k33
@test besselk(-3,3) ≈ true_k33
true_k3m3 = -0.1221703757571835679 - 3.0151549516807985776im
@test_throws DomainError besselk(3,-3)
@test besselk(3,complex(-3)) ≈ true_k3m3
@test besselk(-3,complex(-3)) ≈ true_k3m3
@test_throws OpenSpecFun.AmosException besselk(200,0.01)
# issue #6564
@test besselk(1.0,0.0) == Inf

@test_throws MethodError besselk(1,big(1.0))
@test_throws MethodError besselk(1,complex(big(1.0)))
@test_throws MethodError besselkx(1,big(1.0))
@test_throws MethodError besselkx(1,complex(big(1.0)))


# bessely
y33 = bessely(3,3.)
@test bessely(3,3) == y33
@test bessely(3.,3.) == y33
@test bessely(3,Float32(3.)) ≈ y33
@test bessely(-3,3) ≈ -y33
@test y33 ≈ -0.53854161610503161800
@test_throws DomainError bessely(3,-3)
@test bessely(3,complex(-3)) ≈ 0.53854161610503161800 - 0.61812544451050328724im
@test_throws OpenSpecFun.AmosException bessely(200.5,0.1)
@test_throws DomainError bessely(0.4,-1.0)
@test_throws DomainError bessely(0.4,Float32(-1.0))
@test_throws DomainError bessely(1,Float32(-1.0))
@test_throws DomainError bessely(Cint(3),Float32(-3.))
@test_throws DomainError bessely(Cint(3),Float64(-3.))

@test_throws MethodError bessely(1.2,big(1.0))
@test_throws MethodError bessely(1,complex(big(1.0)))
@test_throws MethodError besselyx(1,big(1.0))
@test_throws MethodError besselyx(1,complex(big(1.0)))


#besselhx
for elty in [Complex64,Complex128]
    z = convert(elty, 1.0 + 1.9im)
    @test besselhx(1.0, 1, z) ≈ convert(elty,-0.5949634147786144 - 0.18451272807835967im)
    @test besselhx(Float32(1.0), 1, z) ≈ convert(elty,-0.5949634147786144 - 0.18451272807835967im)
end

@test_throws MethodError besselh(1,1,big(1.0))
@test_throws MethodError besselh(1,1,complex(big(1.0)))
@test_throws MethodError besselhx(1,1,big(1.0))
@test_throws MethodError besselhx(1,1,complex(big(1.0)))

# issue #6653
for f in (besselj,bessely,besseli,besselk,hankelh1,hankelh2)
    @test f(0,1) ≈ f(0,Complex128(1))
    @test f(0,1) ≈ f(0,Complex64(1))
end

# scaled bessel[ijky] and hankelh[12]
for x in (1.0, 0.0, -1.0), y in (1.0, 0.0, -1.0), nu in (1.0, 0.0, -1.0)
    z = Complex128(x + y * im)
    z == zero(z) || @test hankelh1x(nu, z) ≈ hankelh1(nu, z) * exp(-z * im)
    z == zero(z) || @test hankelh2x(nu, z) ≈ hankelh2(nu, z) * exp(z * im)
    (nu < 0 && z == zero(z)) || @test besselix(nu, z) ≈ besseli(nu, z) * exp(-abs(real(z)))
    (nu < 0 && z == zero(z)) || @test besseljx(nu, z) ≈ besselj(nu, z) * exp(-abs(imag(z)))
    z == zero(z) || @test besselkx(nu, z) ≈ besselk(nu, z) * exp(z)
    z == zero(z) || @test besselyx(nu, z) ≈ bessely(nu, z) * exp(-abs(imag(z)))
end
@test_throws OpenSpecFun.AmosException hankelh1x(1, 0)
@test_throws OpenSpecFun.AmosException hankelh2x(1, 0)
@test_throws OpenSpecFun.AmosException besselix(-1, 0)
@test_throws OpenSpecFun.AmosException besseljx(-1, 0)
@test besselkx(1, 0) == Inf
@test_throws OpenSpecFun.AmosException besselyx(1, 0)
@test_throws DomainError besselix(0.4,-1.0)
@test_throws DomainError besseljx(0.4, -1.0)
@test_throws DomainError besselkx(0.4,-1.0)
@test_throws DomainError besselyx(0.4,-1.0)
