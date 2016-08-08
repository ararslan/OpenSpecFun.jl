using OpenSpecFun
using Base.Test

const SF = OpenSpecFun


# useful test functions for relative error, which differ from isapprox
# in that errc separately looks at the real and imaginay parts
err(z, x) = z == x ? 0.0 : abs(z - x) / abs(x)
errc(z, x) = max(err(real(z),real(x)), err(imag(z),imag(x)))
≅(a,b) = errc(a,b) ≤ 1e-13


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

# SF.beta, SF.lbeta
@testset "beta, lbeta" begin
    @test SF.beta(3/2,7/2) ≈ 5π/128
    @test SF.beta(3,5) ≈ 1/105
    @test SF.lbeta(5,4) ≈ log(SF.beta(5,4))
    @test SF.beta(5,4) ≈ SF.beta(4,5)
    @test SF.beta(-1/2, 3) ≈ SF.beta(-1/2 + 0im, 3 + 0im) ≈ -16/3
    @test SF.lbeta(-1/2, 3) ≈ log(16/3)
    @test SF.beta(Float32(5),Float32(4)) == SF.beta(Float32(4),Float32(5))
    @test SF.beta(3,5) ≈ SF.beta(3+0im,5+0im)
    @test(SF.beta(3.2+0.1im,5.3+0.3im) ≈ exp(SF.lbeta(3.2+0.1im,5.3+0.3im)) ≈
          0.00634645247782269506319336871208405439180447035257028310080 -
          0.00169495384841964531409376316336552555952269360134349446910im)
end

# SF.gamma, SF.lgamma (complex argument)
@testset "gamma, lgamma" begin
    if Base.Math.libm == "libopenlibm"
        @test SF.gamma(Float64[1:25;]) == SF.gamma(1:25)
    else
        @test SF.gamma(Float64[1:25;]) ≈ SF.gamma(1:25)
    end
    for elty in (Float32, Float64)
        @test SF.gamma(convert(elty,1/2)) ≈ convert(elty,sqrt(π))
        @test SF.gamma(convert(elty,-1/2)) ≈ convert(elty,-2sqrt(π))
        @test SF.lgamma(convert(elty,-1/2)) ≈ convert(elty,log(abs(SF.gamma(-1/2))))
    end
    @test SF.lgamma(1.4+3.7im) ≈ -3.7094025330996841898 + 2.4568090502768651184im
    @test SF.lgamma(1.4+3.7im) ≈ log(SF.gamma(1.4+3.7im))
    @test SF.lgamma(-4.2+0im) ≈ SF.lgamma(-4.2)-pi*im
    @test factorial(3.0) == SF.gamma(4.0) == factorial(3)
    for x in (3.2, 2+1im, 3//2, 3.2+0.1im)
        @test factorial(x) == SF.gamma(1+x)
    end
    @test lfact(1) == 0
    @test lfact(2) == SF.lgamma(3)
end

# SF.digamma
@testset "digamma" begin
    for elty in (Float32, Float64)
        @test SF.digamma(convert(elty, 9)) ≈ convert(elty, 2.140641477955609996536345)
        @test SF.digamma(convert(elty, 2.5)) ≈ convert(elty, 0.7031566406452431872257)
        @test SF.digamma(convert(elty, 0.1)) ≈ convert(elty, -10.42375494041107679516822)
        @test SF.digamma(convert(elty, 7e-4)) ≈ convert(elty, -1429.147493371120205005198)
        @test SF.digamma(convert(elty, 7e-5)) ≈ convert(elty, -14286.29138623969227538398)
        @test SF.digamma(convert(elty, 7e-6)) ≈ convert(elty, -142857.7200612932791081972)
        @test SF.digamma(convert(elty, 2e-6)) ≈ convert(elty, -500000.5772123750382073831)
        @test SF.digamma(convert(elty, 1e-6)) ≈ convert(elty, -1000000.577214019968668068)
        @test SF.digamma(convert(elty, 7e-7)) ≈ convert(elty, -1428572.005785942019703646)
        @test SF.digamma(convert(elty, -0.5)) ≈ convert(elty, .03648997397857652055902367)
        @test SF.digamma(convert(elty, -1.1)) ≈ convert(elty,  10.15416395914385769902271)

        @test SF.digamma(convert(elty, 0.1)) ≈ convert(elty, -10.42375494041108)
        @test SF.digamma(convert(elty, 1/2)) ≈ convert(elty, -γ - log(4))
        @test SF.digamma(convert(elty, 1)) ≈ convert(elty, -γ)
        @test SF.digamma(convert(elty, 2)) ≈ convert(elty, 1 - γ)
        @test SF.digamma(convert(elty, 3)) ≈ convert(elty, 3/2 - γ)
        @test SF.digamma(convert(elty, 4)) ≈ convert(elty, 11/6 - γ)
        @test SF.digamma(convert(elty, 5)) ≈ convert(elty, 25/12 - γ)
        @test SF.digamma(convert(elty, 10)) ≈ convert(elty, 7129/2520 - γ)
    end


    for x in -10.2:0.3456:50
        @test 1e-12 > err(SF.digamma(x+0im), SF.digamma(x))
    end

    # Test against Wolfram Alpha
    @test SF.digamma(7+0im) ≅ 1.872784335098467139393487909917597568957840664060076401194232
    @test SF.digamma(7im) ≅ 1.94761433458434866917623737015561385331974500663251349960124 + 1.642224898223468048051567761191050945700191089100087841536im
    @test SF.digamma(-3.2+0.1im) ≅ 4.65022505497781398615943030397508454861261537905047116427511+2.32676364843128349629415011622322040021960602904363963042380im

    @test isa([SF.digamma(x) for x in [1.0]], Vector{Float64})
end

# SF.trigamma
@testset "trigamma" begin
    for elty in (Float32, Float64)
        @test SF.trigamma(convert(elty, 0.1)) ≈ convert(elty, 101.433299150792758817)
        @test SF.trigamma(convert(elty, 1/2)) ≈ convert(elty, π^2/2)
        @test SF.trigamma(convert(elty, 1)) ≈ convert(elty, π^2/6)
        @test SF.trigamma(convert(elty, 2)) ≈ convert(elty, π^2/6 - 1)
        @test SF.trigamma(convert(elty, 3)) ≈ convert(elty, π^2/6 - 5/4)
        @test SF.trigamma(convert(elty, 4)) ≈ convert(elty, π^2/6 - 49/36)
        @test SF.trigamma(convert(elty, 5)) ≈ convert(elty, π^2/6 - 205/144)
        @test SF.trigamma(convert(elty, 10)) ≈ convert(elty, π^2/6 - 9778141/6350400)
    end

    # Test against Wolfram Alpha
    @test SF.trigamma(8+0im) ≅ 0.133137014694031425134546685920401606452509991909746283540546
    @test SF.trigamma(8im) ≅ -0.0078125000000000000029194973110119898029284994355721719150 - 0.12467345030312762782439017882063360876391046513966063947im
    @test SF.trigamma(-3.2+0.1im) ≅ 15.2073506449733631753218003030676132587307964766963426965699+15.7081038855113567966903832015076316497656334265029416039199im

    @test isa([SF.trigamma(x) for x in [1.0]], Vector{Float64})
end

# SF.invdigamma
@testset "invdigamma" begin
    for elty in (Float32, Float64)
        for val in [0.001, 0.01, 0.1, 1.0, 10.0]
            @test abs(SF.invdigamma(SF.digamma(convert(elty, val))) - convert(elty, val)) < 1e-8
        end
    end
    @test abs(SF.invdigamma(2)) == abs(SF.invdigamma(2.))
end

@testset "polygamma" begin
    @test SF.polygamma(20, 7.) ≈ -4.644616027240543262561198814998587152547
    @test SF.polygamma(20, Float16(7.)) ≈ -4.644616027240543262561198814998587152547

    # Test against Wolfram Alpha
    @test SF.polygamma(2, 8.1+0im) ≅ -0.01723882695611191078960494454602091934457319791968308929600
    @test SF.polygamma(30, 8.1+2im) ≅ -2722.8895150799704384107961215752996280795801958784600407589+6935.8508929338093162407666304759101854270641674671634631058im
    @test SF.polygamma(3, 2.1+1im) ≅ 0.00083328137020421819513475400319288216246978855356531898998-0.27776110819632285785222411186352713789967528250214937861im
    @test 1e-11 > err(SF.polygamma(3, -4.2 + 2im),-0.0037752884324358856340054736472407163991189965406070325067-0.018937868838708874282432870292420046797798431078848805822im)
    @test SF.polygamma(13, 5.2 - 2im) ≅ 0.08087519202975913804697004241042171828113370070289754772448-0.2300264043021038366901951197725318713469156789541415899307im
    @test 1e-11 > err(SF.polygamma(123, -47.2 + 0im), 5.7111648667225422758966364116222590509254011308116701029e291)

    @test SF.polygamma(4, -3.1+Inf*im) == SF.polygamma(4, 3.1+Inf*im) == 0
    @test SF.polygamma(4, -0.0) == Inf == -SF.polygamma(4, +0.0)
    @test isa([SF.polygamma(3,x) for x in [1.0]], Vector{Float64})
    @test SF.polygamma(3,5) ≈ SF.polygamma(3,5.)
end

# SF.eta, SF.zeta
@testset "eta, zeta" begin
    @test SF.eta(1) ≈ log(2)
    @test SF.eta(2) ≈ pi^2/12
    @test SF.eta(Float32(2)) ≈ SF.eta(2)
    @test SF.eta(Complex64(2)) ≈ SF.eta(2)
    @test SF.zeta(0) ≈ -0.5
    @test SF.zeta(2) ≈ pi^2/6
    @test SF.zeta(Complex64(2)) ≈ SF.zeta(2)
    @test SF.zeta(4) ≈ pi^4/90
    @test SF.zeta(1,Float16(2.)) ≈ SF.zeta(1,2.)
    @test SF.zeta(1.,Float16(2.)) ≈ SF.zeta(1,2.)
    @test SF.zeta(Float16(1.),Float16(2.)) ≈ SF.zeta(1,2.)
    @test isnan(SF.zeta(NaN))
    @test isnan(SF.zeta(1.0e0))
    @test isnan(SF.zeta(1.0f0))
    @test isnan(SF.zeta(complex(0,Inf)))
    @test isnan(SF.zeta(complex(-Inf,0)))


    @test SF.zeta(4.1+0.3im, -3.2+0.1im) ≅ -281.34474134962502296077659347175501181994490498591796647 + 286.55601240093672668066037366170168712249413003222992205im
    @test SF.zeta(4.1+0.3im, 3.2+0.1im) ≅ 0.0121197525131633219465301571139288562254218365173899270675-0.00687228692565614267981577154948499247518236888933925740902im
    @test SF.zeta(4.1, 3.2+0.1im) ≅ 0.0137637451187986846516125754047084829556100290057521276517-0.00152194599531628234517456529686769063828217532350810111482im
    @test 1e-12 > errc(SF.zeta(1.0001, -4.5e2+3.2im), 10003.765660925877260544923069342257387254716966134385170 - 0.31956240712464746491659767831985629577542514145649468090im)
    @test SF.zeta(3.1,-4.2) ≅ SF.zeta(3.1,-4.2+0im) ≅ 149.7591329008219102939352965761913107060718455168339040295
    @test 1e-15 > errc(SF.zeta(3.1+0im,-4.2), SF.zeta(3.1,-4.2+0im))
    @test SF.zeta(3.1,4.2) ≅ 0.029938344862645948405021260567725078588893266227472565010234
    @test SF.zeta(27, 3.1) ≅ 5.413318813037879056337862215066960774064332961282599376e-14
    @test SF.zeta(27, 2) ≅ 7.4507117898354294919810041706041194547190318825658299932e-9
    @test 1e-12 > err(SF.zeta(27, -105.3), 1.3113726525492708826840989036205762823329453315093955e14)

    @test SF.zeta(4, +0.0) == SF.zeta(4, -0.0) ≅ pi^4 / 90
    @test SF.zeta(5, +0.0) == SF.zeta(5, -0.0) ≅ 1.036927755143369926331365486457034168057080919501912811974
    @test SF.zeta(Inf, 1.) == 1
    @test SF.zeta(Inf, 2.) == 0
    @test isnan(SF.zeta(NaN, 1.))

    @test SF.zeta(2 + 1im, -1.1) ≅ SF.zeta(2 + 1im, -1.1+0im) ≅ -64.580137707692178058665068045847533319237536295165484548 + 73.992688148809018073371913557697318846844796582012921247im

    @test SF.zeta(-3.0, 7.0) ≅ -52919/120
    @test SF.zeta(-3.0, -7.0) ≅ 94081/120
    @test SF.zeta(-3.1, 7.2) ≅ -587.457736596403704429103489502662574345388906240906317350719
    @test SF.zeta(-3.1, -7.2) ≅ 1042.167459863862249173444363794330893294733001752715542569576
    @test SF.zeta(-3.1, 7.0) ≅ -518.431785723446831868686653718848680989961581500352503093748
    @test SF.zeta(-3.1, -7.0) ≅ 935.1284612957581823462429983411337864448020149908884596048161
    @test SF.zeta(-3.1-0.1im, 7.2) ≅ -579.29752287650299181119859268736614017824853865655709516268 - 96.551907752211554484321948972741033127192063648337407683877im
    @test SF.zeta(-3.1-0.1im, -7.2) ≅ 1025.17607931184231774568797674684390615895201417983173984531 + 185.732454778663400767583204948796029540252923367115805842138im
    @test SF.zeta(-3.1-0.1im, 7.2 + 0.1im) ≅ -571.66133526455569807299410569274606007165253039948889085762 - 131.86744836357808604785415199791875369679879576524477540653im
    @test SF.zeta(-3.1-0.1im, -7.2 + 0.1im) ≅ 1035.35760409421020754141207226034191979220047873089445768189 + 130.905870774271320475424492384335798304480814695778053731045im
    @test SF.zeta(-3.1-0.1im, -7.0 + 0.1im) ≅ 929.546530292101383210555114424269079830017210969572819344670 + 113.646687807533854478778193456684618838875194573742062527301im
    @test SF.zeta(-3.1, 7.2 + 0.1im) ≅ -586.61801005507638387063781112254388285799318636946559637115 - 36.148831292706044180986261734913443701649622026758378669700im
    @test SF.zeta(-3.1, -7.2 + 0.1im) ≅ 1041.04241628770682295952302478199641560576378326778432301623 - 55.7154858634145071137760301929537184886497572057171143541058im
    @test SF.zeta(-13.4, 4.1) ≅ -3.860040842156185186414774125656116135638705768861917e6
    @test SF.zeta(3.2, -4) ≅ 2.317164896026427640718298719837102378116771112128525719078
    @test SF.zeta(3.2, 0) ≅ 1.166773370984467020452550350896512026772734054324169010977
    @test SF.zeta(-3.2+0.1im, 0.0) ≅ SF.zeta(-3.2+0.1im, 0.0+0im) ≅ 0.0070547946138977701155565365569392198424378109226519905493 + 0.00076891821792430587745535285452496914239014050562476729610im
    @test SF.zeta(-3.2, 0.0) ≅ SF.zeta(-3.2, 0.0+0im) ≅ 0.007011972077091051091698102914884052994997144629191121056378

    @test 1e-14 > err(SF.eta(1+1e-9), 0.693147180719814213126976796937244130533478392539154928250926)
    @test 1e-14 > err(SF.eta(1+5e-3), 0.693945708117842473436705502427198307157819636785324430166786)
    @test 1e-13 > err(SF.eta(1+7.1e-3), 0.694280602623782381522315484518617968911346216413679911124758)
    @test 1e-13 > err(SF.eta(1+8.1e-3), 0.694439974969407464789106040237272613286958025383030083792151)
    @test 1e-13 > err(SF.eta(1 - 2.1e-3 + 2e-3 * im), 0.69281144248566007063525513903467244218447562492555491581+0.00032001240133205689782368277733081683574922990400416791019im)
    @test 1e-13 > err(SF.eta(1 + 5e-3 + 5e-3 * im), 0.69394652468453741050544512825906295778565788963009705146+0.00079771059614865948716292388790427833787298296229354721960im)
    @test 1e-12 > errc(SF.zeta(1e-3+1e-3im), -0.5009189365276307665899456585255302329444338284981610162-0.0009209468912269622649423786878087494828441941303691216750im)
    @test 1e-13 > errc(SF.zeta(1e-4 + 2e-4im), -0.5000918637469642920007659467492165281457662206388959645-0.0001838278317660822408234942825686513084009527096442173056im)

    # Base Julia issue #7169: (TODO: better accuracy should be possible?)
    @test 1e-9 > errc(SF.zeta(0 + 99.69im), 4.67192766128949471267133846066040655597942700322077493021802+3.89448062985266025394674304029984849370377607524207984092848im)
    @test 1e-12 > errc(SF.zeta(3 + 99.69im), 1.09996958148566565003471336713642736202442134876588828500-0.00948220959478852115901654819402390826992494044787958181148im)
    @test 1e-9 > errc(SF.zeta(-3 + 99.69im), 10332.6267578711852982128675093428012860119184786399673520976+13212.8641740351391796168658602382583730208014957452167440726im)
    @test 1e-13 > errc(SF.zeta(2 + 99.69im, 1.3), 0.41617652544777996034143623540420694985469543821307918291931-0.74199610821536326325073784018327392143031681111201859489991im)
end
