@testset "calm" begin
    a = 5
    function f1()
        a += 1
        a
    end

    function f2()
        5
    end

    cf1 = FinNetValu.calm(f1, 3)
    cf2 = FinNetValu.calm(f2, 4)

    @test [cf2() for i in 1:10] == [5 for i in 1:10]
    @test [cf1() for i in 1:7]  == [6, 6, 6, 7, 7, 7, 8]
end

@testset "bisection" begin
    @test bisect(x -> abs(x) < 10^-5 ? 0 : sign(x), -2, 2) ≈ 0
    @test bisect(x -> abs(x) < 10^-5 ? 0 : sign(x), -3, 2) ≈ 0
    @test bisect(x -> abs(x) < 10^-5 ? 0 : sign(x), -2, 3) ≈ 0
    @test bisect((x,y) -> abs(x+y) < 10^-12 ? 0 : sign(x+y), -6, 6, 2) ≈ -2
end
