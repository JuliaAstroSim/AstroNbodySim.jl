### Static Data Processing Tests
### Tests for find_particle, substract_by_id, diff_by_id on
### StructArrays, Arrays, and distributed Simulations.
###
### Pre-conditions (set up in runtests.jl):
###   - addprocs(4) has been called (procs 1-5 available)
###   - using DataFrames at the top level
###   - @everywhere using AstroNbodySim; @everywhere astro()
###   - @everywhere using Unitful, UnitfulAstro, PhysicalParticles,
###     PhysicalMeshes, PhysicalTrees
###   - exported: find_particle, substract_by_id, diff_by_id
###     (find_particle_local is internal: use AstroNbodySim. prefix)

@testset "StaticDataProcessing" begin

    @testset "find_particle on data" begin
        data = StructArray(Star(uAstro, id = i) for i in 1:5)
        data.Pos[3] = PVector(0.0u"kpc", 1.0u"kpc", 0.0u"kpc")
        data.Pos[5] = PVector(2.0u"kpc", 0.0u"kpc", 0.0u"kpc")

        p3 = find_particle(data, 3)
        @test p3 !== nothing
        @test p3 isa Star
        @test p3.ID == 3
        @test p3.Pos.y ≈ 1.0u"kpc"
        @test p3.Pos.x ≈ 0.0u"kpc"

        p5 = find_particle(data, 5)
        @test p5 !== nothing
        @test p5.Pos.x ≈ 2.0u"kpc"

        @test find_particle(data, 99) === nothing
    end

    @testset "find_particle_local" begin
        data = StructArray(Star(uAstro, id = i) for i in 1:5)
        data.Pos[3] = PVector(0.5u"kpc", 1.0u"kpc", 0.0u"kpc")

        sim = Simulation(
            deepcopy(data);
            pids = [1],
        )

        p3 = AstroNbodySim.find_particle_local(sim, 3)
        @test p3 !== nothing
        @test p3 isa Star
        @test p3.ID == 3
        @test p3.Pos.y ≈ 1.0u"kpc"
        @test p3.Pos.x ≈ 0.5u"kpc"
        @test AstroNbodySim.find_particle_local(sim, 999) === nothing
    end

    @testset "find_particle distributed" begin
        data = StructArray(Star(uAstro, id = i) for i in 1:10)
        data.Mass .= 1.0e9u"Msun"
        data.Pos[7] = PVector(3.0u"kpc", 0.0u"kpc", 0.0u"kpc")

        sim = Simulation(
            deepcopy(data);
            pids = [1, 2, 3, 4],
        )

        p7 = find_particle(sim, 7)
        @test p7 !== nothing
        @test p7 isa Star
        @test p7.ID == 7
        @test p7.Pos.x ≈ 3.0u"kpc"
        @test p7.Pos.y ≈ 0.0u"kpc"
        @test find_particle(sim, 999) === nothing
    end

    @testset "substract_by_id Array" begin
        data = [Star(uAstro, id = i) for i in 1:5]
        @test data isa Array
        @test length(data) == 5

        df = substract_by_id(data, :Acc)
        @test df isa DataFrame
        @test :ID in propertynames(df)
        @test :Acc in propertynames(df)
        @test nrow(df) == 5
        @test df.ID == [1, 2, 3, 4, 5]
        @test df.Acc[1] == data[1].Acc

        # with info = :Pos adds a third column
        df2 = substract_by_id(data, :Acc, info = :Pos)
        @test :Pos in propertynames(df2)
        @test nrow(df2) == 5
        @test df2.Pos[3] == data[3].Pos
    end

    @testset "substract_by_id StructArray" begin
        data = StructArray(Star(uAstro, id = i) for i in 1:5)
        data.Mass[3] = 3.0e10u"Msun"

        df = substract_by_id(data, :Mass)
        @test df isa DataFrame
        @test :ID in propertynames(df)
        @test :Mass in propertynames(df)
        @test nrow(df) == 5
        @test df.ID == [1, 2, 3, 4, 5]
        @test df.Mass[3] ≈ 3.0e10u"Msun"

        df2 = substract_by_id(data, :Mass, info = :Pos)
        @test :Pos in propertynames(df2)
        @test nrow(df2) == 5
    end

    @testset "diff_by_id DataFrame Pos" begin
        data1 = StructArray(Star(uAstro, id = i) for i in 1:5)
        data2 = StructArray(Star(uAstro, id = i) for i in 1:5)
        for i in 1:5
            data2.Pos[i] = PVector(1.0u"kpc", 0.0u"kpc", 0.0u"kpc")
        end

        df1 = substract_by_id(data1, :Pos)
        df2 = substract_by_id(data2, :Pos)
        df = diff_by_id(df1, df2, :Pos)

        @test df isa DataFrame
        @test nrow(df) == 5
        @test df.ID == [1, 2, 3, 4, 5]
        @test :Pos_1 in propertynames(df)
        @test :Pos_2 in propertynames(df)
        @test :diff in propertynames(df)

        norms = [norm(row.diff) for row in eachrow(df)]
        @test all(norms .≈ 1.0u"kpc")
        @test df.diff[3].x ≈ 1.0u"kpc"
    end

    @testset "diff_by_id data overload Mass" begin
        data1 = StructArray(Star(uAstro, id = i) for i in 1:5)
        data1.Mass .= 1.0e10u"Msun"
        data2 = StructArray(Star(uAstro, id = i) for i in 1:5)
        data2.Mass .= 2.0e10u"Msun"

        df = diff_by_id(data1, data2, :Mass)
        @test df isa DataFrame
        @test nrow(df) == 5
        @test :Mass_1 in propertynames(df)
        @test :Mass_2 in propertynames(df)
        @test :diff in propertynames(df)
        @test all(df.diff .≈ 1.0e10u"Msun")
    end

    @testset "diff_by_id missing ID" begin
        data1 = StructArray(Star(uAstro, id = i) for i in 1:5)
        data2 = StructArray(Star(uAstro, id = i) for i in 1:3)
        data1.Mass .= 1.0e9u"Msun"
        data2.Mass .= 1.0e9u"Msun"

        df = diff_by_id(data1, data2, :Mass)
        @test df isa DataFrame
        @test nrow(df) == 3
        @test df.ID == [1, 2, 3]
    end

    @testset "end-to-end Plummer 20" begin
        ic = PlummerStarCluster(NumSamples = 20)
        # NOTE: do NOT pass `MaxRadius = 0.0u"kpc"` here — AstroIC's
        # `plummer.jl:96` enters an infinite loop when MaxRadius <= 0
        # because the resampling condition `while r[i] > MaxRadius` is
        # trivially true for any positive radius. Use the default cutoff
        # (5 * VirialRadius) instead.
        particles = generate(ic)

        ds = Simulation(deepcopy(particles))
        ts = Simulation(deepcopy(particles), GravitySolver = Tree())

        compute_force(ds)
        compute_force(ts)

        dfds = substract_by_id(get_all_data(ds), :Acc, info = :Pos)
        dfts = substract_by_id(get_all_data(ts), :Acc, info = :Pos)
        df = diff_by_id(dfds, dfts, :Acc, info = :Pos)

        @test df isa DataFrame
        @test nrow(df) == 20
        @test :Acc_1 in propertynames(df)
        @test :Acc_2 in propertynames(df)
        @test :diff in propertynames(df)

        # Tree vs DirectSum should be reasonably close
        rel_err = ustrip.(norm.(df.diff)) ./ ustrip.(norm.(df.Acc_1))
        @test all(rel_err .< 1.0)
    end
end
