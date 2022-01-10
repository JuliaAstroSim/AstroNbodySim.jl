function tree_potential_kernel(r2::Number, h::Number, h_inv::Number, mass::Number, G::Number)
    r = sqrt(r2)
    if r >= h
        pot = - mass / r * (1.0 * G)
    else
        u = r * h_inv
        if ustrip(NoUnits, u) < 0.5
            wp = -2.8 + u * u * (5.333333333333 + u * u * (6.4 * u - 9.6))
        else
            wp = -3.2 + 0.066666666667 / u + u * u * (10.666666666667 + u * (-16.0 + u * (9.6 - 2.133333333333 * u)))
        end

        pot = mass * h_inv * wp * (1.0 * G)
    end
    return pot
end

function compute_local_potential_kernel(pos::AbstractPoint, ptype::Collection, aoldtol::Number, tree::Octree, 
                                        TreeOpenAngle::Float64, G::Number, ForceSofteningTable::StaticArray, zeropotpermass::Number)
    treenodes = tree.treenodes
    MaxTreenode = tree.config.MaxTreenode

    pot = zeropotpermass

    h = softlen(ptype, ForceSofteningTable)
    h_inv = 1.0 / h

    no = 1

    ExportFlag = Dict{Int64, Bool}()
    for p in tree.pids
        ExportFlag[p] = false
    end

    target_task = 0
    while no > 0
        if no <= MaxTreenode && treenodes[no].IsAssigned  # single particle / assigned tree leaf
            dp = treenodes[no].MassCenter - pos
            if iszero(norm(dp))      # same particle
                no = treenodes[no].Sibling
                continue
            end
            mass = treenodes[no].Mass
        else
            if no > MaxTreenode  # Pseudo particle
                # Tree force computation needs data in other processors
                # We send the particle data to there and receive the force result later

                target_task = tree.domain.DomainTask[no - MaxTreenode]
                ExportFlag[target_task] = true

                no = tree.NextNodes[no - MaxTreenode]
                continue
            end

            # Tree nodes
            nop = treenodes[no]
            dp = nop.MassCenter - pos
            mass = nop.Mass
        end

        r2 = dp * dp

        # Find the next no
        if treenodes[no].IsAssigned
            no = treenodes[no].Sibling
        else # tree branch

            if TreeOpenAngle > 0.0 # check Barnes-Hut opening criterion
                if nop.SideLength^2 > r2 * TreeOpenAngle^2
                    no = nop.NextNode # yes, open it
                    continue
                end
            else  # check relative opening criterion
                if mass * nop.SideLength^2 > r2^2 * aoldtol
                    no = nop.NextNode
                    continue
                end

                # check in addition whether we lie inside the cell
                if abs(nop.Center.x - pos.x) < 0.6 * nop.SideLength &&
                    abs(nop.Center.y - pos.y) < 0.6 * nop.SideLength &&
                    abs(nop.Center.z - pos.z) < 0.6 * nop.SideLength
                    no = nop.NextNode
                    continue
                end
            end

            no = nop.Sibling # We can use this node, next node is its sibling
        end # if treenodes[no].IsAssigned

        pot += tree_potential_kernel(r2, h, h_inv, mass, G)
    end # while no > 0

    return pot, GravToEvaluate(0, pos, ptype, aoldtol, h), ExportFlag
end

function compute_pseudo_potential_kernel(grav::GravToEvaluate, task_id::Int64, sim::Simulation, 
                                         TreeOpenAngle::Float64, G::Number, zeropotpermass::Number, thread_id::Int, gravCat::Dict)
    tree = sim.simdata.tree
    treenodes = tree.treenodes

    MaxTreenode = tree.config.MaxTreenode

    pos, ptype, aold = extract_grav(grav)

    pot = zeropotpermass

    h = grav.SoftLen
    h_inv = 1.0 / h

    no = 1

    while no > 0
        if no <= MaxTreenode && treenodes[no].IsAssigned  # single particle / assigned tree leaf
            dp = treenodes[no].MassCenter - pos
            mass = treenodes[no].Mass
        else
            if no > MaxTreenode  # Pseudo particle
                no = tree.NextNodes[no - MaxTreenode]
                continue
            end

            # Tree nodes
            nop = treenodes[no]
            dp = nop.MassCenter - pos
            mass = nop.Mass
        end

        r2 = dp * dp

        # Find the next no
        if treenodes[no].IsAssigned
            no = treenodes[no].Sibling
        else # tree branch
            # If it's a top-level node which does not contain local particles
            # continue to do a short-cut
            if (nop.BitFlag & 3) == 1
                no = nop.Sibling
                continue
            end

            if TreeOpenAngle > 0.0 # check Barnes-Hut opening criterion
                if nop.SideLength^2 > r2 * TreeOpenAngle^2
                    no = nop.NextNode
                    continue
                end
            else  # check relative opening criterion
                if mass * nop.SideLength^2 > r2^2 * aold
                    no = nop.NextNode
                    continue
                end

                # check in addition whether we lie inside the cell
                if abs(nop.Center.x - pos.x) < 0.6 * nop.SideLength &&
                    abs(nop.Center.y - pos.y) < 0.6 * nop.SideLength &&
                    abs(nop.Center.z - pos.z) < 0.6 * nop.SideLength
                    no = nop.NextNode
                    continue
                end
            end

            no = nop.Sibling
            if (nop.BitFlag & 1) > 0
                continue
            end
        end # if treenodes[no].IsAssigned

        pot += tree_potential_kernel(r2, h, h_inv, mass, G)
    end # while no > 0

    push!(gravCat[task_id][thread_id], pack_pot_result(grav, pot))
end

function compute_local_potential(sim::Simulation, ::Tree, ::CPU)
    empty!(sim.buffer.gravToEval)
    empty!(sim.buffer.gravToEvalRecv)
    data = sim.simdata.tree.data
    
    PosType = typeof(sim.config.ZeroValues.pos)
    OldAccType = typeof(sim.config.ZeroValues.acc.x)
    LenType = typeof(sim.config.ZeroValues.pos.x)
    ErrTolAcc = sim.simdata.treesimconfig.ErrTolAcc
    zeropotpermass = sim.config.ZeroValues.potpermass

    gravCat = Dict{Int,Any}()
    for p in sim.pids
        sim.buffer.gravToEval[p] = Vector{GravToEvaluate{PosType, OldAccType, LenType}}()
        gravCat[p] = [Vector{GravToEvaluate{PosType, OldAccType, LenType}}() for i in 1:Threads.nthreads()]
    end

    Threads.@threads for k in 1:Threads.nthreads()
        Head, Tail = split_block(length(data), k, Threads.nthreads())
        for i in Head:Tail
            @inbounds pot, result, ExportFlag = compute_local_potential_kernel(
                data.Pos[i],
                data.Collection[i],
                data.OldAcc[i] * ErrTolAcc,
                sim.simdata.tree,
                sim.simdata.treesimconfig.TreeOpenAngle,
                sim.config.constants.G,
                sim.config.grav.ForceSofteningTable,
                zeropotpermass,
            )

            @inbounds data.Potential[i] = pot
            for p in sim.pids
                if ExportFlag[p]
                    #TODO: softening length
                    push!(gravCat[p][k], GravToEvaluate(i, result.Pos, result.Collection, result.OldAccTol, result.SoftLen))
                end
            end
        end
    end

    for p in sim.pids
        sim.buffer.gravToEval[p] = vcat(gravCat[p]...)
    end
end

function compute_pseudo_potential(sim::Simulation, ::Tree, ::CPU)
    empty!(sim.buffer.gravToEval)

    gravCat = Dict{Int,Any}()
    PotType = typeof(sim.config.ZeroValues.potpermass)
    for p in sim.pids
        sim.buffer.potResult[p] = Vector{PotentialResult{PotType}}()
        gravCat[p] = [Vector{PotentialResult{PotType}}() for i in 1:Threads.nthreads()]
    end

    zeropotpermass = sim.config.ZeroValues.potpermass
    for p in keys(sim.buffer.potToEvalRecv)
        grav = sim.buffer.potToEvalRecv[p]
        Threads.@threads for k in 1:Threads.nthreads()
            Head, Tail = split_block(length(grav), k, Threads.nthreads())
            for i in Head:Tail
                @inbounds compute_pseudo_potential_kernel(grav[i], p,
                                                          sim,
                                                          sim.simdata.treesimconfig.TreeOpenAngle,
                                                          sim.config.constants.G,
                                                          zeropotpermass,
                                                          k, gravCat,
                                                          )
            end
        end
    end

    for p in sim.pids
        sim.buffer.potResult[p] = vcat(gravCat[p]...)
    end

    empty!(sim.buffer.potToEvalRecv)
end

function combine_potential(sim::Simulation, ::Tree, ::CPU)
    empty!(sim.buffer.potResult)
    data = sim.simdata.tree.data

    for p in keys(sim.buffer.potResultRecv)
        grav = sim.buffer.potResultRecv[p]
        for i in eachindex(grav)
            index = grav[i].Index
            data.Potential[index] += grav[i].Potential * data.Mass[index]
        end
    end

    # Background potential
    for f in sim.bgpotential
        apply_background_potential(f, data)
    end
    empty!(sim.buffer.potResultRecv)
end

function compute_potential(sim::Simulation, GravSolver::Tree, Device::CPU)
    bcast(sim, compute_local_potential, args = (GravSolver, Device))
    bcast(sim, send_force_buffer)
    bcast(sim, compute_pseudo_potential, args = (GravSolver, Device))
    bcast(sim, send_potential_result_buffer)
    bcast(sim, combine_potential, args = (GravSolver, Device))
end

function compute_local_unit_potential_at_point_kernel(pos::AbstractPoint3D, tree::Octree,
                                                 ErrTolAcc::Float64, TreeOpenAngle::Float64, G::Number, ZeroValues::ZeroValue, h::Number)
    treenodes = tree.treenodes

    MaxTreenode = tree.config.MaxTreenode
    pot = 0.0 * ZeroValues.potpermass

    h_inv = 1.0 / h

    no = 1

    target_task = 0
    while no > 0
        if no <= MaxTreenode && treenodes[no].IsAssigned  # single particle / assigned tree leaf
            dp = treenodes[no].MassCenter - pos
            if iszero(norm(dp))      # same particle
                no = treenodes[no].Sibling
                continue
            end
            mass = treenodes[no].Mass
        else
            if no > MaxTreenode  # Pseudo particle
                # Tree force computation needs data in other processors
                # We send the particle data to there and receive the force result later

                target_task = tree.domain.DomainTask[no - MaxTreenode]

                no = tree.NextNodes[no - MaxTreenode]
                continue
            end

            # Tree nodes
            nop = treenodes[no]
            dp = nop.MassCenter - pos
            mass = nop.Mass
        end

        r2 = dp * dp

        # Find the next no
        if treenodes[no].IsAssigned
            no = treenodes[no].Sibling
        else # tree branch

            if TreeOpenAngle > 0.0 # check Barnes-Hut opening criterion
                if nop.SideLength^2 > r2 * TreeOpenAngle^2
                    no = nop.NextNode # yes, open it
                    continue
                end
            else  # check relative opening criterion

                # check in addition whether we lie inside the cell
                if abs(nop.Center.x - pos.x) < 0.6 * nop.SideLength &&
                    abs(nop.Center.y - pos.y) < 0.6 * nop.SideLength &&
                    abs(nop.Center.z - pos.z) < 0.6 * nop.SideLength
                    no = nop.NextNode
                    continue
                end
            end

            no = nop.Sibling # We can use this node, next node is its sibling
        end # if treenodes[no].IsAssigned

        pot += tree_potential_kernel(r2, h, h_inv, mass, G)
    end # while no > 0

    return pot
end

function compute_local_unit_potential_at_points(sim::Simulation, SoftLength::Number, GravSolver::Tree, Device::CPU)
    pos = sim.buffer.recvbuffer[1]
    sim.buffer.sendbuffer[1] = compute_local_unit_potential_at_point_kernel.(pos, sim.simdata.tree,
                                    sim.simdata.treesimconfig.ErrTolAcc,
                                    sim.simdata.treesimconfig.TreeOpenAngle,
                                    sim.config.constants.G,
                                    sim.config.ZeroValues,
                                    SoftLength)
end

function compute_potential(sim::Simulation, pos::Union{Array{T,1}, T}, SoftLength::Number, GravSolver::Tree, Device::CPU) where T<:AbstractPoint3D
    bcast(sim, :buffer, :recvbuffer, Dict(1 => pos))
    bcast(sim, compute_local_unit_potential_at_points, args = (SoftLength, GravSolver, Device))

    results = gather(sim, :buffer, :sendbuffer)
    pot = results[1][1]
    for p in results[2:end]
        pot += p[1]
    end
    return pot
end

function total_potential(sim::Simulation, GravSolver::Tree, Device::CPU)
    compute_potential(sim, GravSolver, Device)
    return 0.5 * sum(filter(!isnothing, gather(sim.simdata.tree, total_potential, :data)))
end