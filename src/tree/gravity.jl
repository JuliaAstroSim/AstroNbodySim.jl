function tree_force_kernel(r2::Number, h::Number, h_inv::Number, h3_inv::Number, mass::Number, dp::AbstractPoint, G::Number)
    r = sqrt(r2)
    if r >= h
        fac = mass / (r2 * r)
    else
        u = r * h_inv
        if ustrip(NoUnits, u) < 0.5
            fac = mass * h3_inv * (10.666666666667 + u^2 * (32.0 * u - 38.4))
        else
            fac = mass * h3_inv * (21.333333333333 - 48.0 * u + 38.4 * u^2 - 10.666666666667 * u^3 - 0.066666666667 / u^3)
        end
    end
    return dp * fac * (1.0 * G)
end

#TODO @code_warntype
function compute_local_force_kernel(pos::AbstractPoint, ptype::Collection, aoldtol::Number, tree::Octree, 
                                    TreeOpenAngle::Float64, G::Number, ForceSofteningTable::StaticArray)
    treenodes = tree.treenodes
    MaxTreenode = tree.config.MaxTreenode

    acc = PVector(unit(aoldtol))

    ninteractions = 0

    h = softlen(ptype, ForceSofteningTable)
    h_inv = 1.0 / h
    h3_inv = h_inv^3

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

        acc += tree_force_kernel(r2, h, h_inv, h3_inv, mass, dp, G)
        ninteractions += 1
    end # while no > 0

    return acc, ninteractions, GravToEvaluate(0, pos, ptype, aoldtol, h), ExportFlag
end

function compute_pseudo_force_kernel(grav::GravToEvaluate, task_id::Int64, sim::Simulation,
                                     TreeOpenAngle::Float64, G::Number, thread_id::Int, gravCat::Dict)
    tree = sim.simdata.tree
    treenodes = tree.treenodes

    MaxTreenode = tree.config.MaxTreenode

    ninteractions = 0
    
    pos, ptype, aold = extract_grav(grav)

    acc = PVector() * zero(grav.OldAccTol)

    h = grav.SoftLen
    h_inv = 1.0 / h
    h3_inv = h_inv^3

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

        acc += tree_force_kernel(r2, h, h_inv, h3_inv, mass, dp, G)
        ninteractions += 1
    end # while no > 0

    push!(gravCat[task_id][thread_id], pack_grav_result(grav, acc, ninteractions))
end

function compute_local_force(sim::Simulation, ::Tree, ::CPU)
    empty!(sim.buffer.gravToEval)
    empty!(sim.buffer.gravToEvalRecv)

    PosType = typeof(sim.config.ZeroValues.pos)
    OldAccType = typeof(sim.config.ZeroValues.acc.x)
    LenType = typeof(sim.config.ZeroValues.pos.x)
    ErrTolAcc = sim.simdata.treesimconfig.ErrTolAcc

    tree = sim.simdata.tree
    gravCat = Dict{Int,Any}()
    for p in tree.pids
        sim.buffer.gravToEval[p] = Vector{GravToEvaluate{PosType, OldAccType, LenType}}()
        gravCat[p] = [Vector{GravToEvaluate{PosType, OldAccType, LenType}}() for i in 1:Threads.nthreads()]
    end

    data = tree.data
    NumForceUpdate = 0


    Threads.@threads for k in 1:Threads.nthreads()
        Head, Tail = split_block(length(data), k, Threads.nthreads())
        for i in Head:Tail
            if data.Ti_endstep[i] == sim.timeinfo.system_time_int
                acc, ninteractions, result, ExportFlag = compute_local_force_kernel(
                    data.Pos[i],
                    data.Collection[i],
                    data.OldAcc[i] * ErrTolAcc,
                    tree,
                    sim.simdata.treesimconfig.TreeOpenAngle,
                    sim.config.constants.G,
                    sim.config.grav.ForceSofteningTable,
                )

                NumForceUpdate += ninteractions
                data.Acc[i] = acc
                data.GravCost[i] = ninteractions

                for p in tree.pids
                    if ExportFlag[p]
                        #TODO: softening length
                        push!(gravCat[p][k], GravToEvaluate(i, result.Pos, result.Collection, result.OldAccTol, result.SoftLen))
                    end
                end
            end
        end
    end
    
    for p in tree.pids
        sim.buffer.gravToEval[p] = vcat(gravCat[p]...)
    end

    sim.physics.NumForceUpdateSinceLast += NumForceUpdate
end

function compute_pseudo_force(sim::Simulation, ::Tree, ::CPU)
    empty!(sim.buffer.gravToEval)
    
    AccType = typeof(sim.config.ZeroValues.acc)
    gravCat = Dict{Int,Any}()
    for p in sim.simdata.tree.pids
        sim.buffer.gravResult[p] = Vector{GravResult{AccType}}()
        gravCat[p] = [Vector{GravResult{AccType}}() for i in 1:Threads.nthreads()]
    end

    for p in keys(sim.buffer.gravToEvalRecv)
        grav = sim.buffer.gravToEvalRecv[p]

        Threads.@threads for k in 1:Threads.nthreads()
            Head, Tail = split_block(length(grav), k, Threads.nthreads())
            for i in Head:Tail
                compute_pseudo_force_kernel(grav[i], p,
                                            sim,
                                            sim.simdata.treesimconfig.TreeOpenAngle,
                                            sim.config.constants.G,
                                            k, gravCat,
                                            )
            end
        end
    end

    for p in sim.simdata.tree.pids
        sim.buffer.gravResult[p] = vcat(gravCat[p]...)
    end

    empty!(sim.buffer.gravToEvalRecv)
end

function combine_force(sim::Simulation, ::Tree, ::CPU)
    empty!(sim.buffer.gravResult)
    data = get_local_data(sim)
    for p in keys(sim.buffer.gravResultRecv)
        grav = sim.buffer.gravResultRecv[p]
        for i in eachindex(grav)
            index = grav[i].Index
            data.Acc[index] += grav[i].Acc
            data.GravCost[index] += grav[i].NInteractions
        end
    end

    empty!(sim.buffer.gravResultRecv)
end

function compute_force(sim::Simulation, GravSolver::Tree, Device::CPU)
    t_FORCE = time_ns()
    bcast(sim, compute_local_force, args = (GravSolver, Device))
    bcast(sim, send_force_buffer)
    bcast(sim, compute_pseudo_force, args = (GravSolver, Device))
    bcast(sim, send_force_result_buffer)
    bcast(sim, combine_force, args = (GravSolver, Device))
    bcast(sim, postprocessing_force, args = (GravSolver, Device))
    add_timer(sim, FORCE, t_FORCE, time_ns())
end

function compute_local_force_at_point_kernel(pos::AbstractPoint3D, tree::Octree,
                                             TreeOpenAngle::Float64, G::Number, h::Number, zeroacc::AbstractPoint)
    treenodes = tree.treenodes

    MaxTreenode = tree.config.MaxTreenode

    acc = zeroacc

    h_inv = 1.0 / h
    h3_inv = h_inv^3

    no = 1
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
                # Tree force computation needs data in other processors, skip
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

        acc += tree_force_kernel(r2, h, h_inv, h3_inv, mass, dp, G)
    end # while no > 0

    return acc
end

function compute_local_force_at_points(sim::Simulation, SoftLength::Number, GravSolver::Tree, Device::CPU)
    tree = sim.simdata.tree
    pos = tree.recvbuffer[1]
    tree.sendbuffer[1] = compute_local_force_at_point_kernel.(pos, tree,
                                    sim.simdata.treesimconfig.TreeOpenAngle,
                                    sim.config.constants.G,
                                    SoftLength,
                                    sim.config.ZeroValues.acc,
                                    )
end

function compute_force(sim::Simulation, pos::Union{Array{T,1}, T}, SoftLength::Number, GravSolver::Tree, Device::CPU) where T<:AbstractPoint3D
    bcast(sim.simdata.tree, :recvbuffer, Dict(1 => pos))
    bcast(sim, compute_local_force_at_points, args = (SoftLength, GravSolver, Device))

    results = gather(sim.simdata.tree, :sendbuffer)
    acc = results[1][1]
    for a in results[2:end]
        acc += a[1]
    end
    return acc
end