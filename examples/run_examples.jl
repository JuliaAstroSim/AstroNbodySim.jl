function run_examples(;
    examples = [
        "01-binary",
        "02-AutodiffBackground",
        "03-plummer",
        "04-collision",
        "05-TDE-StarCluster",
        
        "07-solarsystem",
    ],
    updatedoc = false,
)
    printstyled("Running examples of AstroNbodySim project", color = :cyan)

    for example in examples
        printstyled("Running example $example...", color = :green)

        if !ispath(joinpath(example, "output"))
            println("Making dir $example/output")
            mkdir(joinpath(example, "output"))
        end

        println("Entering dir $example")
        cd(example)


        run(`julia --color=yes --banner=no $(example).jl`)

        if isfile("output/success")
            printstyled("Example $example done...", color = :green)
        else
            printstyled("Example $example failed...", color = :red)
            break
        end
        cd("../")
    end

    if updatedoc
        #TODO update example pics to doc/
    end
end

run_examples()