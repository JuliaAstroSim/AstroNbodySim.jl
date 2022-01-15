# To generate figures, run codes locally:
# - test
# - examples
# - Benchmark


## Readme


## 01-binary
cp(joinpath(@__DIR__, "../../../../examples/01-binary/output/Binary Circular of every 1 orbit(s).png"),
   joinpath(@__DIR__, "./examples/01-binary/Binary Circular of every 1 orbit(s).png"),
   force=true,
)

cp(joinpath(@__DIR__, "../../../../examples/01-binary/output/Binary Circular with bgforce of every 1 orbit(s).png"),
   joinpath(@__DIR__, "./examples/01-binary/Binary Circular with bgforce of every 1 orbit(s).png"),
   force=true,
)

cp(joinpath(@__DIR__, "../../../../examples/01-binary/output/Elliptic Orbit (const) of every 40 orbit(s).png"),
   joinpath(@__DIR__, "./examples/01-binary/Elliptic Orbit (const) of every 40 orbit(s).png"),
   force=true,
)

cp(joinpath(@__DIR__, "../../../../examples/01-binary/output/Elliptic Orbit (adaptive) of every 40 orbit(s).png"),
   joinpath(@__DIR__, "./examples/01-binary/Elliptic Orbit (adaptive) of every 40 orbit(s).png"),
   force=true,
)

## 02-autodiff-bg
cp(joinpath(@__DIR__, "../../../../examples/02-AutodiffBackground/output/AutodiffBackground.png"),
   joinpath(@__DIR__, "./examples/AutodiffBackground.png"),
   force=true,
)

## 03-Measurements
cp(joinpath(@__DIR__, "../../../../examples/01-binary/output/Uncertainty of elliptic orbit.png"),
   joinpath(@__DIR__, "./examples/01-binary/Uncertainty of elliptic orbit.png"),
   force=true,
)

## 05-plummer
cp(joinpath(@__DIR__, "../../../../examples/03-plummer/output/Plummer-ScaleRadius.png"),
   joinpath(@__DIR__, "./examples/03-plummer/Plummer-ScaleRadius.png"),
   force=true,
)

cp(joinpath(@__DIR__, "../../../../examples/03-plummer/output/Plummer-LagrangianRadii.png"),
   joinpath(@__DIR__, "./examples/03-plummer/Plummer-LagrangianRadii.png"),
   force=true,
)

## 06-GalaxyCollision
cp(joinpath(@__DIR__, "../../../../examples/04-collisions/output/mosaic-collision-DirectSumAdaptiveGPU.png"),
   joinpath(@__DIR__, "./examples/06-GalaxyCollision/mosaic-collision-DirectSumAdaptiveGPU.png"),
   force=true,
)

## 07-TDEcluster
cp(joinpath(@__DIR__, "../../../../examples/05-TDE-StarCluster/output/TDE-elliptic-AccreationHistory.png"),
   joinpath(@__DIR__, "./examples/07-TDEcluster/TDE-elliptic-AccreationHistory.png"),
   force=true,
)

cp(joinpath(@__DIR__, "../../../../examples/05-TDE-StarCluster/output/TDE-elliptic-mosaic.png"),
   joinpath(@__DIR__, "./examples/07-TDEcluster/TDE-elliptic-mosaic.png"),
   force=true,
)
