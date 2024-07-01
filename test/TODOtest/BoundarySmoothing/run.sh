#!/bin/bash

module purge
module load InstantaneousInterface

# First cut the overlapped region between mesh and AFP mesh
mesh_op -op CutOverlappedRegion -i target.ply -ref ref.ply -o output.ply -RayDirection 1 0 0 -box 10 5.9 6.33

# Now cut the teeth like region
mesh_op -op CutTeethlikeFace -i output.ply -o cut.ply

# Now smooth
mesh_op -op curvatureflow -box 10 5.9 6.33 -i cut.ply -o smoothed.ply -pbcOutput true -Decimate true -iteration 25 -lambdadt 0.01 -NumBoundarySmoothing 5

# Calculate the curvature
mesh_op -op jetfit -box 10 5.9 6.33 -i smoothed.ply -fc smoothed_fc.out -neighbors 5 -MongeCoefficient 4 -degree 4 

# Map 3d curvature onto 2d
mesh_op -op ProjectCurvature -i smoothed.ply -ProjectedIndex 1 2 -height 100 -n 1000 1000 -L 5.9 6.33 -FaceCurvatureFile smoothed_fc.out -box 10 5.9 6.33 -RayDirection -1 0 0 -o projected.out

