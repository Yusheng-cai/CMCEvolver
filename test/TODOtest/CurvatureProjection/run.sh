#!/bin/bash


module purge
module load InstantaneousInterface

mesh_op -op ProjectCurvature -i test.ply -o c.out -RayDirection -1 0 0 -ProjectedIndex 1 2 -origin 100 -n 400 400 -L 5.9 6.33 -FaceCurvatureFile faceC.out -FileColumn 1 -box 10 5.9 6.33 -refMesh sbwAFP_cut.ply


