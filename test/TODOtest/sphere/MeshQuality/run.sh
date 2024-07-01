#!/bin/bash

module purge
module load InstantaneousInterface

file=Ice_3200_cut.ply
processfile=Ice_3200_processed.ply
smoothfile=Ice_3200_s.ply
edgeLength=0.1
mesh_op -op MeshCleanup -i ${file} -o ${processfile} -box 10 5.9 6.33 -edgeLengthCutoff ${edgeLength}
mesh_op -op curvatureflow -i ${processfile} -o ${smoothfile} -box 10 5.9 6.33 -iteration 10 -lambdadt 0.1 -pbcOutput true
mesh_op -op jetfit -i ${smoothfile} -box 10 5.9 6.33 -neighbors 4 -MongeCoefficient 4 -degree 4 -o jetfit.out -fc jetfit_fc.out 
mesh_op -op NonPBCFace -i ${smoothfile} -box 10 5.9 6.33 -o nonpbc.out
mesh_op	-op ProjectCurvature -i ${smoothfile} -o projected_jet.out -ProjectedIndex 1 2 -origin 100 -n 200 200 -L 5.9 6.33 -FaceCurvatureFile jetfit_fc.out -box 10 5.9 6.33 -RayDirection -1 0 0

