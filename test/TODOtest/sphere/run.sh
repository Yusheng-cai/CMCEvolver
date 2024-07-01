#!/bin/bash

mesh_op -op curvefit -i sample_sphere.ply -neighbors 2 -o curvefit.out -fc curvefit_face.out
mesh_op -op jetfit -i sample_sphere.ply -neighbors 2 -degree 4 -MongeCoefficient 4 -o jetfit.out -fc jetfit_face.out
mesh_op -op tensorfit -i sample_sphere.ply -o tensorfit.out -fc tensorfit_face.out
mesh_op -op FDMfit -i sample_sphere.ply -o FDMfit.out -fc FDMfit_face.out
