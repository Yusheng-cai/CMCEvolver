#!/bin/bash

module purge
module load InstantaneousInterface

mkdir data
for i in {1..1000}
do
	InstantaneousInterface SampleInputZac.dat
done
