#!/bin/bash

target=build/FEMAddOns

cp -r LICENSE									$target/
cp -r README.md									$target/

cp -r FEMAddOns/DistMesh/COPYING				$target/DistMesh
cp -r FEMAddOns/DistMesh/Kernel					$target/DistMesh

cp -r FEMAddOns/DomainDecomposition/Kernel		$target/DomainDecomposition

