#!/bin/bash

target=build/FEMAddOns

cp -r License.md								$target/
cp -r Readme.md									$target/

cp -r FEMAddOns.m								$target/

cp -r FEMAddOns/DistMesh/COPYING				$target/DistMesh
cp -r FEMAddOns/DistMesh/Kernel					$target/DistMesh

cp -r FEMAddOns/DomainDecomposition/COPYING		$target/DomainDecomposition
cp -r FEMAddOns/DomainDecomposition/Kernel		$target/DomainDecomposition

cp -r FEMAddOns/FEMUtils/COPYING				$target/FEMUtils
cp -r FEMAddOns/FEMUtils/Kernel					$target/FEMUtils

