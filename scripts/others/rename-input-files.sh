#!/bin/bash

# Rename all simulation input files
# To be run in root ecogen folder using ./scripts/others/rename-input-files.sh

find . -type f -name 'mainV5.xml' -exec rename 's/mainV5.xml/main.xml/' {} +
find . -type f -name 'meshV5.xml' -exec rename 's/meshV5.xml/mesh.xml/' {} +
find . -type f -name 'modelV4.xml' -exec rename 's/modelV4.xml/model.xml/' {} +
find . -type f -name 'initialConditionsV4.xml' -exec rename 's/initialConditionsV4.xml/initialConditions.xml/' {} +