#!/bin/csh

# Files required to generate a new model

cp $1/batch.sh     $2/
cp $1/IN_ITS       $2/
cp $1/VADAT        $2/
cp $1/*OUT         $2/
cp $1/GAMMAS       $2/GAMMAS_IN
cp $1/MODEL_SPEC   $2/

cd $2

# Now rename the *OUT files to *IN

out2in
