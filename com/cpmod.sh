#!/bin/csh

# Files required to generate a new model

cp $1/batch.sh           $2/
cp $1/IN_ITS             $2/
cp $1/VADAT              $2/
cp $1/*OUT               $2/
cp $1/GAMMAS             $2/GAMMAS_IN
cp $1/MODEL_SPEC         $2/

#
# Copy GREY_SCL_FAC file across, if it exists
#

#test -e $1/OUT_GREY_SCL_FAC && cp $1/OUT_GREY_SCL_FAC   $2/GREY_SCL_FAC

cd $2

# Now rename the *OUT files to *IN

out2in
