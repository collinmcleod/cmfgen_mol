#!/bin/csh

#
# Use the procedure to copy files required to generate a new model
# For models of normal stars (i.e., not SN).
#

# Test that we will not inadvertantly copy to the same directory
# We also test that the second directory does not already contain batch.sh,
#   possibly indicating model alredy exists in that directory.
 
if($1 == $2)then
  echo "Warning: directories are the same"
  echo "If you wish to continue enter yes"
  set answer=$<
  switch ($answer)
  case [yY][eE][sS]:
    echo "Continuing with the copy command"
    echo "Moving GREY_SCL_FAC_IN to GREY_SCL_FAC_SAVE"
    echo "You may need to rename this file back to GREY_SCL_FAC_IN"
    echo "Continuing with the copy command"
    cp $2/GREY_SCL_FAC_IN $2/GREY_SCL_FAC_SAVE
    breaksw
  default:
    echo "Aborting copy command"
    exit
  endsw
else
  if(-e $2/batch.sh)then
    echo "Warning: batch.sh exists in the new model directory"
    echo "Enter yes if you wish to continue the copy"
    set answer=$<
    switch ($answer)
      case [yY][eE][sS]:
      echo "Continuing with copy command"
      breaksw
    default:
      echo "Aborting copy command"
      exit
    endsw
  endif
endif


cp $1/batch.sh           $2/
cp $1/IN_ITS             $2/
cp $1/VADAT              $2/
cp $1/*OUT               $2/
cp $1/GAMMAS             $2/GAMMAS_IN
cp $1/MODEL_SPEC         $2/

#

cd $2

# Now rename the *OUT files to *IN

out2in
