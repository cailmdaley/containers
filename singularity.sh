IMAGENAME=juliabase
INSTANCENAME=hw4
SIF="$IMAGENAME"_latest

INSTANCE=instance://$NAME

# singularity pull docker://cailmdaley/$IMAGENAME
echo singularity instance start "$IMAGENAME"_latest.sif $INSTANCENAME
singularity instance start "$IMAGENAME"_latest.sif $INSTANCENAME
