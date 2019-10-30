# parse options
while getopts "n:b:p:i:esr" opt; do
  case $opt in
    n)
      IMAGENAME=$OPTARG
      ;;
    b)
      BUILDFROM=$OPTARG
      ;;
    p)
      PULLFROM=$OPTARG
      ;;
    i)
      INSTANCE=$OPTARG
      ;;
    e)
      CMD=exec
      ;;
    r)
      CMD=run
      ;;
    s)
      CMD=shell
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done

export SINGULARITY_DOCKER_USERNAME=cailmdaley
export SINGULARITY_DOCKER_PASSWORD=2851bdd5-cde6-4043-8bc9-5e81ceecde97

# exit cleanly if no options provided
if [[ $OPTIND -eq 1 ]]; then
    echo "No options provided..." >&2
    exit 1
elif [[ $OPTIND -eq 2 && $IMAGENAME ]]; then
    echo "Image name provided, but no command..." >&2
    exit 1
fi


# various singularity commands
if [[ $BUILDFROM ]]; then
    echo "Building writable docker://cailmdaley/$BUILDFROM as $IMAGENAME.sif" >&2
    singularity instance stop $INSTANCE >&2
    chmod -R 777 $IMAGENAME.sif
    rm -rf $IMAGENAME.sif
    singularity build -s $IMAGENAME.sif docker://cailmdaley/$BUILDFROM >&2
fi
if [[ $PULLFROM ]]; then
    echo "Pulling docker://cailmdaley/$PULLFROM as $IMAGENAME.sif" >&2
    singularity pull $IMAGENAME.sif docker://cailmdaley/$PULLFROM >&2
fi


if [[ $INSTANCE ]]; then
    echo "Starting instance://$INSTANCE from $IMAGENAME.sif" >&2
    singularity instance stop $INSTANCE >&2
    singularity instance start -w $IMAGENAME.sif $INSTANCE >&2
fi


if [[ $CMD ]]; then
    singularity $CMD -w $IMAGENAME.sif >&2
fi
