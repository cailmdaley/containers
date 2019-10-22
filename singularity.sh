# parse options
while getopts "n:v:p:sid" opt; do
  case $opt in
    n)
      INSTANCENAME=$OPTARG
      INSTANCE=instance://$INSTANCENAME
      ;;
    v)
      VERSION=$OPTARG
      ;;
    p)
      IMAGE=$OPTARG
      ;;
    s)
      START=1
      ;;
    i)
      INTERACTIVE=1
      ;;
    d)
      STOP=1
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done

# exit cleanly if no options provided
if [[ $OPTIND -eq 1 ]]; then
    echo "No options provided..." >&2
    exit 1
elif [[ $OPTIND -eq 2 && $INSTANCE ]]; then
    echo "Instance name provided, but no command..." >&2
    exit 1
fi

# various singularity commands
if [[ $IMAGE ]]; then
    echo "Saving docker://cailmdaley/$IMAGE as $INSTANCENAME.sif" >&2
    singularity pull $INSTANCENAME.sif docker://cailmdaley/$IMAGE:$VERSION >&2
fi

if [[ $START ]]; then
    echo "Starting $INSTANCENAME.sif" >&2
    singularity instance start "$INSTANCENAME".sif $INSTANCENAME >&2
fi

if [[ $INTERACTIVE ]]; then
    if ! singularity shell $INSTANCE >&2; then
        echo "Couldn't shell $INSTANCE. Running instances:" >&2
        singularity instance list >&2
        exit 1
    fi
elif [[ $STOP ]]; then
    if ! singularity instance stop $INSTANCENAME >&2; then
        echo "Couldn't stop $INSTANCE. Running instances:" >&2
        singularity instance list >&2
        exit 1
    fi
fi
