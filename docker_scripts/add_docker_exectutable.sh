#!/bin/bash

COMMAND=$1
CONTAINER="${2:-sptlab}"
OUTPUT_DIR=$(dirname $0)
FILE=$OUTPUT_DIR/$COMMAND

mkdir -p $OUTPUT_DIR

cat > $FILE <<- EOM
#!/bin/bash

docker exec -it -w \$PWD $CONTAINER $COMMAND "\$@"
EOM

chmod +x $FILE
