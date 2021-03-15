#!/bin/bash

COMMAND=$1
OUTPUT_DIR=$(dirname $0)
FILE=$OUTPUT_DIR/$COMMAND

mkdir -p $OUTPUT_DIR

cat > $FILE <<- EOM
#!/bin/bash

docker exec -it -w \$PWD sptlab $COMMAND "\$@"
EOM

chmod +x $FILE
