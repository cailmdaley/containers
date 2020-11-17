#/bin/bash

[ -z $1 ] && echo "Please provide version tag" && exit
COSCINE_VERSION=$1

docker build \
    -t cailmdaley/coscine:$COSCINE_VERSION \
    -t cailmdaley/coscine:latest \
    coscine

docker build \
    -t cailmdaley/coscine-tweaked:$COSCINE_VERSION \
    -t cailmdaley/coscine-tweaked:latest \
    --build-arg COSCINE_VERSION \
    coscine-tweaked

docker build \
    -t cailmdaley/sptlab:$COSCINE_VERSION \
    -t cailmdaley/sptlab:latest \
    --build-arg COSCINE_VERSION \
    sptlab