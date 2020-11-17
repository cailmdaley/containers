#/bin/bash

[ -z $1 ] && echo "Please provide version tag" && exit
COSCINE_VERSION=$1

docker push cailmdaley/coscine:$COSCINE_VERSION 
docker push cailmdaley/coscine:latest

docker push cailmdaley/coscine-tweaked:$COSCINE_VERSION 
docker push cailmdaley/coscine-tweaked:latest

docker push cailmdaley/sptlab:$COSCINE_VERSION 
docker push cailmdaley/sptlab:latest