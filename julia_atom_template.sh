export JULIA_NUM_THREADS=4

#instance method: probably less good
# uncomment to reset instance
# singularity instance stop julia-atom
# singularity instance start containers/julia.sif julia-atom
singularity exec -w instance://julia-atom julia -O 3 "$@"
