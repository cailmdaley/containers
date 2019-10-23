FROM cailmdaley/juliapython:latest

USER root
RUN apt-get update && apt-get install -y --no-install-recommends \
        libboost-all-dev \
        libfftw3-dev \
        libflac-dev \
        libgsl0-dev \
        libnetcdf-dev \
        python-openssl \
    && rm -rf /var/lib/apt/lists/*

USER $IMAGEUSER
COPY --chown=1000 spt3g_software $HOME/src/spt3g_software
RUN cd $HOME/src/spt3g_software \
    && mkdir build \
    && cd build \
    && cmake -DPYTHON_LIBRARY=~/.pyenv/versions/3.7.3/lib/libpython3.7m.so .. \
    && make \
    && chmod -R 777 $HOME/src/spt3g_software

ENV PORT 8000
RUN mkdir -p /home/$IMAGEUSER/.jupyter \
    && echo "c.NotebookApp.terminado_settings={'shell_command': ['bash']}" \
    >> /home/$IMAGEUSER/.jupyter/jupyter_notebook_config.py \
    && chmod -R 777 home/$IMAGEUSER/.jupyter

CMD $HOME/src/spt3g_software/build/env-shell.sh
# jupyter lab --ip=0.0.0.0 --no-browser --port $PORT
