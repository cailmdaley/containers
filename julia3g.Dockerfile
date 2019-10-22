FROM juliapython

USER root
RUN apt-get update && apt-get install -y --no-install-recommends \
        libboost-all-dev \
        libfftw3-dev \
        libflac-dev \
        libgsl0-dev \
        libnetcdf-dev \
        python-openssl \
    && rm -rf /var/lib/apt/lists/*

RUN mkdir $HOME/src
COPY --chown=$USER spt3g_software $HOME/src/spt3g_software
RUN cd $HOME/src/spt3g_software \
    && mkdir build \
    && cd build \
    && cmake -DPYTHON_LIBRARY=~/.pyenv/versions/3.7.3/lib/libpython3.7m.so .. \
    && make -j 8

USER $USER

ENV PORT 8000
RUN mkdir -p /home/$USER/.jupyter \
    && echo "c.NotebookApp.terminado_settings={'shell_command': ['bash']}" >> /home/$USER/.jupyter/jupyter_notebook_config.py

ENTRYPOINT ["src/spt3g_software/build/env-shell.sh"]
CMD jupyter lab --ip=0.0.0.0 --no-browser --port $PORT
