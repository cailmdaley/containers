ARG COSCINE_VERSION=latest
FROM cailmdaley/coscine:$COSCINE_VERSION


# install nodejs for jupyter lab extensions, respectively
USER root
RUN apt-get update && apt-get install -y --no-install-recommends \
        nodejs \
        npm


### Document Export
################################################################################

# install pandoc:
ENV PANDOC_VERSION=2.11.0.2
RUN wget https://github.com/jgm/pandoc/releases/download/${PANDOC_VERSION}/pandoc-${PANDOC_VERSION}-1-amd64.deb \
    && dpkg -i pandoc-${PANDOC_VERSION}-1-amd64.deb \
    && rm pandoc-${PANDOC_VERSION}-1-amd64.deb
ENV XDG_DATA_HOME=$SRC

RUN chown -R $IMAGE_USER $SRC \
    && fix-permissions $SRC
USER $IMAGE_USER

# install TinyTeX (minimal latex installation) to $TINYTEX_DIR
ENV TINYTEX_DIR=$SRC
RUN wget -qO- "https://yihui.org/tinytex/install-bin-unix.sh"  | sh \
    && $SRC/.TinyTeX/bin/*/tlmgr option sys_bin $SRC/bin \
    && $SRC/.TinyTeX/bin/*/tlmgr path add \ 
    # install useful packages & pandoc pdf dependencies
    && tlmgr install \
        biblatex \
        caption \
        collection-fontsrecommended \
        dvipng \
        epsf \
        epstopdf-pkg \
        fancyhdr \
        grfext \
        lm-math \
        mathtools \
        microtype \
        physics \
        revtex4-1 \
        savesym \
        siunitx \
        textcase \
        ulem \
        unicode-math \
        wrapfig \
        xcolor

# ENV CROSSREF_VERSION=0.3.8.3
# RUN mkdir -p $SRC/pandoc/filters \
#     && curl -L https://github.com/lierdakil/pandoc-crossref/releases/download/v${CROSSREF_VERSION}/pandoc-crossref-Linux-${PANDOC_VERSION}.tar.xz \
#         | tar -xJf - -C $SRC/pandoc/filters pandoc-crossref \
#     && git clone https://github.com/cailmdaley/writing_tools.git $SRC/writing_tools


# # jupyter lab extensions
# COPY --chown=1000 jupyter_notebook_config.py $JUPYTER_CONFIG_DIR/jupyter_notebook_config.py
# # add jupyterlab-git back eventually:
# RUN pip install nbdime jupyter_contrib_nbextensions \
#     && jupyter contrib nbextension install --user \
#     && jupyter labextension install @axlair/jupyterlab_vim
