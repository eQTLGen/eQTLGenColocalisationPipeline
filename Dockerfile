FROM nfcore/base
LABEL authors="Urmo Võsa" \
      description="Docker image containing tools for colocalisation analyses with HyprColoc"

COPY GwasHyprColoc.yml /
RUN apt-get update && apt install -y libgmp-dev && apt install -y build-essential
RUN conda env create -f GwasHyprColoc.yml && conda clean -a
ENV PATH /opt/conda/envs/GwasHyprColoc/bin:$PATH
ENV TAR="/bin/tar"
RUN ln -s /bin/tar /bin/gtar
RUN R -e "devtools::install_github('jrs95/hyprcoloc', build_opts = c('--no-resave-data', '--no-manual', '--no-build-vignettes'), build_vignettes = FALSE)"
