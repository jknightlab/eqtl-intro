## Dockerfile to create a reproducible environment for exercises
FROM bioconductor/release_core
MAINTAINER Peter Humburg <peter.humburg@gmail.com>

## Additional R packages
RUN Rscript -e "biocLite('MatrixEQTL')"
RUN Rscript -e "devtools::install_github('Rsge', user='humburg');devtools::install_github('mePipe', user='jknightlab');"
