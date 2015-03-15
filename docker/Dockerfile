## Dockerfile to create a reproducible environment for exercises
FROM bioconductor/release_core
MAINTAINER Peter Humburg <peter.humburg@gmail.com>

## Additional R packages
RUN Rscript -e "biocLite(c('MatrixEQTL', 'devtools', 'trio'))"
RUN Rscript -e "devtools::install_github('humburg/Rsge');devtools::install_github('jknightlab/mePipe');"

## Data for exercises
ADD data /data
