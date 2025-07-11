## Dockerfile to create a reproducible environment for exercises
FROM bioconductor/bioconductor_docker:latest
LABEL maintainer="Peter Humburg" email="peter.humburg@gmail.com"

## Install additional software packages
RUN apt-get update -y && apt-get install -y haskell-platform 
RUN apt-get install -y texlive && apt-get install -y texlive-latex-extra

## Install pandoc
RUN cabal update && cabal install pandoc

## Additional R packages
RUN Rscript -e "BiocManager::install(c('MatrixEQTL', 'devtools', 'trio', 'dplyr', 'ggplot2', 'tidyr', 'knitr', 'optparse', 'readr'))"
RUN Rscript -e "devtools::install_github('cran/lmSupport');devtools::install_github('humburg/Rsge');devtools::install_github('jknightlab/mePipe')"

## Data for exercises
ADD data /data

