## Dockerfile to create a reproducible environment for exercises
FROM bioconductor/release_core
MAINTAINER Peter Humburg <peter.humburg@gmail.com>

## Install additional software packages
RUN apt-get update -y && apt-get install -y haskell-platform 
RUN apt-get install -y texlive && apt-get install -y texlive-latex-extra

## Install pandoc
RUN cabal update && cabal install pandoc

## Additional R packages
RUN Rscript -e "biocLite(c('MatrixEQTL', 'devtools', 'trio', 'dplyr', 'ggplot2', 'tidyr', 'knitr', 'optparse'))"
RUN Rscript -e "devtools::install_github('humburg/Rsge');devtools::install_github('jknightlab/mePipe');devtools::install_github('hadley/readr')"

## Data for exercises
ADD data /data

