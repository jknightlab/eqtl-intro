# Introduction to eQTL analysis
This repository contains the materials for an introductory 
eQTL analysis course.

## Docker image
The *docker* directory contains the Dockerfile and data 
necessary to build a Docker image that contains all necessary 
components for the exercises that form part of this course.
The resulting Docker image is available from 
[DockerHub](https://registry.hub.docker.com/u/humburg/eqtl-intro/).
`docker pull humburg/eqtl-intro` will pull the latest version
of the image.

### Installing Docker
Using this image requires a working Docker installation.
Installation instructions for [Windows](https://docs.docker.com/installation/windows/),
[Mac](https://docs.docker.com/installation/mac/), 
[Ubuntu](https://docs.docker.com/installation/ubuntulinux/) and a number
of [other platforms](https://docs.docker.com/installation/) are available
through the Docker [documentation](https://docs.docker.com/).

### Using the Docker image
The docker image contains [R](http://www.r-project.org/) and 
[RStudio server](http://www.rstudio.com/) and allows use of the RStudio IDE
as well as the R console. To use the IDE start the Docker container with

```
docker run -p 8787:8787 humburg/eqtl-intro
``` 
To successfully complete the exercises it is necessary to share files
with the local machine. Instructions on how to do this on Windows, Mac and Linux
are available [here](https://github.com/rocker-org/rocker/wiki/Sharing-files-with-host-machine).
