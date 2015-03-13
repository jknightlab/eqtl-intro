# Environment for *Introduction to eQTL Analysis* exercises
This Docker image is build on the Biconductor base image and
is intended to provide a controlled environment for the 
exercises accompanying the course *Introduction to eQTL analysis*.

## Using the Docker image
The docker image contains [R](http://www.r-project.org/) and 
[RStudio server](http://www.rstudio.com/) and allows use of the RStudio IDE
as well as the R console. To use the IDE start the Docker container with

```
docker run -p 8787:8787 humburg/eqtl-intro
``` 
To successfully complete the exercises it is necessary to share files
with the local machine. Instructions on how to do this on Windows, Mac and Linux
are available [here](https://github.com/rocker-org/rocker/wiki/Sharing-files-with-host-machine).
