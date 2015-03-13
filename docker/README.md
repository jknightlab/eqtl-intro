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
RStudio is then available through your web browser by accessing the 
docker host on port 8787. On Linux the docker host is 127.0.0.1 (or localhost) 
by default. So the full URL to RStudio is http://localhost:8787. On Mac or Windows,
running boot2docker, you can determine the docker host with the `boot2docker ip` command.
Log in to RStudio with the username *rstudio* and password *rstudio*.

To successfully complete the exercises it is necessary to share files
with the local machine. Instructions on how to do this on Windows, Mac and Linux
are available [here](https://github.com/rocker-org/rocker/wiki/Sharing-files-with-host-machine).
