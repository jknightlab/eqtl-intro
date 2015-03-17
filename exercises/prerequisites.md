# Prerequisites {-}
## Using the Docker image
All exercises assume the use of the docker container `humburg/eqtl-intro`
to provide the required data as well as the necessary software environment.
This requires a working Docker installation^[installation instructions are available
from the [Docker website](https://docs.docker.com/installation/).].
The docker image can be obtained from DockerHub via

```sh
docker pull humburg/eqtl-intro
```
To run the RStudio server run

```sh
docker run -p 8787:8787 humburg/eqtl-intro
```
RStudio is then accessible at `localhost:8787` or, when using
[boot2docker](http://boot2docker.io/) via the IP address indicated 
by `boot2docker ip`.

## Included data
The image includes a number of simulated and real data sets used for these 
exercises. All data are provided as tab-separated files (typically with a column header).
Files are located in directories below `/data`. All simulated data are located in
`/data/simulated`. Real data can be found in `/data/genotyping`, `/data/expression`
and `/data/annotation` for genotyping, gene expression and annotation data 
respectively.