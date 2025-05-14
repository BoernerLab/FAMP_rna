# Installation and Requirements 

## Docker

The installation of the pipeline is managed through Docker.
Docker can be installed on Linux, macOS, and Windows. We provide pre-built Docker images for all these operating systems. For the best experience, we recommend using <a href="https://docs.docker.com/desktop/" target="_blank"> Docker Desktop</a> to run the pipeline. Below are links to installation guides for Docker Desktop:

- <a href="https://docs.docker.com/desktop/install/linux/" target="_blank">Docker Desktop for Linux</a>
- <a href="https://docs.docker.com/desktop/install/mac-install/" target="_blank">Docker Desktop for macOS</a>
- <a href="https://docs.docker.com/desktop/install/windows-install/" target="_blank">Docker Desktop for Windows</a>

## Our Docker Images

We provide Docker images specifically tailored for machines equipped with graphics cards (GPUs) to enable the execution of molecular dynamics (MD) simulations.

For users without access to GPU-equipped systems, we offer Data_Analysis images. These images allow the post-processing of MD simulations and FRET simulations on standard desktop machines 

A list of available Docker images and their respective system requirements can be found here.

```{dropdown} Here's my dropdown
### Analysis
- Docker 1
- Docker 2
- Docker 3

### GPU Support Linux only?
 
- Docker 1 for GTX ...
- Docker 2 for GTX ....
```


## Running a docker container with an image

````{tab-set}
```{tab-item} Docker Desktop
Wie bekomm ich das Image 
Wie mach ich nen Run 
Was begegnet mir 

```

```{tab-item} Docker CLI
Hier Befehle für Docker CLI einfügen 
```
````

## Installation without Docker

Open a terminal and run the following commands:


- Clone the repository
```bash
git clone https://github.com/felixErichson/FAMP_rna.git
```
- Navigate into the cloned directory
```bash
cd FAMP_rna
```
- Install the package using pip
```bash
pip install .
```