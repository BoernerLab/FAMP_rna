# Installation and Requirements 

## Docker

The installation of the pipeline is managed through Docker. Docker provides a containerized environment that allows for the seamless and consistent deployment of applications across various systems. By using containers, Docker encapsulates all the dependencies, libraries, and settings required for the pipeline into isolated units, ensuring that the pipeline runs uniformly, regardless of the underlying operating system. This reduces compatibility issues and simplifies the setup process for users.

Docker can be installed on Linux, macOS, and Windows. We provide pre-built Docker images for all these operating systems. For the best experience, we recommend using <a href="https://docs.docker.com/desktop/" target="_blank"> Docker Desktop</a> to run the pipeline. Below are links to installation guides for Docker Desktop:

- <a href="https://docs.docker.com/desktop/install/linux/" target="_blank">Docker Desktop for Linux</a>
- <a href="https://docs.docker.com/desktop/install/mac-install/" target="_blank">Docker Desktop for macOS</a>
- <a href="https://docs.docker.com/desktop/install/windows-install/" target="_blank">Docker Desktop for Windows</a>

## Our Docker Images

A Docker Image is a lightweight, standalone, and executable package that includes everything needed to run a piece of software, including the code, runtime, libraries, environment variables, and dependencies. These images ensure consistency and portability, allowing the pipeline to run the same way across different systems.

We provide Docker images specifically tailored for machines equipped with graphics cards (GPUs) to enable the execution of molecular dynamics (MD) simulations.

For users without access to GPU-equipped systems, we offer Data_Analysis images. These images allow the post-processing of MD simulations and FRET simulations on standard machines without the need for a GPU. This makes it possible to use the pipeline on more modest computing resources.

A list of available Docker images and their respective system requirements can be found here.
```{dropdown} Here's my dropdown
jhgfdsdfghjk
```
::::{grid}
:gutter: 2

:::{grid-item}
:outline:
A
:::
:::{grid-item}
:outline:
B
:::
:::{grid-item}
:outline:
C
:::
:::{grid-item}
:outline:
D
:::

::::

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


# OLD

The code can be downloaded under the green button "Code" or cloned via the terminal:
```
git clone https://github.com/felixErichson/FAMP_rna.git
```
<br>
