# Installation Guide

This page provides instructions to install and run the FAMP pipeline. You can either build the required Docker containers yourself using the provided Dockerfiles, or install FAMP directly using `pip` in a local Python environment.

---

## Installation with Docker (recommended)

FAMP includes Dockerfiles to build containerized environments for modeling, simulation, and data analysis. These Dockerfiles allow reproducible and portable setups across systems, including containers with and without GPU support for MD simulations.

```{Note}
We do not publish pre-built Docker images. You must build all containers locally from the provided Dockerfiles due to licensing.  
However, we provide several almost ready-to-use Dockerfile templates designed for straightforward customization with minimal modifications.
```

---

### Installation via Docker Desktop (for inexperienced users)

1. Install Docker Desktop (if not yet installed)
   
   - <a href="https://docs.docker.com/desktop/install/linux/" target="_blank" >Docker Desktop for Linux</a>
   - <a href="https://docs.docker.com/desktop/install/mac-install/" target="_blank">Docker Desktop for macOS</a>
   - <a href="https://docs.docker.com/desktop/install/windows-install/" target="_blank">Docker Desktop for Windows</a>

2. Get FAMP via Git or download the .zip and extract it at the desired location on your system.

3. Download Gromacs (Free) and the Rosetta Software (Free for non-commercial use)

   - <a href="https://manual.gromacs.org/" target="_blank" >Gromacs</a>
   - <a href="https://downloads.rosettacommons.org/software/academic/" target="_blank" >Rosetta Commons</a>
   
3. Edit the desired Dockerfile in the FAMP Folder to match your Gromacs/Rosetta versions and the location of FAMP 

```bash
ENV GROMACS=<Path/to/gromacs.tar.gz>
ENV ROSETTA=<Path/to/rosetta.tar.bz2>
ENV PIPELINE=<Path/to/FAMP_rna-main>
```

For GPU support, make sure to select the appropriate image name on Docker Hub that matches your CUDA version, and update the first line of the Dockerfile accordingly.

- <a href="https://hub.docker.com/r/nvidia/cuda/tags" target="_blank" >Docker Hub CUDA images</a>

```bash
FROM nvidia/cuda:12.2.2-devel-ubuntu22.04
```
    
4. Build the container
   
   - Start the Docker Desktop application on your system.
   - In the left sidebar, go to **Images**, then click on `>_ Terminal` to open an integrated terminal.
   - Use the cd command to move into the folder where the Dockerfiles are located (FAMP/Docker_images). For example:
```bash
cd ~/Downloads/FAMP_rna/Docker_images
```
- Run the following command to build the image. This process may take several minutes depending on your system:
```bash
docker build -t famp_rna:cpu -f <Dockerfile_name> .
```
- After the build finishes, you will see the new image famp_rna:cpu listed under the Images section in Docker Desktop.

5. Run a container

   - Click on Containers
   - Use the docker terminal to run: 
```bash
 docker run -p 8888:8888 -v .\:/src/data/ --rm famp_rna:cpu 
```
| Option                | Description                                                                |
|-----------------------|----------------------------------------------------------------------------|
| `docker run`          | Starts a new container from a Docker image                                 |
| `-p 8888:8889`        | Maps port **8888 in the container** to **port 8889 on your host**          |
| `-v $(pwd):/src/data` | Mounts the **current working directory** into the container at `/src/data` |
| `--rm`                | Automatically removes the container after it stops                         |
| `--gpus=all`          | adds all available GPUs registered to the system to the Container          |
| `famp_rna:cpu`        | Name of the Docker image you built                                         |

```{Note}
If you have built a container with GPU support, remember to include the Docker GPU flag when running it: --gpus=all
```

6. Access the running container
If the container starts successfully, a link to a Jupyter Notebook will be displayed in the terminal output.
You can open this link in your web browser.

Alternatively, if the default settings are used, you can try opening the following link directly in your browser with the Token printed to the Terminal:
http://localhost:8889/TOKEN

---

### Installation via Docker (for experienced users)

1. Download Gromacs (Free) and the Rosetta Software (Free for non-commercial use)

   - <a href="https://manual.gromacs.org/" target="_blank" >Gromacs</a>
   - <a href="https://downloads.rosettacommons.org/software/academic/" target="_blank" >Rosetta Commons</a>
   
2. Edit the desired Dockerfile in the FAMP Folder to match your Gromacs/Rosetta versions and the location of FAMP 

```bash
ENV GROMACS=<Path/to/gromacs.tar.gz>
ENV ROSETTA=<Path/to/rosetta.tar.bz2>
ENV PIPELINE=<Path/to/FAMP_rna-main>
```

3. Build your desired container in your Terminal (Templates can be found in the repository)

```bash
docker build -t famp_rna:cpu -f <Path/to/Dockerfile_name> .
```
- After the build finishes, you will see the new image famp_rna:cpu listed under the docker images.

4. Run your container in your Terminal

```bash
 docker run -p 8888:8888 -v .\:/src/data/ --rm famp_rna:cpu 
```

---

## Installing FAMP without Docker (Linux only)

## Required Dependencies

To use the FAMP pipeline outside of Docker, the following software tools must be **pre-installed** on your system:

- RNAfold (<a href="https://www.tbi.univie.ac.at/RNA/" target="_blank">https://www.tbi.univie.ac.at/RNA/</a>)
- Rosetta (<a href="https://www.rosettacommons.org/software/license-and-download" target="_blank">https://www.rosettacommons.org/software/license-and-download</a>)
- GROMACS (<a href="https://www.gromacs.org" target="_blank">https://www.gromacs.org</a>)

These tools are used for RNA structure prediction (2D and 3D) and molecular dynamics simulations.  
The following sections provide basic installation steps for each dependency on **Linux systems**.  
For detailed information, please consult the official documentation linked above each section.

### GROMACS
```bash
sudo apt update
sudo apt install gromacs
```
Ensure GROMACS is available in your terminal:
```bash
gmx --version
```

If GROMACS was built from source (<a href="https://manual.gromacs.org/documentation/current/install-guide/index.html" target="_blank">installation guide</a>):
Add the GROMACS binary to your PATH:
```bash
echo "source /usr/local/gromacs/bin/GMXRC" >> ~/.bashrc
source ~/.bashrc
```

### RNAfold (from ViennaRNA)
```bash
sudo apt install vienna-rna
```
Test the installation:
```bash
RNAfold --version
```

### Rosetta (requires academic license)
1. Register and download from: <a href="https://rosettacommons.org/software/download/" target="_blank">https://rosettacommons.org/software/download/</a>
2. Unpack and build according to their <a href="https://docs.rosettacommons.org/docs/latest/build_documentation/Build-Documentation" target="_blank">documentation</a> 
3. Add the binary to your PATH:
```bash
echo "export PATH=/path/to/rosetta/main/source/bin:$PATH" >> ~/.bashrc
source ~/.bashrc
```


- Install FAMP_rna using pip
```bash
pip install famp
```