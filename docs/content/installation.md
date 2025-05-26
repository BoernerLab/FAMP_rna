# Installation Guide

This page provides instructions to install and run the FAMP pipeline. You can either build the required Docker containers yourself using the provided Dockerfiles, or install FAMP directly using `pip` in a local Python environment (Linux only).

---

## Installation via Docker (recommended)

FAMP includes Dockerfiles to build containerized environments for modeling, simulation, and data analysis. These Dockerfiles allow reproducible and portable setups across systems, including GPU support for MD simulations.

```{Note}
We do not publish pre-built Docker images. You must build all containers locally from the provided Dockerfiles.
```

### Prerequisites
- Docker Desktop (recommended for inexperienced users)
- Docker Engine (alternative for experienced users or server setups)

### Step-by-step instructions

1. Install Docker Desktop (if not yet installed)
   - <a href="https://docs.docker.com/desktop/install/linux/" target="_blank" >Docker Desktop for Linux</a>
   - <a href="https://docs.docker.com/desktop/install/mac-install/" target="_blank">Docker Desktop for macOS</a>
   - <a href="https://docs.docker.com/desktop/install/windows-install/" target="_blank">Docker Desktop for Windows</a>
---

2. Download FAMP and the Docker files

3. Build the container
   - Start the Docker Desktop application on your system.
   - In the left sidebar, go to **Images**, then click on `>_ Terminal` to open an integrated terminal.
   - Use the cd command to move into the folder where the Dockerfile is located. For example:
```bash
cd ~/Downloads/FAMP_rna/docker_files
```
- Run the following command to build the image. This process may take several minutes, especially on first run:
```bash
docker build -t /famp_rna:gpu famp_placeholder
```
- After the build finishes, you will see the new image famp_rna:gpu listed under the Images section in Docker Desktop.


3. Run a container

   - Click on Containers
   - Use the docker terminal to run: 
```bash
 docker run -p 8888:8888 -v .\:/src/data/ --rm felix/pipeline 
```
| Option                | Description                                                                |
|-----------------------|----------------------------------------------------------------------------|
| `docker run`          | Starts a new container from a Docker image                                 |
| `-p 8888:8889`        | Maps port **8888 in the container** to **port 8889 on your host**          |
| `-v $(pwd):/src/data` | Mounts the **current working directory** into the container at `/src/data` |
| `--rm`                | Automatically removes the container after it stops                         |
| `famp_rna:gpu`        | Name of the Docker image you built                                         |


4. Access the running container
If the container starts successfully, a link to a Jupyter Notebook will be displayed in the terminal output.
You can open this link in your web browser.

Alternatively, if the default settings are used, you can try opening the following link directly in your browser:
http://localhost:8889/TOKEN




## Installing FAMP without Docker (Linux only)

Bla Bla Bla Dependencies bla bla bal 

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
1. Register and download from: [https://rosettacommons.org/software](https://rosettacommons.org/software)
2. Unpack and build according to their <a>documentation</a> 
3. Add the binary to your PATH:
```bash
echo "export PATH=/path/to/rosetta/main/source/bin:$PATH" >> ~/.bashrc
source ~/.bashrc
```


- Install FAMP_rna using pip
```bash
pip install famp
```