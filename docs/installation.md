# Omni2Tree Installation

## Dependencies

This software relies on four external tools: [OMA Standalone](https://omabrowser.org/standalone/), [Rasusa](https://github.com/mbhall88/rasusa?tab=readme-ov-file#install), [czid-dedup](https://github.com/chanzuckerberg/czid-dedup?tab=readme-ov-file#installation), and [Read2Tree](https://github.com/DessimozLab/read2tree/tree/combined?tab=readme-ov-file#installation). It assumes all programs are in your Conda environment or `PATH`.

Below are four ways to install all required dependencies.

## 1. Docker

The repository ships a `Dockerfile` that installs every dependency and the Omni2Tree entry points. No Conda environment or manual tool installation required.

Clone the repository and build the image:

```bash
git clone https://github.com/DanielPAagustinho/omni2tree.git
cd omni2tree
docker build -t omni2tree .
```

Verify the image:

```bash
docker run --rm omni2tree bash -c 'command -v o2t-step1 && command -v read2tree && command -v oma'
```

To run the hRSV demo inside the container, see [Running with Docker](../demo/README.md#running-with-docker).

Two important notes for running the container:
- Use `bash -c '...'`, **not** `bash -lc '...'`. The login-shell flag is not needed; the image sets `PATH` via `ENV` and it is available in every non-login shell.
- The container runs as **root** by default. Do **not** pass `--user`; some tools (OMA, warthog) need write access to directories under `/opt/OMA` on first run.

> **Note:** The default `Dockerfile` does not include `czid-dedup` (read deduplication). To build an image with deduplication support, use `Dockerfile-complete` instead: `docker build -f Dockerfile-complete -t omni2tree .`

## 2. Apptainer / Singularity

For HPC users without Docker or root access, use the prebuilt `omni2tree.sif` image from the GitHub release assets.

Clone the repository, then download the image and checksum:

```bash
wget https://github.com/DanielPAagustinho/omni2tree/releases/download/v0.3.1/omni2tree.sif
wget https://github.com/DanielPAagustinho/omni2tree/releases/download/v0.3.1/omni2tree.sif.sha256
sha256sum -c omni2tree.sif.sha256
```

Verify the image:

```bash
apptainer exec omni2tree.sif bash -c 'command -v o2t-step1 && command -v read2tree && command -v oma'
```

If your system provides `singularity` instead of `apptainer`, replace `apptainer exec` with `singularity exec`.

To run the hRSV demo with Apptainer/Singularity, see [Running with Apptainer / Singularity](../demo/README.md#running-with-apptainer--singularity).

Two important notes for running the image:
- The image runs as the current user, not as root.
- The SIF filesystem is read-only at runtime. Bind writable host directories for outputs and downloaded reads.

## 3. Installation with Conda

[Conda](https://docs.anaconda.com/miniconda/) is a package manager that allows you to install all dependencies quickly and easily.

```bash
conda create -n my_env python=3.10.8 -y
conda activate my_env
conda install -c bioconda rasusa sra-tools entrez-direct -y
```

Notes:

- OMA standalone and czid-dedup are not available via Conda. Please follow the "Installation from source" instructions below.
- Read2Tree must be installed from source (see below). The conda package does not include `--meta` (metagenomic mode), which is required for co-infection analysis.

## 4. Installation from Source

If you prefer to install the tools manually from their source code, use the following commands.

### OMA Standalone

```bash
# Download the latest version; this example uses 2.6.0.
wget -O oma.tgz https://omabrowser.org/standalone/OMA.2.6.0.tgz
tar xvzf oma.tgz
cd OMA.2.6.0

# Choose your install path. If omitted, OMA installs to /usr/local/OMA and may require sudo.
./install.sh /your/install/path

# Make sure the OMA bin folder is in PATH.
echo 'export PATH=$PATH:/your/install/path/OMA/bin' >> ~/.bashrc
source ~/.bashrc
```

### Rasusa

```bash
# When Rasusa is downloaded it is automatically added to your PATH.
curl -sSL rasusa.mbh.sh | sh
```

### czid-dedup

`czid-dedup` requires [rust/cargo](https://www.rust-lang.org/tools/install) for compilation.

```bash
git clone https://github.com/chanzuckerberg/czid-dedup.git
cd czid-dedup
cargo build --release

# Make sure the release directory is in your PATH.
echo 'export PATH=$PATH:your/install/path/czid-dedup/target/release' >> ~/.bashrc
source ~/.bashrc
```

### Read2Tree

```bash
# Create conda environment.
conda create -n r2t python=3.10.8 -y
conda activate r2t

# Get required Python packages.
conda install -c conda-forge biopython numpy Cython ete3 lxml tqdm scipy pyparsing requests natsort pyyaml filelock -y
conda install -c bioconda dendropy pysam -y

# Install required software.
conda install -c bioconda mafft iqtree minimap2 samtools -y

# Clone and install Read2Tree from the combined branch.
git clone -b combined https://github.com/DessimozLab/read2tree.git
cd read2tree
python setup.py install
```

Read2Tree will be placed in the default bin folder of your Conda installation.

The `combined` branch includes both standard and metagenomic (`--meta`) modes.

### SRA Toolkit

```bash
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
tar -xvzf sratoolkit.current-ubuntu64.tar.gz

# Add executable to your path, using your own version. This example uses 3.2.0.
echo 'export PATH="$PATH:/your/install/path/sratoolkit.3.2.0-ubuntu64/bin"' >> ~/.bashrc
source ~/.bashrc
```

### Entrez Direct

```bash
# Get the scripts and download them in an "edirect" folder in the user's home directory.
sh -c "$(wget -q https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh -O -)"
source ~/.bashrc
```

## 5. Verify Dependencies

To verify that all tools are installed and available in your Conda environment or `PATH`, run:

```bash
command -v oma && command -v rasusa && command -v czid-dedup && command -v read2tree && command -v fasterq-dump && command -v esearch
```

## 6. Install Omni2Tree

Clone the repo and run the installer:

```bash
git clone https://github.com/DanielPAagustinho/omni2tree.git
cd omni2tree
./install_omni2tree.sh /your/install/path
```

The installation script creates symlinks to the shell entry points.

- If you omit the install path, symlinks go to `/usr/local/bin` and may require sudo.
- If you use a custom path, make sure it is in your `PATH`.

Check the installation with:

```bash
which o2t-step1 && o2t-step1 --help
which o2t-step2 && o2t-step2 --help
which o2t-step3 && o2t-step3 --help
which o2t-sra   && o2t-sra   --help
```
