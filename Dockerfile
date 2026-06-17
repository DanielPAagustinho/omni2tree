FROM mambaorg/micromamba:1.5.10

USER root

RUN apt-get update && apt-get install -y --no-install-recommends \
    bash \
    bc \
    build-essential \
    ca-certificates \
    curl \
    git \
    make \
    util-linux \
    wget \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/omni2tree

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml

RUN micromamba install -y -n base -f /tmp/environment.yml \
    && micromamba clean -a -y

SHELL ["micromamba", "run", "-n", "base", "/bin/bash", "-c"]

RUN python -m pip install --no-build-isolation "git+https://github.com/DessimozLab/read2tree.git@combined" \
    && python -m pip install "setuptools<81"  # read2tree imports pkg_resources, removed in setuptools>=81

ENV OMA_VERSION="2.7.0"

RUN wget -O /tmp/oma.tgz "https://omabrowser.org/standalone/OMA.${OMA_VERSION}.tgz" \
    && mkdir -p /tmp/oma-src \
    && tar -xzf /tmp/oma.tgz -C /tmp/oma-src --strip-components=1 \
    && /tmp/oma-src/install.sh /opt \
    && cp -a /tmp/oma-src/bin/. "/opt/OMA/OMA.${OMA_VERSION}/bin/" \
    && chmod -R a+rwX /opt/OMA \
    && rm -rf /tmp/oma.tgz /tmp/oma-src

# Build OMA's warthog (GETHOGs bottom-up) venv at build time and fix the
# non-absolute "#!python" shebang that setup.py leaves behind. Without this,
# warthog rebuilds the venv on first use and fails: the micromamba base image
# has no working python venv/virtualenv toolchain on PATH, and the kernel cannot
# resolve a relative shebang interpreter ("required file not found").
RUN OMA_HOME="/opt/OMA/OMA.${OMA_VERSION}" \
    && pip install --no-cache-dir virtualenv \
    && python -m virtualenv "$OMA_HOME/.venv" \
    && "$OMA_HOME/.venv/bin/pip" install --no-cache-dir -r "$OMA_HOME/hog_bottom_up/requirements.txt" \
    && "$OMA_HOME/.venv/bin/pip" install --no-cache-dir "$OMA_HOME/hog_bottom_up" \
    && sed -i "1s|^#!.*|#!$OMA_HOME/.venv/bin/python3|" "$OMA_HOME/.venv/bin/warthogs.py" \
    && "$OMA_HOME/.venv/bin/warthogs.py" --version

ENV PATH="/opt/OMA/OMA.${OMA_VERSION}/bin:/opt/OMA/bin:/opt/conda/bin:${PATH}"

COPY --chown=$MAMBA_USER:$MAMBA_USER . /opt/omni2tree

RUN chmod +x install_omni2tree.sh scripts/*.sh scripts/*.py scripts/*.R utils/*.py view/omni2treeview.py \
    && ./install_omni2tree.sh /usr/local/bin \
    && command -v o2t-step1 \
    && command -v read2tree \
    && command -v oma \
    && command -v oma-status

WORKDIR /work

CMD ["bash"]
