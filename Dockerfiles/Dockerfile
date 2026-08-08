FROM ubuntu:24.04

LABEL maintainer="Michal Kadlof, Michal Lazniewski <m.lazniewski@pzh.gov.pl>"

ENV DEBIAN_FRONTEND=noninteractive
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

# ---------------------------
# System packages + python 3.11
# ---------------------------
RUN apt update && apt install -y software-properties-common && \
    add-apt-repository ppa:deadsnakes/ppa && \
    apt update && apt upgrade -y && \
    apt install -y --no-install-recommends \
    build-essential pkgconf \
    python3.11 python3.11-distutils \
    curl unzip git gcc wget \
    zlib1g-dev libbz2-dev libcurl4-gnutls-dev libssl-dev liblzma-dev \
    python3-pycurl mash seqtk roary bison flex default-jre cmake npm seqkit vim jq 

# mafft (fixed veriosn)
RUN wget https://mafft.cbrc.jp/alignment/software/mafft_7.520-1_amd64.deb && \
    dpkg -i mafft_7.520-1_amd64.deb && \
    rm -f /mafft_7.520-1_amd64.deb


# ---------------------------
# Update system Python
# ---------------------------
RUN update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.11 1 && \
    curl -sS https://bootstrap.pypa.io/get-pip.py | python3.11 && \
    ln -s /usr/local/bin/pip /usr/bin/pip3

# ---------------------------
# Install all Python packages globally
# ---------------------------
RUN pip install ete3==3.1.3 \
                nextstrain-augur==30.0.0 \
                click==8.1.8 \
                requests==2.32.3 \
                Biopython==1.85 \
                networkx==3.6.1 \
                pandas==2.3.0 \
                plotly==6.3.0

WORKDIR /usr/bin
RUN ln -s python3 python

# ---------------------------
# Bacterial: RAxML-NG pinned specific version
# ---------------------------
WORKDIR /opt/docker
WORKDIR /opt/docker

RUN git clone https://github.com/amkozlov/raxml-ng && \
    cd raxml-ng && \
    git checkout 805318cef87bd5d67064efa299b5d1cf948367fd && \
    git submodule update --init --recursive && \
    mkdir build && cd build && \
    cmake .. && \
    make -j 4 && \
    make install

# ---------------------------
# Bacterial: IQtree2 v. 2.4.0
# ---------------------------
RUN curl -L -o iqtree-2.4.0-Linux-intel.tar.gz https://github.com/iqtree/iqtree2/releases/download/v2.4.0/iqtree-2.4.0-Linux-intel.tar.gz && \
    tar -zxf iqtree-2.4.0-Linux-intel.tar.gz && \
    cp iqtree-2.4.0-Linux-intel/bin/* /usr/bin && \
    cd /opt/docker

# ---------------------------
# Bacterial: Auspice pinned specific version
# ---------------------------
RUN git clone https://github.com/nextstrain/auspice.git && \
    cd auspice && \
    git checkout a9b74fe && \
    npm ci && \
    # npm audit fix --force && \
    npm install --global . && \
    cd /opt/docker

ENV HOST=0.0.0.0

# ---------------------------
# Copy scripts and configs
# Bacterial and viral scripts paths are not identical for now so the same script must be copied sometimes to two different locations
# ---------------------------
RUN mkdir -p /opt/docker/custom_scripts /opt/docker/config /opt/scripts


# Bacterial scripts/config
COPY bin/* /opt/docker/custom_scripts/
COPY config/* /opt/docker/config/

# Viral scripts/config
COPY bin/* /opt/scripts/
COPY data/auspice_config.json /etc/auspice/auspice_config.json
COPY data/microreact_config_bacteria.microreact /opt/docker/config/microreact_config_bacteria.microreact
COPY VERSION /opt/docker/

CMD ["/bin/bash"]

