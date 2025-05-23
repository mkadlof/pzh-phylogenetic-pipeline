FROM ubuntu:24.04

LABEL maintainer="Michal Kadlof <mkadlof@pzh.gov.pl>"

ENV DEBIAN_FRONTEND=noninteractive
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

RUN apt update -y && apt install --no-install-recommends -y \
    python3 \
    python3-venv \
    mafft \
    iqtree

RUN python3 -m venv /opt/venv

ENV PATH="/opt/venv/bin:$PATH:/opt/scripts:$PATH"

RUN pip install nextstrain-augur \
        Biopython

COPY bin/* /opt/scripts/
COPY data/auspice_config.json /etc/auspice/auspice_config.json