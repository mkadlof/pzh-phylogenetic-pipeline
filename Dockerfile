FROM ubuntu:24.04

LABEL maintainer="Michal Kadlof <mkadlof@pzh.gov.pl>"

ENV DEBIAN_FRONTEND=noninteractive
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

RUN apt update -y && apt install --no-install-recommends -y \
    python3 \
    python3-venv

RUN python3 -m venv /opt/venv

ENV PATH="/opt/venv/bin:$PATH"

RUN pip install nextstrain-augur

RUN apt install --no-install-recommends -y \
    mafft