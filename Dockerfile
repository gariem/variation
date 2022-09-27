FROM ubuntu:20.04

ARG DEBIAN_FRONTEND=noninteractive

ENV TZ=Europe/London
ENV LD_LIBRARY_PATH=/usr/local/lib

ENV DISPLAY=:99
ENV DISPLAY_CONFIGURATION=1024x768x24

# Update
RUN apt-get update && apt-get install -y build-essential git zlib1g-dev wget \ 
    libbz2-dev liblzma-dev libcurl4-gnutls-dev autoconf libncurses5-dev \
    libncursesw5-dev libssl-dev libxml-xpath-perl libjson-perl bedtools \
    xvfb xorg x11-utils bc curl


# htslib, bcftools, and samtools
RUN git clone -b 1.13 https://github.com/samtools/htslib.git && \
    cd htslib && git submodule update --init --recursive && make -j `nproc` && make install && \
    cd .. && rm -rf htslib

RUN git clone -b 1.13 https://github.com/samtools/bcftools.git && \
    cd bcftools && autoheader && autoconf && ./configure && make -j `nproc` && make install && \
    cd .. && rm -rf bcftools

RUN git clone -b 1.13 https://github.com/samtools/samtools.git && \
    cd samtools && autoheader && autoconf && ./configure && make -j `nproc` && make install && \
    cd .. && rm -rf samtools

# Survivor
RUN git clone https://github.com/gariem/SURVIVOR.git && cd SURVIVOR/Debug && make -j `nproc` && \
    ln /SURVIVOR/Debug/SURVIVOR /sbin/SURVIVOR

# Minigraph
RUN git clone https://github.com/lh3/minigraph.git && cd minigraph && make -j `nproc` && \
    ln /minigraph/minigraph /sbin/minigraph

# Minimap2
RUN curl -L https://github.com/lh3/minimap2/releases/download/v2.24/minimap2-2.24_x64-linux.tar.bz2 | tar -jxvf - && \
    ln ./minimap2-2.24_x64-linux/minimap2 /sbin/minimap2

# User mouse
RUN useradd -ms /bin/bash mouse
USER mouse
WORKDIR /home/mouse

# IGV
RUN wget https://data.broadinstitute.org/igv/projects/downloads/2.12/IGV_Linux_2.12.2_WithJava.zip -O IGV_Linux_2.12.2_WithJava.zip && \
    unzip IGV_Linux_2.12.2_WithJava.zip

# Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py37_4.11.0-Linux-x86_64.sh -O /home/mouse/miniconda.sh
RUN yes | bash /home/mouse/miniconda.sh -b -p /home/mouse/miniconda
ENV PATH="$PATH:/home/mouse/miniconda/condabin:/home/mouse/miniconda/bin"

RUN conda install -y -c bioconda scipy matplotlib pandas numpy pbsv
RUN conda install -y -c bioconda bedops
RUN conda install -y -c hcc smrtlink-tools

# PATH
ENV PATH="$PATH:/home/mouse/miniconda/condabin:/home/mouse/miniconda/bin:/home/mouse/IGV_Linux_2.12.2/jdk-11/bin"

CMD ["bash"]