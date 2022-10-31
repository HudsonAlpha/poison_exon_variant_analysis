# Docker Image based on general-purpose docker image developed for CSER Consotium by RRG & SAF.
FROM us.gcr.io/broad-dsp-gcr-public/terra-jupyter-base:1.0.11
# At the time of analysis, terra-jupyter-base:1.0.0 was used as the base image.
# terra-jupyter-base:1.0.0 currently does not compile due to dependency changes, so updated 10/31/2022 make this dockerfile functional JMJL.
# Redundant package installs removed 10/31/2022 JMJL.

ENV DEBIAN_FRONTEND noninteractive

USER root

# Update packages
RUN apt-get update
# Install dependencies
RUN apt-get update && apt-get install -y \
	autoconf \
	automake \
	build-essential \
	curl \
	git \
	libbz2-dev \
	libcurl4 \
	libcurl4-openssl-dev \
	libgsl0-dev \
	liblzma-dev \
	libncurses5-dev \
	libncursesw5-dev \
	libperl-dev \
	libssl-dev \
	libtool \
	libz-dev \
	make \
	parallel \
	tabix \
	zlib1g \
	zlib1g-dev \
&& rm -rf /var/lib/apt/lists/*

WORKDIR /tmp
ARG htsversion=1.9
RUN curl -L https://github.com/samtools/htslib/releases/download/${htsversion}/htslib-${htsversion}.tar.bz2 | tar xj && \
    (cd htslib-${htsversion} && ./configure --enable-plugins --with-plugin-path='$(libexecdir)/htslib:/usr/libexec/htslib' && make install) && \
    ldconfig && \
    curl -L https://github.com/samtools/samtools/releases/download/${htsversion}/samtools-${htsversion}.tar.bz2 | tar xj && \
    (cd samtools-${htsversion} && ./configure --with-htslib=system && make install) && \
    curl -L https://github.com/samtools/bcftools/releases/download/${htsversion}/bcftools-${htsversion}.tar.bz2 | tar xj && \
    (cd bcftools-${htsversion} && ./configure --enable-libgsl --enable-perl-filters --with-htslib=system && make install) && \
    git clone --depth 1 https://github.com/samtools/htslib-plugins.git && \
    (cd htslib-plugins && make PLUGINS='hfile_cip.so hfile_mmap.so' install)


# Install vcftools
ENV ZIP=vcftools-0.1.15.tar.gz
ENV URL=https://github.com/vcftools/vcftools/releases/download/v0.1.15/
ENV FOLDER=vcftools-0.1.15
ENV DST=/tmp

RUN wget $URL/$ZIP -O $DST/$ZIP && \
  tar xvf $DST/$ZIP -C $DST && \
  rm $DST/$ZIP && \
  cd $DST/$FOLDER && \
  ./configure && \
  make && \
  make install && \
  cd / && \
  rm -rf $DST/$FOLDER

RUN wget https://github.com/vcftools/vcftools/blob/master/examples/floats.vcf && \
    vcftools --vcf floats.vcf && rm floats.vcf

###


# Install Bedtools
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools-2.30.0.tar.gz -O /tmp/bedtools.tar.gz && \
	tar zxvf /tmp/bedtools.tar.gz -C /opt/ && rm /tmp/bedtools.tar.gz && \
	cd /opt/bedtools2 && \
	make


ENV PATH="/opt/bedtools2/bin/:${PATH}"

# Update user paths
RUN export PATH="$PATH:/usr/bin/bcftools-1.9"
RUN export PATH="$PATH:/usr/bin/samtools-1.9"
RUN export PATH="$PATH:/usr/bin/htslib-1.9"
RUN export PATH="$PATH:/usr/bin/bedtools-2.30.0"


USER $USER
