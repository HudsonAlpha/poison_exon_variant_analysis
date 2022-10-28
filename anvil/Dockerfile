##### Docker Image
FROM us.gcr.io/broad-dsp-gcr-public/terra-jupyter-base:1.0.11

################################
# Custom packages RRG 8/1/21
# SFelker Edits 11/17/21: uncomment out bedtools load and paths
################################

#docker build script

ENV DEBIAN_FRONTEND noninteractive

USER root

RUN apt-get update
RUN apt-get install -y automake
RUN apt-get install -y autoconf
RUN apt-get install -y libtool

# Update packages

RUN apt-get install make && apt-get install libbz2-dev && apt-get install zlib1g-dev && apt-get install libncurses5-dev && apt-get install libncursesw5-dev && apt-get install liblzma-dev && apt-get install -y libcurl4 && apt install zlib1g && apt install zlib1g-dev && apt-get install libz-dev


# Install vcftools
RUN apt-get update && apt-get install -y \
	build-essential \
	curl \
	git \
	libbz2-dev \
	libcurl4-openssl-dev \
	libgsl0-dev \
	liblzma-dev \
	libncurses5-dev \
	libperl-dev \
	libssl-dev \
	parallel \
	tabix \
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
    git clone --depth 1 git://github.com/samtools/htslib-plugins && \
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

#Update user paths
RUN export PATH="$PATH:/usr/bin/bcftools-1.9"
RUN export PATH="$PATH:/usr/bin/samtools-1.9"
RUN export PATH="$PATH:/usr/bin/htslib-1.9"
RUN export PATH="$PATH:/usr/bin/bedtools-2.30.0"

#RUN source ~/.profile

USER $USER
