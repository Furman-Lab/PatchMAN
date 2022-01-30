FROM python:3.5

ENV PYTHONUNBUFFERED definitely
RUN mkdir -p /patchman/app


###### Rosetta ######
ENV ROSETTA_BIN=/rosetta_src_2019.14.60699_bundle/main/source/bin
ENV ROSETTA_DB=/rosetta_src_2019.14.60699_bundle/main/database

# Get Rosetta credentials from build-arg
ARG ROSETTA_USER
ARG ROSETTA_PASS

# Download and Compile Rosetta
RUN mkdir -p /rosetta
WORKDIR /rosetta
RUN curl -f -o rosetta.tar.gz  -u ${ROSETTA_USER}:${ROSETTA_PASS} https://www.rosettacommons.org/downloads/academic/2019/wk14/rosetta_src_2019.14.60699_bundle.tgz --keepalive-time 2 -H 'Expect:'
RUN tar --extract --file rosetta.tar.gz
RUN rm rosetta.tar.gz

WORKDIR /rosetta/rosetta_src_2019.14.60699_bundle/main/source


###### MASTER ######
# Download MASTER
WORKDIR /patchman
RUN curl -f -o  master-v1.6.tar.gz https://grigoryanlab.org/master/master-v1.6.tar.gz
RUN tar xvzf master-v1.6.tar.gz
RUN rm master-v1.6.tar.gz
RUN cp master-v1.6/src/Match.cpp master-v1.6/src/Match_old.cpp

# patch MASTER
RUN head -106 master-v1.6/src/Match_old.cpp > master-v1.6/src/Match.cpp
RUN echo "double **R=((Match*)(&m))->getRotation();"  >> master-v1.6/src/Match.cpp
RUN echo "double *T=((Match*)(&m))->getTranslation();"  >> master-v1.6/src/Match.cpp
RUN echo ""  >> master-v1.6/src/Match.cpp
RUN echo "os << \" T: \" << T[0]    << \" \" << T[1]    << \" \" << T[2]    << \" \";"  >> master-v1.6/src/Match.cpp
RUN echo "os << \" U: \" << R[0][0] << \" \" << R[0][1] << \" \" << R[0][2] << \" \""   >> master-v1.6/src/Match.cpp 
RUN echo "             << R[1][0] << \" \" << R[1][1] << \" \" << R[1][2] << \" \""     >> master-v1.6/src/Match.cpp
RUN echo "             << R[2][0] << \" \" << R[2][1] << \" \" << R[2][2] << \" ===\" ;" >> master-v1.6/src/Match.cpp
RUN tail -n +107 master-v1.6/src/Match_old.cpp >> master-v1.6/src/Match.cpp

# Get mslib for MASTER
WORKDIR /patchman/master-v1.6/
RUN curl -f -o msl-static-Linux-x86-64_1.2.2.7.tar.gz https://grigoryanlab.org/msl/msl-static-Linux-x86-64_1.2.2.7.tar.gz
RUN tar --extract --file msl-static-Linux-x86-64_1.2.2.7.tar.gz

# Compile MASTER
RUN make all
# STOPSIGNAL 0


##### PyRosetta #####

# Get PyRosetta credentials from build-arg
ARG PYROSETTA_USER
ARG PYROSETTA_PASS

# Download and install PyRosetta
WORKDIR /patchman
ARG PYROSETTA_LINK=https://graylab.jhu.edu/download/PyRosetta4/archive/release/PyRosetta4.Release.python35.linux.wheel/pyrosetta-2020.25+release.d2d9f90-cp35-cp35m-linux_x86_64.whl
RUN curl -f -o pyrosetta-2020.25+release.d2d9f90-cp35-cp35m-linux_x86_64.whl -u ${PYROSETTA_USER}:${PYROSETTA_PASS} $PYROSETTA_LINK --keepalive-time 2 -H 'Expect:'
RUN pip install pyrosetta-2020.25+release.d2d9f90-cp35-cp35m-linux_x86_64.whl
RUN rm pyrosetta-2020.25+release.d2d9f90-cp35-cp35m-linux_x86_64.whl

COPY *py /patchman/app/
RUN mv /patchman/master-v1.6/ /master-v1.6


##### OpenMPI #####
ENV OMPI_DIR=/opt/ompi
ENV OMPI_VERSION=2.1.6
ENV OMPI_URL="https://download.open-mpi.org/release/open-mpi/v2.1/openmpi-$OMPI_VERSION.tar.bz2"
RUN mkdir -p /tmp/ompi
RUN mkdir -p /opt

# Download OpenMPI
WORKDIR /tmp/ompi 
RUN wget -O openmpi-$OMPI_VERSION.tar.bz2 $OMPI_URL 
RUN tar -xjf openmpi-$OMPI_VERSION.tar.bz2

# Compile and install
WORKDIR /tmp/ompi/openmpi-$OMPI_VERSION 
RUN ./configure --prefix=$OMPI_DIR 
RUN make install

# Set env variables so we can compile our application
ENV PATH=$OMPI_DIR/bin:$PATH
ENV LD_LIBRARY_PATH=$OMPI_DIR/lib:$LD_LIBRARY_PATH
ENV MANPATH=$OMPI_DIR/share/man:$MANPATH

# Compile Rosetta
WORKDIR /rosetta/rosetta_src_2019.14.60699_bundle/main/source
RUN ./scons.py -j 4 extras=mpi bin

WORKDIR /patchman
ENTRYPOINT []
