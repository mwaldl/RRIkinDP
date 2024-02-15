FROM continuumio/miniconda3

RUN apt update
RUN apt install --yes build-essential cmake

RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge
RUN conda config --set channel_priority strict

RUN conda install mamba -n base -c conda-forge
RUN mamba install viennarna intarna

RUN cd /
RUN git clone https://github.com/amrayn/easyloggingpp.git /easyloggingpp
RUN mkdir /easyloggingpp/build
RUN cd /easyloggingpp/build; cmake -Dtest=OFF ../
RUN cd /easyloggingpp/build; make
RUN cd /easyloggingpp/build; make install

ENV LD_LIBRARY_PATH="/opt/conda/lib"
COPY src/RRIkinDP/RRIkinDP.cpp RRIkinDP.cpp
RUN g++ -std=c++1y  RRIkinDP.cpp -o RRIkinDP  -I"`conda info --base`/include" -L"`conda info --base`/lib" -lboost_regex -lboost_program_options -lboost_filesystem -lboost_system -lIntaRNA -fopenmp -lboost_regex -lRNA -leasylogging -fpermissive

