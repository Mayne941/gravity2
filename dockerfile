FROM continuumio/miniconda3

ENV CACHE_TTL=15
ENV CACHE_MAX_SIZE=1024

EXPOSE 80

COPY ./app /workspace/app
COPY ./cli /workspace/cli
COPY ./install_deps_local.sh /workspace/
COPY ./requirements.txt /workspace/

WORKDIR /workspace
RUN mkdir /logs

RUN conda create -y -n gravity python=3.10
RUN echo "source activate gravity" > ~/.bashrc
ENV PATH /opt/conda/envs/gravity/bin:$PATH

RUN pip install --upgrade pip
RUN pip install -r requirements.txt --prefer-binary
RUN bash install_deps_local.sh
