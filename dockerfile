FROM python:3.9-slim-buster

EXPOSE 80

COPY ./app /workspace/app
COPY ./data /workspace/data
COPY ./ /workspace

WORKDIR /workspace
RUN apt-get -y update
RUN apt install -y g++ make cmake git wget
RUN bash install-reqs.sh

RUN python3 -m pip install --upgrade pip
RUN pip install -r requirements.txt --prefer-binary

ENV HHLIB="~/programs/hh-suite"
ENV PATH="$PATH:$HHLIB/bin:$HHLIB/scripts"

CMD ["uvicorn", "app.api:app", "--host", "0.0.0.0", "--port", "80"]
