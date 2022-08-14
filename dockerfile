FROM python:3.9-slim-buster

EXPOSE 80

COPY ./app /workspace/app 
#COPY other files/folders # RM <
WORKDIR /workspace
RUN mkdir logs etc

RUN apt install hmmer
RUN apt install mcl
# RM < Automate BLAST install
# RM < automate hhsuite https://github.com/soedinglab/hh-suite/wiki#installation-of-the-hhsuite-and-its-databases
# RM < Automate muscle install https://drive5.com/muscle5/manual/install.html
RUN apt-get update
RUN apt-get install --upgrade pip
RUN pip install -r requirements.txt --prefer-binary

CMD ["uvicorn", "app.api:app", "--host", "0.0.0.0", "--port", "80"]