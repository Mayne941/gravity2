FROM python:3.9-slim-buster

EXPOSE 80

COPY ./app /workspace/app 
COPY ./data /workspace/data
COPY ./ /workspace
WORKDIR /workspace

# RUN apt install hmmer
# RUN apt install mcl
# RM < Automate BLAST install
# RM < automate hhsuite https://github.com/soedinglab/hh-suite/wiki#installation-of-the-hhsuite-and-its-databases
# RM < Automate muscle install https://drive5.com/muscle5/manual/install.html
RUN apt update
RUN python3 -m pip install --upgrade pip
RUN pip install -r requirements.txt --prefer-binary

CMD ["uvicorn", "app.api:app", "--host", "0.0.0.0", "--port", "80"]