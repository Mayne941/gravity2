FROM python:3.9-slim-buster

EXPOSE 80

COPY ./app /workspace/app 
#COPY other files/folders
WORKDIR /workspace
RUN mkdir logs etc

#RUN sh file with cmd line tools 
RUN apt-get update
RUN apt-get install --upgrade pip
RUN pip install -r requirements.txt --prefer-binary

CMD ["uvicorn", "app.api:app", "--host", "0.0.0.0", "--port", "80"]