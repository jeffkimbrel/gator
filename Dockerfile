# For more information, please refer to https://aka.ms/vscode-docker-python
FROM ubuntu:22.04

# Keeps Python from generating .pyc files in the container
ENV PYTHONDONTWRITEBYTECODE=1

# Turns off buffering for easier container logging
ENV PYTHONUNBUFFERED=1

# Install pip requirements (could come from pip freeze > requirements.txt)
COPY requirements.txt .

RUN apt-get update
RUN apt-get -y install python3.10 python3-pip python3-setuptools python3-dev
RUN apt-get -y install ncbi-blast+
RUN apt-get -y install git

RUN pip3 install pandas
RUN pip3 install hmmer
RUN pip3 install tqdm
RUN pip3 install git+https://github.com/jeffkimbrel/jakomics.git

# 
WORKDIR /app
COPY gator.py /app/
COPY metadata.py /app/
COPY pathway.py /app/
COPY gator_db.xlsx /app/

# Creates a non-root user with an explicit UID and adds permission to access the /app folder
# For more info, please refer to https://aka.ms/vscode-docker-python-configure-containers
RUN adduser -u 5678 --disabled-password --gecos "" appuser && chown -R appuser /app
USER appuser

# During debugging, this entry point will be overridden. For more information, please refer to https://aka.ms/vscode-docker-python-debug
ENTRYPOINT ["python3", "gator.py"]
