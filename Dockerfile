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

RUN pip3 install pandas
RUN pip3 install hmmer
RUN pip3 install tqdm

RUN apt-get -y install git
RUN pip3 install git+https://github.com/jeffkimbrel/jakomics.git

# 
WORKDIR /app
COPY gator.py /app/
COPY metadata.py /app/
COPY pathway.py /app/
COPY gator_db.xlsx /app/
COPY gator.faa /app/

# kofamscan
RUN apt-get -y install parallel 
RUN apt-get -y install wget 
RUN apt-get -y install ruby
RUN wget https://www.genome.jp/ftp/tools/kofam_scan/kofam_scan-1.3.0.tar.gz && \
      tar xvzf kofam_scan-1.3.0.tar.gz && \
      rm -f kofam_scan-1.3.0.tar.gz && \
      cp -r kofam_scan-1.3.0/* /app && \
      rm -rf kofam_scan-1.3.0

# Creates a non-root user with an explicit UID and adds permission to access the /app folder
# For more info, please refer to https://aka.ms/vscode-docker-python-configure-containers
RUN adduser -u 5678 --disabled-password --gecos "" appuser && chown -R appuser /app
USER appuser

ENV PATH="$PATH:/app/"

# During debugging, this entry point will be overridden. For more information, please refer to https://aka.ms/vscode-docker-python-debug
ENTRYPOINT ["python3", "gator.py", "--docker"]
