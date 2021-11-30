FROM python:3.8-buster

COPY requirements.txt .
RUN pip install -r requirements.txt

RUN mkdir /root/tools
COPY antenna_count_reads.py /root/tools/ 
COPY antenna_tag_reads.py /root/tools/ 
