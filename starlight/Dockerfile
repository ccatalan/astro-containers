FROM ubuntu:focal
MAINTAINER Cristina Catalan <cataln.cris@gmail.com>

WORKDIR /home/starlight/STARLIGHTv04

COPY STARLIGHTv04/* /home/starlight/STARLIGHTv04/

RUN ls
RUN pwd

RUN echo "deb [trusted=yes] http://old-releases.ubuntu.com/ubuntu/ hardy universe" >> /etc/apt/sources.list
RUN echo "deb-src [trusted=yes] http://old-releases.ubuntu.com/ubuntu/ hardy universe" >> /etc/apt/sources.list
RUN echo "deb [trusted=yes] http://old-releases.ubuntu.com/ubuntu/ hardy-updates universe" >> /etc/apt/sources.list
RUN echo "deb-src [trusted=yes] http://old-releases.ubuntu.com/ubuntu/ hardy-updates universe" >> /etc/apt/sources.list

RUN apt update -y
RUN apt install -y g77
RUN chmod +x /home/starlight/STARLIGHTv04/StarlightChains_v04.amd64_g77-3.4.6-r1_static.exe
ENTRYPOINT ["/home/starlight/STARLIGHTv04/StarlightChains_v04.amd64_g77-3.4.6-r1_static.exe"]
