FROM ubuntu:focal
MAINTAINER Cristina Catalan <cataln.cris@gmail.com>

# Steckmap initial data
ENV HOME /home/steckmap
COPY sources /home/steckmap/Yorick/
WORKDIR /home/steckmap/Yorick/

# Install main steps
RUN apt update -y

# Configure TZDATA
ENV TZ=Europe/Madrid
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# Misc dependencies
RUN apt install -y wget python3
RUN apt-get install -y libtool build-essential autoconf automake

# Yorick install
RUN apt-get install -y yorick

# For the X server
RUN apt-get  install -y xvfb x11vnc x11-xkb-utils xfonts-100dpi xfonts-75dpi xfonts-scalable xfonts-cyrillic x11-apps
ADD xvfb_init /etc/init.d/xvfb
RUN chmod a+x /etc/init.d/xvfb
ADD xvfb_daemon_run /usr/bin/xvfb-daemon-run
RUN chmod a+x /usr/bin/xvfb-daemon-run
ENV DISPLAY :99

ADD start_steckmap /usr/bin/start_steckmap
RUN chmod a+x /usr/bin/start_steckmap

ENTRYPOINT ["start_steckmap"]
