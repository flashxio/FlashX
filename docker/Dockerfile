# FlashX Docker image with ssh port forwarding and general ubuntu hackery


FROM ubuntu:16.04
MAINTAINER Alexander Niculescu <al3xander.niculescu@gmail.com>, Da Zheng <zhengda1936@gmail.com>

RUN apt-get update && apt-get install -y openssh-server
RUN mkdir /var/run/sshd
RUN echo 'root:screencast' | chpasswd
RUN sed -i 's/PermitRootLogin without-password/PermitRootLogin yes/' /etc/ssh/sshd_config

# SSH login fix. Otherwise user is kicked off after login
RUN sed 's@session\s*required\s*pam_loginuid.so@session optional pam_loginuid.so@g' -i /etc/pam.d/sshd


###FLASHX CONF COMBINED FROM DOCKERFILE &&FLASHX QUICKSTART###
#https://github.com/icoming/FlashX/wiki/FlashX-Quick-Start-Guide
#https://github.com/wking/dockerfile

RUN apt-get update
RUN apt-get install -y git cmake g++
RUN apt-get install -y libboost-dev libboost-system-dev libboost-filesystem-dev libnuma-dev libaio-dev libhwloc-dev libatlas-base-dev zlib1g-dev
#wget is for trilinos
RUN apt-get install wget

RUN git clone https://github.com/flashxio/FlashX.git

WORKDIR /FlashX
RUN mkdir build
WORKDIR build
RUN cmake ..
RUN make -j4


####Install and compile R
#https://www.digitalocean.com/community/tutorials/how-to-set-up-r-on-ubuntu-14-04
RUN apt-get -y install r-base

#run R >> intstall igraph install.packages("igraph")
RUN R -e "install.packages('Rcpp', repos = 'http://cran.rstudio.com/')"
RUN R -e "install.packages('RSpectra', repos = 'http://cran.rstudio.com/')"
RUN R -e "install.packages('igraph', repos = 'http://cran.rstudio.com/')"

WORKDIR /FlashX
RUN ./install_FlashR.sh
RUN ./install_FlashGraphR.sh

#check to see if it's there ^^^?

####R finished####


###FLASHX CONF END ###


ENV NOTVISIBLE "in users profile"
RUN echo "export VISIBLE=now" >> /etc/profile

EXPOSE 22
CMD ["/usr/sbin/sshd", "-D"]
