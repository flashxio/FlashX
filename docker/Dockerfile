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

RUN echo "deb http://cran.rstudio.com/bin/linux/ubuntu xenial/" >> /etc/apt/sources.list
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
RUN apt-get update
RUN apt-get install -y git cmake g++
RUN apt-get install -y libboost-dev libboost-system-dev libboost-filesystem-dev libaio-dev libatlas-base-dev zlib1g-dev
RUN apt-get install -y numactl xfsprogs
RUN if [ `numactl -H | grep "nodes" | awk '{print $2}'` -gt 1 ]; then apt-get install -y libnuma-dev libhwloc-dev; fi
#wget is for trilinos
RUN apt-get install -y wget vim

RUN git clone https://github.com/flashxio/FlashX.git
RUN git clone https://github.com/flashxio/FlashR.git
RUN git clone https://github.com/flashxio/FlashGraphR.git

WORKDIR /FlashX
RUN mkdir build
WORKDIR build
RUN cmake ..
RUN make


####Install and compile R
#https://www.digitalocean.com/community/tutorials/how-to-set-up-r-on-ubuntu-14-04
RUN apt-get update
RUN apt-get install -y r-base

# Install Dependency
RUN apt-get install -y libcurl4-openssl-dev libssl-dev
RUN R -e "install.packages(c('repr', 'IRdisplay', 'evaluate', 'crayon', 'pbdZMQ', 'devtools', 'uuid', 'digest'), repos = 'http://cran.rstudio.com/')"
RUN R -e "install.packages('Rcpp', repos = 'http://cran.rstudio.com/')"
RUN R -e "install.packages('RSpectra', repos = 'http://cran.rstudio.com/')"
RUN R -e "require(devtools); install_version('igraph', version='1.0.1', repos='https://cran.rstudio.org/')"

# Install R packages in FlashX
WORKDIR /FlashX
RUN ln -sf ../FlashR FlashR
RUN ln -sf ../FlashGraphR FlashGraphR
RUN ./install_FlashR.sh
RUN ./install_FlashGraphR.sh
RUN R -e "install.packages('lbfgs', repos = 'http://cran.rstudio.com/')"
RUN R -e "devtools::install_github('flashxio/FlashR-learn')"

#check to see if it's there ^^^?

# Install Python packages.
RUN apt-get install -y python-pip python-numpy python-scipy cython
RUN pip install --upgrade pip
RUN pip install jupyter
RUN R -e "devtools::install_github('IRkernel/IRkernel')"
RUN R -e "IRkernel::installspec()"
RUN R -e "IRkernel::installspec(user = FALSE)"

RUN git clone https://github.com/flashxio/FlashPy.git
WORKDIR /FlashX/FlashPy
RUN python setup_sep.py install
WORKDIR /FlashX

RUN git clone https://github.com/flashxio/scikit-learn.git
RUN mv scikit-learn/fplearn/ /usr/local/lib/python2.7/dist-packages/

###FLASHX CONF END ###


ENV NOTVISIBLE "in users profile"
RUN echo "export VISIBLE=now" >> /etc/profile

EXPOSE 22
CMD ["/usr/sbin/sshd", "-D"]

EXPOSE 8888
CMD ["jupyter", "notebook", "--ip=0.0.0.0", "--no-browser", "--allow-root"]
