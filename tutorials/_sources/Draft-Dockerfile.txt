.. code:: {.sh}

    docker run -itP ubuntu:14.04

    sudo apt-get update
    sudo apt-get update
    sudo apt-get install -y git cmake g++
    sudo apt-get install -y libboost-dev libboost-system-dev libboost-filesystem-dev libnuma-dev libaio-dev
    sudo apt-get install -y libstxxl-dev zlib1g-dev
    git clone https://github.com/icoming/FlashGraph.git

    sudo apt-get install wget

    wget http://trilinos.csbsju.edu/download/files/trilinos-12.0.1-Source.tar.bz2
    tar -jxvf trilinos-12.0.1-Source.tar.bz2
    rm trilinos-12.0.1-Source.tar.bz2

    sudo apt-get install -y gfortran libatlas-dev liblapack-dev

    cd trilinos-12.0.1-Source/

    echo ' #!/bin/sh
    EXTRA_ARGS=$@
    SOURCE_BASE=..
    cmake \
        -D CMAKE_BUILD_TYPE:STRING=RELEASE \
        -D BUILD_SHARED_LIBS:BOOL=ON \
        -D Trilinos_ENABLE_TESTS=OFF \
        $EXTRA_ARGS \
        ${SOURCE_BASE}
    ' > do-configure


    chmod u+x do-configure

    mkdir build
    cd build/
    ../do-configure -DTrilinos_ENABLE_Anasazi=ON
    make -j32
    make install

    cd /FlashGraph

    git checkout dev

    apt-get install -y libhwloc-dev
    mkdir build; cd build; cmake ../; make
