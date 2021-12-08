FROM ubuntu:18.04
RUN apt update && apt-get -y install bc wget build-essential zlib1g-dev libbz2-dev liblzma-dev libncurses5-dev libncursesw5-dev libcurl4-openssl-dev
RUN mkdir -p /opt/tools
WORKDIR /opt/tools
RUN wget "https://downloads.sourceforge.net/project/samtools/samtools/1.14/samtools-1.14.tar.bz2?ts=gAAAAABhsM7LzTZl5jRULHl2VHSsINSj7l8mYwUSjyXhNBDANG6xG5BAAn24ET5RNgAZStxyww77QgASvHlcQ7A0BaI34VEvkQ%3D%3D&r=https%3A%2F%2Fsourceforge.net%2Fprojects%2Fsamtools%2Ffiles%2Fsamtools%2F1.14%2Fsamtools-1.14.tar.bz2%2Fdownload" -O samtools-1.14.tar.bz2 && tar xvjf samtools-1.14.tar.bz2
WORKDIR /opt/tools/samtools-1.14
RUN ./configure --prefix=/opt/tools/samtools/ --without-curses
RUN make & make install
ENV PATH=/opt/tools/samtools/bin/:$PATH
