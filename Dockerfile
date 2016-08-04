# This is a comment
FROM ubuntu:14.04
MAINTAINER David Raffelt <draffelt@gmail.com>
RUN apt-get update && apt-get install -y git g++ python libeigen3-dev zlib1g-dev wget bsdtar
RUN wget -qO- http://github.com/MRtrix3/mrtrix3/archive/master.zip | bsdtar -xvf- && cd mrtrix3-master && python configure -nogui -verbose && NUMBER_OF_PROCESSORS=2 python build && find . -type d -exec chmod 755 {} + && find scripts -maxdepth 1 -type f -exec chmod 775 {} +

RUN wget -O- http://neuro.debian.net/lists/trusty.us-ca.full | tee /etc/apt/sources.list.d/neurodebian.sources.list
RUN apt-key adv --recv-keys --keyserver hkp://pgp.mit.edu:80 0xA5D32F012649A5A9 && apt-get update && DEBIAN_FRONTEND=noninteractive apt-get -y install fsl-5.0-eddy-nonfree ants

ENV FSLDIR=/usr/share/fsl/5.0
ENV FSLOUTPUTTYPE=NIFTI_GZ
ENV PATH=/usr/lib/fsl/5.0:$PATH
ENV PATH=/usr/lib/ants/:$PATH
ENV PATH=/mrtrix3-master/release/bin:$PATH
ENV PATH=/mrtrix3-master/scripts:$PATH
ENV FSLMULTIFILEQUIT=TRUE
ENV LD_LIBRARY_PATH=/usr/lib/fsl/5.0
ENV PYTHONPATH=/mrtrix3-master/scripts:$PYTHONPATH

RUN mkdir /bids_input
RUN mkdir /output
COPY run.py /code/run.py

ENTRYPOINT ["/code/run.py"]
