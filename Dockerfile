FROM bids/base_validator

# MAINTAINER David Raffelt <draffelt@gmail.com>

ARG DEBIAN_FRONTEND="noninteractive"

RUN apt-get update -qq && \
    apt-get install -qq --no-install-recommends -y git python libeigen3-dev zlib1g-dev wget bsdtar software-properties-common

RUN add-apt-repository ppa:ubuntu-toolchain-r/test
RUN apt-get update && apt-get install -y g++-5

RUN wget --progress=dot:giga -O- http://neuro.debian.net/lists/trusty.us-ca.full | tee /etc/apt/sources.list.d/neurodebian.sources.list
RUN apt-key adv --recv-keys --keyserver hkp://pgp.mit.edu:80 0xA5D32F012649A5A9 && \
    apt-get update  -qq && \
    apt-get install -qq y --no-install-recommends fsl-5.0-eddy-nonfree ants

ENV CXX=/usr/bin/g++-5

RUN git clone https://github.com/MRtrix3/mrtrix3.git mrtrix3 && \
    cd mrtrix3 && \
    git checkout tag_0.3.16 && \
    python configure -nogui
RUN if [ "$CIRCLECI" = "true" ]; then cd mrtrix3 && NUMBER_OF_PROCESSORS=1 python build; else cd mrtrix3 && python build; fi

ENV FSLDIR=/usr/share/fsl/5.0
ENV FSLOUTPUTTYPE=NIFTI_GZ
ENV PATH=/usr/lib/fsl/5.0:$PATH
ENV PATH=/usr/lib/ants/:$PATH
ENV PATH=/mrtrix3/release/bin:$PATH
ENV PATH=/mrtrix3/scripts:$PATH
ENV FSLMULTIFILEQUIT=TRUE
ENV LD_LIBRARY_PATH=/usr/lib/fsl/5.0
ENV PYTHONPATH=/mrtrix3/scripts:$PYTHONPATH

RUN mkdir /bids_input && \
    mkdir /output
COPY run.py /code/run.py

COPY version /version

ENTRYPOINT ["/code/run.py"]
