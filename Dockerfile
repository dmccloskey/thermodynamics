FROM dmccloskey/python3cobrapy
USER root

# install openbabel
RUN apt-get -y update && \
    apt-get install -y \
	swig \
	libeigen3-dev \
    cmake \
    make && \
    # apt-get clean && \
    # apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
	cd /usr/local \
	&& wget https://sourceforge.net/projects/openbabel/files/openbabel/2.4.1/openbabel-2.4.1.tar.gz	\
	&& tar -zxvf openbabel-2.4.1.tar.gz \
	&& mkdir build \
	&& cd build \
	&& cmake ../openbabel-2.4.1 -DPYTHON_BINDINGS=ON \
	&& make -j4 \
	&& make install
# # Cannot use openbabel/pybel
# # https://github.com/openbabel/openbabel/issues/368
# # No module named 'DLFCN'

# install PTVS
EXPOSE 3000
RUN pip3 install --no-cache-dir \
		ptvsd==3.0.0 \
	&&pip3 install --upgrade

# # Custom modules
# ENV CC_VERSION feature/dgf
# ENV CC_REPOSITORY https://github.com/dmccloskey/component-contribution.git
# RUN cd /usr/local/ && \
# 	#install io_utilities
# 	git clone ${CC_REPOSITORY} && \
# 	cd /usr/local/component-contribution/ && \
# 	git checkout ${CC_VERSION} && \
# 	python3 setup.py install && \
# 	cd /usr/local/

USER user