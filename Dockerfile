FROM dmccloskey/docker-openms:mrm_trgroup
USER root

# install PTVS
EXPOSE 3000
RUN pip3 install --no-cache-dir \
		ptvsd==3.0.0 \
	&&pip3 install --upgrade

# Custom modules
ENV IOUTILITIES_VERSION master
ENV IOUTILITIES_REPOSITORY https://github.com/dmccloskey/io_utilities.git
ENV BASE_VERSION master
ENV BASE_REPOSITORY https://github.com/dmccloskey/SBaaS_base.git
RUN cd /usr/local/ && \
	#install io_utilities
	git clone ${IOUTILITIES_REPOSITORY} && \
	cd /usr/local/io_utilities/ && \
	git checkout ${IOUTILITIES_VERSION} && \
	python3 setup.py install && \
	cd /usr/local/ && \
	#install SBaaS_base
	git clone ${BASE_REPOSITORY} && \
	cd /usr/local/SBaaS_base/ && \
	git checkout ${BASE_VERSION} && \
	python3 setup.py install

USER user