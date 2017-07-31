FROM dmccloskey/python3cobrapy
USER root

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