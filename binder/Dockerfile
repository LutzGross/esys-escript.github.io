FROM aellery/esys.escript.binder

LABEL maintainer="a.ellery@uq.edu.au"

# set up the user
ARG NB_USER=jovyan
ARG NB_UID=1000
ENV USER ${NB_USER}
ENV NB_UID ${NB_UID}
ENV HOME /home/${NB_USER}

RUN adduser --disabled-password \
    --gecos "Default user" \
    --uid ${NB_UID} \
    ${NB_USER}

COPY . ${HOME}
USER root
RUN chown -R ${NB_UID} ${HOME}
USER ${NB_USER}

ENV PYTHONPATH=/usr/esys/

CMD ["/bin/bash"]
