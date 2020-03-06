FROM tensorflow/tensorflow:latest-py3-jupyter

COPY requirements.txt /tmp/
COPY requirements_git.txt /tmp/

RUN apt-get update && apt-get -y install \
	gcc \
    graphviz \
    git

RUN pip3 install -r /tmp/requirements.txt
RUN pip3 install -r /tmp/requirements_git.txt

EXPOSE 8888
CMD ["jupyter","lab","--allow-root","--ip","0.0.0.0"]
