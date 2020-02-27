IMAGE := kolmogorv:1.0
CONTAINER := kolmogorov-container

THIS_FILE := $(realpath $(lastword $(MAKEFILE_LIST)))
THIS_FILE_DIR := $(shell dirname $(THIS_FILE))

build: 
	DOCKER_BUILDKIT=1 docker build --tag $(IMAGE)  \
	    -f /$(THIS_FILE_DIR)/Dockerfile \
        .

run:
	docker run --rm \
		-it \
        -p 8888:8888 \
	    -v $(THIS_FILE_DIR)/../../code:/tf \
	    -v $(THIS_FILE_DIR)/../../data:/tf \
	    --name $(CONTAINER) \
	    $(IMAGE)

