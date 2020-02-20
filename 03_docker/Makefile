IMAGE = metropolis:1.0
CONTAINER = metropolis-container


THIS_FILE := $(realpath $(lastword $(MAKEFILE_LIST)))
THIS_FILE_DIR := $(shell dirname $(THIS_FILE))

build: 
	DOCKER_BUILDKIT=1 docker build --file $(THIS_FILE_DIR)/Dockerfile \
    --tag $(IMAGE)  \
    .

run:
	docker run --rm \
		-itd \
		-p 7772:22 \
	    -v $(THIS_FILE_DIR):/home \
	    --name $(CONTAINER) \
	    $(IMAGE)
