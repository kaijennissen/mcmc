build_all:
	-make -C docker/julia build
	-make -C docker/jupyter build
	-make -C docker/python build
	-make -C docker/rstudio build

run_all:
	-make -C docker/julia run
	-make -C docker/jupyter run
	-make -C docker/python run
	-make -C docker/rstudio run

rm_all:
	-make -C docker/julia rm 
	-make -C docker/jupyter rm
	-make -C docker/python rm
	-make -C docker/rstudio rm

