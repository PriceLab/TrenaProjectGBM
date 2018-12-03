all:  docs install

docs:
	R -e "devtools::document()"
build:
	(cd ..; R CMD build --no-build-vignettes TrenaProjectGBM)

install:
	(cd ..; R CMD INSTALL TrenaProjectGBM)

check:
	(cd ..; R CMD check `ls -t TrenaProjectGBM) | head -1`)

tests:
	for x in inst/unitTests/test_*.R; do echo $$x; R -f $$x; done
