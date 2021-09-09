#
# Copyright 2013-2021 INRIA
#

fn = note

all:
	latex $(fn).tex
	latex $(fn).tex
	dvips $(fn)
	ps2pdf $(fn).ps

open:
	open $(fn).pdf

clean:
	rm -rf $(fn).aux
	rm -rf $(fn).log
	rm -rf $(fn).out
	rm -rf $(fn).dvi
	rm -rf $(fn).pfg
	rm -rf $(fn).bbl
	rm -rf $(fn).blg
	rm -rf $(fn).toc
	rm -rf $(fn).ps
#	rm -rf $(fn).pdf

# dummy targets
.PHONY: clean
