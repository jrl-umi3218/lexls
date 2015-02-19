HOT_START=hot_start

all: html hot_start

html:
	doxygen Doxyfile

hot_start:
	latex $(HOT_START).tex
	latex $(HOT_START).tex
	dvips -o $(HOT_START).ps $(HOT_START).dvi
	ps2pdf -sPAPERSIZE=a4 $(HOT_START).ps $(HOT_START).pdf

clean:
	rm -rf html 
	rm -rf $(HOT_START).aux 
	rm -rf $(HOT_START).log
	rm -rf $(HOT_START).dvi
	rm -rf $(HOT_START).out
	rm -rf $(HOT_START).ps
	rm -rf $(HOT_START).pdf

# dummy targets
.PHONY: clean
