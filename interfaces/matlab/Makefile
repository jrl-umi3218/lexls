build:
#	matlab -nosplash -nodesktop -nojvm -r 'make; quit' | tail +10
	octave --no-history --silent make.m

clean:
	rm -f *.o
	rm -f *.mexmaci*

forceclean: clean
	rm -f *.mex
