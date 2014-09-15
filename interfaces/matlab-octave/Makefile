info:
	@echo "Please use 'make matlab' or 'make octave'."

matlab:
	matlab -nosplash -nodesktop -nojvm -r 'make-matlab; quit' | tail +10

octave:
	octave --no-history --silent make-octave.m

clean:
	rm -f *.o
	rm -f *.mexmaci*

forceclean: clean
	rm -f *.mex
