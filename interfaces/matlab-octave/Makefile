#
# Copyright 2013-2021 INRIA
#

info:
	@echo "Please use 'make matlab' or 'make octave'."

matlab:
	matlab -nosplash -nodesktop -nojvm -r 'make_matlab; quit' | tail +11

octave:
	octave --no-history --silent make_octave.m

clean:
	rm -f *.o
	rm -f *.mexmaci*
	rm -f *.mex
