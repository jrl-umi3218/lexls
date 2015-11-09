include const_def.mk

#------------------------------------------------t
FILE = test
#-------------------------------------------------

test:
	$(CPP) $(CFLAGS) -o a $(FILE).cpp -I$(LEXLS_INCLUDE) $(EIGEN_INCLUDES) $(EIGEN_FLAGS)

clean:
	rm -f a

# dummy targets
.PHONY: clean
