test: currents_opt_v2.f90 mathtool.o test_iso.f90
	gfortran -o $@ $^

%.o: %.f90
	gfortran -o $@ -c $^

.PHONY: clean

clean:
	$(RM) *.mod *.o test
