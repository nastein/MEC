test: currents_opt_EW_iso_sep.f90 mathtool.o test_iso.f90
	gfortran -o $@ $^

%.o: %.f90
	gfortran -o $@ -c $^

.PHONY: clean

clean:
	$(RM) *.mod *.o test
