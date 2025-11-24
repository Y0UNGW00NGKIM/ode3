# ode 3

Run make in ode 3 directory. Will produce executable named vterm in src directory.

For energy conservation study, use ./vterm -E -n N where N is number of RK4 steps. For terminal velocity study, run ./vterm.
It will generate vterm.pdf. 

Energy conservation test (no air resistance)
  nsteps             = 100
  t_max              = 5 s
  E_min              = 200 J
  E_max              = 200 J
  relative variation = 8.2423e-15

Energy conservation test (no air resistance)
  nsteps             = 500
  t_max              = 5 s
  E_min              = 200 J
  E_max              = 200 J
  relative variation = 7.84439e-14

Energy conservation test (no air resistance)
  nsteps             = 1000
  t_max              = 5 s
  E_min              = 200 J
  E_max              = 200 J
  relative variation = 3.48166e-14

Energy conservation test (no air resistance)
  nsteps             = 5000
  t_max              = 5 s
  E_min              = 200 J
  E_max              = 200 J
  relative variation = 3.28697e-13

  Relative variation is (Emax - Emin)/Eavg, so considering that Emin/Emax are very close + relative variation is extremely small, this shows that energy is essentially conserved
  in vacuum case.

For each mass in the scan, 
(phys56xx) [bab6cw src]$ ./vterm
m = 0.001 kg, v_t (numeric)  = 0.313209 m/s, v_t (analytic) = 0.313209 m/s
m = 0.257385 kg, v_t (numeric)  = 5.02488 m/s, v_t (analytic) = 5.02488 m/s
m = 0.513769 kg, v_t (numeric)  = 7.09935 m/s, v_t (analytic) = 7.09935 m/s
m = 0.770154 kg, v_t (numeric)  = 8.69207 m/s, v_t (analytic) = 8.69207 m/s
m = 1.02654 kg, v_t (numeric)  = 10.0351 m/s, v_t (analytic) = 10.0351 m/s
m = 1.28292 kg, v_t (numeric)  = 11.2185 m/s, v_t (analytic) = 11.2185 m/s
m = 1.53931 kg, v_t (numeric)  = 12.2885 m/s, v_t (analytic) = 12.2885 m/s
m = 1.79569 kg, v_t (numeric)  = 13.2724 m/s, v_t (analytic) = 13.2724 m/s
m = 2.05208 kg, v_t (numeric)  = 14.1883 m/s, v_t (analytic) = 14.1883 m/s
m = 2.30846 kg, v_t (numeric)  = 15.0486 m/s, v_t (analytic) = 15.0486 m/s
m = 2.56485 kg, v_t (numeric)  = 15.8623 m/s, v_t (analytic) = 15.8623 m/s
m = 2.82123 kg, v_t (numeric)  = 16.6362 m/s, v_t (analytic) = 16.6362 m/s
m = 3.07762 kg, v_t (numeric)  = 17.3757 m/s, v_t (analytic) = 17.3757 m/s
m = 3.334 kg, v_t (numeric)  = 18.0849 m/s, v_t (analytic) = 18.0849 m/s
m = 3.59038 kg, v_t (numeric)  = 18.7674 m/s, v_t (analytic) = 18.7674 m/s
m = 3.84677 kg, v_t (numeric)  = 19.426 m/s, v_t (analytic) = 19.426 m/s
m = 4.10315 kg, v_t (numeric)  = 20.0629 m/s, v_t (analytic) = 20.0629 m/s
m = 4.35954 kg, v_t (numeric)  = 20.6802 m/s, v_t (analytic) = 20.6802 m/s
m = 4.61592 kg, v_t (numeric)  = 21.2796 m/s, v_t (analytic) = 21.2796 m/s
m = 4.87231 kg, v_t (numeric)  = 21.8626 m/s, v_t (analytic) = 21.8626 m/s
m = 5.12869 kg, v_t (numeric)  = 22.4304 m/s, v_t (analytic) = 22.4304 m/s
m = 5.38508 kg, v_t (numeric)  = 22.9843 m/s, v_t (analytic) = 22.9843 m/s
m = 5.64146 kg, v_t (numeric)  = 23.525 m/s, v_t (analytic) = 23.525 m/s
m = 5.89785 kg, v_t (numeric)  = 24.0537 m/s, v_t (analytic) = 24.0537 m/s
m = 6.15423 kg, v_t (numeric)  = 24.5709 m/s, v_t (analytic) = 24.5709 m/s
m = 6.41062 kg, v_t (numeric)  = 25.0775 m/s, v_t (analytic) = 25.0775 m/s
m = 6.667 kg, v_t (numeric)  = 25.5741 m/s, v_t (analytic) = 25.5741 m/s
m = 6.92338 kg, v_t (numeric)  = 26.0611 m/s, v_t (analytic) = 26.0612 m/s
m = 7.17977 kg, v_t (numeric)  = 26.5393 m/s, v_t (analytic) = 26.5393 m/s
m = 7.43615 kg, v_t (numeric)  = 27.009 m/s, v_t (analytic) = 27.009 m/s
m = 7.69254 kg, v_t (numeric)  = 27.4706 m/s, v_t (analytic) = 27.4707 m/s
m = 7.94892 kg, v_t (numeric)  = 27.9247 m/s, v_t (analytic) = 27.9247 m/s
m = 8.20531 kg, v_t (numeric)  = 28.3714 m/s, v_t (analytic) = 28.3715 m/s
m = 8.46169 kg, v_t (numeric)  = 28.8112 m/s, v_t (analytic) = 28.8113 m/s
m = 8.71808 kg, v_t (numeric)  = 29.2445 m/s, v_t (analytic) = 29.2445 m/s
m = 8.97446 kg, v_t (numeric)  = 29.6713 m/s, v_t (analytic) = 29.6714 m/s
m = 9.23085 kg, v_t (numeric)  = 30.0922 m/s, v_t (analytic) = 30.0923 m/s
m = 9.48723 kg, v_t (numeric)  = 30.5072 m/s, v_t (analytic) = 30.5073 m/s
m = 9.74362 kg, v_t (numeric)  = 30.9166 m/s, v_t (analytic) = 30.9168 m/s
m = 10 kg, v_t (numeric)  = 31.3207 m/s, v_t (analytic) = 31.3209 m/s
I compared numerical and analytical terminal speed values across the given range, and also printed out the relative error vs mass graph.
It can be clearly seen from the data printed by the code that they agree with maximum discrepancy being +-0.0002, so I would say that they are reasonably accurate.
