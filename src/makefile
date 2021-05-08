TDFEM.exe:*.f90
	ifort -g -traceback  -mkl -O3 -o TDFEM.exe CALCULATE11.f90 DELTAT.f90 GET_PP.f90 GRID.f90 HECHENG.f90 YOUHECHENG.f90 TDFEM.f90 

clean:
	 rm -rf *.exe *.mod  ooc*.* nohup* E* YD* *.DAT
