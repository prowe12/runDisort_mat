# get them all
everything: disort_driver_mat

# This is a disort driver that can handle multiple pmom layers
disort_driver_mat: disort_driver_mat.f disort.f getmom.f errpack.f linpak.f rdi1mach.f bdref.f
	gfortran -O3 -o disort_driver_mat disort_driver_mat.f disort.f getmom.f errpack.f linpak.f rdi1mach.f bdref.f

