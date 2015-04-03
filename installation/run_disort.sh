rm disort.exe
export OutputFile="my_output.txt"
date > $OutputFile
cat disotest.f disort.f errpack.f linpak.f rdi1mach.f > code.f
# compile into disort.exe executable
gfortran -O3 code.f -o disort.exe >> $OutputFile
# make disort.exe executable
chmod u+x ./disort.exe
# "./" ensures local disort.exe is executed
./disort.exe >> $OutputFile
