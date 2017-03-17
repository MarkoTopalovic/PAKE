echo OFF
setlocal

set BLACSdir="..\..\..\blacs\SRC\MPI"
set blacs_pinfo=%BLACSdir%\blacs_pinfo_.c

copy /Y %blacs_pinfo% %BLACSdir%\blacs_pinfo_MainInF77.c
copy /Y %blacs_pinfo% %BLACSdir%\blacs_pinfo_MainInF77_CallFromC.c
copy /Y %blacs_pinfo% %BLACSdir%\blacs_pinfo_MainInC.c
copy /Y %blacs_pinfo% %BLACSdir%\blacs_pinfo_MainInC_CallFromC.c

endlocal


