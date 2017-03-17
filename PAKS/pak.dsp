# Microsoft Developer Studio Project File - Name="pak" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=pak - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "pak.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "pak.mak" CFG="pak - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "pak - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "pak - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "pak - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
# ADD BASE F90 /compile_only /nologo /warn:nofileopt
# ADD F90 /compile_only /nologo /optimize:0 /warn:nofileopt
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib /nologo /subsystem:console /machine:I386

!ELSEIF  "$(CFG)" == "pak - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
# ADD BASE F90 /check:bounds /compile_only /dbglibs /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD F90 /compile_only /dbglibs /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt
# SUBTRACT F90 /check:bounds
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib /nologo /subsystem:console /incremental:no /debug /machine:I386 /pdbtype:sept

!ENDIF 

# Begin Target

# Name "pak - Win32 Release"
# Name "pak - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;f90;for;f;fpp"
# Begin Source File

SOURCE=.\pak0.for
DEP_F90_PAK0_=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak00.for
DEP_F90_PAK00=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak01.for
DEP_F90_PAK01=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak02.for
DEP_F90_PAK02=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak021.for
DEP_F90_PAK021=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak03.for
DEP_F90_PAK03=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak031.for
# End Source File
# Begin Source File

SOURCE=.\pak03g.for
# End Source File
# Begin Source File

SOURCE=.\pak04.for
# End Source File
# Begin Source File

SOURCE=.\pak05.for
DEP_F90_PAK05=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak051.for
DEP_F90_PAK051=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak06.for
DEP_F90_PAK06=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak061.for
DEP_F90_PAK061=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak061a.for
# End Source File
# Begin Source File

SOURCE=.\pak062.for
DEP_F90_PAK062=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak062a.for
# End Source File
# Begin Source File

SOURCE=.\pak07.for
DEP_F90_PAK07=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak071.for
DEP_F90_PAK071=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak08.for
DEP_F90_PAK08=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak081.for
DEP_F90_PAK081=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak09.for
DEP_F90_PAK09=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak10.for
DEP_F90_PAK10=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak101.for
DEP_F90_PAK101=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak102.for
DEP_F90_PAK102=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak103.for
# End Source File
# Begin Source File

SOURCE=.\pak105.for
DEP_F90_PAK105=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak11.for
DEP_F90_PAK11=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak111.for
DEP_F90_PAK111=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak112.for
DEP_F90_PAK112=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak113.for
DEP_F90_PAK113=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak12.for
DEP_F90_PAK12=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak121.for
DEP_F90_PAK121=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak122.for
DEP_F90_PAK122=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak123.for
DEP_F90_PAK123=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak13.for
DEP_F90_PAK13=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak1405.for
DEP_F90_PAK14=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak1411.for
DEP_F90_PAK141=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak1424.for
DEP_F90_PAK142=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak1430.for
DEP_F90_PAK143=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak15.for
DEP_F90_PAK15=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak21.for
DEP_F90_PAK21=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak22.for
DEP_F90_PAK22=\
	".\paka.inc"\
	".\pakxfem.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak221.for
DEP_F90_PAK221=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak222.for
DEP_F90_PAK222=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak223.for
DEP_F90_PAK223=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak22b1.for
DEP_F90_PAK22B=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak23.for
DEP_F90_PAK23=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak231.for
# End Source File
# Begin Source File

SOURCE=.\pak2405.for
DEP_F90_PAK24=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak2406.for
DEP_F90_PAK240=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak2407.for
DEP_F90_PAK2407=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak2408.for
DEP_F90_PAK2408=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak2409.for
DEP_F90_PAK2409=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak2410.for
DEP_F90_PAK241=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak2416.for
DEP_F90_PAK2416=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak2419.for
DEP_F90_PAK2419=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak2420.for
DEP_F90_PAK242=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak2421.for
DEP_F90_PAK2421=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak2422.for
DEP_F90_PAK2422=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak2425.for
DEP_F90_PAK2425=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak2426.for
DEP_F90_PAK2426=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak2428.for
DEP_F90_PAK2428=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak2430.for
DEP_F90_PAK243=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak2431.for
DEP_F90_PAK2431=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak2440.for
DEP_F90_PAK244=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak2451.for
DEP_F90_PAK245=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak2461.for
DEP_F90_PAK246=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak25.for
DEP_F90_PAK25=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak26.for
# End Source File
# Begin Source File

SOURCE=.\pak31.for
DEP_F90_PAK31=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak32.for
DEP_F90_PAK32=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak32b1.for
DEP_F90_PAK32B=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak33.for
DEP_F90_PAK33=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak3405.for
DEP_F90_PAK34=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak3406.for
DEP_F90_PAK340=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak3406k.for
DEP_F90_PAK3406=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak3407.for
DEP_F90_PAK3407=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak3409.for
DEP_F90_PAK3409=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak3410.for
DEP_F90_PAK341=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak3413.for
DEP_F90_PAK3413=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak3416.for
DEP_F90_PAK3416=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak3419.for
DEP_F90_PAK3419=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak3420.for
DEP_F90_PAK342=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak3421.for
DEP_F90_PAK3421=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak3422.for
DEP_F90_PAK3422=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak3423.for
DEP_F90_PAK3423=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak3424.for
DEP_F90_PAK3424=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak3425z.for
DEP_F90_PAK3425=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak3426.for
DEP_F90_PAK3426=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak3427.for
DEP_F90_PAK3427=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak3428.for
DEP_F90_PAK3428=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak3429.for
DEP_F90_PAK3429=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak3430.for
DEP_F90_PAK343=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak3431.for
DEP_F90_PAK3431=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak3441.for
DEP_F90_PAK344=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak3461.for
DEP_F90_PAK346=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak35.for
DEP_F90_PAK35=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak36.for
# End Source File
# Begin Source File

SOURCE=.\pak41.for
DEP_F90_PAK41=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak42.for
DEP_F90_PAK42=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak422.for
DEP_F90_PAK422=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak43.for
DEP_F90_PAK43=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak44.for
DEP_F90_PAK44=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak45.for
DEP_F90_PAK45=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak51.for
DEP_F90_PAK51=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak52.for
DEP_F90_PAK52=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak53.for
DEP_F90_PAK53=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak61.for
DEP_F90_PAK61=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak62.for
DEP_F90_PAK62=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak621.for
# End Source File
# Begin Source File

SOURCE=.\pak622.for
# End Source File
# Begin Source File

SOURCE=.\pak63.for
DEP_F90_PAK63=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak65.for
DEP_F90_PAK65=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak66.for
# End Source File
# Begin Source File

SOURCE=.\pak71.for
DEP_F90_PAK71=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak72.for
DEP_F90_PAK72=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak73.for
DEP_F90_PAK73=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak81.for
DEP_F90_PAK81=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak82.for
DEP_F90_PAK82=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak82b.for
DEP_F90_PAK82B=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak82f.for
DEP_F90_PAK82F=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak82i.for
DEP_F90_PAK82I=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak82p.for
# End Source File
# Begin Source File

SOURCE=.\pak82t.for
# End Source File
# Begin Source File

SOURCE=.\pak83.for
DEP_F90_PAK83=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak8405.for
DEP_F90_PAK84=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak8406.for
DEP_F90_PAK840=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak8410.for
DEP_F90_PAK841=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak8416.for
DEP_F90_PAK8416=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak8419.for
DEP_F90_PAK8419=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak8420.for
DEP_F90_PAK842=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak8421.for
DEP_F90_PAK8421=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak8432.for
DEP_F90_PAK843=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak8433.for
DEP_F90_PAK8433=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak8434.for
DEP_F90_PAK8434=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak8461.for
DEP_F90_PAK846=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak85.for
DEP_F90_PAK85=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak86.for
# End Source File
# Begin Source File

SOURCE=.\pak91.for
DEP_F90_PAK91=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak92.for
DEP_F90_PAK92=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak920.for
DEP_F90_PAK920=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak921.for
DEP_F90_PAK921=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak922.for
DEP_F90_PAK922=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak923a.for
DEP_F90_PAK923=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak93.for
DEP_F90_PAK93=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak931.for
DEP_F90_PAK931=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak932.for
DEP_F90_PAK932=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\PAK932p.for
# End Source File
# Begin Source File

SOURCE=.\pak933.for
DEP_F90_PAK933=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak935.for
DEP_F90_PAK935=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak9406.for
DEP_F90_PAK94=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pak9414.for
DEP_F90_PAK941=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pakf10.for
# End Source File
# Begin Source File

SOURCE=.\pakf11.for
# End Source File
# Begin Source File

SOURCE=.\pakf12.for
# End Source File
# Begin Source File

SOURCE=.\pakf2.for
# End Source File
# Begin Source File

SOURCE=.\pakf3.for
# End Source File
# Begin Source File

SOURCE=.\pakf4.for
# End Source File
# Begin Source File

SOURCE=.\pakf5.for
DEP_F90_PAKF5=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pakf6.for
# End Source File
# Begin Source File

SOURCE=.\pakf7.for
# End Source File
# Begin Source File

SOURCE=.\pakf8.for
# End Source File
# Begin Source File

SOURCE=.\pakf9.for
# End Source File
# Begin Source File

SOURCE=.\pakj.for
DEP_F90_PAKJ_=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pakj2d.for
# End Source File
# Begin Source File

SOURCE=.\pakj2dx.for
DEP_F90_PAKJ2=\
	".\pakxfem.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pakj3d.for
DEP_F90_PAKJ3=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pakjcir2d.for
# End Source File
# Begin Source File

SOURCE=.\pakjcir2dx.for
DEP_F90_PAKJC=\
	".\pakxfem.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pakjcir3d.for
# End Source File
# Begin Source File

SOURCE=.\pakjcor.for
# End Source File
# Begin Source File

SOURCE=.\pakjfc.for
# End Source File
# Begin Source File

SOURCE=.\pakjfront.for
# End Source File
# Begin Source File

SOURCE=.\pakjsec.for
# End Source File
# Begin Source File

SOURCE=.\pakkt.for
DEP_F90_PAKKT=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pakn.for
DEP_F90_PAKN_=\
	".\paka.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\pakn8.for
# End Source File
# Begin Source File

SOURCE=.\paku.for
# End Source File
# Begin Source File

SOURCE=.\pakx.for
DEP_F90_PAKX_=\
	".\pakxfem.inc"\
	
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl;fi;fd"
# End Group
# Begin Group "Resource Files"

# PROP Default_Filter "ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe"
# End Group
# End Target
# End Project
