# Microsoft Developer Studio Project File - Name="encoder" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=encoder - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "encoder.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "encoder.mak" CFG="encoder - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "encoder - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "encoder - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "encoder - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD BASE RSC /l 0x1009 /d "NDEBUG"
# ADD RSC /l 0x1009 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 3gplib.lib er-libisomedia.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386 /out:"..\..\encoder.exe" /libpath:"../../c-code/3gplib/Release/"

!ELSEIF  "$(CFG)" == "encoder - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /FR /YX /FD /GZ /c
# ADD BASE RSC /l 0x1009 /d "_DEBUG"
# ADD RSC /l 0x1009 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 3gplib.lib er-libisomedia.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /out:"..\..\encoder.exe" /pdbtype:sept /libpath:"../../c-code/3gplib/Debug/"

!ENDIF 

# Begin Target

# Name "encoder - Win32 Release"
# Name "encoder - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE="..\..\c-code\encoder\avq_cod.c"
# End Source File
# Begin Source File

SOURCE="..\..\c-code\encoder\c_stereo_x.c"
# End Source File
# Begin Source File

SOURCE="..\..\c-code\encoder\cod_ace.c"
# End Source File
# Begin Source File

SOURCE="..\..\c-code\encoder\cod_cp_state.c"
# End Source File
# Begin Source File

SOURCE="..\..\c-code\encoder\cod_hf.c"
# End Source File
# Begin Source File

SOURCE="..\..\c-code\encoder\cod_hi_stereo.c"
# End Source File
# Begin Source File

SOURCE="..\..\c-code\encoder\cod_lf.c"
# End Source File
# Begin Source File

SOURCE="..\..\c-code\encoder\cod_lf_b.c"
# End Source File
# Begin Source File

SOURCE="..\..\c-code\encoder\cod_main.c"
# End Source File
# Begin Source File

SOURCE="..\..\c-code\encoder\cod_tcx.c"
# End Source File
# Begin Source File

SOURCE="..\..\c-code\encoder\cod_tcx_stereo.c"
# End Source File
# Begin Source File

SOURCE="..\..\c-code\encoder\enc_prm.c"
# End Source File
# Begin Source File

SOURCE="..\..\c-code\encoder\enc_wbplus.c"
# End Source File
# Begin Source File

SOURCE="..\..\c-code\encoder\get_gain.c"
# End Source File
# Begin Source File

SOURCE="..\..\c-code\encoder\lag_wind.c"
# End Source File
# Begin Source File

SOURCE="..\..\c-code\encoder\nclass.c"
# End Source File
# Begin Source File

SOURCE="..\..\c-code\encoder\q_gain2p.c"
# End Source File
# Begin Source File

SOURCE="..\..\c-code\encoder\re8_cod.c"
# End Source File
# Begin Source File

SOURCE="..\..\c-code\encoder\segsnr.c"
# End Source File
# Begin Source File

SOURCE="..\..\c-code\encoder\wb_vad.c"
# End Source File
# End Group
# End Target
# End Project
