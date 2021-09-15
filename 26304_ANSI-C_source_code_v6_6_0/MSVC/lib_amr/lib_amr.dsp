# Microsoft Developer Studio Project File - Name="lib_amr" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=lib_amr - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "lib_amr.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "lib_amr.mak" CFG="lib_amr - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "lib_amr - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "lib_amr - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "lib_amr - Win32 Release"

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
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD BASE RSC /l 0x1009 /d "NDEBUG"
# ADD RSC /l 0x1009 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "lib_amr - Win32 Debug"

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
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD BASE RSC /l 0x1009 /d "_DEBUG"
# ADD RSC /l 0x1009 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ENDIF 

# Begin Target

# Name "lib_amr - Win32 Release"
# Name "lib_amr - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE="..\..\c-code\lib_amr\dec_acelp.c"
# End Source File
# Begin Source File

SOURCE="..\..\c-code\lib_amr\dec_dtx.c"
# End Source File
# Begin Source File

SOURCE="..\..\c-code\lib_amr\dec_gain.c"
# End Source File
# Begin Source File

SOURCE="..\..\c-code\lib_amr\dec_if.c"
# End Source File
# Begin Source File

SOURCE="..\..\c-code\lib_amr\dec_lpc.c"
# End Source File
# Begin Source File

SOURCE="..\..\c-code\lib_amr\dec_main.c"
# End Source File
# Begin Source File

SOURCE="..\..\c-code\lib_amr\dec_rom.c"
# End Source File
# Begin Source File

SOURCE="..\..\c-code\lib_amr\dec_util.c"
# End Source File
# Begin Source File

SOURCE="..\..\c-code\lib_amr\enc_acelp.c"
# End Source File
# Begin Source File

SOURCE="..\..\c-code\lib_amr\enc_dtx.c"
# End Source File
# Begin Source File

SOURCE="..\..\c-code\lib_amr\enc_gain.c"
# End Source File
# Begin Source File

SOURCE="..\..\c-code\lib_amr\enc_if.c"
# End Source File
# Begin Source File

SOURCE="..\..\c-code\lib_amr\enc_lpc.c"
# End Source File
# Begin Source File

SOURCE="..\..\c-code\lib_amr\enc_main.c"
# End Source File
# Begin Source File

SOURCE="..\..\c-code\lib_amr\enc_rom.c"
# End Source File
# Begin Source File

SOURCE="..\..\c-code\lib_amr\enc_util.c"
# End Source File
# Begin Source File

SOURCE="..\..\c-code\lib_amr\if_rom.c"
# End Source File
# End Group
# End Target
# End Project
