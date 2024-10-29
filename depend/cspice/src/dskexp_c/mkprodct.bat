rem 
rem    mkcdskexp.bat
rem 
rem    Creates dskexp.exe for MS Visual C++ and moves it to the
rem    appropriate Toolkit directory.
rem 
rem
rem    Version 1.1.0  19-OCT-2003 (BVS)
rem
rem       added -DNON_ANSI_STDIO compile option.
rem
rem    Version 1.0.0  29-DEC-1998 (NJB) 
rem


set cl= /c /O2 -D_COMPLEX_DEFINED -DMSDOS -DNON_ANSI_STDIO

copy dskexp.pgm main.c

for %%f in (*.c) do cl %%f 

dir /b *.obj > temp.lst

link -lib /out:dskexp.lib  @temp.lst

copy main.x dskexp.c

cl dskexp.c

link /STACK:16000000 dskexp.obj dskexp.lib   ..\..\lib\csupport.lib ..\..\lib\cspice.lib

move dskexp.exe  ..\..\exe

del *.obj
del dskexp.lib
del temp.lst

