@echo off
REM ============================================================
REM  build.bat  --  Compila THS MC con gfortran MSYS2 ucrt64
REM  Uso dalla cartella ths_mc\:
REM    build.bat          (release: -O3)
REM    build.bat debug    (debug:   -O0 -fcheck=all)
REM    build.bat clean    (rimuove obj/ e ths_mc.exe)
REM ============================================================

set MSYS2=C:\msys64
set PATH=%MSYS2%\ucrt64\bin;%MSYS2%\usr\bin;%PATH%
set TMP=%LOCALAPPDATA%\Temp
set TEMP=%LOCALAPPDATA%\Temp

if "%1"=="clean" (
    echo Pulizia...
    make clean
    goto :end
)

if "%1"=="debug" (
    echo Compilazione DEBUG...
    make debug
    goto :end
)

echo Compilazione RELEASE...
make

:end
if %ERRORLEVEL%==0 (
    echo.
    echo === Build OK: ths_mc.exe ===
) else (
    echo.
    echo === ERRORE di compilazione ===
)
