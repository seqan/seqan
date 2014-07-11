rem rmdir /Q /S dist
del /s dist\*.py
del /s dist\*.pyc
del /s dist\*.pyd
C:\Python27\python.exe setup.py py2exe
rmdir /Q /S build
pause