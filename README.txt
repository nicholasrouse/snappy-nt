Oct-18-2021

Known issues:
If one tries to run

$ sage --python
>>> import snappynt.QuaternionAlgebraNF

one will get a mysterious import error. The solution is to either run sage rather than
sage's python, or first import sage.all. See https://groups.google.com/g/sage-support/c/3lXInDgHJe4/m/SdOpuuOkqPUJ
for more details.
