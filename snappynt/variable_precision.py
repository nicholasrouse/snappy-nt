"""
Module for allowing functions that depends on some precision parameters to be varied.
I.e. if there is a function that needs to have precision variables specified at call
time, this module allows us to call the function many times for varying precision. This
was born out functions that can take a long time to run if you call them with high
precision even if the answer can be found much more quickly with lower precision, but
trying a lot of smaller precisions by hand is not ideal.

I don't really know how I want to do this just yet, so the module is blank and
untracked for now.
"""

    