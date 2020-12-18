"""
This module provides functions that decode a string consisting of multiple distinct JSON
objects (and hence is not strictly speaking a JSON string itself). This is a totally
general procedure, and the logic is a little involved, so it gets its own module.
"""

def multiple_object_decoder_generator(s, decoder=ManifoldAP_Encoder):
    """
    For a string s that looks like multiple JSON objects, this returns a generator
    for each the decoded objects. In our applications, this will be the decoded
    ManifoldAP objects, hence the decoder default. Of course, there is nothing special
    here about the ManifoldAP_Encoder, so one could use this function with any json
    decoder.

    The trick we use is to pair braces. I.e. one object ends if and only if the number
    of left braces (or curly braces or whatever one wants to call them; this symbol: { )
    is equal to the number of right braces. We alternatively could look for something
    like "}{", but this could get messed up by whitespace. We could alternatively put
    in some marker to delineate between objects, but I think this is the most robust
    though perhaps not the cleanest or fastest way.
    """
    pass

def inside_string_positions(s, chars=["{", "}"]):
    """
    This takes some string s returns a list of indices that have characters that are
    between paired double quotation marks. The point is that for the type of strings
    we're considering (i.e. those that have multiple JSON entries), we don't want to
    parse anything that appears as some a string there. The particular issue we need to
    sidestep with multiple object decoding is that if we're pairing braces we might have
    some string like
        '{"name" : "bill}"}'
    If we were naively to pair braces, we would pair braces to find out where the object
    ends, we would think it ends in the wrong place, so we need to ignore things inside
    the paired double quotes. For JSON objects, by the way, only double quotes are
    allowed for strings.

    If I bothered to learn regular expressions, this could perhaps be simplified.
    """
    starting_index = 0
    while "\"" in s[starting_index:]:
        first_quote_mark_index = s.find("\"")
        second_quote_mark_index = s.find("\"", first_quote_mark_index+1)
