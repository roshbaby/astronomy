from math import fmod

########
# Positive fmod
########
def pfmod(x,y):
    result = fmod(x,y)
    if result < 0: result += y
    return result
