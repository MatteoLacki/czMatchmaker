from dplython import DelayFunction

@DelayFunction
def Round_and_Multi(x, times = 100):
    return [ int(round(a*times)) for a in x ]

@DelayFunction
def Round(x):
    return [ round(a) for a in x ]

@DelayFunction
def crossprod( x, y):
    return [ xx*yy for xx,yy in zip(x,y) ]
