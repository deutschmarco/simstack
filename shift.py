def shift_twod(seq, x, y):
    from numpy import roll
    out = roll(roll(seq, int(x), axis = 1), int(y), axis = 0)
    return out 

def shift(seq, x):
    from numpy import roll
    out = roll(seq, int(x))
    return out 
