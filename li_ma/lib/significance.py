import math

def eq_5(n_on, n_off, alpha):
    num = n_on - alpha * n_off
    den = math.sqrt(n_on + alpha**2 * n_off)
    return num / den

def eq_9(n_on, n_off, alpha):
    num = n_on - alpha * n_off
    den = math.sqrt(alpha * (n_on + n_off))
    return num / den

def eq_17(n_on, n_off, alpha):
    """
    with n_on = 0, n_off > 0, alpha > 0 we have:
        * f = 0, so
        * first as n * log(n) with n = 0, that is undefined, but limit = 0 for n -> 0
        * with n_on = 0, log(g) is only dependent from 1+alpha
        * so log(g) > 0 always
        * second = n_off * log(g)
    """
    if n_on <= 0 or n_off <= 0 or alpha <= 0:
        raise ValueError("Invalid parameter. Need: n_on > 0, n_off > 0, alpha > 0")
    fc = 1 / alpha * (1 + alpha)
    fb = n_on / (n_on + n_off)
    f  = fc * fb
    gc = 1 + alpha
    gb = n_off / (n_on + n_off)
    g  = gc * gb
    first  = n_on  * math.log(f)
    second = n_off * math.log(g)
    fullb  = first + second
    return math.sqrt(2) * math.sqrt(fullb)

