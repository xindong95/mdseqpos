import math

def qnorm(p, upper, mu = None, sigma2 = None):
    
    if p < 0 or p > 1:
        raise Exception("Illegal argument %f for qnorm(p)." % p)
    
    split = 0.42
    a0 = 2.50662823884
    a1 = -18.61500062529
    a2 = 41.39119773534
    a3 = -25.44106049637
    b1 = -8.47351093090
    b2 = 23.08336743743
    b3 = -21.06224101826
    b4 = 3.13082909833
    c0 = -2.78718931138
    c1 = -2.29796479134
    c2 = 4.85014127135
    c3 = 2.32121276858
    d1 = 3.54388924762
    d2 = 1.63706781897
    q = p - 0.5
    
    r = 0.0
    ppnd = 0.0
    
    if abs(q) <= split:
        r = q * q
        ppnd = q * (((a3 * r + a2) * r + a1) * r + a0) / ((((b4 * r + b3) * r + b2) * r + b1) * r + 1)
    else:
        r = p
        if q > 0:
            r = 1 - p
        
        if r > 0:
            r = math.sqrt(- math.log(r))
            ppnd = (((c3 * r + c2) * r + c1) * r + c0) / ((d2 * r + d1) * r + 1)
            
            if q < 0:
                ppnd = - ppnd
        else:
            ppnd = 0
            
    if upper:
        ppnd = - ppnd
    
    if mu != None and sigma2 != None:
        return ppnd * math.sqrt(sigma2) + mu
    else:
        return ppnd

def pnorm(z, upper, mu = 0.0, sigma2 = 1.0):
    
    z = (z - mu) / math.sqrt(sigma2)
    
    ltone = 7.0
    utzero = 18.66
    con = 1.28
    a1 = 0.398942280444
    a2 = 0.399903438504
    a3 = 5.75885480458
    a4 = 29.8213557808
    a5 = 2.62433121679
    a6 = 48.6959930692
    a7 = 5.92885724438
    b1 = 0.398942280385
    b2 = 3.8052e-8
    b3 = 1.00000615302
    b4 = 3.98064794e-4
    b5 = 1.986153813664
    b6 = 0.151679116635
    b7 = 5.29330324926
    b8 = 4.8385912808
    b9 = 15.1508972451
    b10 = 0.742380924027
    b11 = 30.789933034
    b12 = 3.99019417011

    y = 0.0
    alnorm = 0.0
    
    if z < 0:
        upper = not upper
        z = - z
    
    if z <= ltone or upper and z <= utzero:
        y = 0.5 * z * z
        if z > con:
            alnorm = b1 * math.exp(- y) / (z - b2 + b3 / (z + b4 + b5 / (z - b6 + b7 / (z + b8 - b9 / (z + b10 + b11 / (z + b12))))))
        else:
            alnorm = 0.5 - z * (a1 - a2 * y / (y + a3 - a4 / (y + a5 + a6 / (y + a7))))
    else:
        alnorm = 0.0
    
    if not upper:
        alnorm = 1.0 - alnorm
    
    return alnorm