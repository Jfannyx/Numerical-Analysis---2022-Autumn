import gmpy2
import math
import decimal

## Q1 Summation of series
def zeta(x, precision):
    precision = gmpy2.mpfr(precision)
    i = 1
    N = gmpy2.root(gmpy2.div(1,gmpy2.mul(x - 1, precision / 2)), int(x - 1))
    p = int(gmpy2.rint_ceil(gmpy2.log10(gmpy2.div(N, gmpy2.mul(2, precision / 2)))))
    ksi = gmpy2.mpfr(0)

    while i < N:
        ksi = ksi + round(gmpy2.mpfr(1 / i ** x), p)
        i = i + 1
    ksi = ksi + round(gmpy2.mpfr(1 / i ** x), p)

    return round(ksi, p)

## Q2 Root of the equation (Newton's method)

def f2(x, precision):
    x = gmpy2.mpfr(x)
    precision = gmpy2.mpfr(precision)
    i = 1
    ksi = gmpy2.mpfr(0)
    dksi = gmpy2.mpfr(0)
    while gmpy2.div(1, i ** x) > precision:
        ksi = ksi + gmpy2.mpfr(1 / i ** x)
        dksi = dksi + gmpy2.mpfr((1 / i) ** x * gmpy2.log (1 / i))
        i = i + 1
    ksi = ksi + gmpy2.mpfr(1 / i ** x)
    dksi = dksi + gmpy2.mpfr((1 / i) ** x * gmpy2.log (1 / i))
    return ksi, dksi

def solveEquation(a, precision):
    if a > 1.1:
        x = 2
    else:
        x = 4
    while 1:
        ksi, dksi = f2(x, precision ** 2)
        x = x -  (ksi - a) / dksi
        if gmpy2.mul(f2(x - precision, precision ** 2)[0] - a, f2(x + precision, precision ** 2)[0] - a) <= 0:
            break
    x = decimal.Decimal(str(x), decimal.getcontext())
    return x.__round__(math.floor(-math.log10(precision * 2)))

## Q3 Numerical solution of definite integrals
def rightError(x, Delta):
    b = gmpy2.root(gmpy2.div(2,gmpy2.mul(x - 2, Delta)), math.floor(x - 2))
    return b

def interval3(x, b, Delta):
    df2 = gmpy2.div((x - 1) ** (x - 1), gmpy2.exp(x - 1))
    h = gmpy2.root(gmpy2.div(gmpy2.div(gmpy2.mul(6, Delta), b), df2),2)
    return h


def integral(x, precision):
    b = rightError(x, gmpy2.div(gmpy2.root(precision, 2),2))
    h = interval3(x, b, gmpy2.div(gmpy2.root(precision, 2),2))
    xk = gmpy2.mpfr(0)
    x = gmpy2.mpfr(3.5)
    s = gmpy2.div(xk ** (x - 1), gmpy2.exp(xk))
    xk = xk + h
    while(xk < b):
        s = s + 2 * xk ** (x - 1) / gmpy2.exp(xk)
        xk = xk + h
    s = (s + b ** (x - 1) / (gmpy2.exp(b) - 1)) * h / 2 * zeta(x, precision)
    s = decimal.Decimal(str(s), decimal.getcontext())
    return s.__round__(math.floor(-math.log10(precision * 2)))

## Q4 Numerical solution of differential equations

def interval4(x, precision):
    m = gmpy2.rint_floor(gmpy2.log10(gmpy2.root(gmpy2.div(precision, 2), 2)))
    h = gmpy2.exp10(m)
    return h

def solveODE(x0, y0, x, precision):
    '''
    x0,y0: the start point
    x: the end
    '''
    xn = x0
    yn = y0
    h = interval4(x, precision / 2)
    p = math.ceil(-math.log10(precision) + math.log10(2) + 1) 
    while (xn < x - h):
        fyn = zeta(yn, precision / 2)
        yp = round(yn + h * fyn, p) 
        yn1 = round(yn + h / 2 * (fyn + zeta(yp, precision / 2)),p)
        xn = xn + h
        yn = yn1
    yn1 = decimal.Decimal(str(yn1), decimal.getcontext())
    return yn1.__round__(math.floor(-math.log10(precision * 2)))

if __name__ == "__main__":
    gmpy2.get_context().precision = 100
    precision1 = gmpy2.mpfr(0.5 * 10 ** (-28))
    precision2 = gmpy2.mpfr(0.5 * (10 ** (-4)))

    # print("zeta(2) = " + str(zeta(2, precision2)))
    # print("pi = " + "{0:.20Df}".format(round(gmpy2.root(gmpy2.mul(9450, zeta(8, precision1)), 8), 20)))
    # print("solution of the zeta(x) = a is " + str(solveEquation(1.5, precision2)))
    # print("the result of integral is " + str(integral(3.5, precision2)))
    # print("y' = zeta(y), y(2) = 2, then y(4) = " + str(solveODE(2, 2, 4, precision2)))