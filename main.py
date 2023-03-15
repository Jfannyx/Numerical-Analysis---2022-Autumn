from algorithm import zeta, solveEquation, integral, solveODE
import gmpy2

if __name__ == "__main__":
    gmpy2.get_context().precision = 100
    precision1 = gmpy2.mpfr(0.5 * 10 ** (-28))
    precision2 = gmpy2.mpfr(0.5 * (10 ** (-4)))

    print("zeta(2) = " + str(zeta(2, precision2)))
    print("pi = " + "{0:.20Df}".format(round(gmpy2.root(gmpy2.mul(9450, zeta(8, precision1)), 8), 20)))
    print("solution of the zeta(x) = a is " + str(solveEquation(1.5, precision2)))
    print("the result of integral is " + str(integral(3.5, precision2)))
    print("y' = zeta(y), y(2) = 2, then y(4) = " + str(solveODE(2, 2, 4, precision2)))

    # print(zeta(2, precision2))
    # print(solveEquation(1.5, precision2))
    # print(integral(3.5, precision2)
    # print(solveODE(2, 2, 4, precision2))
