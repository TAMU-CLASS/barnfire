"""
Andrew Till
Fall 2012

Handy plotting functions for Python
"""

#stdlib
import math

#TPL: Numpy and Matplotlib
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

#Return color option
def get_simple_colors(it_cur, it_tot=0):
    if it_tot == 0:
        it_tot = it_cur
    col = ['b', 'r', 'g', 'c', 'm', 'k', 'y']
    if it_tot < len(col):
        return col[it_cur]
    else:
        return str(float(it_cur)/it_tot)

#Return color from a colormap (cm)
def get_colors(size, colorMap='Paired', offsetLower=0, offsetUpper=0):
    colorLocs = np.linspace(0+offsetLower, 1-offsetUpper, size)
    cm = getattr(mpl.cm, colorMap)
    return cm(colorLocs)

#Return linestyle option
def get_linestyle(it_cur):
    lin = ['-','--','-.',':']
    return lin[it_cur % len(lin)]

#Return marker option
def get_marker(it_cur):
    mark = ['D', 's', '^', 'v', 'd', 'h', '*', 'o', 'p', 'H', '<', '>', r'$\clubsuit$', r'$\spadesuit$', '8', '1', '3', '2', '4', '+', 'x', ',']
    return mark[it_cur % len(mark)]

#Get best-fit line of form y = a*x^order
def get_order_line(xs, ys, order):
    xs = np.array(xs, dtype=np.float)
    ys = np.array(ys, dtype=np.float)
    xsOrder = np.power(xs, order)
    coeff = np.dot(ys, xsOrder) / np.dot(xsOrder, xsOrder)
    return coeff * xsOrder

# Get least-squares fit of form y = a*x^b
def get_fit_log(xs, ys):
    xs = np.array(xs, dtype=np.float)
    ys = np.array(ys, dtype=np.float)
    mask = np.logical_and(xs > 0, ys > 0)
    xsLog = np.log(xs[mask])
    ysLog = np.log(ys[mask])
    b, aLog = np.polyfit(xsLog, ysLog, 1)
    a = np.exp(aLog)
    return a * np.power(xs, b)

#Get best-fit line of form y = a*x^order, but use log L2 fit
def get_order_line_log(xs, ys, order):
    xs = np.array(xs, dtype=np.float)
    ys = np.array(ys, dtype=np.float)
    mask = np.logical_and(xs > 0, ys > 0)
    xsLog = np.log(xs[mask])
    ysLog = np.log(ys[mask])
    coeffLog = np.mean(ysLog - order * xsLog)
    coeff = np.exp(coeffLog)
    return coeff * np.power(xs, order)

#Plot stairs
def get_stairs(x, y):
    """
    Given x and y, return points necessary for a step plot
    x is of size n+1 and y is of size n.

    :Parameters:
     - `x`: [x_0,...,x_n] list of x coordinates as numpy array
     - `y`: [y_0,...,y_{n-1}] list of y coordinates as numpy array
    :rtype: tuple of lists
    :return: (xnew,ynew) where
        xnew = (x_0,x_1,x_2,...,x_n,x_n)
        ynew = (y_0,y_0,y_1,...,y_{n-1})
    :requries: numpy
    :author: based on pytrix.py of Google Code
    """
    nobs = len(x)
    assert len(x) == (len(y)+1), "x must be one larger than y in length."
    xnew = np.repeat(x, 2)[1: -1]
    ynew = np.repeat(y, 2)
    return xnew, ynew
