"""
Original Paper's Simulation
Taking Normal Probability Distribution

Exactly same variables are defined
"""
from math import e
from random import uniform
from cvxopt import solvers

Omega = range(30)
s1_low = 0
s1_high = 8
s2_low = 10
s2_high = 50
delta = 0.9
PI = []


def c(w):
    """
    A simple approximation to cost of intrusion
    We take intrusion as highly undesirable, that's why exponential
    :param w: State.
    :return:
    """
    alpha = 10
    beta = 0.05
    return alpha * (e ^ (w * beta) - 1)


def xi1(w, s1, s2):
    """

    :param w: Current State
    :param s1: Intruder's strategy
    :param s2: IDS's strategy
    :return: Cost for above variables
    """
    c1_high = 8
    c1_low = 0

    if s1 == s1_high:
        return -c(w) + c1_high
    else:
        return -c(w) + c1_low


def xi2(w, s1, s2):
    """

    :param w: Current State
    :param s1: Intruder's strategy
    :param s2: IDS's strategy
    :return: Cost for above variables
    """
    c2_high = 10
    c2_low = 2

    if s2 == s2_high:
        return -c(w) + c2_high
    else:
        return -c(w) + c2_low


def p(s1, s2):
    """
    PD - IDS detection rate
    Taking Constant for now, average(pd) = 0.8

    What's T?
    Take data transmission speed = 250 KBps
    Take malicious data packet size = 5 Kb
    So time taken = T = 5/250

    :param s1: Intruder's strategy
    :param s2: IDS's strategy
    :return: Probability of Defense
    """
    p_low = 0.8
    p_high = 0.8
    t = 5 / 250
    pd = uniform(p_low, p_high)

    return min([pd * t * s1 / s2, 1])


def _pi(_w, w, _p):
    W = len(Omega)
    if _w == w - 1:
        return w / W * _p
    elif _w == w:
        return w / W * (1 - _p) + (1 - w / W) * _p
    elif _w == w + 1:
        return (1 - w / W) * _p
    else:
        return 0


def main():
    solvers

if __name__ == '__main__':
    main()
