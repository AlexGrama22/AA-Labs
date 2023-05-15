import time
import matplotlib.pyplot as plt
import numpy as np
from mpmath import mp


# Algorithm 1: Chudnovsky Algorithm
def chudnovsky_algorithm(n):
    mp.dps = n + 1
    C = 426880 * mp.sqrt(10005)
    M = 1
    L = 13591409
    X = 1
    K = 6
    S = L
    for _ in range(1, n):
        M = (K ** 3 - 16 * K) * M // K ** 3
        L += 545140134
        X *= -262537412640768000
        S += (M * L) // X
        K += 12
    return C / S


# Algorithm 2: Nilakantha Series
def nilakantha_series(n):
    pi = mp.mpf(3)
    sign = mp.mpf(-1)
    term = mp.mpf(2)
    for i in range(1, n):
        pi += sign * (4 / (term * (term + 1) * (term + 2)))
        sign *= -1
        term += 2
    return pi


# Algorithm 3: Bailey–Borwein–Plouffe (BBP) Formula
def bbp_formula(n):
    mp.dps = n + 1
    pi = 0
    for k in range(n):
        pi += (1 / (16 ** k)) * ((4 / (8 * k + 1)) - (2 / (8 * k + 4)) - (1 / (8 * k + 5)) - (1 / (8 * k + 6)))
    return pi


def analyze_algorithms():
    decimal_places = [100, 1000, 10000, 100000]
    algorithms = [chudnovsky_algorithm, nilakantha_series, bbp_formula]
    timings = [[] for _ in range(len(algorithms))]  # List of lists to store timings for each algorithm

    for n in decimal_places:
        for i, algorithm in enumerate(algorithms):
            start_time = time.time()
            pi = algorithm(n)
            end_time = time.time()
            timings[i].append(end_time - start_time)
            print(f"{algorithm.__name__} ({n} decimal places): {pi}, time: {end_time - start_time:.6f} seconds")

    # Plotting the results on the same graph
    plt.figure()
    for i, algorithm in enumerate(algorithms):
        plt.plot(decimal_places, timings[i], marker='o', label=algorithm.__name__)

    plt.xlabel('Decimal Places')
    plt.ylabel('Execution Time (seconds)')
    plt.title('Performance Comparison of Pi Calculation Algorithms')
    plt.legend()
    plt.show()


if __name__ == '__main__':
    analyze_algorithms()
