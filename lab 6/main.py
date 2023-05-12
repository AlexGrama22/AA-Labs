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


# Algorithm 2: Leibniz Formula
def leibniz_formula(n):
    pi = 0
    for i in range(n):
        pi += ((-1) ** i) / (2 * i + 1)
    return 4 * pi


# Algorithm 3: Bailey–Borwein–Plouffe (BBP) Formula
def bbp_formula(n):
    mp.dps = n + 1
    pi = 0
    for k in range(n):
        pi += (1 / (16 ** k)) * ((4 / (8 * k + 1)) - (2 / (8 * k + 4)) - (1 / (8 * k + 5)) - (1 / (8 * k + 6)))
    return pi


def analyze_algorithms():
    n = int(input("Enter the number of decimal places (n): "))
    algorithms = [chudnovsky_algorithm, leibniz_formula, bbp_formula]
    timings = []
    for algorithm in algorithms:
        start_time = time.time()
        pi = algorithm(n)
        end_time = time.time()
        timings.append(end_time - start_time)
        print(f"{algorithm.__name__}: {pi}, time: {end_time - start_time:.6f} seconds")

    plt.bar(np.arange(len(algorithms)), timings, align='center', alpha=0.5)
    plt.xticks(np.arange(len(algorithms)), [algorithm.__name__ for algorithm in algorithms])
    plt.ylabel('Execution Time (seconds)')
    plt.title('Performance Comparison of Pi Calculation Algorithms')
    plt.show()


if __name__ == '__main__':s
    analyze_algorithms()
