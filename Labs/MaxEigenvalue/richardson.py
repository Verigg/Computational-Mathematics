from numpy.linalg import norm
import numpy as np


def richardson(a, b, x0):
    k = 0
    omega = 2 / (min(np.linalg.eigvals(a)) + max(np.linalg.eigvals(a)))
    epsilon = 0.001
    while k <= 500:
        x_next = x0 - omega * (np.dot(a, x0) - b)
        if norm(x0 - x_next) < epsilon:
            return x0, k
        x0 = x_next
        k += 1
    return x0, k


def main():
    a = np.array([[2.0, 1.0, 0.0],
                  [1.0, 2.0, 1.0],
                  [0.0, 1.0, 3.0]])

    b = np.array([2.0, 1.0, 4.0])

    x0 = np.zeros(3)
    x, k = richardson(a, b, x0)

    print("Calculated Solution vector:", x)
    print("True Solution vector:", np.linalg.solve(a,b))
    print("Obtained in", k, "steps.")

    a = np.array([[6.22, 1.42, -1.72, 1.91],
                  [1.42, 5.33, 1.11, -1.82],
                  [-1.72, 1.11, 5.24, 1.42],
                  [1.91, -1.82, 1.42,  6.55]])
    b = np.array([7.53, 6.06, 8.05, 8.06])

    x0 = np.zeros(4)
    x, k = richardson(a, b, x0)

    print("Solution vector:", x)
    print("True Solution vector:", np.linalg.solve(a,b))
    print("Obtained in", k, "steps.")

if __name__ == "__main__":
    main()


