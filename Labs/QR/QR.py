import numpy as np


def inner_product(vec1, vec2):
    return np.dot(vec1.T, vec2)


def project_on(vec, u):
    inner_prod = inner_product(u, u)
    return inner_product(u, vec) * u / inner_prod


def qr_decompose(matrix):
    assert (
        np.linalg.det(matrix) != 0.0
    )

    N = matrix.shape[0]

    q = np.zeros_like(matrix)

    for i in range(N):
        q[:, i] = matrix[:, i]

        if i > 0:
            projection = np.zeros_like(matrix[:, i])
            for j in range(i):
                projection += project_on(matrix[:, i],q[:, j])
            q[:, i] -= projection

    q = np.divide(q, np.linalg.norm(q, axis=0))

    r = np.zeros((N, N))

    for i in range(N):
        for j in range(N):
            if i <= j:
                r[i, j] = inner_product(q[:, i], matrix[:, j])

    return q, r

def solve_linear_system(A, b): # решаем системы с помошью QR разложения
    Q, R = qr_decompose(A) # Находим Q,R
    y = np.dot(Q.T, b.reshape(-1, 1)) # находим Y, y=Q^t*b
    x = np.linalg.solve(R, y) # находим X
    return x.flatten()


def test_accuracy(result, expected):
    return np.allclose(result, expected, atol=1e-4)


def do_tests():
    A1 = np.array([[2.1, -4.5, -2.0], [3.0, 2.5, 4.3], [-6.0, 3.5, 2.5]])
    b1 = np.array([19.07, 3.21, -18.25])
    res1 = np.array([1.34025, -4.75798, 2.5771])
    test_accuracy(solve_linear_system(A1, b1), res1)

    A2 = np.array([[5.0, -1.0, 5.0], [-3.0, 6.0, 2.0], [10.0, -7.0, 0.0]])
    b2 = np.array([3.2, 5.4, -1.2])
    res2 = np.array([0.7297, 1.2138, 0.1531])
    test_accuracy(solve_linear_system(A2, b2), res2)

    A3 = np.array([[5.0, 2.0, 3.0], [1.0, 6.0, 1.0], [3.0, -4.0, -2.0]])
    b3 = np.array([3.0, 5.0, 8.0])
    res3 = np.array([2.0, 1.0, -3.0])
    test_accuracy(solve_linear_system(A3, b3), res3)

    A4 = np.array([[1.0, 2.0, 1.0, 4.0], [2.0, 0.0, 4.0, 3.0], [4.0, 2.0, 2.0, 1.0], [-3.0, 1.0, 3.0, 2.0]])
    b4 = np.array([13.0, 28.0, 20.0, 6.0])
    res4 = np.array([3.0, -1.0, 4.0, 2.0])
    test_accuracy(solve_linear_system(A4, b4), res4)

    A5 = np.array([[2.0, 1.0, 3.0], [11.0, 7.0, 5.0], [9.0, 8.0, 4.0]])
    b5 = np.array([1.0, -6.0, -5.0])
    res5 = np.array([-1.0, 0.0, 1.0])
    test_accuracy(solve_linear_system(A5, b5), res5)


do_tests()


def solve_systems():
    A1 = np.array([[1., 2., 3.], [4., 6., 7.], [8., 9., 0.]])
    b1 = np.array([6., 12., 24.])
    res1 = np.array([-11.538, 12.923, -2.769])
    x1 = solve_linear_system(A1, b1)
    test_accuracy(x1, res1)
    print("A =")
    print(A1)
    print("b =", b1)
    print("x =", x1)
    print("x* =", res1)
    print("x*-x", res1-x1, "\n")

    A2 = np.array([[6.03, 13., -17.], [13., 29.03, -38.], [-17., -38., 50.03]])
    b2 = np.array([2.0909, 4.1509, -5.1191])
    res2 = np.array([1.03, 1.03, 1.03])
    x2 = solve_linear_system(A2, b2)
    test_accuracy(x2, res2)
    print("A =")
    print(A2)
    print("b =", b2)
    print("x =", x2)
    print("x*-x", res2 - x2, "\n")

    A3 = np.array([[2., 0., 1.], [0., 1., -1.], [1., 1., 1.]])
    b3 = np.array([3., 0., 3.])
    x3 = solve_linear_system(A1, b1)
    print("A =")
    print(A3)
    print("b =", b3)
    print("x =", x3)


def main():
    solve_systems()


if __name__ == "__main__":
    main()