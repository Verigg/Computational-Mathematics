import numpy as np
from tabulate import tabulate

def givens_rot(n, ij, theta):
    rot_mat = np.identity(n)
    for i in range(n):
        for j in range(n):
            if (i == ij[0] and j == ij[0]) or (i == ij[1] and j == ij[1]):
                rot_mat[i][j] = np.cos(theta)
            elif i == ij[0] and j == ij[1]:
                rot_mat[i][j] = np.sin(theta)
            elif i == ij[1] and j == ij[0]:
                rot_mat[i][j] = -np.sin(theta)
            else:
                rot_mat[i][j] = rot_mat[i][j]
    return rot_mat


def eigenvalues(a_mat, epsilon):
    def eigenvalues_recur(k, d_mat):
        nonlocal a_mat
        i, j = np.unravel_index(np.argmax(np.abs(d_mat - np.diag(np.diag(d_mat)))), (n, n))
        dij = d_mat[i, j]
        dii = d_mat[i, i]
        djj = d_mat[j, j]

        if np.abs(dij) < epsilon or k > 1000000:
            return k, np.diag(d_mat)

        theta = np.arctan(2 * dij / (dii - djj)) / 2 if dii != djj else np.sign(dij) * np.pi / 4 # Находим угол поворота
        s_mat = givens_rot(n, (i, j), theta) # Составим матрицу вращения
        return eigenvalues_recur(k + 1, s_mat @ d_mat @ s_mat.T)

    n = a_mat.shape[0]
    return eigenvalues_recur(0, np.array(a_mat))


def generate_symmetric_matrix(rng, dim):
    temp = np.random.uniform(low=rng[0], high=rng[1], size=(dim, dim))
    return temp @ temp.T


def main():
    salt = 421
    np.random.seed(salt)
    rng = (-10, 10)

    epsilon = [1e-3, 1e-7]


    for size in (3,5,7,100):
        a = generate_symmetric_matrix(rng, size)
        k, eigenvals = eigenvalues(a, epsilon[0])
        real_eigenvals = np.linalg.eig(a)[0]
        print("Matrix:")
        print(tabulate(a, tablefmt="grid"))
        print("Epsilon = ", epsilon[0])
        print("Done in", k, "steps")
        print("Calculated eigenvalues", sorted(eigenvals, reverse=True))
        print("Real eigenvalues", sorted(real_eigenvals, reverse=True))

        k, eigenvals = eigenvalues(a, epsilon[1])
        real_eigenvals = np.linalg.eig(a)[0]
        print("Epsilon = ", epsilon[1])
        print("Done in", k, "steps")
        print("Calculated eigenvalues", sorted(eigenvals, reverse=True))
        print("Real eigenvalues", sorted(real_eigenvals, reverse=True))
        print("\n")



if __name__ == "__main__":
    main()