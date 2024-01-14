import numpy as np
from tabulate import tabulate


def eigenvalues(a, epsilon, step=0, v_cur=0):
    if step == 0:
        v_cur=np.ones((a.shape[0], 1)) # Начальное приближение
    v_next = np.matmul(a, v_cur) / np.linalg.norm(np.matmul(a, v_cur)) # Находим значение вектора и осуществляем нормировку
    if np.linalg.norm(v_next - v_cur) < epsilon or step > 1e6: # Проверяем условие |𝜆𝑘+1 − 𝜆𝑘| ≤ ε.
        return (np.matmul(np.matmul(v_cur.T, a), v_cur) / np.matmul(v_cur.transpose(), v_cur)).item(), step
    else:
        return eigenvalues(a, epsilon,  step + 1, v_next)


def generate_symmetric_matrix(rng, dim):
    temp = np.random.uniform(low=rng[0], high=rng[1], size=(dim, dim))
    return temp @ temp.T


def main():
    salt = 4221
    np.random.seed(salt)
    rng = (0, 10)

    epsilon = 1e-4

    for size in (3,5,7):
        a = generate_symmetric_matrix(rng, size)

        print("Matrix:")
        print(tabulate(a, tablefmt="grid"))
        print("Max. eigenvalue")
        print(eigenvalues(a, epsilon)[0])
        print("Steps")
        print(eigenvalues(a, epsilon)[1])
        print("All eigenvalues")
        print(np.linalg.eig(a)[0])
        print()


if __name__ == "__main__":
    main()