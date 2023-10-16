import numpy as np

# Лабораторная №5 метод монотонной прогонки

def is_valid(matrix):
    if matrix.shape[0] != matrix.shape[1] - 1: # Проверяем является ли матрица квадратной
        return False

    for i, row in enumerate(matrix):
        if i == 0:  # выбираем индексы по которым должны быть не нулевые значения в стоке
            js = [i, i + 1]
        elif i == matrix.shape[0] - 1:
            js = [i - 1, i]
        else:
            js = [i - 1, i, i + 1]

        for j, x_ij in enumerate(row):
            if j < matrix.shape[0] and j not in js and x_ij != 0: # Проверяем на отсутствие нулей в нужных индексах
                return False

    return True

class coeff:
    def __init__(self, A, B, C, F):
        self.A = A
        self.B = B
        self.C = C
        self.F = F

def get_coeffs(matrix):
    coeffs = np.zeros((matrix.shape[0], 4))

    for i, row in enumerate(matrix): # собираем коэффициенты учитывая первую и последнюю строку
        coeffs[i, 1] = row[i]  # b
        if i > 0:
            coeffs[i, 0] = row[i - 1]  # a
        if i < matrix.shape[0] - 1:
            coeffs[i, 2] = row[i + 1]  # c
        coeffs[i, 3] = row[matrix.shape[1] - 1]  # d

    A = coeffs[:, 0]
    B = coeffs[:, 1]
    C = coeffs[:, 2]
    F = coeffs[:, 3]

    return coeff(A, B, C, F)


def monotonic_method(matrix):
    assert is_valid(matrix)

    coeffs = get_coeffs(matrix)
    n = matrix.shape[0]

    A_ = np.zeros(n)
    B_ = np.zeros(n)

    A_[0] = coeffs.F[0] / coeffs.B[0]
    B_[0] = -coeffs.C[0] / coeffs.B[0]
    for i in range(1, n):
        A_[i] = (coeffs.F[i] - coeffs.A[i] * A_[i - 1]) / (coeffs.A[i] * B_[i - 1] + coeffs.B[i])
        B_[i] = -coeffs.C[i] / (coeffs.A[i] * B_[i - 1] + coeffs.B[i])

    X = np.zeros(n)

    X[n - 1] = A_[n - 1]
    for i in range(n - 2, -1, -1):
        X[i] = A_[i] + B_[i] * X[i + 1]

    return X

# проверка
matrix = np.array([[5,2,0,0,10],
                    [3,9,4,0,21],
                    [0,6,8,7,32],
                    [0,0,3,5,15]])

x_array =  monotonic_method(matrix)
print(matrix)
print(f"x_array:{x_array}")
# проверка ответов
print("5x₁ + 2x₂ = ", x_array[0] * 5 + x_array[1] * 2)
print("3x₁ + 9x₂ + 4x₃ = ", x_array[0] * 3 + x_array[1] * 9 + x_array[2] * 4)
print("6x₂ + 8x₃ + 7x₄ = ", x_array[1] * 6 + x_array[2] * 8 + x_array[3] * 7)
print("3x₃ + 5x₄ = ", x_array[2] * 3 + x_array[3] * 5)