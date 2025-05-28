
#include "task2_2.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>

namespace {
    std::vector<std::vector<double>> multiplyMatrices(const std::vector<std::vector<double>>& A,
        const std::vector<std::vector<double>>& B) {
        const size_t n = A.size();
        const size_t m = B[0].size();
        const size_t p = B.size();

        std::vector<std::vector<double>> result(n, std::vector<double>(m, 0.0));

        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < m; ++j) {
                double sum = 0.0;
                for (size_t k = 0; k < p; ++k) {
                    sum += A[i][k] * B[k][j];
                }
                result[i][j] = sum;
            }
        }
        return result;
    }

    std::vector<std::vector<double>> inverseMatrix(const std::vector<std::vector<double>>& A) {
        size_t n = A.size();
        std::vector<std::vector<double>> aug(n, std::vector<double>(2 * n, 0.0));
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j)
                aug[i][j] = A[i][j];
            for (size_t j = n; j < 2 * n; ++j)
                aug[i][j] = (j - n == i) ? 1.0 : 0.0;
        }
        for (size_t i = 0; i < n; ++i) {
            double pivot = aug[i][i];
            if (std::fabs(pivot) < 1e-9) {
                bool swapped = false;
                for (size_t j = i + 1; j < n; ++j) {
                    if (std::fabs(aug[j][i]) > 1e-9) {
                        std::swap(aug[i], aug[j]);
                        pivot = aug[i][i];
                        swapped = true;
                        break;
                    }
                }
                if (!swapped) {
                    std::cerr << "Ошибка: матрица необратима!\n";
                    return std::vector<std::vector<double>>();
                }
            }
            for (size_t j = 0; j < 2 * n; ++j)
                aug[i][j] /= pivot;
            for (size_t k = 0; k < n; ++k) {
                if (k == i)
                    continue;
                double factor = aug[k][i];
                for (size_t j = 0; j < 2 * n; ++j)
                    aug[k][j] -= factor * aug[i][j];
            }
        }
        std::vector<std::vector<double>> inv(n, std::vector<double>(n, 0.0));
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j)
                inv[i][j] = aug[i][j + n];
        }
        return inv;
    }

    void printMatrix(const std::vector<std::vector<double>>& matrix) {
        for (size_t i = 0; i < matrix.size(); ++i) {
            for (size_t j = 0; j < matrix[i].size(); ++j) {
                std::cout << matrix[i][j] << " ";
            }
            std::cout << "\n";
        }
    }
}

std::vector<std::vector<double>> matrixPower(std::vector<std::vector<double>> matrix, int power) {
    if (power < 0) {
        std::vector<std::vector<double>> inv = inverseMatrix(matrix);
        if (inv.empty()) {
            return inv;
        }
        return matrixPower(inv, -power);
    }

    if (power == 0) {
        const size_t n = matrix.size();
        std::vector<std::vector<double>> identity(n, std::vector<double>(n, 0.0));
        for (size_t i = 0; i < n; ++i) {
            identity[i][i] = 1.0;
        }
        return identity;
    }

    std::vector<std::vector<double>> result = matrix;
    for (int i = 1; i < power; ++i) {
        result = multiplyMatrices(result, matrix);
    }
    return result;
}

void runTask2_2() {
    std::cout << "=== Задача 2.2: Возведение матрицы в степень (2D вектор) ===\n";
    size_t rows, cols;
    std::cout << "Введите размеры матрицы (строки и столбцы): ";
    std::cin >> rows >> cols;

    std::vector<std::vector<double>> matrix(rows, std::vector<double>(cols));

    std::cout << "Введите элементы матрицы построчно:\n";
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            std::cin >> matrix[i][j];
        }
    }

    int power;
    std::cout << "Введите степень, в которую нужно возвести матрицу: ";
    std::cin >> power;

    if (rows != cols && power != 1) {
        std::cout << "Ошибка: матрица не квадратная, возведение в степень (кроме 1) невозможно!\n";
        return;
    }

    std::vector<std::vector<double>> result = matrixPower(matrix, power);

    if (!result.empty()) {
        std::cout << "Результат возведения в степень " << power << ":\n";
        printMatrix(result);
    }
}

std::vector<double> multiplyMatrices1D(const std::vector<double>& A, const std::vector<double>& B, int size) {
    const int matrix_dim = size;
    const size_t total_elements = static_cast<size_t>(matrix_dim * matrix_dim);

    if (A.size() != total_elements || B.size() != total_elements) {
        std::cerr << "Ошибка: Неверный размер матриц для умножения (1D представление).\n";
        return {};
    }

    std::vector<double> result(total_elements, 0.0);

    for (int i = 0; i < matrix_dim; ++i) {
        for (int j = 0; j < matrix_dim; ++j) {
            for (int k = 0; k < matrix_dim; ++k) {
                result[i * matrix_dim + j] += A[i * matrix_dim + k] * B[k * matrix_dim + j];
            }
        }
    }
    return result;
}

std::vector<double> inverseMatrix1D(const std::vector<double>& A, int n) {
    const int matrix_dim = n;
    const size_t total_elements = static_cast<size_t>(matrix_dim * matrix_dim);
    const int augmented_row_size = matrix_dim * 2;

    if (A.size() != total_elements) {
        std::cerr << "Ошибка: Размер матрицы не соответствует 'n' для 1D представления.\n";
        return {};
    }

    std::vector<double> augmented(matrix_dim * augmented_row_size, 0.0);
    for (int i = 0; i < matrix_dim; i++) {
        for (int j = 0; j < matrix_dim; j++) {
            augmented[i * augmented_row_size + j] = A[i * matrix_dim + j];
        }
        augmented[i * augmented_row_size + (matrix_dim + i)] = 1.0;
    }

    for (int i = 0; i < matrix_dim; i++) {
        double pivot = augmented[i * augmented_row_size + i];
        if (std::fabs(pivot) < 1e-9) {
            bool swapped = false;
            for (int k = i + 1; k < matrix_dim; ++k) {
                if (std::fabs(augmented[k * augmented_row_size + i]) > 1e-9) {
                    for (int j = 0; j < augmented_row_size; ++j) {
                        std::swap(augmented[i * augmented_row_size + j], augmented[k * augmented_row_size + j]);
                    }
                    pivot = augmented[i * augmented_row_size + i];
                    swapped = true;
                    break;
                }
            }
            if (!swapped) {
                std::cerr << "Ошибка: матрица является необратимой (сингулярной)!\n";
                return {};
            }
        }

        for (int j = 0; j < augmented_row_size; j++) {
            augmented[i * augmented_row_size + j] /= pivot;
        }

        for (int k = 0; k < matrix_dim; k++) {
            if (k != i) {
                double factor = augmented[k * augmented_row_size + i];
                for (int j = 0; j < augmented_row_size; j++) {
                    augmented[k * augmented_row_size + j] -= factor * augmented[i * augmented_row_size + j];
                }
            }
        }
    }

    std::vector<double> inverse(total_elements);
    for (int i = 0; i < matrix_dim; i++) {
        for (int j = 0; j < matrix_dim; j++) {
            inverse[i * matrix_dim + j] = augmented[i * augmented_row_size + (j + matrix_dim)];
        }
    }
    return inverse;
}

std::vector<double> powarray(const std::vector<double> array, int power, int size) {
    const int matrix_dim = size;
    const size_t total_elements = static_cast<size_t>(matrix_dim * matrix_dim);

    if (array.size() != total_elements) {
        std::cerr << "Ошибка: Размер предоставленного массива не соответствует квадратной матрице " << matrix_dim << "x" << matrix_dim << ".\n";
        return {};
    }

    if (power == 0) {
        std::vector<double> identity(total_elements, 0.0);
        for (int i = 0; i < matrix_dim; ++i) {
            identity[i * matrix_dim + i] = 1.0;
        }
        return identity;
    }

    if (power < 0) {
        std::vector<double> inv = inverseMatrix1D(array, matrix_dim);
        if (inv.empty()) {
            return {};
        }
        return powarray(inv, -power, matrix_dim);
    }

    std::vector<double> result = array;
    for (int i = 1; i < power; ++i) {
        result = multiplyMatrices1D(result, array, matrix_dim);
    }
    return result;
}

void printMatrix1(const std::vector<double>& M, int n) {
    const int matrix_dim = n;
    const size_t total_elements = static_cast<size_t>(matrix_dim * matrix_dim);

    if (M.size() != total_elements) {
        std::cerr << "Ошибка: Размер массива не соответствует матрице " << matrix_dim << "x" << matrix_dim << " для вывода.\n";
        return;
    }
    for (int i = 0; i < matrix_dim; i++) {
        for (int j = 0; j < matrix_dim; j++) {
            std::cout << M[i * matrix_dim + j] << "\t";
        }
        std::cout << std::endl;
    }
}

void runTask2_21() {
    std::cout << "=== Задача 2.21: Возведение матрицы в степень (1D вектор) ===\n";
    int size_int;
    std::cout << "Введите размерность квадратной матрицы (n для n x n): ";
    std::cin >> size_int;

    if (size_int <= 0) {
        std::cerr << "Ошибка: Размерность матрицы должна быть положительным числом.\n";
        return;
    }
    const int matrix_dim = size_int;
    const size_t total_elements = static_cast<size_t>(matrix_dim * matrix_dim);

    std::vector<double> array(total_elements);

    std::cout << "Введите элементы матрицы построчно (всего " << total_elements << " элементов):\n";
    for (size_t i = 0; i < total_elements; ++i) {
        std::cin >> array[i];
    }

    std::cout << "Введенная матрица:\n";
    printMatrix1(array, matrix_dim);

    int power;
    std::cout << "Введите степень, в которую нужно возвести матрицу: ";
    std::cin >> power;

    std::vector<double> result = powarray(array, power, matrix_dim);

    if (!result.empty()) {
        std::cout << "Результат возведения в степень " << power << ":\n";
        printMatrix1(result, matrix_dim);
    }
    else {
        std::cout << "Операция возведения в степень не может быть выполнена.\n";
    }
}