
#include "task2_2.h"
#include <iostream>
#include <vector>
#include <cmath>  // ��� std::fabs

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

    // ����� ������� 
    std::vector<std::vector<double>> inverseMatrix(const std::vector<std::vector<double>>& A) {
        size_t n = A.size();
        // ������ ����������� ������� [A | I]
        std::vector<std::vector<double>> aug(n, std::vector<double>(2 * n, 0.0));
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j)
                aug[i][j] = A[i][j];
            for (size_t j = n; j < 2 * n; ++j)
                aug[i][j] = (j - n == i) ? 1.0 : 0.0;
        }
        // ������ ��� ��������������
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
                    std::cerr << "������: ������� ����������!\n";
                    return std::vector<std::vector<double>>();
                }
            }
            // ������������ ������
            for (size_t j = 0; j < 2 * n; ++j)
                aug[i][j] /= pivot;
            // ��������� ������� i � ��������� �������
            for (size_t k = 0; k < n; ++k) {
                if (k == i)
                    continue;
                double factor = aug[k][i];
                for (size_t j = 0; j < 2 * n; ++j)
                    aug[k][j] -= factor * aug[i][j];
            }
        }
        // ���������� �������� ������� �� �����������
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
        // ��� ������������� ������� ������� ��������� �������� �������
        std::vector<std::vector<double>> inv = inverseMatrix(matrix);
        if (inv.empty()) {
            // ���� ������� ����������, ������� inverseMatrix ��� ������ ��������� �� ������
            return inv;
        }
        // �������� �������� ������� � ������������� ������� (-power)
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
    std::cout << "=== ������ 2.2: ���������� ������� � ������� ===\n";
    size_t
        rows, cols;
    std::cout << "������� ������� ������� (������ � �������): ";
    std::cin >> rows >> cols;

    std::vector<std::vector<double>> matrix(rows, std::vector<double>(cols));

    std::cout << "������� �������� ������� ���������:\n";
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            std::cin >> matrix[i][j];
        }
    }

    int power;
    std::cout << "������� �������, � ������� ����� �������� �������: ";
    std::cin >> power;

    if (rows != cols && power != 1) {
        std::cout << "������: ������� �� ����������, ���������� � ������� ����������!\n";
        return;
    }

    std::vector<std::vector<double>> result = matrixPower(matrix, power);

    if (!result.empty()) {
        std::cout << "��������� ���������� � ������� " << power << ":\n";
        printMatrix(result);
    }
}

std::vector<double> powerplus(const std::vector<double> array, int power, int size)
{
    std::vector<double> res(pow(size, 2));

    if (power == 1)
        return array;
    for (int pwr = 0;pwr < power - 1;pwr++) {
        for (int i = 0; i < pow(size, 2); i++) {
            for (int j = 0;j < size;j++) {
                res[i] += array[i - i % size + j] * array[i % size + j * size];
            }
        }
    }
    return res;
}

std::vector<double> inverseMatrix(const std::vector<double>& A, int n) {
    std::vector<double> augmented(n * 2 * n, 0);

    // ��������� ����������� ������� (A | I)
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            augmented[i * 2 * n + j] = A[i * n + j];
        }
        augmented[i * 2 * n + (n + i)] = 1;
    }

    // ������ ��� ������ ������-�������
    for (int i = 0; i < n; i++) {
        if (augmented[i * 2 * n + i] == 0) {
            std::cerr << "������� �������� �����������!" << std::endl;
            exit(1);
        }

        double diagElement = augmented[i * 2 * n + i];
        for (int j = 0; j < 2 * n; j++) {
            augmented[i * 2 * n + j] /= diagElement;
        }

        for (int k = 0; k < n; k++) {
            if (k != i) {
                double factor = augmented[k * 2 * n + i];
                for (int j = 0; j < 2 * n; j++) {
                    augmented[k * 2 * n + j] -= factor * augmented[i * 2 * n + j];
                }
            }
        }
    }

    // ��������� �������� �������
    std::vector<double> inverse(n * n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            inverse[i * n + j] = augmented[i * 2 * n + (j + n)];
        }
    }
    return inverse;
}

std::vector<double> powarray(const std::vector<double> array, int power, int size) {
    std::vector<double> res(array.size());
    if (power == 0) {
        for (int i = 0;i < array.size();i++) {
            res[i] = 1.0;
        }
        return res;
    }

    if (power > 0) {
        for (int i = 0; i < array.size(); i++) {
            res[i] = 0;
        }
        return powerplus(array, power, size);
    }

    res = inverseMatrix(array, size);
    return powerplus(res, -1 * power, size);

}



void printMatrix1(const std::vector<double>& M, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            std::cout << M[i * n + j] << "\t";
        }
        std::cout << std::endl;
    }
}

void runTask2_21() {
    std::cout << "=== ������ 2.2: ���������� ������� � ������� ===\n";
    size_t size;
    std::cout << "������� ����������� �������: ";
    std::cin >> size;

    std::vector<double> array(pow(size, 2));

    std::cout << "������� �������� ������� ���������:\n";
    for (size_t i = 0; i < pow(size, 2); ++i) {
        std::cin >> array[i];
    }

    for (int i = 0;i < size;++i) {
        for (int j = 0;j < size;++j) {
            std::cout << array[i * size + j] << ' ';
        }
        std::cout << std::endl;
    }

    int power;
    std::cout << "������� �������, � ������� ����� �������� �������: ";
    std::cin >> power;

    /*if (rows != cols && power != 1) {                     � ����� ���� ����� ��� ������� ������� size * size
std::cout << "������: ������� �� ����������, ���������� � ������� ����������!\n";
        return;
    }*/

    std::vector<double> result = powarray(array, power, size);

    if (!result.empty()) {
        std::cout << "��������� ���������� � ������� " << power << ":\n";
        printMatrix1(result, size);
    }
}