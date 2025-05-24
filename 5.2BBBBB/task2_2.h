#pragma once
#include <vector>
#include <iostream>
#include <cmath>

std::vector<double> matrixPower(const std::vector<double>& matrix, int size, int power);
std::vector<double> inverseMatrix(const std::vector<double>& matrix, int size);
std::vector<double> multiplyMatrices(const std::vector<double>& A, const std::vector<double>& B, int size);
void printMatrix(const std::vector<double>& matrix, int size);
void runTask2_2();
void runTask2_21();