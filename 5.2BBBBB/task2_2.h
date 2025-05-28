#pragma once
#include <vector>
#include <iostream>
#include <cmath>

std::vector<std::vector<double>> matrixPower(std::vector<std::vector<double>> matrix, int power);
void runTask2_2();

std::vector<double> powarray(const std::vector<double> array, int power, int size);
std::vector<double> inverseMatrix1D(const std::vector<double>& matrix, int size);
std::vector<double> multiplyMatrices1D(const std::vector<double>& A, const std::vector<double>& B, int size);
void printMatrix1(const std::vector<double>& matrix, int size);
void runTask2_21();