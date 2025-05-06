#include "task2_1.h"
#include <iostream>
#include <vector>

std::vector<int> findRepeatedNumbers(const std::vector<int>& numbers) {
    std::vector<int> repeatedNumbers;

    if (numbers.empty()) {
        return repeatedNumbers;
    }

    // Находим максимальное число
    int max_num = numbers[0];
    for (size_t i = 1; i < numbers.size(); ++i) {
        if (numbers[i] > max_num) {
            max_num = numbers[i];
        }
    }

    std::vector<int> count(max_num + 1, 0);

    // Первый проход: подсчёт количества
    for (size_t i = 0; i < numbers.size(); ++i) {
        count[numbers[i]]++;
    }

    // Второй проход: сбор повторяющихся чисел
    for (int num = 0; num <= max_num; ++num) {
        if (count[num] > 1) {
            repeatedNumbers.push_back(num);
        }
    }

    return repeatedNumbers;
}

void runTask2_1() {
    std::cout << "=== Задача 2.1: Поиск повторяющихся чисел ===\n";
    std::cout << "Введите числа через пробел (закончите ввод символом конца строки): ";

    std::vector<int> numbers;
    int num;
    while (std::cin >> num) {
        numbers.push_back(num);
    }

    std::cin.clear();
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    std::vector<int> repeated = findRepeatedNumbers(numbers);

    if (repeated.empty()) {
        std::cout << "Все числа уникальны, повторяющихся нет.\n";
    }
    else {
        std::cout << "Числа, встречающиеся более одного раза: ";
        for (size_t i = 0; i < repeated.size(); ++i) {
            std::cout << repeated[i] << " ";
        }
        std::cout << "\n";
    }
}