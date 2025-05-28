#include <iostream>
#include <clocale> 
#include "task2_2.h" // Для задач 2.2 и 2.21
#include "task2_1.h" // Для задачи 2.1

int main() {
    setlocale(LC_ALL, "ru");
    int choice;
    do {
        std::cout << "\nВыберите задачу:\n";
        std::cout << "1. Задача 2.1 (поиск повторяющихся чисел)\n"; 
        std::cout << "2. Задача 2.2 (матрица в степени n, 2D вектор)\n";
        std::cout << "3. Задача 2.21 (матрица в степени n, 1D вектор)\n";
        std::cout << "0. Выход\n";
        std::cout << "Введите номер задачи: ";
        std::cin >> choice;

        switch (choice) {
        case 1:
            runTask2_1();
            break;
        case 2:
            runTask2_2();
            break;
        case 3:
            runTask2_21();
            break;
        case 0:
            std::cout << "Выход из программы.\n";
            break;
        default:
            std::cout << "Ошибка: неверный номер задачи!\n";
        }
    } while (choice != 0);

    return 0;
}