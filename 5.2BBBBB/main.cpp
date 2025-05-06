#include "task2_1.h"
#include "task2_2.h"
#include <iostream>

int main() {
    setlocale(LC_ALL, "ru");
    int choice;
    do {
        std::cout << "\n�������� ������:\n";
        std::cout << "1. ������ 2.1 (������������� �����)\n";
        std::cout << "2. ������ 2.2 (������� � ������� n)\n";
        std::cout << "0. �����\n";
        std::cout << "������� ����� ������: ";
        std::cin >> choice;

        switch (choice) {
        case 1:
            runTask2_1();
            break;
        case 2:
            runTask2_2();
            break;
        case 0:
            std::cout << "����� �� ���������.\n";
            break;
        default:
            std::cout << "������: �������� ����� ������!\n";
        }
    } while (choice != 0);

    return 0;
}