#include <iostream>
#include <clocale> 
#include "task2_2.h" // ��� ����� 2.2 � 2.21
#include "task2_1.h" // ��� ������ 2.1

int main() {
    setlocale(LC_ALL, "ru");
    int choice;
    do {
        std::cout << "\n�������� ������:\n";
        std::cout << "1. ������ 2.1 (����� ������������� �����)\n"; 
        std::cout << "2. ������ 2.2 (������� � ������� n, 2D ������)\n";
        std::cout << "3. ������ 2.21 (������� � ������� n, 1D ������)\n";
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
        case 3:
            runTask2_21();
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