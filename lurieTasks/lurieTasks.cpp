#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>

#include <Windows.h>
#include <ctime>
#include <cmath>

#define PI 3.14159
#define g  9.81

using namespace std;

struct PipeModel {
    double p0, pl;
    double L, D, delta, d, thickness, eps;
    double z0, zl;
    double ro, visc;
    double Q;
    double speed;

    PipeModel()
    {
        // p0 = ???
        pl = 6e5;
        L = 8e4;
        D = 0.720;
        thickness = 0.010;
        delta = 15e-6;
        z0 = 50;
        zl = 100;
        ro = 870;
        visc = 15e-6;
        Q = 3500.0 / 3600.0;

        // Внутренний диаметр трубопровода
        d = D - 2 * thickness;
        cout << "Внутренний диаметр: " << d << endl;

        // Относительная жквивалетная шероховатость
        eps = delta / d;
        cout << "Относительная жквивалетная шероховатость: " << eps << endl;

        // Скорость потока
        speed = 4 * Q / (PI * pow(d, 2));
        cout << "Скорость потока: " << speed << endl;
    }
};

/// @brief Класс солвера для решения задач
class PipeSolver
{
public:
    PipeSolver(double& lym) 
        : lymbda{lym}
    {

    }

    void QP(PipeModel& pipe)
    {
        pipe.p0 = pipe.pl / (pipe.ro * g) + pipe.zl - pipe.z0 + lymbda * pipe.L / pipe.d * pow(pipe.speed, 2) / 2 / g;
        pipe.p0 *= pipe.ro * g;
    }

private:
    double lymbda;
};

int main()
{
    int time_count = clock();

    // Установка кодировки консоли
    setlocale(LC_ALL, "Russian");

    // Структура для трубопровода
    PipeModel pipeData;

    //________________________Начало решения___________________
    // Число Рейнольдса
    double Re = pipeData.speed * pipeData.d / pipeData.visc;
    cout << "Число Рейнольдса: " << Re << endl;
    
    // Формула Альтшуля
    double lymbda = 0.11 * pow((pipeData.eps + 68 / Re), 0.25);

    PipeSolver solv(lymbda);

    solv.QP(pipeData);
    cout << "Давление в начале трубопровода: " << pipeData.p0 << endl;

    // Вывод затраченного времени
    printf("Затраченное время: %i ms\n", time_count);

    return 0;
}
