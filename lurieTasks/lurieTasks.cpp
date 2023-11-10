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
        p0 = 0;
        pl = 6e+5;
        L = 8e+4;
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

        // Относительная эквивалетная шероховатость
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
    /*PipeSolver(double& lym) 
        : lymbda{lym}
    {

    }*/

    double find_Re(PipeModel& pipe)
    {
        return pipe.speed * pipe.d / pipe.visc;
    }

    double find_lymbda(float Re, float eps) {
        if (Re <= 2300)
            return 64 / Re;
        else if (2300 < Re < 10e+3)
        {
            double gamma = 1 - exp(-0.002 * (Re - 2300));
            return 64 / Re * (1 - gamma) + 0.3164 / pow(Re, 0.25) * gamma;
        }
        else
        {
            if ((10e+3 < Re < 27 / pow(eps, 1.143)) && (Re < 10e+4))
                return 0.3164 / pow(Re, 1 / 4);
            else if (Re > (500 / eps))
                return 0.11 * pow(eps, 0.25);
            else
                return 0.11 * pow(eps + 68 / Re, 0.25);
        }
    }

    void QP(PipeModel& pipe)
    {
        Re = find_Re(pipe);
        lymbda = find_lymbda(Re, pipe.eps);
        pipe.p0 = pipe.pl / (pipe.ro * g) + pipe.zl - pipe.z0 + lymbda * pipe.L / pipe.d * pow(pipe.speed, 2) / (2 * g) * (pipe.ro * g);
    }

    void PP(PipeModel& pipe)
    {
        double lym_b, Re;
        int itr_stop = 0;
        double a = 2 * g * pipe.d / pipe.L * ((pipe.p0 - pipe.pl) / (pipe.ro * g) + pipe.z0 - pipe.zl);
       
        pipe.p0 = 5e+6;
        pipe.pl = 8e+5;
        lymbda = 0.02;
        do
        {
            itr_stop++;
            if (itr_stop > 1000)
            {
                cout << "Расчёт неудался!" << endl;
                break;
            }
                
            lym_b = lymbda;
            pipe.speed = sqrt(a / lymbda);
            Re = pipe.speed * pipe.d / pipe.visc;
            //lymbda = 0.11 * pow(pipe.eps + 68 / Re, 0.25);
            find_lymbda(Re, pipe.eps);

        } while (abs(lymbda - lym_b) > 0.0005);

        cout << "\n\n\nКоличество итераций: " << itr_stop << endl;
        cout << "Lymbda: " << lymbda << endl;

        pipe.Q = PI * pow(pipe.d, 2) * pipe.speed / 4;
    }

private:
    double lymbda, Re;
};

int main()
{
    int time_count = clock();

    // Установка кодировки консоли
    setlocale(LC_ALL, "Russian");

    // Структура для трубопровода
    PipeModel pipeData;

    //________________________Начало решения___________________
    //// Число Рейнольдса
    //double Re = pipeData.speed * pipeData.d / pipeData.visc;
    //cout << "Число Рейнольдса: " << Re << endl;
    //
    //// Формула Альтшуля
    //double lymbda = 0.11 * pow((pipeData.eps + 68 / Re), 0.25);

    PipeSolver solv;

    solv.QP(pipeData);
    cout << "Давление в начале трубопровода: " << pipeData.p0 << endl;

    solv.PP(pipeData);
    cout << "Решение задачи PP: \nQ = " << pipeData.Q << " м3 / с\n" << endl;
    cout << "Решение задачи PP: \nQ = " << pipeData.Q * 3600 << " м3 / ч\n" << endl;

    // Вывод затраченного времени
    printf("Затраченное время: %i ms\n", time_count);

    return 0;
}
