#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>

#include <iomanip>
#include <fixed/fixed.h>
#include <pde_solvers/pde_solvers.h>

#include <Windows.h>
#include <ctime>
#include <cmath>

#define PI 3.14
#define g  9.81

using namespace std;

/// @brief Вывод в файл
/// @param press 
/// @param  
void writeFun(vector<double>& press, double& dx) {
    ofstream my_file;
    int profCount = press.size();
    my_file.open("res.csv");
    my_file << "time,x,Pressure" << endl;

    for (int i = 0; i < profCount; i++)
    {
        my_file << 0 << "," << i * dx;
        my_file << "," << press[i] << endl;
    }

    my_file.close();
}

/// @brief Структура для параметров трубопровода
struct PipeModel {
    double p0, pl;
    double L, D, delta, d, thickness, eps;
    double z0, zl;
    double ro, visc;
    double Q;
    double speed;

    /// @brief Инициализация значения параметров
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

        d = D - 2 * thickness; // Внутренний диаметр трубопровода

        eps = delta / d; // Относительная эквивалетная шероховатость

        speed = 4 * Q / (PI * pow(d, 2)); // Скорость потока
    }
};

class sample_system : public fixed_system_t<1>
{
    using fixed_system_t<1>::var_type;
    PipeModel& pipe;
    double Re, lymbda;

public:
    sample_system(PipeModel& pipeM)
        : pipe{ pipeM }
    {}
    // Задание функции невязок
    var_type residuals(const var_type& x) {

        Re = x * pipe.d / pipe.visc;
        lymbda = hydraulic_resistance_isaev(Re, pipe.eps);
        return
        {
            lymbda * (pipe.L * pow(x, 2)) / (2 * pipe.d * g) + (pipe.pl / (pipe.ro * g) + pipe.zl) - (pipe.p0 / (pipe.ro * g) + pipe.z0)
        };
    }
};

/// @brief Класс солвера для решения задач
class PipeSolver
{
public:
    /// @brief Расчёт числа Рейнольдса
    /// @param pipe Структура с параметрами трубопровода
    /// @return Возвращает число Рейнольдса
    double find_Re(PipeModel& pipe)
    {
        return pipe.speed * pipe.d / pipe.visc;
    }

    double diff(PipeModel& pipe, double dx, double dz)
    {
        double tau = lymbda / 8 * pipe.ro * pow(pipe.speed, 2);
        return -4 / pipe.d * tau - pipe.ro * g * dz / dx;
    }

    void euler(PipeModel& pipe, vector<double>& p, double h)
    {
        lymbda = hydraulic_resistance_isaev(find_Re(pipe), pipe.eps);
        double dz = (pipe.z0 - pipe.zl) / (p.size() - 1);
        
        for (int i = p.size() - 2; i >= 0; i--)
        {
            p[i] = p[i + 1] - h * diff(pipe, h, dz);
        }
    }

    /// @brief Решение первой задачи
    /// @param pipe Структура с параметрами трубопровода
    void QP(PipeModel& pipe)
    {
        Re = find_Re(pipe); // Расчёт числа Рейнольдса
        lymbda = hydraulic_resistance_isaev(Re, pipe.eps); // Расчёт коэффициента лямбда
        pipe.p0 = (pipe.pl / (pipe.ro * g) + pipe.z0 - pipe.zl + lymbda * pipe.L / pipe.d * pow(pipe.speed, 2) / (2 * g)) * (pipe.ro * g);
        PP_Newton(pipe);
    }

    /// @brief Решение второй задачи
    /// @param pipe Структура с параметрами трубопровода
    void PP(PipeModel& pipe)
    {
        double lym_b;
        int itr_stop = 0; // Ограничитель итераций

        pipe.p0 = 5e+6;
        pipe.pl = 8e+5;

        double a = 2 * g * pipe.d / pipe.L * ((pipe.p0 - pipe.pl) / (pipe.ro * g) + pipe.z0 - pipe.zl); // Правая часть уравнения
       
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
            lymbda = hydraulic_resistance_isaev(Re, pipe.eps);

        } while (abs(lymbda - lym_b) > 0.0005);

        /*cout << "\n\n\nКоличество итераций: " << itr_stop << endl;
        cout << "Lymbda: " << lymbda << endl;*/

        pipe.Q = PI * pow(pipe.d, 2) * pipe.speed / 4;
    }

    void PP_Newton(PipeModel& pipe)
    {
        /*pipe.p0 = 5e+6;
        pipe.pl = 8e+5;*/

        sample_system test(pipe);
        // Задание настроек решателя по умолчанию
        fixed_solver_parameters_t<1, 0> parameters;
        // Создание структуры для записи результатов расчета
        fixed_solver_result_t<1> result;
        // Решение системы нелинейныйх уравнений <2> с помощью решателя Ньютона - Рафсона
        // { 0, 0 } - Начальное приближение
        fixed_newton_raphson<1>::solve_dense(test, { 1 }, parameters, &result);

        cout << "\nРешение Ньютона: " << result.argument << endl;
        cout << "Решение Ньютона в м3/с: " << result.argument * PI * pow(pipe.d, 2) / 4 << endl;
        cout << "Решение Ньютона в м3/ч: " << result.argument * PI * pow(pipe.d, 2) / 4 * 3600 << endl << endl;
        

    }

private:
    double lymbda, Re;
};

int main()
{
    int time_count = clock();

    // Установка кодировки консоли
    setlocale(LC_ALL, "Russian");

    PipeModel pipeData; // Структура для трубопровода

    PipeSolver solv; // Cолвер

    solv.QP(pipeData); // Решение первой задачи
    cout << "Решение задачи QP: \np0 =  " << pipeData.p0 << endl;

    double dx = 8;
    int count = (int)(pipeData.L / dx + 1);
    vector<double> press_profile(count, pipeData.pl);
    solv.euler(pipeData, press_profile, dx);
    cout << "\nРешение задачи QP методом Эйлера: \np0 = " << press_profile[0] << " Па" << endl;
    writeFun(press_profile, dx);

    solv.PP(pipeData);
    cout << "\nРешение задачи PP: \nQ = " << pipeData.Q << " м3/с" " или в м3/ч: Q = " << pipeData.Q * 3600 << " м3/ч\n";

    //solv.PP_Newton(pipeData);

    // Вывод затраченного времени
    printf("\nЗатраченное время: %i ms\n", time_count);

    // Построение графика
    system("py charts.py");

    return 0;
}
