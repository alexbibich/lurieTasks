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

/// @brief Класс с функцией невязок для решения методом Ньютона
class solver_Newton : public fixed_system_t<1>
{
    PipeModel& pipe;

public:
    solver_Newton(PipeModel& pipeM)
        : pipe{ pipeM }
    {}
    
    double residuals(const double& v) {

        double Re = v * pipe.d / pipe.visc;
        double lambda = hydraulic_resistance_isaev(Re, pipe.eps);
        double H0 = (pipe.p0 / (pipe.ro * g) + pipe.z0);
        double HL = (pipe.pl / (pipe.ro * g) + pipe.zl);
        double dH = lambda * (pipe.L * pow(v, 2)) / (2 * pipe.d * g);
        double result = dH + HL - H0; // Задание функции невязок
        return result;
    }
};

class solver_Newton_mix_pp_qp : public fixed_system_t<1>
{
    PipeModel& pipe;

public:
    solver_Newton_mix_pp_qp(PipeModel& pipeM)
        : pipe{ pipeM }
    {}

    /// @brief Рассчёт производной
    /// @param pipe 
    double diff(PipeModel& pipe, double dx, double dz)
    {
        double lambda = hydraulic_resistance_isaev(pipe.speed * pipe.d / pipe.visc, pipe.eps);
        double tau = lambda / 8 * pipe.ro * pow(pipe.speed, 2);
        return -4 / pipe.d * tau - pipe.ro * g * dz / dx;
    }

    /// @brief Метод Эйлера для нахождения давления в начале участка
    /// @param pipe Структура с параметрами трубопровода
    /// @param p Профиль давления
    /// @param h Шаг по координате
    void euler(PipeModel& pipe, vector<double>& p, double h, bool direct = 1)
    {
        
        double dz = (pipe.zl - pipe.z0) / (p.size() - 1);
        double diff_p = diff(pipe, h, dz);

        if (!direct) {
            for (int i = p.size() - 2; i >= 0; i--)
            {
                p[i] = p[i + 1] - h * diff_p;
            }
        }
        else {
            for (int i = 1; i < p.size(); i++)
            {
                p[i] = p[i - 1] + h * diff_p;
            }
        }
    }

    double residuals(const double& v) {
        pipe.speed = v;
        double dx = 8;
        int count = (int)(pipe.L / dx + 1);
        vector<double> press_prof(count, pipe.pl);
        euler(pipe, press_prof, dx, false);
        double result = press_prof[0] - pipe.p0; // Задание функции невязок
        return result;
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

    /// @brief Решение первой задачи по формуле
    /// @param pipe Структура с параметрами трубопровода
    void QP(PipeModel& pipe)
    {
        Re = find_Re(pipe); // Расчёт числа Рейнольдса
        lambda = hydraulic_resistance_isaev(Re, pipe.eps); // Расчёт коэффициента лямбда
        pipe.p0 = (pipe.pl / (pipe.ro * g) + pipe.zl - pipe.z0 + lambda * pipe.L / pipe.d * pow(pipe.speed, 2) / (2 * g)) * (pipe.ro * g);
        cout << "Решение задачи QP: \np0 =  " << pipe.p0 << endl;
    }

    /// @brief Решение второй задачи методом простой итерации
    /// @param pipe Структура с параметрами трубопровода
    void PP(PipeModel& pipe)
    {
        double lym_b;
        int itr_stop = 0; // Ограничитель итераций

        pipe.p0 = 5e+6;
        pipe.pl = 8e+5;

        double a = 2 * g * pipe.d / pipe.L * ((pipe.p0 - pipe.pl) / (pipe.ro * g) + pipe.z0 - pipe.zl); // Правая часть уравнения
       
        lambda = 0.02;
        do
        {
            itr_stop++;
            if (itr_stop > 1000)
            {
                cout << "Расчёт неудался!" << endl;
                break;
            }
                
            lym_b = lambda;
            pipe.speed = sqrt(a / lambda);
            Re = pipe.speed * pipe.d / pipe.visc;
            lambda = hydraulic_resistance_isaev(Re, pipe.eps);

        } while (abs(lambda - lym_b) > 0.00005);

        /*cout << "\n\n\nКоличество итераций: " << itr_stop << endl;
        cout << "lambda: " << lambda << endl;*/

        pipe.Q = PI * pow(pipe.d, 2) * pipe.speed / 4;
    }

    /// @brief Вспомогательный расчёт для Эйлера
    /// @param pipe Структура с параметрами трубопровода
    /// @param dx
    /// @param dz 
    /// @return 
    double diff(PipeModel& pipe, double dx, double dz)
    {
        double tau = lambda / 8 * pipe.ro * pow(pipe.speed, 2);
        return -4 / pipe.d * tau - pipe.ro * g * dz / dx;
    }

    /// @brief Метод Эйлера для нахождения давления в начале участка
    /// @param pipe Структура с параметрами трубопровода
    /// @param p Профиль давления
    /// @param h Шаг по координате
    void euler(PipeModel& pipe, vector<double>& p, double h, bool direct = true)
    {
        lambda = hydraulic_resistance_isaev(find_Re(pipe), pipe.eps);
        double dz = (pipe.zl - pipe.z0) / (p.size() - 1);
        double diff_p = diff(pipe, h, dz);

        if (direct) {
            for (int i = p.size() - 2; i >= 0; i--)
            {
                p[i] = p[i + 1] - h * diff_p;
            }
        }
        else {
            for (int i = 1; i < p.size(); i++)
            {
                p[i] = p[i - 1] + h * diff_p;
            }
        }
    }

    /// @brief Расчёт скорости методом Ньютона
    /// @param pipe Структура с параметрами трубопровода
    void PP_Newton(PipeModel& pipe)
    {
        solver_Newton test(pipe);
        // Задание настроек решателя по умолчанию
        fixed_solver_parameters_t<1, 0> parameters;
        // Создание структуры для записи результатов расчета
        fixed_solver_result_t<1> result;
        // Решение системы нелинейныйх уравнений <1> с помощью решателя Ньютона - Рафсона
        // { 0, 0 } - Начальное приближение
        fixed_newton_raphson<1>::solve_dense(test, { 1 }, parameters, &result);

        cout << "\nРешение Ньютона: " << result.argument << endl;
        cout << "Решение Ньютона в м3/с: " << result.argument * PI * pow(pipe.d, 2) / 4 << endl;
        cout << "Решение Ньютона в м3/ч: " << result.argument * PI * pow(pipe.d, 2) / 4 * 3600 << endl << endl;
    }

    void PP_QP_mix(PipeModel& pipe)
    {
        solver_Newton_mix_pp_qp test(pipe);
        // Задание настроек решателя по умолчанию
        fixed_solver_parameters_t<1, 0> parameters;
        // Создание структуры для записи результатов расчета
        fixed_solver_result_t<1> result;
        // Решение системы нелинейныйх уравнений <1> с помощью решателя Ньютона - Рафсона
        // { 0, 0 } - Начальное приближение
        fixed_newton_raphson<1>::solve_dense(test, { 1 }, parameters, &result);

        cout << "\nРешение Ньютона поверх Эйлера: " << result.argument << endl;
        cout << "Решение Ньютона поверх Эйлера в м3/с: " << result.argument * PI * pow(pipe.d, 2) / 4 << endl;
        cout << "Решение Ньютона поверх Эйлера в м3/ч: " << result.argument * PI * pow(pipe.d, 2) / 4 * 3600 << endl << endl;

    }

private:
    double lambda, Re;
};

int main()
{
    int time_count = clock();

    // Установка кодировки консоли
    setlocale(LC_ALL, "Russian");

    // Структура для трубопровода
    PipeModel pipeData; 

    // объекта Cолвер
    PipeSolver solv; 

    // Решение первой задачи
    solv.QP(pipeData); 
    
    double dx = 8;
    int count = (int)(pipeData.L / dx + 1);
    vector<double> press_profile(count, pipeData.pl);
    // Решение первой задачи методом Эйлера
    solv.euler(pipeData, press_profile, dx);
    cout << "\nРешение задачи QP методом Эйлера: \np0 = " << press_profile[0] << " Па" << endl;
    writeFun(press_profile, dx);

    // Решение второй задачи методом итераций
    solv.PP(pipeData); 
    cout << "\nРешение задачи PP: \nQ = " << pipeData.Q << " м3/с" " или в м3/ч: Q = " << pipeData.Q * 3600 << " м3/ч\n";

    pipeData.p0 = 5e+6;
    pipeData.pl = 8e+5;
    // Решение второй задачи методом Ньютона
    solv.PP_Newton(pipeData); 

    // Решение второй задачи методом Ньютона поверх Эйлера
    solv.PP_QP_mix(pipeData);

    // Вывод затраченного времени
    printf("\nЗатраченное время: %i ms\n", time_count);

    // Построение графика
    system("py charts.py");

    return 0;
}
