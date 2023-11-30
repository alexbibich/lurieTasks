#pragma once

#include <fstream>
#include <vector>
#include <cstdlib>

#include <iomanip>
#include <fixed/fixed.h>
#include <pde_solvers/pde_solvers.h>

#include <cmath>

#define PI 3.14
#define g  9.81

using namespace std;

/// @brief Класс для решения задач QP
class QP_tasks_solver
{
public:
	/// @brief Решение задачи QP по формуле из учебника Лурье 
	/// @param p Известное давление
	/// @param Q Известный расход
	/// @param direction Означает направление расчёта
	/// @return Возвращает неизвестное давление
	double QP_formula_solve(double p, double Q, bool direction = true)
	{
		// Расчёт относительной шероховатости
		double eps = pipe.wall.getCompressionRatio();           
		// Расчёт скорости
		double speed = 4 * Q / (PI * pow(pipe.wall.diameter, 2));          
		// Расчёт числа Рейнольдса
		double Re = speed * pipe.wall.diameter / oil.viscosity.nominal_viscosity; 
		// Расчёт коэффициента лямбда
		double lambda = hydraulic_resistance_isaev(Re, eps); 
		// Длина трубопровода
		double L = pipe.profile.coordinates.back() - pipe.profile.coordinates.front();
		// Перепал высот
		double dz = pipe.profile.heights.back() - pipe.profile.heights.front();
		// Потеря напора на трение
		double dH = lambda * L * pow(speed, 2) / (pipe.wall.diameter * 2 * g); 
		// rho * g
		double rho_g = oil.density.nominal_density * g;
		// Расчёт неизвестного давления в зависимости от выбранного направления
		double p_ans = direction ? (p / rho_g - dz - dH) * rho_g : (p / rho_g + dz + dH) * rho_g;

		return p_ans;
	}

	/// @brief Решение задачи QP методом Эйлера
	/// @param p Известное давление
	/// @param speed Известная скорость
	/// @param direction Направление расчёта профиля
	/// @return Профиль давления
	vector<double> QP_Euler_solver(double p, double speed, bool direction = true)
	{
		// Расчёт относительной шероховатости
		double eps = pipe.wall.getCompressionRatio();
		// Расчёт числа Рейнольдса
		double Re = speed * pipe.wall.diameter / oil.viscosity.nominal_viscosity;
		// Расчёт коэффициента лямбда
		double lambda = hydraulic_resistance_isaev(Re, eps);
		// Плотность (для краткости записи)
		double rho = oil.density.nominal_density;
		// Перепад высот между двумя ближайшими точками (принимаем что перепад высот постоянен по всему трубопроводу)
		double dz = pipe.profile.heights[1] - pipe.profile.heights[0];
		// Шаг по координате
		double dx = pipe.profile.coordinates[1] - pipe.profile.coordinates[0];
		// Количество точек профиля
		size_t num_dots = pipe.profile.getPointCount();
		// Профиль давления
		vector<double> press_prof(num_dots, p);
		// Расчёт тау и производной 
		double tau = lambda / 8 * rho * pow(speed, 2);
		double diff = -4 / pipe.wall.diameter * tau - rho * g * dz / dx;

		for (size_t i = 1; i < num_dots; i++)
		{
			size_t index = direction ? i : num_dots - 1 - i;
			if (direction)
				press_prof[index] = press_prof[index - 1] + dx * diff;
			else
				press_prof[index] = press_prof[index + 1] - dx * diff;
		}

		return press_prof;
	}
	/// @brief Конструктор класса
	/// @param pipe_m Ссылка на структуру трубы
	/// @param oil_m Ссылка на структуру нефти
	QP_tasks_solver(const pipe_properties_t& pipe_m, const oil_parameters_t& oil_m)
		: pipe{ pipe_m }, oil{ oil_m }
	{
	}

protected:
	const pipe_properties_t& pipe;
	const oil_parameters_t& oil;

};

/// @brief Класс с функцией невязок для решения методом Ньютона
class solver_Newton : public fixed_system_t<1>
{
public:
	solver_Newton(const pipe_properties_t& pipe_t, const oil_parameters_t& oil_t, double p0, double pl, double epsilon)
		: pipe{ pipe_t }, oil{ oil_t }, p{ p0, pl }, eps{ epsilon }
	{
	}

	double residuals(const double& v) {
		// Расчёт числа Рейнольдса
		double Re = v * pipe.wall.diameter / oil.viscosity.nominal_viscosity;
		double lambda = hydraulic_resistance_isaev(Re, eps);
		// Расчёт длины трубопровода
		double L = pipe.profile.coordinates.back() - pipe.profile.coordinates.front();

		double H0 = p[0] / (oil.density.nominal_density * g) + pipe.profile.heights.front();
		double HL = p[1] / (oil.density.nominal_density * g) + pipe.profile.heights.back();
		double dH = lambda * (L * pow(v, 2)) / (2 * pipe.wall.diameter * g);
		double result = dH + HL - H0; // Задание функции невязок
		return result;
	}
protected:
	const pipe_properties_t& pipe;
	const oil_parameters_t& oil;
	double p[2];
	double eps;
};

/// @brief Солвер для Ньютона поверх Эйлера
class solver_Newton_Euler : public fixed_system_t<1>
{
public:
	double residuals(const double& v) {
		QP_tasks_solver euler_solver(pipe, oil);
		vector<double> press_prof = euler_solver.QP_Euler_solver(p[1], v, dir);
		// Задание функции невязок
		double result = press_prof[0] - p[0]; 
		return result;
	}

	/// @brief Конструктор солвера
	/// @param pipe_t Ссылка на структуру трубы
	/// @param oil_t Ссылка на структуру нефти
	/// @param p0 Давление в начале трубопровода
	/// @param pl Давление в конце трубопровода
	/// @param direction Направление расчёта
	solver_Newton_Euler(const pipe_properties_t& pipe_t, const oil_parameters_t& oil_t, double p0, double pl, bool direction = false)
		: pipe{ pipe_t }, oil{ oil_t }, p{ p0, pl }, dir{ direction }
	{
	}

protected:
	const pipe_properties_t& pipe;
	const oil_parameters_t& oil;
	double p[2];
	bool dir;
};

/// @brief Солвер для решения задач PP
class PP_tasks_solver
{
public:
	/// @brief Решение задачи PP методом простых итераций
	/// @param p0 Давление в начале трубопровода
	/// @param pl Давление в конце трубопровода
	/// @return Возвращает расход
	double PP_Iteration_solve(double p0, double pl)
	{
		double lym_b, lambda, Re, speed;
		// Расчёт относительной шероховатости
		double eps = pipe.wall.getCompressionRatio(); 
		// Ограничитель итераций
		size_t itr_stop = 0; 
		size_t itr_max = 1000;
		
		// Для компактности записей
		double L = pipe.profile.coordinates.back() - pipe.profile.coordinates.front();
		double d = pipe.wall.diameter;
		double visc = oil.viscosity.nominal_viscosity;
		double rho = oil.density.nominal_density;
		double z0 = pipe.profile.heights.front();
		double zl = pipe.profile.heights.back();

		// Правая часть уравнения
		double a = 2 * g * d / L * ((p0 - pl) / (rho * g) + z0 - zl); 

		lambda = 0.02;
		do
		{
			itr_stop++;
			if (itr_stop > itr_max)
				throw std::runtime_error("Reached maximum number of iterations");

			lym_b = lambda;
			speed = sqrt(a / lambda);
			Re = speed * d / visc;
			lambda = hydraulic_resistance_isaev(Re, eps);

		} while (abs(lambda - lym_b) > 0.000005);

		double Q = PI * pow(d, 2) * speed / 4;

		return Q;
	}

	/// @brief Решение задачи РР методом Ньютона
	/// @param p0 Давление в начале трубопровода
	/// @param pl Давление в конце трубопровода
	/// @return Возвращает расход
	double PP_Newton_solve(double p0, double pl)
	{
		double eps = pipe.wall.getCompressionRatio();
		solver_Newton test(pipe, oil, p0, pl, eps);
		// Задание настроек решателя по умолчанию
		fixed_solver_parameters_t<1, 0> parameters;
		// Создание структуры для записи результатов расчета
		fixed_solver_result_t<1> result;
		// Решение системы нелинейныйх уравнений <1> с помощью решателя Ньютона - Рафсона
		// { 0, 0 } - Начальное приближение
		fixed_newton_raphson<1>::solve_dense(test, { 1 }, parameters, &result);

		double Q = result.argument * PI * pow(pipe.wall.diameter, 2) / 4;

		return Q;
	}

	/// @brief Решение задачи PP Ньютоном поверх Эйлера
	/// @param p0 Давление в начале трубопровода
	/// @param pl Давление в конце трубопровода
	/// @return Возвращает расход
	double PP_Newton_Euler_solve(double p0, double pl)
	{
		solver_Newton_Euler solver(pipe, oil, p0, pl);
		// Задание настроек решателя по умолчанию
		fixed_solver_parameters_t<1, 0> parameters;
		// Создание структуры для записи результатов расчета
		fixed_solver_result_t<1> result;
		// Решение системы нелинейныйх уравнений <1> с помощью решателя Ньютона - Рафсона
		// { 0, 0 } - Начальное приближение
		fixed_newton_raphson<1>::solve_dense(solver, { 1 }, parameters, &result);

		double Q = result.argument * PI * pow(pipe.wall.diameter, 2) / 4;

		return Q;
	}

	PP_tasks_solver(const pipe_properties_t& pipe_m, const oil_parameters_t& oil_m)
		: pipe{ pipe_m }, oil{ oil_m }
	{
	}

protected:
	const pipe_properties_t& pipe;
	const oil_parameters_t& oil;
};



