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
typedef function<double(const double& v)> residual_func_t;
typedef vector<double> profile_t;


/// @brief Расчёт коэффициента лямбда
/// @param pipe Структура трубы
/// @param oil Структура нефти
/// @param speed Скорость в трубопроводе
/// @return 
double getLambda(const pipe_properties_t& pipe, const oil_parameters_t& oil, double speed)
{
	// Расчёт относительной шероховатости
	double eps = pipe.wall.relativeRoughness();
	// Расчёт числа Рейнольдса
	double Re = speed * pipe.wall.diameter / oil.viscosity.nominal_viscosity;
	// Расчёт коэффициента лямбда
	return hydraulic_resistance_isaev(Re, eps);
}

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
		// Расчёт скорости
		double speed = 4 * Q / (PI * pow(pipe.wall.diameter, 2));     
		// Расчёт коэффициента лямбда
		double lambda = getLambda(pipe, oil, speed);
		// Длина трубопровода
		double L = pipe.profile.getLength();
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
	profile_t QP_Euler_solver(profile_t& press_prof, double speed, bool direction = true)
	{
		//Расчёт коэффициента лямбда
		double lambda = getLambda(pipe, oil, speed);
		// Плотность (для краткости записи)
		double rho = oil.density.nominal_density;
		// Перепад высот между двумя ближайшими точками (принимаем что перепад высот постоянен по всему трубопроводу)
		double dz = pipe.profile.heights[1] - pipe.profile.heights[0];
		// Шаг по координате
		double dx = pipe.profile.coordinates[1] - pipe.profile.coordinates[0];
		// Количество точек профиля
		size_t num_dots = pipe.profile.getPointCount();
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

/// @brief Решение задачи PP методом простых итераций
/// @param pipe Ссылка на структуру трубы
/// @param oil Ссылка на структуру нефти
/// @param p0 Давление в начале трубопровода
/// @param pl Давление в конце трубопровода
/// @return Возвращает расход
double PP_Iteration_solve(const pipe_properties_t& pipe, const oil_parameters_t& oil, double p0, double pl)
{
	double lym_b, lambda, Re, speed;
	// Расчёт относительной шероховатости
	double eps = pipe.wall.relativeRoughness();
	// Ограничитель итераций
	size_t itr_stop = 0;
	size_t itr_max = 1000;

	// Для компактности записей
	double L = pipe.profile.getLength();
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


/// @brief Класс с функцией невязок для решения методом Ньютона-Рафсона
class solver_Newton : public fixed_system_t<1>
{
public:
	solver_Newton(const residual_func_t& res_fun)
		: res_func{ res_fun }
	{
	}

	double residuals(const double& v) {
		double result = res_func(v);
		return result;
	}

	/// @brief Решение задачи РР методом Ньютона
	/// @param p0 Давление в начале трубопровода
	/// @param pl Давление в конце трубопровода
	/// @return Возвращает расход
	double solve()
	{
		fixed_solver_parameters_t<1, 0> parameters;
		// Создание структуры для записи результатов расчета
		fixed_solver_result_t<1> result;
		// Решение системы нелинейныйх уравнений <1> с помощью решателя Ньютона - Рафсона
		// { 0, 0 } - Начальное приближение
		fixed_newton_raphson<1>::solve_dense(*this, { 1 }, parameters, &result);

		double speed = result.argument;

		return speed;
	}
protected:
	const residual_func_t& res_func;
};

/// @brief Солвер для решения задачи PP методом Ньютона-Рафсона
class PP_solver
{
public:
	PP_solver(const pipe_properties_t& pipe, const oil_parameters_t& oil, double p0, double pl)
		: pipe{ pipe }, oil{ oil }, p0{ p0 }, pl{ pl }
	{}

	/// @brief Решает задачу PP методом Ньютона по формуле
	/// @return Возвращает расход
	double solve_newton()
	{
		residual_func_t res_fun =
			[this](const double& v)
			{
				// Расчёт коэффициента лямбда
				double lambda = getLambda(pipe, oil, v);
				// Расчёт длины трубопровода
				double L = pipe.profile.getLength();

				double H0 = p0 / (oil.density.nominal_density * g) + pipe.profile.heights.front();
				double HL = pl / (oil.density.nominal_density * g) + pipe.profile.heights.back();
				double dH = lambda * (L * pow(v, 2)) / (2 * pipe.wall.diameter * g);
				return dH + HL - H0; // Задание функции невязок
			};

		solver_Newton Newton_solver(res_fun);
		double Q = Newton_solver.solve() * PI * pow(pipe.wall.diameter, 2) / 4;

		return Q;
	}

	/// @brief Решает задачу PP методом Ньютона поверх Эйлера
	/// @return Возвращает расход 
	double solve_newton_euler()
	{
		residual_func_t res_fun =
			[this](const double& v)
			{
				// Количество точек профиля
				size_t dots_count = pipe.profile.getPointCount();
				// Профиль давлений
				profile_t press_profile(dots_count, pl);
				// Объект для решения QP методом Эйлера
				QP_tasks_solver solver(pipe, oil);
				// Решение задачи QP методом Эйлера
				solver.QP_Euler_solver(press_profile, v, false);
				// Функция невязок
				return press_profile.front() - p0;
			};

		solver_Newton Newton_solver(res_fun);
		double Q = Newton_solver.solve() * PI * pow(pipe.wall.diameter, 2) / 4;

		return Q;
	}

protected:
	const pipe_properties_t& pipe;
	const oil_parameters_t& oil;
	double p0; 
	double pl;
};

