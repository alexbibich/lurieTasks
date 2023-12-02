#pragma once

/// @brief Вывод ответа с пояснением
/// @param lable Пояснение
/// @param answer Ответ
/// @param filename Имя файла
void write_ans(string lable, auto answer, string filename="output/answers.txt")
{
    ofstream answ;
    answ.open(filename, ios::app);
    answ << lable << answer << endl;
    answ.close();
}

/// @brief Запись профиля давления в файл
/// @param press Профиль давления
/// @param dx Шаг по координате
void write_profile(profile_t& press, double& dx, string filename= "output/res.csv") {
    ofstream press_file;
    size_t profCount = press.size();
    press_file.open(filename);
    press_file << "time,x,Pressure" << endl;

    for (int i = 0; i < profCount; i++)
    {
        press_file << 0 << "," << i * dx;
        press_file << "," << press[i] << endl;
    }

    press_file.close();
}

/// @brief Инициализация данных в структурах
/// @param pipe Ссылка на структуру трубы
/// @param oil Ссылка на структуру нефти
void init_cond(pipe_properties_t& pipe, oil_parameters_t& oil)
{
    double L = 8e+4;
    double x0 = 0;
    double xl = 8e4;
    double D = 0.720;
    double thickness = 0.010;
    double delta = 15e-6;
    double z0 = 50;
    double zl = 100;
    double ro = 870;
    double visc = 15e-6;
    double p_capacity = 10e6;
    size_t n = 1000;

    pipe.profile = PipeProfile::create(n, x0, xl, z0, zl, p_capacity);
    pipe.wall.wallThickness = thickness;
    pipe.wall.diameter = D - 2 * pipe.wall.wallThickness;
    pipe.wall.equivalent_roughness = delta;

    oil.density.nominal_density = ro;
    oil.viscosity.nominal_viscosity = visc;
}

/// @brief Решение задачи 1 (QP) по формуле
TEST(QP_task, Formula)
{
    // Модель трубопровода
    pipe_properties_t pipe_prop;
    // Модель нефти
    oil_parameters_t oil;

    init_cond(pipe_prop, oil);

    double pl = 6e5;
    double Q = 3500.0 / 3600.0;
    bool direction = false;

    QP_tasks_solver solver(pipe_prop, oil);

    double p0 = solver.QP_formula_solve(pl, Q, direction);
    write_ans("Решение задачи QP по формуле, Па: ", p0);
}

/// @brief Решение задачи 2 (PP) по методом простой итерации
TEST(PP_task, Iteration)
{
    // Модель трубопровода
    pipe_properties_t pipe_prop;
    // Модель нефти
    oil_parameters_t oil;

    init_cond(pipe_prop, oil);

    double p0 = 5e6;
    double pl = 0.8e6;
    double Q = PP_Iteration_solve(pipe_prop, oil, p0, pl);
    write_ans("Решение задачи PP методом простых итераций, м3/c: ", Q);
    write_ans("Решение задачи PP методом простых итераций, м3/ч: ", Q * 3600);
}

/// @brief Решение задачи 1 (QP) по методом Эйлера
TEST(QP_task, Euler)
{
    // Модель трубопровода
    pipe_properties_t pipe_prop;
    // Модель нефти
    oil_parameters_t oil;

    init_cond(pipe_prop, oil);

    double pl = 6e5;
    double Q = 3500.0 / 3600.0;
    double speed = 4 * Q / (PI * pow(pipe_prop.wall.diameter, 2));
    bool direction = false;
    double dx = pipe_prop.profile.coordinates[1] - pipe_prop.profile.coordinates[0];
    size_t dots_count = pipe_prop.profile.getPointCount();

    profile_t press_profile(dots_count, pl);

    QP_tasks_solver solver(pipe_prop, oil);

    solver.QP_Euler_solver(press_profile, speed, direction);
    write_ans("Решение задачи QP методом Эйлера, Па: ", press_profile.front());
    write_profile(press_profile, dx);
}

/// @brief Решение задачи 2 (PP) по методом Ньютона
TEST(PP_task, Newton)
{
    // Модель трубопровода
    pipe_properties_t pipe_prop;
    // Модель нефти
    oil_parameters_t oil;

    init_cond(pipe_prop, oil);

    double p0 = 5e6;
    double pl = 0.8e6;

    // Задание функции невязок
    residual_func_t res_fun =
        [pipe_prop, oil, p0, pl](const double& v)
        {
            // Расчёт коэффициента лямбда
            double lambda = getLambda(pipe_prop, oil, v);
            // Расчёт длины трубопровода
            double L = pipe_prop.profile.getLength();

            double H0 = p0 / (oil.density.nominal_density * g) + pipe_prop.profile.heights.front();
            double HL = pl / (oil.density.nominal_density * g) + pipe_prop.profile.heights.back();
            double dH = lambda * (L * pow(v, 2)) / (2 * pipe_prop.wall.diameter * g);
            return dH + HL - H0; // Задание функции невязок
        };

    solver_Newton Newton_solver(res_fun);

    double Q = Newton_solver.solve() * PI * pow(pipe_prop.wall.diameter, 2) / 4;
    write_ans("Решение задачи PP методом Ньютона, м3/c: ", Q);
    write_ans("Решение задачи PP методом Ньютона, м3/ч: ", Q * 3600);

}

/// @brief Решение задачи 2 (PP) по методом Ньютона поверх Эйлера
TEST(PP_task, Newton_Euler)
{
    // Модель трубопровода
    pipe_properties_t pipe_prop;
    // Модель нефти
    oil_parameters_t oil;

    init_cond(pipe_prop, oil);

    double p0 = 5e6;
    double pl = 0.8e6;

    // Задание функции невязок
    residual_func_t res_fun =
        [pipe_prop, oil, p0, pl](const double& v)
        {
            size_t dots_count = pipe_prop.profile.getPointCount();
            profile_t press_profile(dots_count, pl);

            QP_tasks_solver solver(pipe_prop, oil);

            solver.QP_Euler_solver(press_profile, v, false);
            
            return press_profile[0] - p0;
        };

    solver_Newton Newton_solver(res_fun);

    double Q = Newton_solver.solve() * PI * pow(pipe_prop.wall.diameter, 2) / 4;
    write_ans("Решение задачи PP методом Ньютона поверх Эйлера, м3/c: ", Q);
    write_ans("Решение задачи PP методом Ньютона поверх Эйлера, м3/ч: ", Q * 3600);
}