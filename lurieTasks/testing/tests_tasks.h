#pragma once

/// @brief ����� ������ � ����������
/// @param lable ���������
/// @param answer �����
/// @param filename ��� �����
void write_ans(string lable, auto answer, string filename="answers.txt")
{
    ofstream answ;
    answ.open(filename, ios::app);
    answ << lable << answer << endl;
    answ.close();
}

/// @brief ������ ������� �������� � ����
/// @param press ������� ��������
/// @param dx ��� �� ����������
void write_profile(vector<double>& press, double& dx) {
    ofstream press_file;
    size_t profCount = press.size();
    press_file.open("res.csv");
    press_file << "time,x,Pressure" << endl;

    for (int i = 0; i < profCount; i++)
    {
        press_file << 0 << "," << i * dx;
        press_file << "," << press[i] << endl;
    }

    press_file.close();
}

/// @brief ������������� ������ � ����������
/// @param pipe ������ �� ��������� �����
/// @param oil ������ �� ��������� �����
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

/// @brief ������� ������ 1 (QP) �� �������
TEST(QP_task, Formula)
{
    // ������ ������������
    pipe_properties_t pipe_prop;
    // ������ �����
    oil_parameters_t oil;

    init_cond(pipe_prop, oil);

    double pl = 6e5;
    double Q = 3500.0 / 3600.0;
    bool direction = false;

    QP_tasks_solver solver(pipe_prop, oil);

    double p0 = solver.QP_formula_solve(pl, Q, direction);
    write_ans("������� ������ QP �� �������, ��: ", p0);
}

/// @brief ������� ������ 2 (PP) �� ������� ������� ��������
TEST(PP_task, Iteration)
{
    // ������ ������������
    pipe_properties_t pipe_prop;
    // ������ �����
    oil_parameters_t oil;

    init_cond(pipe_prop, oil);

    double p0 = 5e6;
    double pl = 0.8e6;

    PP_tasks_solver solver(pipe_prop, oil);

    double Q = solver.PP_Iteration_solve(p0, pl);
    write_ans("������� ������ PP ������� ������� ��������, �3/c: ", Q);
    write_ans("������� ������ PP ������� ������� ��������, �3/�: ", Q * 3600);
}

/// @brief ������� ������ 1 (QP) �� ������� ������
TEST(QP_task, Euler)
{
    // ������ ������������
    pipe_properties_t pipe_prop;
    // ������ �����
    oil_parameters_t oil;

    init_cond(pipe_prop, oil);

    double pl = 6e5;
    double Q = 3500.0 / 3600.0;
    double speed = 4 * Q / (PI * pow(pipe_prop.wall.diameter, 2));
    bool direction = false;
    double dx = pipe_prop.profile.coordinates[1] - pipe_prop.profile.coordinates[0];

    QP_tasks_solver solver(pipe_prop, oil);

    vector<double> press_profile = solver.QP_Euler_solver(pl, speed, direction);
    write_ans("������� ������ QP ������� ������, ��: ", press_profile.front());
    write_profile(press_profile, dx);
}

/// @brief ������� ������ 2 (PP) �� ������� �������
TEST(PP_task, Newton)
{
    // ������ ������������
    pipe_properties_t pipe_prop;
    // ������ �����
    oil_parameters_t oil;

    init_cond(pipe_prop, oil);

    double p0 = 5e6;
    double pl = 0.8e6;

    PP_tasks_solver solver(pipe_prop, oil);

    double Q = solver.PP_Newton_solve(p0, pl);
    write_ans("������� ������ PP ������� �������, �3/c: ", Q);
    write_ans("������� ������ PP ������� �������, �3/�: ", Q * 3600);

}

/// @brief ������� ������ 2 (PP) �� ������� ������� ������ ������
TEST(PP_task, Newton_Euler)
{
    // ������ ������������
    pipe_properties_t pipe_prop;
    // ������ �����
    oil_parameters_t oil;

    init_cond(pipe_prop, oil);

    double p0 = 5e6;
    double pl = 0.8e6;

    PP_tasks_solver solver(pipe_prop, oil);

    double Q = solver.PP_Newton_Euler_solve(p0, pl);
    write_ans("������� ������ PP ������� ������� ������ ������, �3/c: ", Q);
    write_ans("������� ������ PP ������� ������� ������ ������, �3/�: ", Q * 3600);
}