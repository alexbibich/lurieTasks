﻿#define GTEST_BREAK_ON_FAILURE 1
#define GTEST_CATCH_EXCEPTIONS 0
#define GTEST_HAS_SEH 0
#define _VARIADIC_MAX 10 /* for gtest */
#include "gtest/gtest.h"

#include <fixed/fixed.h>
#include <pde_solvers/pde_solvers.h>
#include "solvers/solver_test.h"

#include <iostream>
#include <fstream>
#include <filesystem>

/// @brief Возвращает тестовую строку в формате TestBundle.TestName
inline std::string get_test_string() {
    auto test_info = ::testing::UnitTest::GetInstance()->current_test_info();
    auto test_string = std::string(test_info->test_case_name()) + "." + std::string(test_info->name());
    return test_string;
}

/// @brief Создает путь для теста в папке "<pde_solvers_root>/testing_out/"
/// Возвращает путь к созданной папке
inline std::string prepare_test_folder()
{
    std::string path = std::string("../testing_out/") + get_test_string() + "/";
    std::filesystem::create_directories(path);
    return path;
}

/// @brief Очищает содержимое файла
/// @param filename Название файла
void clear_file(string filename = "output/answers.txt")
{
    ofstream ofs;
    ofs.open(filename, ofstream::out | ofstream::trunc);
    ofs.close();
}

#include "testing/tests_tasks.h"

int main(int argc, char** argv)
{
    clear_file(); // Очищайм содержимое файла
    ::testing::InitGoogleTest(&argc, argv);
#ifdef _WIN32
    std::wcout.imbue(std::locale("rus_rus.866"));
#endif
    int res = RUN_ALL_TESTS();
    return res;
}

