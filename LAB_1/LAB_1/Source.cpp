// Лабораторная работа №1. Методы решения систем линейных алгебраических уравнений
// Вариант 13
//   3,8*x1 + 6,7*x2 - 1,2*x3 =  5,2
//   6,4*x1 + 1,3*x2 - 2,7*x3 =  3,8
//   2,4*x1 - 4,5*x2 + 3,5*x3 = -0,6
//
// Задание 1: метод Гаусса с выбором главного элемента по всей матрице
// Задание 2: метод итераций (точность 1e-5)
//
// Формат входного файла MY_SLAU.txt (десятичный разделитель - ТОЧКА):
//   n
//   a11 a12 ... a1n b1
//   ...
//   an1 an2 ... ann bn

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <locale>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <utility>

#include <windows.h>

using Vec = std::vector<double>;
using Mat = std::vector<Vec>;

const double ZERO_EPS = 1e-10;  // всё, что по модулю меньше, считаем нулём
const double ITER_EPS = 1e-5;   // требуемая точность метода итераций
const int    MAX_ITER = 10000;  // защита от бесконечного цикла

// ----------------------------------------------------------------------------
//                                  ВВОД / ВЫВОД
// ----------------------------------------------------------------------------

// Печать расширенной матрицы системы [A | b]
void PrintMatrix(const Mat& A, const Vec& b) {
    int n = (int)b.size();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            std::cout << std::setw(12) << A[i][j] << ' ';
        }
        std::cout << "| " << std::setw(12) << b[i] << '\n';
    }
}

void PrintVector(const std::string& name, const Vec& x, int precision = 8) {
    std::cout << name << ":\n";
    for (size_t i = 0; i < x.size(); i++) {
        std::cout << "  x" << i + 1 << " = " << std::fixed
            << std::setprecision(precision) << x[i] << '\n';
    }
    std::cout.unsetf(std::ios::fixed);
}

// Чтение системы из файла (Ваш вариант считывания, обобщённый на любое n)
bool ReadSystem(const std::string& filename, int& matrix_size, Mat& A, Vec& b) {
    std::ifstream SLAU_file(filename);
    if (!SLAU_file) {
        std::cout << "Ошибка: не удалось открыть файл " << filename << '\n';
        return false;
    }
    SLAU_file.imbue(std::locale::classic());

    std::vector<double> numbers_left;   // коэффициенты при неизвестных (по строкам)
    std::vector<double> numbers_right;  // правые части

    if (!(SLAU_file >> matrix_size) || matrix_size <= 0) {
        std::cout << "Ошибка: в начале файла должен быть размер системы n > 0\n";
        return false;
    }

    double number;
    int ctr = 1;
    while (SLAU_file >> number) {
        if (ctr % (matrix_size + 1) == 0) {   // каждое (n+1)-е число - правая часть
            numbers_right.push_back(number);
        }
        else {
            numbers_left.push_back(number);
        }
        ctr++;
    }

    if ((int)numbers_left.size() != matrix_size * matrix_size ||
        (int)numbers_right.size() != matrix_size) {
        std::cout << "Ошибка: в файле неверное количество чисел "
            "(проверьте, что используется точка, а не запятая)\n";
        return false;
    }

    // Переносим плоский массив в двумерную матрицу
    A.assign(matrix_size, Vec(matrix_size));
    for (int i = 0; i < matrix_size; i++)
        for (int j = 0; j < matrix_size; j++)
            A[i][j] = numbers_left[i * matrix_size + j];
    b = numbers_right;
    return true;
}

// Максимальная по модулю невязка  max|A*x - b|
double Residual(const Mat& A, const Vec& b, const Vec& x) {
    double r = 0;
    for (size_t i = 0; i < b.size(); i++) {
        double s = -b[i];
        for (size_t j = 0; j < b.size(); j++) s += A[i][j] * x[j];
        r = std::max(r, std::fabs(s));
    }
    return r;
}

// ----------------------------------------------------------------------------
//             ЗАДАНИЕ 1. МЕТОД ГАУССА С ВЫБОРОМ ГЛАВНОГО ЭЛЕМЕНТА
// ----------------------------------------------------------------------------

enum class GaussStatus { UNIQUE, INFINITE, NO_SOLUTION };

struct GaussResult {
    GaussStatus status;
    Vec x;        // решение (при INFINITE - одно из частных, свободные неизвестные = 0)
    int rank;     // найденный ранг матрицы
    double det;   // определитель
};

// Матрицы A и b передаются по значению, чтобы исходные данные не портились
GaussResult GaussFullPivot(Mat A, Vec b) {
    int n = (int)b.size();
    std::vector<int> col(n);              // col[k] - номер исходной неизвестной в k-м столбце
    std::iota(col.begin(), col.end(), 0);

    int rank = n;
    double det = 1.0;

    // ---- Прямой ход ----
    for (int k = 0; k < n; k++) {
        // Поиск главного элемента - максимального по модулю в подматрице [k..n-1]x[k..n-1]
        double mx = 0;
        int pi = k, pj = k;
        for (int i = k; i < n; i++)
            for (int j = k; j < n; j++)
                if (std::fabs(A[i][j]) > mx) { mx = std::fabs(A[i][j]); pi = i; pj = j; }

        // Главный элемент равен нулю => деление на него невозможно.
        // Вместо аварийного завершения фиксируем ранг и выходим из цикла.
        if (mx < ZERO_EPS) {
            rank = k;
            det = 0.0;
            break;
        }

        // Перестановка строк
        if (pi != k) {
            std::swap(A[pi], A[k]);
            std::swap(b[pi], b[k]);
            det = -det;
        }
        // Перестановка столбцов (запоминаем, какая неизвестная куда ушла)
        if (pj != k) {
            for (int i = 0; i < n; i++) std::swap(A[i][pj], A[i][k]);
            std::swap(col[pj], col[k]);
            det = -det;
        }

        det *= A[k][k];

        // Исключение x_k из строк ниже (здесь A[k][k] != 0 гарантированно)
        for (int i = k + 1; i < n; i++) {
            double m = A[i][k] / A[k][k];
            for (int j = k; j < n; j++) A[i][j] -= m * A[k][j];
            b[i] -= m * b[k];
        }
    }

    GaussResult res;
    res.rank = rank;
    res.det = det;
    res.x.assign(n, 0.0);

    // ---- Вырожденный случай: проверка совместности ----
    if (rank < n) {
        for (int i = rank; i < n; i++) {
            if (std::fabs(b[i]) > ZERO_EPS) {      // 0 = b_i != 0
                res.status = GaussStatus::NO_SOLUTION;
                return res;
            }
        }
        res.status = GaussStatus::INFINITE;        // 0 = 0 - уравнение лишнее
    }
    else {
        res.status = GaussStatus::UNIQUE;
    }

    // ---- Обратный ход (свободные неизвестные, если они есть, равны 0) ----
    Vec y(n, 0.0);   // решение в порядке переставленных столбцов
    for (int i = rank - 1; i >= 0; i--) {
        double s = b[i];
        for (int j = i + 1; j < rank; j++) s -= A[i][j] * y[j];
        y[i] = s / A[i][i];
    }
    for (int k = 0; k < n; k++) res.x[col[k]] = y[k];   // возвращаем исходный порядок

    return res;
}

// ----------------------------------------------------------------------------
//                       ЗАДАНИЕ 2. МЕТОД ИТЕРАЦИЙ
// ----------------------------------------------------------------------------

// Норма матрицы, подчинённая октаэдрической норме вектора ||x||_1 = sum|x_i|
// (максимальная сумма модулей по СТОЛБЦАМ)
double NormOctahedral(const Mat& B) {
    int n = (int)B.size();
    double r = 0;
    for (int j = 0; j < n; j++) {
        double s = 0;
        for (int i = 0; i < n; i++) s += std::fabs(B[i][j]);
        r = std::max(r, s);
    }
    return r;
}

// Норма матрицы, подчинённая кубической норме вектора ||x||_inf = max|x_i|
// (максимальная сумма модулей по СТРОКАМ)
double NormCubic(const Mat& B) {
    int n = (int)B.size();
    double r = 0;
    for (int i = 0; i < n; i++) {
        double s = 0;
        for (int j = 0; j < n; j++) s += std::fabs(B[i][j]);
        r = std::max(r, s);
    }
    return r;
}

// Приведение системы к виду, удобному для итераций (диагональное преобладание).
// Исходная система для варианта 13 НЕ имеет диагонального преобладания,
// поэтому заменяем уравнения их линейными комбинациями:
//   новое 1-е = 2-е исходное                         (|6,4| > |1,3| + |-2,7|)
//   новое 2-е = 1-е исходное                         (|6,7| > |3,8| + |-1,2|)
//   новое 3-е = 1-е - 2-е + 3-е исходные             (|5,0| > |-0,2| + |0,9|)
// Матрица T задаёт, из каких исходных уравнений (с какими коэффициентами) собрано
// каждое новое. det T = -1 != 0, значит новая система эквивалентна исходной.
void PrepareForIterations(const Mat& A, const Vec& b, Mat& A2, Vec& b2) {
    int n = (int)b.size();
    Mat T(n, Vec(n, 0.0));
    if (n == 3) {
        T = { { 0,  1, 0 },
              { 1,  0, 0 },
              { 1, -1, 1 } };
    }
    else {
        for (int i = 0; i < n; i++) T[i][i] = 1.0;   // для другого n - без преобразования
        std::cout << "(Преобразование T задано для n = 3, для другого n система не меняется)\n";
    }

    A2.assign(n, Vec(n, 0.0));
    b2.assign(n, 0.0);
    for (int i = 0; i < n; i++)
        for (int k = 0; k < n; k++) {
            b2[i] += T[i][k] * b[k];
            for (int j = 0; j < n; j++) A2[i][j] += T[i][k] * A[k][j];
        }
}

void SolveByIterations(const Mat& A, const Vec& b, const Vec& x_gauss) {
    int n = (int)b.size();

    // 1. Приводим систему к виду, удобному для итераций
    Mat A2; Vec b2;
    PrepareForIterations(A, b, A2, b2);
    std::cout << "Система, подготовленная для итераций:\n";
    PrintMatrix(A2, b2);

    // 2. Строим x = B*x + c   (x_i = (b_i - sum_{j!=i} a_ij x_j) / a_ii)
    Mat B(n, Vec(n, 0.0));
    Vec c(n, 0.0);
    for (int i = 0; i < n; i++) {
        if (std::fabs(A2[i][i]) < ZERO_EPS) {
            std::cout << "Диагональный элемент a" << i + 1 << i + 1
                << " равен нулю - метод итераций неприменим.\n";
            return;
        }
        for (int j = 0; j < n; j++)
            if (j != i) B[i][j] = -A2[i][j] / A2[i][i];
        c[i] = b2[i] / A2[i][i];
    }

    std::cout << "\nМатрица B и вектор c (x = B*x + c):\n";
    PrintMatrix(B, c);

    // 3. Проверка условия сходимости: достаточно, чтобы хотя бы одна норма была < 1
    double q1 = NormOctahedral(B);   // ||B||_1   (по столбцам)
    double qinf = NormCubic(B);        // ||B||_inf (по строкам)
    std::cout << std::fixed << std::setprecision(6)
        << "\n||B||_1   (октаэдрическая) = " << q1
        << "\n||B||_inf (кубическая)     = " << qinf << '\n';
    std::cout.unsetf(std::ios::fixed);

    double q = std::min(q1, qinf);
    if (q >= 1.0) {
        std::cout << "Обе нормы >= 1: достаточное условие сходимости не выполнено, "
            "итерации не запускаем.\n";
        return;
    }
    bool useCubic = (qinf <= q1);   // в какой норме оцениваем погрешность
    std::cout << "Условие сходимости выполнено (q = " << q << " < 1), метод сходится.\n\n";

    // 4. Итерации. Два одномерных массива: предыдущая и текущая итерации
    Vec x_prev = c;        // начальное приближение x^(0) = c
    Vec x_curr(n, 0.0);

    std::cout << "Итерации (точность " << ITER_EPS << "):\n";
    int it = 0;
    bool converged = false;
    while (it < MAX_ITER) {
        it++;
        for (int i = 0; i < n; i++) {
            double s = c[i];
            for (int j = 0; j < n; j++) s += B[i][j] * x_prev[j];
            x_curr[i] = s;
        }

        // ||x^(k) - x^(k-1)|| в той же норме, что и q
        double d = 0;
        for (int i = 0; i < n; i++) {
            double diff = std::fabs(x_curr[i] - x_prev[i]);
            d = useCubic ? std::max(d, diff) : d + diff;
        }

        std::cout << std::setw(3) << it << ": ";
        std::cout << std::fixed << std::setprecision(7);
        for (int i = 0; i < n; i++) std::cout << x_curr[i] << ' ';
        std::cout << "  |dx| = " << std::scientific << std::setprecision(2) << d << '\n';
        std::cout.unsetf(std::ios::floatfield);

        x_prev = x_curr;   // текущая итерация становится предыдущей

        // Апостериорная оценка: ||x^(k) - x*|| <= q/(1-q) * ||x^(k) - x^(k-1)||
        if (q / (1.0 - q) * d < ITER_EPS) { converged = true; break; }
    }

    if (!converged) {
        std::cout << "Достигнут предел числа итераций.\n";
    }
    std::cout << "\nКоличество итераций: " << it << '\n';
    PrintVector("Решение методом итераций", x_curr, 6);

    double maxdiff = 0;
    for (int i = 0; i < n; i++) maxdiff = std::max(maxdiff, std::fabs(x_curr[i] - x_gauss[i]));
    std::cout << std::scientific << std::setprecision(2)
        << "max|x_итер - x_Гаусс| = " << maxdiff << '\n';
    std::cout << "Невязка max|Ax - b| = " << Residual(A, b, x_curr) << '\n';
    std::cout.unsetf(std::ios::floatfield);
}

// ----------------------------------------------------------------------------
//                                    MAIN
// ----------------------------------------------------------------------------

int main() {
    SetConsoleOutputCP(CP_UTF8);   // чтобы русский текст отображался в консоли Windows

    int matrix_size = 0;
    Mat A;
    Vec b;
    if (!ReadSystem("MY_SLAU.txt", matrix_size, A, b)) return 1;

    std::cout << "Исходная система (n = " << matrix_size << "):\n";
    PrintMatrix(A, b);

    // ---------------- Задание 1 ----------------
    std::cout << "\n===== Задание 1. Метод Гаусса с выбором главного элемента =====\n";
    GaussResult g = GaussFullPivot(A, b);
    std::cout << "Ранг матрицы: " << g.rank << ", определитель: " << g.det << '\n';

    switch (g.status) {
    case GaussStatus::UNIQUE:
        PrintVector("Единственное решение", g.x, 8);
        std::cout << std::scientific << std::setprecision(2)
            << "Невязка max|Ax - b| = " << Residual(A, b, g.x) << '\n';
        std::cout.unsetf(std::ios::floatfield);
        break;
    case GaussStatus::INFINITE:
        std::cout << "Определитель = 0, система совместна: решений бесконечно много.\n";
        PrintVector("Одно из частных решений (свободные неизвестные = 0)", g.x, 8);
        break;
    case GaussStatus::NO_SOLUTION:
        std::cout << "Определитель = 0, система несовместна: решений нет.\n";
        break;
    }

    // ---------------- Задание 2 ----------------
    std::cout << "\n===== Задание 2. Метод итераций =====\n";
    if (g.status != GaussStatus::UNIQUE) {
        std::cout << "Матрица вырождена - метод итераций не применяется.\n";
        return 0;
    }
    SolveByIterations(A, b, g.x);

    return 0;
}