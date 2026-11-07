#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <algorithm>
#include <limits>

constexpr double EPSILON = 1e-6;
constexpr int MAX_ITERATIONS = 100000;

// Вариант №4:
// x^3 - 0.1x^2 + 0.4x - 1.5 = 0

double F(double x) {
    return x * x * x - 0.1 * x * x + 0.4 * x - 1.5;
}

// Первая производная
double dF(double x) {
    return 3.0 * x * x - 0.2 * x + 0.4;
}

// Точка пересечения хорды с осью Ox
double chordPoint(double left, double right) {
    const double fLeft = F(left);
    const double fRight = F(right);

    return left - fLeft * (right - left) / (fRight - fLeft);
}

struct Result {
    std::string method;
    std::string start;
    double root;
    int iterations;
    double errorBound;
};

// Проверка наличия корня на отрезке по смене знака
bool hasRootBracket(double left, double right) {
    return F(left) * F(right) <= 0.0;
}

// 1. Аналитическое отделение корня
void separateRoot() {
    std::cout << "===== 1. ОТДЕЛЕНИЕ КОРНЯ =====\n";

    std::cout << "f(x) = x^3 - 0.1x^2 + 0.4x - 1.5\n";
    std::cout << "f'(x) = 3x^2 - 0.2x + 0.4\n";

    // Дискриминант производной отрицателен,
    // а коэффициент при x^2 положителен.
    const double discriminant = (-0.2) * (-0.2) - 4.0 * 3.0 * 0.4;

    std::cout << "D(f') = " << discriminant << " < 0\n";
    std::cout << "Следовательно, f'(x) > 0 для любого x.\n";
    std::cout << "Функция строго возрастает и имеет не более одного корня.\n\n";

    std::cout << "f(1.0) = " << F(1.0) << " < 0\n";
    std::cout << "f(1.1) = " << F(1.1) << " > 0\n";
    std::cout << "Следовательно, единственный корень находится на [1.0; 1.1].\n\n";
}

// 2. Метод половинного деления
Result bisectionMethod(double left, double right) {
    int iterations = 0;

    if (!hasRootBracket(left, right)) {
        return { "Половинного деления", "[1.0; 1.1]", std::numeric_limits<double>::quiet_NaN(), 0, std::numeric_limits<double>::quiet_NaN() };
    }

    // Половина длины отрезка является оценкой погрешности середины.
    while ((right - left) / 2.0 > EPSILON && iterations < MAX_ITERATIONS) {
        const double mid = left + (right - left) / 2.0;
        const double fMid = F(mid);

        if (fMid == 0.0) {
            left = right = mid;
            ++iterations;
            break;
        }

        if (std::signbit(fMid) == std::signbit(F(left))) {
            left = mid;
        }
        else {
            right = mid;
        }

        ++iterations;
    }

    const double root = left + (right - left) / 2.0;
    const double error = (right - left) / 2.0;

    return { "Половинного деления", "[1.0; 1.1]", root, iterations, error };
}

// 3. Метод хорд (с сохранением отрезка смены знака)
Result chordMethod(double left, double right) {
    const double minDerivative = dF(left);

    int iterations = 0;
    double x = chordPoint(left, right);
    double error = std::abs(F(x)) / minDerivative;

    while (error > EPSILON && iterations < MAX_ITERATIONS) {
        if (std::signbit(F(x)) == std::signbit(F(left))) {
            left = x;
        }
        else {
            right = x;
        }

        const double next = chordPoint(left, right);
        ++iterations;

        // Защита от остановки из-за точности типа double.
        if (next == x) {
            break;
        }

        x = next;
        error = std::abs(F(x)) / minDerivative;
    }

    error = std::abs(F(x)) / minDerivative;

    return { "Хорд", "[1.0; 1.1]", x, iterations, error };
}

// 4. Метод касательных (Ньютона)
Result newtonMethod(double x) {
    const double minDerivative = dF(1.0);
    int iterations = 0;
    double error = std::abs(F(x)) / minDerivative;

    while (error > EPSILON && iterations < MAX_ITERATIONS) {
        const double derivative = dF(x);

        if (std::abs(derivative) < 1e-14) {
            break;
        }

        const double next = x - F(x) / derivative;
        ++iterations;

        if (next == x) {
            break;
        }

        x = next;
        error = std::abs(F(x)) / minDerivative;
    }

    error = std::abs(F(x)) / minDerivative;

    return { "Касательных (Ньютона)", "x0 = 1.1", x, iterations, error };
}

// 5. Комбинированный метод хорд и касательных
Result combinedMethod(double left, double right) {
    int iterations = 0;

    if (!hasRootBracket(left, right)) {
        return { "Хорд и касательных", "[1.0; 1.1]", std::numeric_limits<double>::quiet_NaN(), 0, std::numeric_limits<double>::quiet_NaN() };
    }

    // На данном отрезке f'(x) > 0 и f''(x) > 0:
    // хорда уточняет левую границу, касательная — правую.
    while ((right - left) / 2.0 > EPSILON && iterations < MAX_ITERATIONS) {
        const double newLeft = chordPoint(left, right);
        const double newRight = right - F(right) / dF(right);

        left = newLeft;
        right = newRight;
        ++iterations;
    }

    const double root = left + (right - left) / 2.0;
    const double error = (right - left) / 2.0;

    return { "Хорд и касательных", "[1.0; 1.1]", root, iterations, error };
}

// 6. Метод простых итераций
Result iterationMethod(double left, double right) {
    // m и M — минимальное и максимальное значения f'(x)
    // на отрезке [left; right]. Здесь f''(x) > 0, поэтому f'(x) возрастает.
    const double m = dF(left);
    const double M = dF(right);

    // Подбираем lambda так, чтобы коэффициент сжатия был минимальным.
    const double lambda = 2.0 / (m + M);

    const double q = std::max(
        std::abs(1.0 - lambda * m),
        std::abs(1.0 - lambda * M)
    );

    // Проверяем, что phi(x) = x - lambda*f(x) отображает отрезок в себя.
    const auto phi = [lambda](double value) {
        return value - lambda * F(value);
        };

    double phiMin = std::min(phi(left), phi(right));
    double phiMax = std::max(phi(left), phi(right));

    // Внутренняя критическая точка phi: phi'(x) = 0, то есть f'(x) = 1/lambda.
    const double discriminant = 0.04 - 12.0 * (0.4 - 1.0 / lambda);
    if (discriminant >= 0.0) {
        const double criticalPoint = (0.2 + std::sqrt(discriminant)) / 6.0;
        if (criticalPoint >= left && criticalPoint <= right) {
            phiMax = std::max(phiMax, phi(criticalPoint));
            phiMin = std::min(phiMin, phi(criticalPoint));
        }
    }

    const bool mapsIntervalToItself = phiMin >= left && phiMax <= right;
    if (q >= 1.0 || !mapsIntervalToItself) {
        return { "Итераций", "x0 = 1.05", std::numeric_limits<double>::quiet_NaN(), 0, std::numeric_limits<double>::quiet_NaN() };
    }

    // x(n+1) = x(n) - lambda * f(x(n))
    double x = (left + right) / 2.0;
    double error = std::numeric_limits<double>::infinity();
    int iterations = 0;

    while (error > EPSILON && iterations < MAX_ITERATIONS) {
        const double next = phi(x);

        // Апостериорная оценка погрешности метода итераций.
        error = q / (1.0 - q) * std::abs(next - x);

        x = next;
        ++iterations;
    }

    return { "Итераций", "x0 = 1.05", x, iterations, error };
}

void printResult(const Result& result) {
    std::cout << result.method << " | "
        << result.start << " | "
        << result.root << " | "
        << result.iterations << " | "
        << result.errorBound << '\n';
}

int main() {
	setlocale(LC_ALL, "ru-RU.UTF-8");

    separateRoot();

    const double left = 1.0;
    const double right = 1.1;

    // Проверка условий сходимости метода простых итераций.
    const double m = dF(left);
    const double M = dF(right);
    const double lambda = 2.0 / (m + M);
    const double q = std::max(
        std::abs(1.0 - lambda * m),
        std::abs(1.0 - lambda * M)
    );
    const auto phi = [lambda](double x) {
        return x - lambda * F(x);
        };
    double phiMin = std::min(phi(left), phi(right));
    double phiMax = std::max(phi(left), phi(right));
    const double phiDiscriminant = 0.04 - 12.0 * (0.4 - 1.0 / lambda);
    if (phiDiscriminant >= 0.0) {
        const double criticalPoint = (0.2 + std::sqrt(phiDiscriminant)) / 6.0;
        if (criticalPoint >= left && criticalPoint <= right) {
            phiMin = std::min(phiMin, phi(criticalPoint));
            phiMax = std::max(phiMax, phi(criticalPoint));
        }
    }

    std::cout << "===== УСЛОВИЯ СХОДИМОСТИ =====\n";
    std::cout << "На [1.0; 1.1]: f'(x) > 0, f''(x) = 6x - 0.2 > 0.\n";
    std::cout << "Для метода Ньютона выбрано x0 = 1.1, поскольку f(x0) * f''(x0) > 0.\n";
    std::cout << "Для метода итераций: m = " << m
        << ", M = " << M
        << ", lambda = " << lambda
        << ", q = " << q << " < 1.\n";
    std::cout << "phi([1.0; 1.1]) лежит примерно в [" << phiMin
        << "; " << phiMax << "], то есть внутри исходного отрезка.\n\n";

    const Result results[] = {
        bisectionMethod(left, right),
        chordMethod(left, right),
        newtonMethod(right),
        combinedMethod(left, right),
        iterationMethod(left, right)
    };

    std::cout << "===== РЕЗУЛЬТАТЫ МЕТОДОВ =====\n";
    std::cout << std::fixed << std::setprecision(10);
    std::cout << "Метод | Начальные данные | Корень xn | Итерации | Оценка ошибки\n";
    std::cout << std::string(95, '-') << '\n';

    for (const Result& result : results) {
        printResult(result);
        std::cout << "   f(xn) = " << F(result.root) << "\n";
    }

    return 0;
}
