#include <iostream>
#include <cmath>
#include <vector>
#include <stdexcept>
using namespace std;

// ---- Вспомогательные функции для работы с векторами и матрицами ----

// Евклидова норма вектора
double norm(const vector<double>& v) {
    double s = 0.0;
    for (double val : v) s += val * val;
    return sqrt(s);
}

// Умножение матрицы на вектор
vector<double> matVecMul(const vector<vector<double>>& A, const vector<double>& v) {
    size_t n = A.size();
    vector<double> res(n, 0.0);
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j)
            res[i] += A[i][j] * v[j];
    return res;
}

// Сложение вектора со скаляром, умноженным на вектор
vector<double> vecAddScaled(const vector<double>& x, double t, vector<double> d) {
    vector<double> y = x;
    for (size_t i = 0; i < y.size(); ++i) y[i] += t * d[i];
    return y;
}

// ---- Функция f(x) и её обёртка для одномерного поиска ----

// Значение функции в точке x
double f_value(const vector<double>& x) {
    return 2 * x[0] * x[0] - 4 * x[0] + x[1] * x[1] - 8 * x[1] + 3;
}

// Обёртка: f(X + t * D) для использования в swann и gold
double f_wrap(vector<double> X, double t, vector<double> D) {
    vector<double> x = vecAddScaled(X, t, D);
    return f_value(x);
}

// ---- Градиент и матрица Гессе (аналитические) ----

vector<double> gradient(const vector<double>& x) {
    return {4 * x[0] - 4, 2 * x[1] - 8};
}

vector<vector<double>> hessian(const vector<double>& x) {
    // H = [[4, 0], [0, 2]]
    return {{4.0, 0.0}, {0.0, 2.0}};
}

// Обратная матрица для 2x2
vector<vector<double>> inverse2x2(const vector<vector<double>>& H) {
    double det = H[0][0] * H[1][1] - H[0][1] * H[1][0];
    if (fabs(det) < 1e-15) throw runtime_error("Matrix is singular");
    vector<vector<double>> inv(2, vector<double>(2));
    inv[0][0] = H[1][1] / det;
    inv[0][1] = -H[0][1] / det;
    inv[1][0] = -H[1][0] / det;
    inv[1][1] = H[0][0] / det;
    return inv;
}

// Проверка положительной определённости матрицы (критерий Сильвестра для 2x2)
bool isPositiveDefinite(const vector<vector<double>>& H) {
    return (H[0][0] > 1e-12) && (H[0][0] * H[1][1] - H[0][1] * H[1][0] > 1e-12);
}

// ---- Одномерный поиск (метод Свенна + золотое сечение) ----
// Скопированы из main.cpp и адаптированы к сигнатуре f_wrap

pair<double,double> swann(double(*f)(vector<double>, double, vector<double>),
                          vector<double> X, vector<double> Di, double x0,
                          const double &t, bool print = false) {
    if(t <= 0) throw invalid_argument("t must be non-negative!");
    int k = 0;
    if(print) cout << "(1.0) x.0 = " << x0 << "; t = " << t << "; k = " << k << ";" << endl;
    if(print) cout << "(2.0) f(x.0-t) = " << f(X, x0-t, Di)
                   << "; f(x.0) = " << f(X, x0, Di)
                   << "; f(x.0+t) = " << f(X, x0+t, Di) << ";" << endl;
    if(f(X, x0-t, Di) >= f(X, x0, Di) && f(X, x0, Di) <= f(X, x0+t, Di)) {
        pair<double,double> ab = make_pair(x0-t, x0+t);
        if(print) cout << "(3.0) -> [a.0; b.0] = [" << x0-t << "; " << x0+t << "]" << endl;
        return ab;
    }
    if(f(X, x0-t, Di) <= f(X, x0, Di) && f(X, x0, Di) >= f(X, x0+t, Di))
        throw invalid_argument("Function is not unimodal!");
    if(print) cout << "(3.0) termination condition not met;" << endl;

    double delta;
    pair<double,double> ab;
    double x1;
    if(f(X, x0-t, Di) >= f(X, x0, Di) && f(X, x0, Di) >= f(X, x0+t, Di)) {
        delta = t;
        ab.first = x0;
        x1 = x0 + t;
        k = 1;
        if(print) cout << "(4.0) delta = " << t << "; a0 = " << x0 << "; x.1 = " << x0+t << "; k = 1;" << endl;
    }
    if(f(X, x0-t, Di) <= f(X, x0, Di) && f(X, x0, Di) <= f(X, x0+t, Di)) {
        delta = -t;
        ab.second = x0;
        x1 = x0 - t;
        k = 1;
        if(print) cout << "(4.0) delta = " << -t << "; b0 = " << x0 << "; x.1 = " << x0-t << "; k = 1;" << endl;
    }

    bool end = false;
    int iter = 0;
    do {
        double x2 = x1 + pow(2, k) * delta;
        if(print) cout << "(5." << iter << ") x." << k+1 << " = " << x2 << ";" << endl;
        if(f(X, x2, Di) < f(X, x1, Di) && delta == t) {
            ab.first = x1;
            k++;
            if(print) cout << "(6." << iter << ") a.0 = " << ab.first << "; k = " << k << ";" << endl;
        }
        if(f(X, x2, Di) < f(X, x1, Di) && delta == -t) {
            ab.second = x1;
            k++;
            if(print) cout << "(6." << iter << ") b.0 = " << ab.second << "; k = " << k << ";" << endl;
        }
        if(f(X, x2, Di) >= f(X, x1, Di)) {
            end = true;
            if(delta == t) ab.second = x2;
            if(delta == -t) ab.first = x2;
            if(print) cout << "(6." << iter << ") [a0; b0] = [" << ab.first << "; " << ab.second << "]." << endl;
        }
        x0 = x1;
        x1 = x2;
        iter++;
    } while(!end);
    return ab;
}

double gold(double(*f)(vector<double>, double, vector<double>),
            vector<double> X, vector<double> Di, pair<double,double> ab,
            double l, bool print = false) {
    if(l <= 0) throw invalid_argument("l must be non-negative!");
    if(print) cout << "(1.0) L0 = [" << ab.first << "; " << ab.second << "]; l = " << l << ";" << endl;
    double k = 0;
    if(print) cout << "(2.0) k = 0;" << endl;
    const double g = (3 - sqrt(5)) / 2;
    double y0 = ab.first + g * (ab.second - ab.first);
    double z0 = ab.first + ab.second - y0;
    double y1, z1, answ;
    if(print) cout << "(3.0) y.0 = " << y0 << "; z.0 = " << z0 << ";" << endl;

    int iter = 0;
    bool end = false;
    do {
        if(print) cout << "(4." << iter << ") f(y." << k << ") = " << f(X, y0, Di)
                       << "; f(z." << k << ") = " << f(X, z0, Di) << ";" << endl;
        if(f(X, y0, Di) <= f(X, z0, Di)) {
            ab.second = z0;
            y1 = ab.first + ab.second - y0;
            z1 = y0;
            if(print) cout << "(5." << iter << ") f(y) <= f(z) -> L0 = [" << ab.first << "; " << ab.second
                           << "]; y." << k+1 << " = " << y1 << "; z." << k+1 << " = " << z1 << ";" << endl;
        } else {
            ab.first = y0;
            y1 = z0;
            z1 = ab.first + ab.second - z0;
            if(print) cout << "(5." << iter << ") f(y) > f(z) -> L0 = [" << ab.first << "; " << ab.second
                           << "]; y." << k+1 << " = " << y1 << "; z." << k+1 << " = " << z1 << ";" << endl;
        }
        double delta = abs(ab.first - ab.second);
        if(print) cout << "(6." << iter << ") delta = " << delta << "; ";
        if(delta <= l) {
            end = true;
            answ = (ab.first + ab.second) / 2;
            if(print) cout << " x* = " << answ << "." << endl;
        } else if(print) cout << endl;
        k++;
        iter++;
        y0 = y1;
        z0 = z1;
    } while(!end);
    return answ;
}

// ---- Основной алгоритм метода Ньютона ----
vector<double> newton_method(double(*f)(vector<double>, double, vector<double>),
                             vector<double> x0,
                             double eps1, double eps2, int M) {
    // Шаг 1: задаём начальные параметры
    int k = 0;
    vector<double> x = x0;
    bool prev_ok = false;   // флаг для двукратного выполнения условий

    while (true) {
        // Шаг 3: вычислить градиент в x^k
        vector<double> grad = gradient(x);
        double grad_norm = norm(grad);

        // Шаг 4: проверка ||grad|| <= eps1
        if (grad_norm <= eps1) {
            return x;
        }

        // Шаг 5: проверка k >= M
        if (k >= M) {
            return x;
        }

        // Шаг 6: вычислить матрицу Гессе H(x^k)
        vector<vector<double>> H = hessian(x);

        // Шаг 7: вычислить обратную матрицу H^{-1}
        vector<vector<double>> H_inv;
        try {
            H_inv = inverse2x2(H);
        } catch (...) {
            H_inv = vector<vector<double>>(2, vector<double>(2, 0.0));
        }

        // Шаг 8: проверить H^{-1} > 0 (эквивалентно положительной определённости H)
        bool posDef = isPositiveDefinite(H);

        vector<double> d(2);
        double t;

        if (posDef) {
            // Шаг 9: d^k = -H^{-1} * grad
            d = matVecMul(H_inv, grad);
            for (int i = 0; i < 2; ++i) d[i] = -d[i];
            t = 1.0;  // шаг по Ньютону
        } else {
            // Шаг 10 (б): d^k = -grad, выбрать t_k из условия уменьшения функции
            d = grad;
            for (int i = 0; i < 2; ++i) d[i] = -d[i];

            // Одномерный поиск (точная минимизация)
            double step0 = 0.1;
            double tol = eps1 / 10.0;
            auto interval = swann(f, x, d, 0.0, step0, false);
            t = gold(f, x, d, interval, tol, false);
        }

        // Шаг 10 (продолжение): x^{k+1} = x^k + t_k * d^k
        vector<double> x_new = vecAddScaled(x, t, d);

        // Шаг 11: проверка условий останова (двукратное выполнение)
        double dx_norm = norm(vecAddScaled(x_new, -1.0, x));
        double df = fabs(f_value(x_new) - f_value(x));
        bool current_ok = (dx_norm < eps2) && (df < eps2);

        if (current_ok && prev_ok) {
            return x_new;
        }

        prev_ok = current_ok;
        x = x_new;
        k++;
    }
}


int main() {
    vector<double> x0 = {0.0, 0.0};

    double eps1 = 1e-6;
    double eps2 = 1e-6;
    int M = 100;
    vector<double> result = newton_method(f_wrap, x0, eps1, eps2, M);

    cout << result[0] << " " << result[1]<<endl;

    return 0;
}
