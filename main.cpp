#include <iostream>
#include <cmath>
#include <vector>
#include <stdexcept>
using namespace std;

//  одномерный поиск
pair<double,double> swann(double(*f)(vector<double>, double, vector<double>),
                          vector<double> X, vector<double> Di, double x0,
                          const double &t, bool print = true) {
    if(t <= 0) throw invalid_argument("t must be non-negative!");
    int k = 0;
    if(print) cout << "(1.0) x.0 = " << x0 << "; t = " << t << "; k = " << k << ";" << endl;
    if(print) cout << "(2.0) f(x.0-t) = " << f(X, x0-t, Di)
                   << "; f(x.0) = " << f(X, x0, Di)
                   << "; f(x.0+t) = " << f(X, x0+t, Di) << ";" << endl;
    if(f(X, x0-t, Di) >= f(X, x0, Di) and f(X, x0, Di) <= f(X, x0+t, Di)) {
        pair<double,double> ab = make_pair(x0-t, x0+t);
        if(print) cout << "(3.0) -> [a.0; b.0] = [" << x0-t << "; " << x0+t << "]" << endl;
        return ab;
    }
    if(f(X, x0-t, Di) <= f(X, x0, Di) and f(X, x0, Di) >= f(X, x0+t, Di))
        throw invalid_argument("Function is not unimodal!");
    if(print) cout << "(3.0) termination condition not met;" << endl;

    double delta;
    pair<double,double> ab;
    double x1;
    if(f(X, x0-t, Di) >= f(X, x0, Di) and f(X, x0, Di) >= f(X, x0+t, Di)) {
        delta = t;
        ab.first = x0;
        x1 = x0 + t;
        k = 1;
        if(print) cout << "(4.0) delta = " << t << "; a0 = " << x0 << "; x.1 = " << x0+t << "; k = 1;" << endl;
    }
    if(f(X, x0-t, Di) <= f(X, x0, Di) and f(X, x0, Di) <= f(X, x0+t, Di)) {
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
        if(f(X, x2, Di) < f(X, x1, Di) and delta == t) {
            ab.first = x1;
            k++;
            if(print) cout << "(6." << iter << ") a.0 = " << ab.first << "; k = " << k << ";" << endl;
        }
        if(f(X, x2, Di) < f(X, x1, Di) and delta == -t) {
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
            double l, bool print = true) {
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

//определитель
double determinant(vector<vector<double>>& mat) {
    int n = mat.size();
    if(n == 0) return 1.0;
    for(auto& row : mat) if(row.size() != (size_t)n) return 0.0;
    vector<vector<double>> a = mat;
    double det = 1.0;
    for(int i = 0; i < n; ++i) {
        int pivot = i;
        double max_val = fabs(a[i][i]);
        for(int j = i+1; j < n; ++j)
            if(fabs(a[j][i]) > max_val) { max_val = fabs(a[j][i]); pivot = j; }
        if(max_val < 1e-12) return 0.0;
        if(pivot != i) { swap(a[i], a[pivot]); det = -det; }
        det *= a[i][i];
        for(int j = i+1; j < n; ++j) {
            double factor = a[j][i] / a[i][i];
            for(int k = i+1; k < n; ++k) a[j][k] -= factor * a[i][k];
        }
    }
    return det;
}

//сравнение векторов
bool vec_equal(const vector<double>& a, const vector<double>& b, double eps) {
    double norm = 0.0;
    for(size_t i = 0; i < a.size(); ++i) norm += (a[i]-b[i])*(a[i]-b[i]);
    return sqrt(norm) < eps;
}

// (2x1^2 -4x1 + x2^2 -8x2 +3)
double f(vector<double> x, double ti, vector<double> di) {
    return 2*pow(x[0]+ti*di[0],2) - 4*(x[0]+ti*di[0]) +
           pow(x[1]+ti*di[1],2) - 8*(x[1]+ti*di[1]) + 3;
}


//  метод сопряжённых направлений Пауэлла
vector<double> powell(double(*f)(vector<double>, double, vector<double>),
                      vector<double> x0, double e = 0.1) {
    int n = x0.size();
    // Шаг 1: задать начальные направления поиска (единичные векторы)
    vector<vector<double>> d(n+1, vector<double>(n, 0.0));
    for(int i = 1; i <= n; ++i) d[i][i-1] = 1.0;   // d1..dn
    d[0] = d[n];                 // d0 = dn
    double eps1 = e / 10.0;      // точность одномерного поиска
    vector<double> y0 = x0;      // y^0 = x^0
    int k = 0;                   // номер итерации
    vector<double> x_prev;       // для проверки |x^{k+1} - x^k| < eps

    while(true) {
        // Шаг 2: циклический поиск по направлениям d0, d1, ..., dn
        vector<vector<double>> y(1, y0);   // y[0] = y^0
        int i = 0;
        // поиск по d0, d1, ..., d_{n-1} (всего n шагов)
        while(i <= n-1) {
            auto interval = swann(f, y[i], d[i], 0.0, 0.1, false);
            double ti = gold(f, y[i], d[i], interval, eps1, false);
            vector<double> y_next(n);
            for(int j = 0; j < n; ++j)
                y_next[j] = y[i][j] + ti * d[i][j];
            y.push_back(y_next);
            // Шаг 3: проверка после i = n-1
            if(i == n-1) {
                if(vec_equal(y[n], y[0], e)) {
                    return y[n];          // y^n == y^0 -> минимум
                } else {
                    i++;                  // i = n, переходим к поиску по d_n
                    continue;
                }
            }
            i++;
        }
        // поиск по d_n (i == n)
        auto interval = swann(f, y[n], d[n], 0.0, 0.1, false);
        double tn = gold(f, y[n], d[n], interval, eps1, false);
        vector<double> y_next(n);
        for(int j = 0; j < n; ++j)
            y_next[j] = y[n][j] + tn * d[n][j];
        y.push_back(y_next);               // y[n+1]

        // Шаг 3 (продолжение): проверка y^{n+1} == y^n (а не y^1!)
        if(vec_equal(y[n+1], y[n], e))
            return y[n+1];

        // Шаг 4: обновление направления
        vector<double> x_new = y[n+1];
        // Проверка |x^{k+1} - x^k| < eps
        if(k > 0) {
            if(vec_equal(x_new, x_prev, e))
                return x_new;
        }

        // новое направление: d0 = dn = y^{n+1} - y^1
        vector<double> new_dir(n);
        for(int j = 0; j < n; ++j)
            new_dir[j] = y[n+1][j] - y[1][j];

        vector<vector<double>> old_d = d;   // сохраняем старые направления

        d[0] = new_dir;
        d[n] = new_dir;
        // сдвиг d1..d_{n-1} = старые d2..dn
        for(int i = 1; i <= n-1; ++i)
            d[i] = old_d[i+1];

        // проверка линейной независимости новой системы (d1..dn)
        vector<vector<double>> mat(n, vector<double>(n));
        for(int i = 1; i <= n; ++i)
            mat[i-1] = d[i];
        if(fabs(determinant(mat)) < e) {
            // ранг < n -> оставляем старые направления
            d = old_d;
        }

        // подготовка к следующей итерации
        x_prev = x_new;
        y0 = x_new;
        k++;
    }
}

int main() {
    vector<double> x = {0, 0};
    auto r = powell(f, x, 0.0001);
    cout << r[0] << " " << r[1] << endl;
    return 0;
}
