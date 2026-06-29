#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>
using namespace std;

//I (целевая функция)
double f1(vector<double> x){
    return pow(x[0],2)+pow(x[1],2)+x[2];
}

//II (целевая функция)
double f2(vector<double> x){
    return pow(x[0],2)+pow(x[1],2);
}

//I (вспомогательная функция)
double F1(vector<double> x, double r){
    return pow(x[0],2)+pow(x[1],2)+x[2]+ //f(x)
            (r/2)* //(r^k)/2
            ((pow(x[0]+x[1]+x[2]-4,2)+pow(2*x[0]-3*x[1]-12,2))+ //сумма из квадратов g
            0); //сумма из квадратов g+
}

//II (вспомогательная функция)
double F2(vector<double> x, double r){
    return pow(x[0],2)+pow(x[1],2)+ //f(x)
            (r/2)* //(r^k)/2
            (0+ //сумма из квадратов g
            (pow(fmax(0,-x[0]-x[1]+4),2)+pow(fmax(0,x[0]+2*x[1]-8),2))); //сумма из квадратов g+
}

//обёртка для вспомогательных функций для шага 3
double wrapF(double(*F)(vector<double>, double), vector<double> x, double r, double ti, vector<double> di){
    vector<double> x_new;
    int s = x.size();
    for(int i = 0; i < s; i++)
        x_new.push_back(x[i]+ti*di[i]);

    return F(x_new, r);
}

//Дальше кусок из task2_1------------(переработанный для шага 3 с учётом обёртки для вспомогательной функции)
//  одномерный поиск
pair<double,double> swann(double(*F)(vector<double>, double), double r,
                          vector<double> X, vector<double> Di, double x0,
                          const double &t, bool print = true) {
    if(t <= 0) throw invalid_argument("t must be non-negative!");
    int k_ = 0;
    if(print) cout << "(1.0) x.0 = " << x0 << "; t = " << t << "; k = " << k_ << ";" << endl;
    if(print) cout << "(2.0) f(x.0-t) = " << wrapF(F, X, r, x0-t, Di)
                   << "; f(x.0) = " << wrapF(F, X, r, x0, Di)
                   << "; f(x.0+t) = " << wrapF(F, X, r, x0+t, Di) << ";" << endl;
    if(wrapF(F, X, r, x0-t, Di) >= wrapF(F, X, r, x0, Di) and wrapF(F, X, r, x0, Di) <= wrapF(F, X, r, x0+t, Di)) {
        pair<double,double> ab = make_pair(x0-t, x0+t);
        if(print) cout << "(3.0) -> [a.0; b.0] = [" << x0-t << "; " << x0+t << "]" << endl;
        return ab;
    }
    if(wrapF(F, X, r, x0-t, Di) <= wrapF(F, X, r, x0, Di) and wrapF(F, X, r, x0, Di) >= wrapF(F, X, r, x0+t, Di))
        throw invalid_argument("Function is not unimodal!");
    if(print) cout << "(3.0) termination condition not met;" << endl;

    double delta;
    pair<double,double> ab;
    double x1;
    if(wrapF(F, X, r, x0-t, Di) >= wrapF(F, X, r, x0, Di) and wrapF(F, X, r, x0, Di) >= wrapF(F, X, r, x0+t, Di)) {
        delta = t;
        ab.first = x0;
        x1 = x0 + t;
        k_ = 1;
        if(print) cout << "(4.0) delta = " << t << "; a0 = " << x0 << "; x.1 = " << x0+t << "; k = 1;" << endl;
    }
    if(wrapF(F, X, r, x0-t, Di) <= wrapF(F, X, r, x0, Di) and wrapF(F, X, r, x0, Di) <= wrapF(F, X, r, x0+t, Di)) {
        delta = -t;
        ab.second = x0;
        x1 = x0 - t;
        k_ = 1;
        if(print) cout << "(4.0) delta = " << -t << "; b0 = " << x0 << "; x.1 = " << x0-t << "; k = 1;" << endl;
    }

    bool end = false;
    int iter = 0;
    do {
        double x2 = x1 + pow(2, k_) * delta;
        if(print) cout << "(5." << iter << ") x." << k_+1 << " = " << x2 << ";" << endl;
        if(wrapF(F, X, r, x2, Di) < wrapF(F, X, r, x1, Di) and delta == t) {
            ab.first = x1;
            k_++;
            if(print) cout << "(6." << iter << ") a.0 = " << ab.first << "; k = " << k_ << ";" << endl;
        }
        if(wrapF(F, X, r, x2, Di) < wrapF(F, X, r, x1, Di) and delta == -t) {
            ab.second = x1;
            k_++;
            if(print) cout << "(6." << iter << ") b.0 = " << ab.second << "; k = " << k_ << ";" << endl;
        }
        if(wrapF(F, X, r, x2, Di) >= wrapF(F, X, r, x1, Di)) {
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

double gold(double(*F)(vector<double>, double), double r,
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
        if(print) cout << "(4." << iter << ") f(y." << k << ") = " << wrapF(F, X, r, y0, Di)
                       << "; f(z." << k << ") = " << wrapF(F, X, r, z0, Di) << ";" << endl;
        if(wrapF(F, X, r, y0, Di) <= wrapF(F, X, r, z0, Di)) {
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

//  метод сопряжённых направлений Пауэлла
vector<double> powell(double(*F)(vector<double>, double), double r,
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
            auto interval = swann(F, r, y[i], d[i], 0.0, 0.1, false);
            double ti = gold(F, r, y[i], d[i], interval, eps1, false);
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
        auto interval = swann(F, r, y[n], d[n], 0.0, 0.1, false);
        double tn = gold(F, r, y[n], d[n], interval, eps1, false);
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
//Выше кусок из task2_1------------(переработанный для шага 3 с учётом обёртки для вспомогательной функции)

//Метод Штрафов
vector<double> penalty(double(*F)(vector<double>, double), double(*f)(vector<double>), vector<double> x0, double r, double C, double e){
    //Шаг 1
    int k = 0;

    do{
    //Шаг 2 (Вспомогательная функция уже задана и передана)

    //Шаг 3
    vector<double> x_Answ = powell(F, r, x0, e);
    double PENALTY = F(x_Answ,r) - f(x_Answ);
    //Шаг 4
    //a)
    if(PENALTY<=e)
        return x_Answ;

    //б)
    r = C*r;
    x0 = x_Answ;
    k=k+1;
    } while(true);
}


int main()
{
    vector<double> x01(3,0);
    auto x = penalty(F1, f1, x01, 1, 10, 0.00001);
    cout << x[0]<<" "<<x[1]<<" "<<x[2]<<endl;
    cout<<f1(x)<<endl<<endl;

    vector<double> x02(2,0);
    x = penalty(F2, f2, x02, 1, 10, 0.00001);
    cout << x[0]<<" "<<x[1]<< endl;
    cout<<f2(x)<<endl<<endl;
    return 0;
}
