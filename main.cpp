#include <iostream>
#include <math.h>
using namespace std;

pair<double,double> swann(double x0, const double &t, bool print = true){
    if(t<=0)
        throw invalid_argument("t must be non-negative!");

    //1 пункт
    int k = 0;
    if(print)
        cout<<"(1.0) x.0 = "<<x0<<"; t = "<<t<<"; k = "<<k<<";"<<endl;

    //2 пункт
    if(print)
        cout<<"(2.0) f(x.0-t) = "<<(10*(x0-t)*log(x0-t) - pow(x0-t,2)/2.0)<<
              "; f(x.0) = "<<(10*(x0)*log(x0) - pow(x0,2)/2.0)<<
              "; f(x.0+t) = "<<(10*(x0+t)*log(x0+t) - pow(x0+t,2)/2.0)<<";"<<endl;

    //3 пункт
    if((10*(x0-t)*log(x0-t) - pow(x0-t,2)/2.0)>=(10*(x0)*log(x0) - pow(x0,2)/2.0) and (10*(x0)*log(x0) - pow(x0,2)/2.0)<=(10*(x0+t)*log(x0+t) - pow(x0+t,2)/2.0)){
        pair<double,double> ab = make_pair(x0-t, x0+t);
        if(print)
            cout<<"(3.0) f(x.0-t)>=f(x.0)<=f(x.0+t) -> [a.0; b.0] = ["<<x0-t<<"; "<<x0+t<<"]"<<endl;
        return ab;}
    if((10*(x0-t)*log(x0-t) - pow(x0-t,2)/2.0)>=(10*(x0)*log(x0) - pow(x0,2)/2.0) and (10*(x0)*log(x0) - pow(x0,2)/2.0)<=(10*(x0+t)*log(x0+t) - pow(x0+t,2)/2.0))
        throw invalid_argument("The function is not unimodal, so it is recommended to set a different starting point!");
    if(print)
        cout<<"(3.0) the termination condition is not met;"<<endl;

    //4 пункт
    double delta;
    pair<double,double> ab = make_pair(0,0);
    double x1;
    if((10*(x0-t)*log(x0-t) - pow(x0-t,2)/2.0)>=(10*(x0)*log(x0) - pow(x0,2)/2.0) and (10*(x0)*log(x0) - pow(x0,2)/2.0)>=(10*(x0+t)*log(x0+t) - pow(x0+t,2)/2.0)){
        delta = t;
        ab.first = x0;
        x1 = x0+t;
        k = 1;
        if(print)
            cout<<"(4.0) f(x.0-t)>=f(x.0)>=f(x.0+t) -> delta = "<<t<<"; a0 = "<<x0<<"; x.1 = "<<x0+t<<"; k = "<<1<<";"<<endl;
    }
    if((10*(x0-t)*log(x0-t) - pow(x0-t,2)/2.0)<=(10*(x0)*log(x0) - pow(x0,2)/2.0) and (10*(x0)*log(x0) - pow(x0,2)/2.0)<=(10*(x0+t)*log(x0+t) - pow(x0+t,2)/2.0)){
        delta = -t;
        ab.second = x0;
        x1 = x0 - t;
        k = 1;
        if(print)
            cout<<"(4.0) f(x.0-t)<=f(x.0)<=f(x.0+t) -> delta = "<<-t<<"; b0 = "<<x0<<"; x.1 = "<<x0-t<<"; k = "<<1<<";"<<endl;
    }

    //5 и 6 пункты
    bool end = false;
    int iterations = 0;
    do{
        //5 пункт
        double x2 = x1 + pow(2,k)*delta;
        if(print)
            cout<<"(5."<<iterations<<") x."<<k+1<<" = "<<x2<<";"<<endl;

        //6 пункт
        if((10*(x2)*log(x2) - pow(x2,2)/2.0)<(10*(x1)*log(x1) - pow(x1,2)/2.0) and delta == t){
            ab.first = x1;
            k = k+1;
            if(print)
                cout<<"(6."<<iterations<<") a.0 = "<<x0<<"; k = "<<k<<";"<<endl;
        }
        if((10*(x2)*log(x2) - pow(x2,2)/2.0)<(10*(x1)*log(x1) - pow(x1,2)/2.0) and delta == -t){
            ab.second = x0;
            k = k+1;
            if(print)
                cout<<"(6."<<iterations<<") b.0 = "<<x0<<"; k = "<<k<<";"<<endl;
        }
        if((10*(x2)*log(x2) - pow(x2,2)/2.0)>=(10*(x1)*log(x1) - pow(x1,2)/2.0)){
            end = true;
            if(delta == t)
                ab.second = x0;
            if(delta == -t)
                ab.first = x2;
            if(print)
                cout<<"(6."<<iterations<<") [a0; b0] = ["<<ab.first<<"; "<<ab.second<<"]."<<endl;
        }
        x0 = x1;
        x1 = x2;
        iterations++;
    }while(!end);

    return ab;
}

double gold(pair<double,double> ab, double l, bool print = true){
    if(l<=0)
        throw invalid_argument("l must be non-negative!");

    //1 пункт
    if(print)
        cout<<"(1.0) L0 = ["<<ab.first<<"; "<<ab.second<<"]; l = "<<l<<";"<<endl;

    //2 пункт
    double k = 0;
    if(print)
        cout<<"(2.0) k = 0;"<<endl;

    //3 пункт
    const double g = (3-sqrt(5))/2;
    double y0 = ab.first + g*(ab.second - ab.first);
    double z0 = ab.first + ab.second - y0;
    double y1;
    double z1;
    double answ;
    if(print)
        cout<<"(3.0) y.0 = "<<y0<<"; z.0 = "<<z0<<";"<<endl;

    //4-6 пункты
    int iterations = 0;
    bool end = false;
    do{
        //4 пункт
        if(print)
            cout<<"(4."<<iterations<<") f(y."<<k<<") = "<<(10*(y0)*log(y0) - pow(y0,2)/2.0)<<
                  "; f(z."<<k<<") = "<<(10*(z0)*log(z0) - pow(z0,2)/2.0)<<";"<<endl;

        //5 пункт
        if((10*(y0)*log(y0) - pow(y0,2)/2.0)<=(10*(z0)*log(z0) - pow(z0,2)/2.0)){
            ab.first = ab.first;
            ab.second = z0;
            y1 = ab.first + ab.second - y0;
            z1 = y0;
            if(print)
                cout<<"(5."<<iterations<<") f("<<y0<<") <= f("<<z0<<") -> L0 = ["<<ab.first<<"; "<<ab.second<<"]; y."<<k+1<<" = "<<y1<<"; z."<<k+1<<" = "<<z1<<";"<<endl;
        }
        if((10*(y0)*log(y0) - pow(y0,2)/2.0)>(10*(z0)*log(z0) - pow(z0,2)/2.0)){
            ab.first = y0;
            ab.second = ab.second;
            y1 = z0;
            z1 = ab.first + ab.second - z0;
            if(print)
                cout<<"(5."<<iterations<<") f("<<y0<<") > f("<<z0<<") -> L0 = ["<<ab.first<<"; "<<ab.second<<"]; y."<<k+1<<" = "<<y1<<"; z."<<k+1<<" = "<<z1<<";"<<endl;
        }

        //6 пункт
        double delta = abs(ab.first - ab.second);
        if(print)
            cout<<"(6."<<iterations<<") delta = "<<delta<<";";
        if(delta <= l){
            end = true;
            answ = (ab.first + ab.second)/2;
            if(print)
                cout<<" x* = "<<answ<<"."<<endl;
        } else cout<<endl;
        k = k + 1;
        iterations++;
        y0 = y1;
        z0 = z1;
    }while(!end);

    return answ;
}

double quadInt(double x1, double dx, double e1, double e2, bool print = true){
    if(dx<=0)
        throw invalid_argument("dx must be non-negative!");
    if(e1<=0)
        throw invalid_argument("e1 must be non-negative!");
    if(e2<=0)
        throw invalid_argument("e2 must be non-negative!");

    //1 пункт
    if(print)
        cout<<"(1.0) x1 = "<<x1<<"; dx = "<<dx<<"; e1 = "<<e1<<"; e2 = "<<e2<<";"<<endl;

    //2-8 пункты
    bool end = false;
    bool pass2_5 = false;
    bool pass8 = false;
    double answ, x2, x3, f1, f2, f3, Fmin, Xmin, x_, f_;
    int iterantion = 0;
    do {
        if(!pass2_5){
        //2 пункт
        x2 = x1+dx;
        if(print)
            cout<<"(2."<<iterantion<<") x2 = "<<x2<<";"<<endl;
        }

        if(!pass2_5){
        //3 пункт
        f1 = 10*(x1)*log(x1) - pow(x1,2)/2.0;
        f2 = 10*(x2)*log(x2) - pow(x2,2)/2.0;
        if(print)
            cout<<"(3."<<iterantion<<") f1 = "<<f1<<"; f2 = "<<f2<<";"<<endl;
        }

        if(!pass2_5){
        //4 пункт
        if(f1>f2)
            x3 = x1 + 2*dx;
        else
            x3 = x1 - dx;
        if(print)
            cout<<"(4."<<iterantion<<") f1"<<(f1>f2?">":"<=")<<"f2 -> x3 = "<<x3<<";"<<endl;
        }

        if(!pass2_5){
        //5 пункт
        f3 = 10*(x3)*log(x3) - pow(x3,2)/2.0;
        if(print)
            cout<<"(5."<<iterantion<<") f3 = "<<f3<<";"<<endl;
        }

        //6 пункт
        pass2_5 = false;
        Fmin = min(f1,min(f2,f3));
        Xmin = Fmin==f1?x1:Fmin==f2?x2:x3;
        if(print)
            cout<<"(6."<<iterantion<<") Fmin = "<<Fmin<<"; Xmin = "<<Xmin<<";"<<endl;

        //7 пункт
        pass8 = false;
        double temp_chislitel = (x2-x3)*f1+(x3-x1)*f2+(x1-x2)*f3;
        if(temp_chislitel == 0){
            pass8 = true;
            x1 = Xmin;
        } else{
            x_ = 0.5 * ((pow(x2,2)-pow(x3,2))*f1+(pow(x3,2)-pow(x1,2))*f2+(pow(x1,2)-pow(x2,2))*f3)/temp_chislitel;
            f_ = 10*(x_)*log(x_) - pow(x_,2)/2.0;
        }
        if(!pass8)
            cout<<"(7."<<iterantion<<") x_ = "<<x_<<"; f_ = "<<f_<<";"<<endl;
        else
            cout<<"(7."<<iterantion<<") x1 = "<<Xmin<<";"<<endl;

        //8 пункт
        if(!pass8) {
            if(abs((Fmin - f_)/f_) < e1 and abs((Xmin - x_)/x_) < e2) {
                answ = x_;
                end = true;
            } else {
                double left = min(x1, x3);
                double right = max(x1, x3);
                if(x_ >= left and x_ <= right) {
                    double xx[4] = {x1, x2, x3, x_};
                    double ff[4] = {f1, f2, f3, f_};
                    int imin = 0;
                    for(int i = 1; i < 4; ++i) if(ff[i] < ff[imin]) imin = i;
                    double xl = xx[imin], xr = xx[imin];
                    for(int i = 0; i < 4; ++i) {
                        if(i == imin) continue;
                        if(xx[i] < xx[imin] and (xl == xx[imin] or xx[i] > xl)) xl = xx[i];
                        if(xx[i] > xx[imin] and (xr == xx[imin] or xx[i] < xr)) xr = xx[i];
                    }
                    x1 = xl; f1 = 10*x1*log(x1) - x1*x1/2.0;
                    x2 = xx[imin]; f2 = ff[imin];
                    x3 = xr; f3 = 10*x3*log(x3) - x3*x3/2.0;
                    pass2_5 = true;
                } else {
                    x1 = x_;
                    pass2_5 = false;
                }
            }
            if(print) cout << "(8." << iterantion << ") x1 = " << x1 << "; x2 = " << x2 << "; x3 = " << x3 << ";\n";
        }
    iterantion++;
    } while (!end);

    return answ;
}


int main()
{
    cout<<"Swann's Method: "<<endl;
    auto ab = swann(0.5,0.1);
    cout<<"["<<ab.first<<"; "<<ab.second<<"]"<<endl;
    cout<<endl<<"The Golden Ratio: "<<endl;
    double answ1 = gold(ab,0.05);
    cout<<"f("<<answ1<<") = "<<10*(answ1)*log(answ1) - pow(answ1,2)/2.0<<endl;
    cout<<endl<<"Quadratic Interpolation Method: "<<endl;
    double answ2 = quadInt(0.5,0.2,0.1,0.1);
    cout<<"f("<<answ2<<") = "<<10*(answ2)*log(answ2) - pow(answ2,2)/2.0<<endl;
    return 0;
}
