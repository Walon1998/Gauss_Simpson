#include <cmath>
#include <iostream>
#include <functional>
#include <vector>

using namespace std;

double GaussLegendre5(const std::function<double(double)> &f, double a, double b) {
    const std::vector<double> c = {-0.90617984593, -0.53846931010, 0.0, 0.53846931010, 0.90617984593};
    const std::vector<double> w = {0.23692688505, 0.47862867049, 0.56888888888, 0.47862867049, 0.23692688505};

    double result = 0;
    double factor = (b - a) / 2;
    for (int i = 0; i < 5; ++i) {

        double sum1 = ((b - a) / 2) * c[i];
        double sum2 = ((b + a) / 2);
        double fvalue = f(sum1 + sum2) * w[i];
        result += fvalue;

    }


    return factor * result;
}

double CompositeSimpson(const std::function<double(double)> &f, const std::vector<double> &x) {
//    int n = x.size();
//    double result = 0;
//    double a = x[0];
//    double b = x[n - 1];
//
//    double h = (b - a) / (2 * n);
//
//    double sum1 = 0;
//    double sum2 = 0;
//
//    for (int i = 0; i < n; ++i) {
//        sum1 += f(a + 2 * h * i);
//        sum2 += f(a + 2 * h * (2 * i - 1));
//    }
//
//
//    result = (h / 3) * (f(a) + f(b) + sum1 + sum2);
//
//    return result;


    long n = x.size();
    double result = 0;
    double a = x[0];
    double b = x[n - 1];
//
//    double h = (b - a) / (2 * n);
//
    double sum1 = (1.0 / 6.0) * (x[1] - x[0]) * f(a);
    double sum2 = (1.0 / 6.0) * (x[n - 1] - x[n - 2]) * f(b);

    double sum3 = 0;
    double sum4 = 0;



    for (int i = 1; i < n - 1; ++i) {
        sum3 += (1.0 / 6.0) * (x[i + 1] - x[i - 1]) * f(x[i]);
    }
    for (int i = 1; i < n; ++i) {
        sum4 += (2.0 / 3.0) * (x[i] - x[i - 1]) * f((1.0 / 2.0) * (x[i] + x[i - 1]));
    }


    result += sum1;
    result += sum2;
    result += sum3;
    result += sum4;
    return result;
}

std::vector<double> LinSpace(int n, double a, double b) {
    std::vector<double> x(n);
    double d = (b - a) / (n - 1);
    for (int i = 0; i < n; ++i) {
        x[i] = i * d + a;
    }
    return x;
}

double f(double x) {
    return std::sqrt(x);
}

double F(double x) {
    double y = std::sqrt(x);
    return 2.0 / 3.0 * y * y * y;
}

int main() {
    int n = 5;
    double a = 0.0;
    double b = 1.0;

    std::cout << "Gauss-Legendre: " << GaussLegendre5(f, a, b) << std::endl;
    std::cout << "Simpson: " << CompositeSimpson(f, LinSpace(n, a, b)) << std::endl;
    std::cout << "Exact value: " << F(b) - F(a) << std::endl;

    return 0;
}

