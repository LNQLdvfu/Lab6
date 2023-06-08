#include <iostream>
#include <vector>

using namespace std;

void X(vector<float>&a,float b, float c, int n);
void Fx(vector<float>&a, vector<float>& b);
float FPO(float x) { return 2*x + M_PI * sin(x);}
float FPN(float x) { return 2*x + M_PI * sin(x);}
float Cubic_Spline(vector<float>&x, vector<float>&y, float fbo, float fpn, int x_x, float cal_x);
int main() {
    const float a1{0.1};
    const float b1{0.6};
    int n[]{3,5, 10, 20, 30, 40, 50, 60, 70, 80, 90,100};

    float calculate_x{};
    cout<<"Input x:" ;
    cin>>calculate_x;
    for(auto & i : n){
        int i_th = static_cast<int>((calculate_x - a1)/((b1-a1)/i));
        vector<float>x(i+1);
        vector<float>y(i+1);
        X(x,a1,b1,i);
        Fx(x,y);
        cout<<Cubic_Spline(x,y, FPO(x[0]), FPO(x[i]), i_th, calculate_x)<<endl;
    }

}
void X(vector<float>&a,float b, float c, int n){
    float h = (c-b)/n;
    for(int i{}; i<=n; i++){
        a[i] = b + i*h;
    }
}
void Fx(vector<float>&a, vector<float>& b) {
    for(int i{}; i <= b.size(); i++) {
        b[i] = powf(a[i],2) - cos(M_PI * a[i]);
    }
}
float  Cubic_Spline(vector<float>&x, vector<float>&y, float fbo, float fpn, int x_x,float cal_x)
{
    vector<float>h(y.size());
    for(int i{}; i< y.size(); i ++){
        h[i] = x[i+1] - x[i] ;
    }
    vector<float>alpha(y.size());
     alpha[0] = 3*(y[1] - y[0]) / h[0] -3*fbo;
     alpha[y.size()] = 3 * fpn -(y[y.size()] -y[y.size() -1]) /h[y.size() -1];
    for(int i{1};i<y.size();i++)
    {
        alpha[i] = 3 / h[i] * ( y[i+1] - y[i]) -3/ h[i -1] *(y[i] - y[i-1]);
    }

    vector<float>l(y.size());
    vector<float>u(y.size());
    vector<float>z(y.size());
    l[0] = 2*h[0];
    u[0] = 0.5;
    z[0] = alpha[0] / l[0];

    for(int i{1} ;i < y.size(); i++) {
        l[i] = 2 * (x[i+1] - x[i -1]) - h[i-1]*u[i-1];
        u[i] = h[i]/ l[i];
        z[i] = (alpha[i] - h[i-1]*z[i-1])/l[i];
    }
    vector<float>b(y.size());
    vector<float>c(y.size());
    vector<float>d(y.size());

    l[y.size()] = h[y.size() -1]* (2 - u[y.size()] -1);
    z[y.size()] = (alpha[y.size()] -h[y.size() -1]*z[y.size() -1]) / l[y.size()];
    c[y.size()] =z[y.size()];

    for(unsigned long i{y.size() -1}; i<y.size() && i >= 0;i--){
        c[i] = z[i] - u[i]*c[i+1];
        b[i] = (y[i+1] -y[i]) / h[i] -h[i]*(c[i+1] +2 * c[i]) / 3;
        d[i] = (c[i+1] - c[i])/(3*h[i]);
    }
    return(y[x_x] + b[x_x]*(cal_x -x[x_x])+c[x_x]* powf(cal_x -x[x_x],2) +d[x_x]* powf(cal_x-x[x_x] ,3));
}
