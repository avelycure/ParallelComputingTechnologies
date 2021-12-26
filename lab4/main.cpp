#include <iostream>
#include <math.h>
#include <vector>

const double G = 6.67e-11;

double generate_random(double a=0, double b=1)
{
    return ((double) rand() / (RAND_MAX)) * (b - a) + a;
}

struct point
{
    double x;
    double y;
    double z;

    point()
    {
        x = generate_random(-100, 100);
        y = generate_random(-100, 100);
        z = generate_random(-100, 100);
    }
    point(double x_, double y_, double z_)
    {
        x = x_;
        y = y_;
        z = z_;
    }
};

double point_module(point p1)
{
    return sqrt(p1.x*p1.x + p1.y*p1.y + p1.z+p1.z);
}

point point_substraction(point p1, point p2)
{
    return {p1.x-p2.x, p1.y-p2.y, p1.z-p2.z};
}

point point_sum(point p1, point p2)
{
    return {p1.x+p2.x, p1.y+p2.y, p1.z+p2.z};
}

point point_by_double(point p, double d)
{
    return {p.x * d, p.y * d, p.z * d};
}

struct left_ // material point
{
    point r;
    point v;

    double m;

    left_()
    {
        r = point();
        v = point();
    }
};

struct right_ // abstract velocity-acceleration object
{
    point v;
    point a;

    right_()
    {
        v = point();
        a = point();
    }

    void null_a()
    {
        a.x = 0;
        a.y = 0;
        a.z = 0;
    }
};

class NBody
{
public:
    std::vector<left_> bodys;
    double t;
    unsigned int iter;
    bool test_flag;

    NBody(size_t, double, unsigned int, bool);
    ~NBody();

    void step();
    void solve();
};

NBody::NBody(size_t n = 4, double t_ = 0.001, unsigned int iter_ = 100, bool test_flag_=false)
{
    t = t_;
    iter = iter_;
    test_flag = test_flag_;

    if(test_flag)
    {
        std::cout << "test\n";

        bodys.resize(4);
        // 8810324116.227   1.0  0.0  0.0   0.0  0.9  0.0
        bodys[0].m = 8810324116.227;
        bodys[0].r.x = 1;
        bodys[0].r.y = 0;
        bodys[0].r.z = 0;
        bodys[0].v.x = 0;
        bodys[0].v.y = 0.9;
        bodys[0].v.z = 0;

        bodys[1].m = 8810324116.227;
        bodys[1].r.x = 0;
        bodys[1].r.y = 1;
        bodys[1].r.z = 0;
        bodys[1].v.x = -0.9;
        bodys[1].v.y = 0;
        bodys[1].v.z = 0;

        bodys[2].m = 8810324116.227;
        bodys[2].r.x = -1;
        bodys[2].r.y = 0;
        bodys[2].r.z = 0;
        bodys[2].v.x = 0;
        bodys[2].v.y = -0.9;
        bodys[2].v.z = 0;

        bodys[3].m = 8810324116.227;
        bodys[3].r.x = 0;
        bodys[3].r.y = -1;
        bodys[3].r.z = 0;
        bodys[3].v.x = 0.9;
        bodys[3].v.y = 0;
        bodys[3].v.z = 0;
    }
    else
    {
        std::cout << "no test\n";

        bodys.resize(n);
        for(int i = 0; i < n; ++i)
        {
            bodys[n] = left_();
        }
    }
}
NBody::~NBody()
{}

void NBody::step()
{
    right_ temp;
    point res;
    double multiplier;

    point dr, dv;
    //calculate a
    for(int i = 0; i < bodys.size(); ++i)
    {
        
        temp.null_a();

        for(int j = 0; j < bodys.size(); ++j)
        {
            temp.null_a();
            if(i != j)
            {
                res = point_substraction(bodys[i].r, bodys[j].r);
                multiplier = G * bodys[j].m / pow(point_module(res), 3);
                res = point_by_double(res, multiplier);
                temp.a = point_substraction(temp.a, res);
            }
            dr = point_by_double(bodys[i].v, t);
            dv = point_by_double(temp.a, t);
            bodys[i].r = point_sum(bodys[i].r, dr);
            bodys[i].v = point_sum(bodys[i].v, dv);
            if(test_flag && i == 0)
            {
                std::cout << bodys[i].r.x << " " << bodys[i].r.y << " " << bodys[i].r.z << "\n";
            }
        }
    }
}

void NBody::solve()
{
    for(int i = 0; i < iter; ++i)
    {
        step();
    }
}

int main()
{
    srand(42);
    NBody sas(4, 0.1, 100, true);
    sas.solve();

    return 0;
}