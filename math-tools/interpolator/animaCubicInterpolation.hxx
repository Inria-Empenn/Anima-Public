#pragma once
#include "animaCubicInterpolation.h"

namespace anima
{

template <class T>
T
Cubic (T x0, T x1, T x2, T x3, T y0, T y1, T y2, T y3, T x)
{
    double p0 = y1;
    double p1 = y2;

    double m0 = (y2-y0)/(x2-x0);
    double m1 = (y3-y1)/(x3-x1);

    double t = (x - x1)/(x2 - x1);
    double t2 = t*t;
    double t3 = t2*t;
    double h00 = 2*t3 -3*t2 +1;
    double h10 = t3 -2*t2 +t;
    double h01 = -2*t3 + 3*t2;
    double h11 = t3 - t2;

    return h00*p0 + h10*m0*(x2 - x1) + h01*p1 + h11*m1*(x2 - x1);
}

template <class T>
void
InverseCubicInterpolator(std::vector <T> &inputVect, std::vector <T> &outputVect, T step)
{
    double x0, x1, x2, x3, y0, y1, y2, y3;
    int i=0;
    outputVect.clear();
    outputVect.push_back(0);
    double value = step + inputVect.front();

    while (value < inputVect.back())
    {
        while(inputVect[i] <= value)
        {
            i++;
        }
        x1 = i-1;
        x2 = i;
        x3 = i+1;

        y1 = inputVect[i-1];
        y2 = inputVect[i];

        if (i == 1)
        {
            x0 = x1;;
            y0 = y1;
        }
        else
        {
            y0 = inputVect[i-2];
            x0 = i-2;
        }

        if (i == inputVect.size() - 1)
        {
            y3 = y2 ;
            x3 = x2;
        }
        else
        {
            y3 = inputVect[i+1];
            x3 = i+1;
        }

        outputVect.push_back(Cubic(-y3, -y2, -y1, -y0, x3, x2, x1, x0, -value));
        value += step;
    }
    outputVect.push_back(inputVect.size() - 1);
}

template <class T>
void
CubicInterpolator(std::vector <T> &transfVect, std::vector <T> &scale, std::vector <T> &outputVect, T LengthLine)
{
    double x0, x1, x2, x3, y0, y1, y2, y3;
    int i = 0;
    int it = 0;
    outputVect.clear();
    outputVect.push_back(transfVect[0]);
    i++;
    while(i < LengthLine -1)
    {
        while (scale[it] < i)
        {
            it++;
        }

        x1 = scale[it-1];
        x2 = scale[it];
        y1 = transfVect[it-1];
        y2 = transfVect[it];

        if(it ==1)
        {
            x0 = x1 ;
            y0 = y1;
        }
        else
        {
            x0 = scale[it-2];
            y0 = transfVect[it-2];
        }

        if (it == scale.size() -1)
        {
            x3 = x2;
            y3 = y2;
        }
        else
        {
            x3 = scale[it+1];
            y3 = transfVect[it+1];
        }

        outputVect.push_back(Cubic<double>(x0, x1, x2, x3, y0, y1, y2, y3, i));
        i++;
    }

    outputVect.push_back(transfVect.back());
}

} // end namespace anima
