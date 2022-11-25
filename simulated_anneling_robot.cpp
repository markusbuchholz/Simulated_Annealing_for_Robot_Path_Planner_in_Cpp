// Markus Buchholz
// g++ simulated_anneling_robot.cpp -o t -I/usr/include/python3.8 -lpython3.8

#include <iostream>
#include <vector>
#include <tuple>
#include <math.h>
#include <algorithm>
#include <random>

#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

//--------Path Planner---------------------------------------------------------------------------

float x1min = 0.0;
float x1max = 50.0;
float x2min = 0.0;
float x2max = 50.0;

float obsX = 25.0;
float obsY = 25.0;
float obsR = 5.0;

float goalX = 45.0;
float goalY = 45.0;

float startX = 2.0;
float startY = 2.0;

float factorKp = 10.0;
float deltaK = 10.0;

//---------Simulated Anneling--------------------------------------------------------------------------

int N = 40;                        // number of evolutions
float NA = 0.0;                    // number of accepted solutions
float p1 = 0.9;                    // start probablitty to accept soltution
float p2 = 0.0001;                 // end probablitty to accept soltution
float temp1 = -1.0 / std::log(p1); // initial temperatoure
float temp2 = -1.0 / std::log(p2); // final temperature

//---------------------------------------------------------------------------------------------

struct Pos
{

    float x;
    float y;
};

//--------------------------------------------------------------------------------

float generateRandom()
{

    std::random_device engine;
    std::uniform_real_distribution<float> distrib(0.0, 1.0);
    return distrib(engine);
}

//--------------------------------------------------------------------------------

float generateRandomX()
{

    std::random_device engine;
    std::uniform_real_distribution<float> distrib(-1.0, 1.0);
    return distrib(engine);
}

//--------------------------------------------------------------------------------

float valueGenerator(float low, float high)
{

    return low + generateRandom() * (high - low);
}

//--------------------------------------------------------------------------------

std::vector<Pos> initPosXY()
{

    std::vector<Pos> pos;
    pos.push_back({startX, startY});

    auto y = [](float x)
    {
        float m = (goalY - startY) / (goalX - startX);

        return m * (x - startX);
    };

    for (float ii = startX + 4.0; ii < goalX; ii = ii + 4.0)
    {

        pos.push_back({ii, y(ii)});
    }

    pos.push_back({goalX, goalY});
    return pos;
}

//--------------------------------------------------------------------------------

float functionX(std::vector<Pos> pos)
{

    // std::vector<float> funcValue;

    float factor1 = 0.0;
    float factor2 = 0.0;

    for (int ii = 0; ii < pos.size() - 1; ii++)
    {

        factor1 = factor1 + std::sqrt(std::pow(pos[ii].x - pos[ii + 1].x, 2) + std::pow(pos[ii].y - pos[ii + 1].y, 2));
        factor2 = factor2 + obsR * obsR + (std::pow(pos[ii].x - pos[ii + 1].x, 2) + std::pow(pos[ii].y - pos[ii + 1].y, 2)) - (pos[ii].x * pos[ii + 1].y + pos[ii + 1].x * pos[ii].y);

        // funcValue.push_back(factor1 + factorKp * factor2);
    }

    return factor1 + factorKp * factor2;
}

//--------------------------------------------------------------------------------
float funcX(Pos a, Pos b)
{

    float factor1 = 0.0;
    float factor2 = 0.0;

    factor1 = factor1 + std::sqrt(std::pow(a.x - b.x, 2) + std::pow(a.y - b.y, 2));
    factor2 = factor2 + obsR * obsR + (std::pow(a.x - b.x, 2) + std::pow(a.y - b.y, 2)) - (a.x * b.y + b.x * a.y);

    return factor1 + factorKp * factor2;
}

//--------------------------------------------------------------------------------

std::vector<Pos> positionUpdate(std::vector<Pos> actualPos)
{
    std::vector<Pos> updatedPos;

    for (int ii = 0; ii < actualPos.size(); ii++)
    {
        Pos Pnew;

        Pnew.x = actualPos[ii].x + generateRandomX() * deltaK;
        Pnew.y = actualPos[ii].y + generateRandomX() * deltaK;

        if (Pnew.x < startX)
        {
            Pnew.x = startX;
        }

        if (Pnew.x > goalX)
        {
            Pnew.x = goalX;
        }

        if (Pnew.y < startY)
        {
            Pnew.y = startY;
        }

        if (Pnew.y > goalY)
        {
            Pnew.y = goalY;
        }
        updatedPos.push_back(Pnew);
    }
    return updatedPos;
}

//--------------------------------------------------------------------------------
std::vector<Pos> runSA()
{

    std::vector<Pos> positionsCurrent = initPosXY();
    float funcValueCurrent = functionX(positionsCurrent);
    float DelataE_avg = 0.0;
    float temp = temp1;
    bool accept = false;
    NA = NA + 1.0;

    for (int ii = 0; ii < N; ii++)
    {

        std::vector<Pos> newPoses= positionUpdate(positionsCurrent);
        //std::vector<Pos> memory = positionsCurrent;
        //positionsCurrent = newPos1;
        float funcValueNew = functionX(newPoses);

        float DelatE = std::abs(funcValueNew - funcValueCurrent);
        if (funcValueNew > funcValueCurrent)
        {
            if (ii == 0)
            {
                DelataE_avg = DelatE;
            }

            float accaptance_prob = std::exp(-DelatE / (DelataE_avg * temp));
            if (generateRandom() < accaptance_prob)
            {
                accept = true;
            }
            else
            {
                accept = false;
            }
        }
        else
        {
            accept = true;
        }

        if (accept == true)
        {
            positionsCurrent = newPoses;    
            NA = NA + 1;
            DelataE_avg = (DelataE_avg * (NA - 1.0) + DelatE) / NA;
        }

        float DELTA_T = std::pow(temp2 / temp1, (1.0 / (N - 1.0))); // step
        temp = temp - DELTA_T;
    }

    return positionsCurrent;
}

//--------------------------------------------------------------------------------------------

std::tuple<std::vector<float>, std::vector<float>> gen_circle(float a, float b, float r)
{

    std::vector<float> xX;
    std::vector<float> yY;

    for (float dt = -M_PI; dt < M_PI; dt += 0.01)
    {

        xX.push_back(a + r * std::cos(dt));
        yY.push_back(b + r * std::sin(dt));
    }
    return std::make_tuple(xX, yY);
}

//----------------------------------------------------------------------------------------------

void plot2D(std::vector<float> xX, std::vector<float> yY)
{
    std::sort(xX.begin(), xX.end());
    std::sort(yY.begin(), yY.end());

    std::tuple<std::vector<float>, std::vector<float>> circle = gen_circle(obsX, obsY, obsR);

    std::vector<float> xObs = std::get<0>(circle);
    std::vector<float> yObs = std::get<1>(circle);

    plt::plot(xX, yY);
    plt::plot(xObs, yObs);
    plt::xlabel("X");
    plt::ylabel("Y");
    plt::show();
}

//---------------------------------------------------------------------------------------------

void plot2DX(std::tuple<std::vector<float>, std::vector<float>> data)
{

    std::vector<float> xX = std::get<0>(data);
    std::vector<float> yY = std::get<1>(data);

    plt::plot(xX, yY);
    plt::show();
}

//---------------------------------------------------------------------------------------------
int main()
{
    std::vector<Pos> path = runSA();

    std::vector<float> xX;
    std::vector<float> yY;

    for (auto &ii : path)
    {
        xX.push_back(ii.x);
        yY.push_back(ii.y);

        std::cout << ii.x << " ," << ii.y << "\n";
    }

    plot2D(xX, yY);
}
