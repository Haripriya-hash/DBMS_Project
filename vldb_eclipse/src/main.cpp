// cvldb.cpp : Defines the entry point for the application.
//

#include <string>
#include <vector>
#include <fstream>
#include <math.h>
#include <iostream>
#include "time_point.h"
#include "main.h"

using namespace std;

static std::string filename = "data.txt";

void Tokenize(const string& str, const string& delimiter, vector<string>& tokens)
{
    int start = 0;
    int end = str.find(delimiter);
    while (end != -1)
    {
        tokens.push_back(str.substr(start, end - start));
        start = end + delimiter.size();
        end = str.find(delimiter, start);
    }
    string e = str.substr(start, end - start);
    if (!e.empty())
    {
        tokens.push_back(e);
    }
}

TimeSeries ReadData(const string& filename, int index, const string& delimiter)
{
    TimeSeries timeSeries;

    ifstream file(filename);
    string line;
    while (getline(file, line))
    {
        vector<string> tokens;
        Tokenize(line, delimiter, tokens);

        long timestamp = stol(tokens.at(0));
        double value = stod(tokens.at(index));
        TimePoint tp(timestamp, value);
        timeSeries.push_back(tp);
    }

    return timeSeries;
}

void ReadLabel(const string& filename, int index, const string& delimiter, vector<bool>& labelList)
{
    ifstream file(filename);
    string line;
    while (getline(file, line))
    {
        vector<string> tokens;
        Tokenize(line, delimiter, tokens);

        bool value = false;
        if (strncmp(tokens.at(index).c_str(), "true", 4) == 0) {
        	value = true;
        }
        labelList.push_back(value);
    }
}

double CalcRMS(TimeSeries truthSeries, TimeSeries resultSeries, vector<bool> labelList)
{
    double cost = 0;
    int len = truthSeries.size();
    int labelNum = 0;

    for (int i = 0; i < len; i++)
    {
        if (labelList.at(i))
        {
            labelNum++;
            continue;
        }

        double delta = resultSeries.at(i).value - truthSeries.at(i).value;
        cost += delta * delta;
    }

    cost /= (len - labelNum);
    return sqrt(cost);
}

void TestAlgorithm1()
{
    cout << "Test Algorithm 1" << endl;

    string delimiter(",");
    const TimeSeries& dirtySeries = ReadData(filename, 1, delimiter);
    const TimeSeries& labelSeries = ReadData(filename, 2, delimiter);
    const TimeSeries& truthSeries = ReadData(filename, 3, delimiter);
    vector<bool> labelList;
    ReadLabel(filename, 4, delimiter, labelList);

    double rmsDirty = CalcRMS(truthSeries, dirtySeries, labelList);
    cout << "Dirty RMS error is " << rmsDirty << endl;

    int p = 3;
    double delta = 0.1;
    int maxIterations = 100000;

    const TimeSeries& resultSeries = IMR_Alorithm(dirtySeries, labelSeries, labelList, p, delta, maxIterations);

    double rms = CalcRMS(truthSeries, resultSeries, labelList);
    cout << "RMS error is " << rms << endl;
}

void TestAlgorithm2()
{
    cout << endl << "Test Algorithm 2" << endl;

    string delimiter(",");
    const TimeSeries& dirtySeries = ReadData(filename, 1, delimiter);
    const TimeSeries& labelSeries = ReadData(filename, 2, delimiter);
    const TimeSeries& truthSeries = ReadData(filename, 3, delimiter);
    vector<bool> labelList;
    ReadLabel(filename, 4, delimiter, labelList);

    double rmsDirty = CalcRMS(truthSeries, dirtySeries, labelList);
    cout << "Dirty RMS error is " << rmsDirty << endl;

    int p = 3;
    double delta = 0.1;
    int maxIterations = 100000;

    const TimeSeries& resultSeries = IMRIC_Algorithm(dirtySeries, labelSeries, labelList, p, delta, maxIterations);

    double rms = CalcRMS(truthSeries, resultSeries, labelList);
    cout << "RMS error is " << rms << endl;
}

int main()
{
    TestAlgorithm1();
    TestAlgorithm2();
    return 0;
}
