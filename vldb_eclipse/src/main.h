#pragma once
#include "time_point.h"

TimeSeries IMR_Alorithm(const TimeSeries& dirtySeries, const TimeSeries& labelSeries,
    const std::vector<bool>& labelList, int p, double delta, int maxIterations);

TimeSeries IMRIC_Algorithm(const TimeSeries& dirtySeries, const TimeSeries& labelSeries,
    const std::vector<bool>& labelList, int p, double delta, int maxIterations);