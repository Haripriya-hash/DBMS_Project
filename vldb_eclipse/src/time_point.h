#pragma once

#include <iostream>
#include <vector>
#include <cstdint>

class TimePoint
{
public:
    TimePoint()
    {
        this->timestamp = 0;
        this->value = 0;
    }

    TimePoint(int64_t timestamp, double value)
    {
        this->timestamp = timestamp;
        this->value = value;
    }

    void toString() const
    {
        std::cout << this->timestamp << ", " << this->value << std::endl;
    }

public:
    int64_t timestamp;
    double value;
};

typedef std::vector<TimePoint> TimeSeries;
