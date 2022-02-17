#include <vector>
#include <cstdint>
#include <algorithm>
#include "time_point.h"
#include "vldb.h"

void print_time_series(const TimeSeries& timeSeries)
{
#ifdef _DEBUG
    std::cout << "TimeSeries: size=" << timeSeries.size();
    for (int i = 0; i < std::min(10, (int)timeSeries.size()); i++)
    {
        timeSeries.at(i).toString();
    }
#endif
}

void print_matrix(const_matrix mat)
{
#ifdef _DEBUG
    std::cout << "Matrix:" << std::endl;
    std::cout << mat << std::endl;
#endif
}

void print_time_point(const TimePoint& timepoint)
{
#ifdef _DEBUG
    timepoint.toString();
#endif
}

TimeSeries IMR_Alorithm(const TimeSeries& dirtySeries, const TimeSeries& labelSeries,
    const std::vector<bool>& labelList, int p, double delta, int maxIterations)
{
    Vldb vldb(labelList, p, delta, maxIterations);

    print_time_series(dirtySeries);
    print_time_series(labelSeries);

    int size = dirtySeries.size();
    int rows = size - p;

    // form z
    std::vector<double> ary_diff;
    for (int i = 0; i < size; i++)
    {
        double diff = labelSeries.at(i).value - dirtySeries.at(i).value;
        ary_diff.push_back(diff);
    }

    // build x,y for params estimation
    matrix xMatrix(rows, p);
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < p; j++)
        {
            xMatrix(i, j) = ary_diff.at(p + i - j - 1);
        }
    }
    print_matrix(xMatrix);

    matrix yMatrix(rows, 1);
    for (int i = 0; i < rows; i++)
    {
        yMatrix(i, 0) = (double)ary_diff.at(p + i);
    }
    print_matrix(yMatrix);

    // begin iteration
    int iterations = 1;
    do 
    {
        const_matrix olsmatrix = vldb.initOLSParams(xMatrix, yMatrix);
        const_matrix combinedMatrix = vldb.combine(olsmatrix, xMatrix);

        int index = vldb.repair(combinedMatrix, yMatrix);
        if (index < 0)
            break;

        double combval= combinedMatrix(index, 0);
        yMatrix(index, 0) = combval;

        // update x
        for (int j = 0; j < p; j++)
        {
            int i = index + 1 + j; // p+i-j-1 \Leftrightarrow p+i = index+p
            if (i < 0)
                continue;

            if (i >= rows)
                break;

            xMatrix(i, j) = combval;
        }
    } while (++iterations <= vldb.maxIterations);

    std::cout << "Stop after " << iterations << " iterations " << std::endl;

    // form result series
    TimeSeries resultSeries;
    for (int i = 0; i < size; i++)
    {
        double modify(0.0);
        uint64_t timestamp = labelSeries.at(i).timestamp;
        if (!labelList.at(i))
        {
            modify = dirtySeries.at(i).value + yMatrix(i - p, 0);
        }
        else
        {
            modify = labelSeries.at(i).value;
        }
        TimePoint timepoint(timestamp, modify);
        print_time_point(timepoint);
        resultSeries.push_back(timepoint);
    }

    print_time_series(resultSeries);
    return resultSeries;
}

void initMatrix(Vldb& vldb, matrix& A, matrix& B, std::vector<double>& aryDiff, int size)
{
    // aii
    for (int i = 1; i <= vldb.p; ++i)
    {
        double val = 0;
        for (int l = vldb.p - i; l <= size - i; ++l)
        {
            val += aryDiff.at(l) * aryDiff.at(l);
        }
        A(i - 1, i - 1) = val;
    }
    print_matrix(A);

    // aij=aji
    for (int i = 1; i <= vldb.p; ++i)
    {
        for (int j = i + 1; j <= vldb.p; ++j)
        {
            int u = j - i;
            double val = 0;
            for (int l = vldb.p - i; l < size - i; ++l)
            {
                val += aryDiff.at(l) * aryDiff.at(l - u);
            }
            A(j - 1, i - 1) = val;
            A(i - 1, j - 1) = val;
        }
    }
    print_matrix(A);

    // bi
    for (int i = 1; i <= vldb.p; ++i)
    {
        int u = i;
        double val = 0;
        for (int l = vldb.p; l < size; ++l)
        {
            val += aryDiff.at(l) * aryDiff.at(l - u);
        }
        B(i - 1, 0) = val;
    }
    print_matrix(B);
}

void update(Vldb& vldb, int index, double preVal, double val, matrix& A, matrix& B, std::vector<double>& aryDiff, int size)
{
    int n = size - 1;
    int zPos = index + vldb.p;
    aryDiff[zPos] = val;

    // A:p*p, B:p*1
    double aiiVal = val * val - preVal * preVal;
    double aijVal = val - preVal;
    double addVal = 0;

    // the border will minus 1 since the index n already minus 1
    for (int i = 1; i <= vldb.p; ++i)
    {
        if (zPos < (vldb.p + 1 - i) - 1 || zPos > n - i)
            continue;
        A(i - 1, i - 1) += aiiVal;
    }
    print_matrix(A);

    // aij
    for (int i = 1; i <= vldb.p; ++i)
    {
        for (int j = i + 1; j <= vldb.p; ++j)
        {
            if (zPos < (vldb.p + 1 - j) - 1 || zPos > n - i)
            {
                addVal = 0;
            }
            else if (zPos > n - j && zPos <= n - i)
            {
                addVal = aijVal * aryDiff[zPos - j + i];
            }
            else if (zPos >= (vldb.p + 1 - j) - 1 && zPos < (vldb.p + 1 - i) - 1)
            {
                addVal = aijVal * aryDiff[zPos + j - i];
            }
            else
            {
                addVal = aijVal * (aryDiff[zPos - j + i] + aryDiff[zPos + j - i]);
            }
            A(j - 1, i - 1) += addVal;
            A(i - 1, j - 1) += addVal;
        }
    }
    print_matrix(A);

    // bi
    for (int i = 1; i <= vldb.p; ++i)
    {
        if (zPos < (vldb.p + 1 - i) - 1)
        {
            addVal = 0;
        }
        else if (zPos >= (vldb.p + 1 - i) - 1 && zPos < (vldb.p + 1 + i) - 1)
        {
            addVal = aijVal * aryDiff[zPos + i];
        }
        else if (zPos > n - i)
        {
            addVal = aijVal * aryDiff[zPos - i];
        }
        else
        {
            addVal = aijVal * (aryDiff[zPos - i] + aryDiff[zPos + i]);
        }
        B(i - 1, 0) += addVal;
    }
    print_matrix(B);
}

void Compute(Vldb& vldb, matrix& xMatrix, matrix& yMatrix, std::vector<double>& aryDiff, int size)
{
    int rows = yMatrix.size1();

    // begin iteration
    matrix aMatrix(vldb.p, vldb.p), bMatrix(vldb.p, 1);
    initMatrix(vldb, aMatrix, bMatrix, aryDiff, size);

    int iteration = 0;
    do
    {
        const_matrix icmatrix = vldb.initICParams(aMatrix, bMatrix);
        const_matrix yhatMatrix = vldb.combine(icmatrix, xMatrix);

        int index = vldb.repair(yhatMatrix, yMatrix);
        if (index < 0)
            break;

        double preVal = yMatrix(index, 0);
        double val = yhatMatrix(index, 0);
        yMatrix(index, 0) = val;

        // update x
        for (int j = 0; j < vldb.p; ++j)
        {
            int i = index + j + 1; // p+i-j-1 \Leftrightarrow p+i = index+p
            if (i < 0)
                continue;

            if (i >= rows)
                break;

            xMatrix(i, j) = val;
        }

        int zPos = index + vldb.p;
        aryDiff[zPos] = val;
        update(vldb, index, preVal, val, aMatrix, bMatrix, aryDiff, size);

    } while (++iteration <= vldb.maxIterations);

    std::cout << "Stop after " << iteration << " iterations" << std::endl;
}

void incrase_1(Vldb& vldb, matrix xMatrix, matrix yMatrix, std::vector<double>& aryDiff, int size) {
    int rowNum = size - vldb.p;

    // initial alpha and beta
    double alpha = 0, beta = 0;

    for (int i = 1; i < size - 1; i++)
    {
        double itmval = aryDiff.at(i);
        beta += itmval * aryDiff.at(i - 1);
        alpha += itmval * itmval;
    }
    beta += aryDiff.at(size - 1) * aryDiff.at(size - 2);
    alpha += aryDiff.at(0) * aryDiff.at(0);

    int iterations = 1;
    matrix phiArray(1, 1);
    do
    {
        phiArray(0, 0) = beta / alpha;
        
        matrix phi(phiArray);
        const_matrix yhatMatrix = vldb.combine(phi, xMatrix);

        int index = vldb.repair(yhatMatrix, yMatrix);
        if (index == -1)
            break;

        double preVal = yMatrix(index, 0);
        double val = yhatMatrix(index, 0);
        yMatrix(index, 0) = val;

        // update x
        for (int j = 0; j < vldb.p; ++j)
        {
            int i = index + 1 + j;
            if (i >= rowNum)
                break;
            if (i < 0)
                continue;

            xMatrix(i, j) = val;
        }

        // update alpha
        int zPos = index + vldb.p;
        aryDiff[zPos] = val;
        if (zPos <= size - 2) {
            alpha = alpha - preVal * preVal + val * val;
        }

        // update beta
        double zPosPreVal = aryDiff.at(zPos - 1);
        double zPosNextVal = zPos < size - 1 ? aryDiff.at(zPos + 1) : 0;
        beta = beta + (val - preVal) * (zPosPreVal + zPosNextVal);

    } while (++iterations <= vldb.maxIterations);

    std::cout << "Stop after " << iterations << " iterations" << std::endl;
}

void incrase_2(Vldb& vldb, matrix& xMatrix, matrix& yMatrix, std::vector<double>& aryDiff, int size) {
    int rowNum = size - vldb.p;

    // initial alpha, beta, gamma
    double alpha = 0, beta = 0, gamma = 0;

    // [3,n-2] -> [2,size-3]
    for (int i = 2; i < size - 2; ++i)
    {
        alpha += aryDiff.at(i) * aryDiff.at(i);
        beta += aryDiff.at(i) * aryDiff[i - 1];
        gamma += aryDiff.at(i) * aryDiff[i - 2];
    }

    double alpha1 = alpha + aryDiff.at(0) * aryDiff.at(0);
    double alpha2 = alpha + aryDiff.at(size-2) * aryDiff.at(size-2);
    double beta2 = beta + aryDiff.at(1) * aryDiff.at(0);
    double beta3 = beta + aryDiff.at(size-1) * aryDiff.at(size-2);
    double gamma3 = gamma + aryDiff.at(size-2) * aryDiff.at(size-4) + aryDiff.at(size-1) * aryDiff.at(size-3);
    double det = alpha2 * alpha1 - beta2 * beta2;
    
    double alphachange = 0;
    matrix phiArray(2, 1);
    int iterations = 1;
    do
    {
        phiArray(0, 0) = (beta3 * alpha2 - gamma3 * beta2) / det;
        phiArray(1, 0) = (-beta3 * beta2 + gamma3 * alpha1) / det;

        matrix phi(phiArray);
        // Matrix phi = OLSInitParams(xMatrix, yMatrix);
        const_matrix yhatMatrix = vldb.combine(phi, xMatrix);

        int index = vldb.repair(yhatMatrix, yMatrix);
        if (index == -1)
            break;

        double preVal = yMatrix(index, 0);
        double val = yhatMatrix(index, 0);
        yMatrix(index, 0) = val;

        for (int j = 0; j < vldb.p; ++j)
        {
            int i = index + 1 + j;
            if (i >= rowNum)
                break;
            if (i < 0)
                continue;

            xMatrix(i, j) = val;
        }

        // update alpha
        int zPos = index + vldb.p;
        aryDiff[zPos] = val;
        alphachange = -preVal * preVal + val * val;
        if (zPos <= size - 2 && zPos >= 2) {
            alpha2 += alphachange;
        }
        if (zPos <= size - 3 && zPos >= 1) {
            alpha1 += alphachange;
        }

        // update beta
        double zPosPreValb2 = aryDiff[zPos - 1];
        double zPosPreValb3 = zPos >= 2 ? aryDiff[zPos - 1] : 0;
        double zPosNextValb2 = zPos <= size - 2 ? aryDiff[zPos + 1] : 0;
        double zPosNextValb3 = zPos <= size - 2 ? aryDiff[zPos + 1] : 0; // not to
                                                                    // exceed size

        beta2 = beta2 + (val - preVal) * (zPosPreValb2 + zPosNextValb2);
        beta3 = beta3 + (val - preVal) * (zPosPreValb3 + zPosNextValb3);
        // update gamma
        double zPosPreValg3 = zPos >= 2 ? aryDiff[zPos - 2] : 0;
        double zPosNextValg3 = zPos <= size - 3 ? aryDiff[zPos + 2] : 0;

        gamma3 = gamma3 + (val - preVal) * (zPosPreValg3 + zPosNextValg3);

        // update det
        det = alpha2 * alpha1 - beta2 * beta2;

    } while (++iterations <= vldb.maxIterations);

    std::cout << "Stop after " << iterations << " iterations" << std::endl;
}

void incrase_3(Vldb& vldb, matrix& xMatrix, matrix& yMatrix, std::vector<double>& aryDiff, int size)
{
    int p = vldb.p;
    int rowNum = size - p;

    // initial alpha, beta, gamma
    double alpha = 0, beta = 0, gamma = 0, zeta = 0;

    // [4,n-3] -> [3,size-4]
    for (int i = 3; i <= size - 4; ++i)
    {
        alpha += aryDiff.at(i) * aryDiff.at(i);
        beta += aryDiff.at(i) * aryDiff[i - 1];
        gamma += aryDiff.at(i) * aryDiff[i - 2];
        zeta += aryDiff.at(i) * aryDiff[i - 3];
    }

    double alpha1 = alpha + aryDiff.at(2) * aryDiff.at(2) + aryDiff.at(1) * aryDiff.at(1) + aryDiff.at(0) * aryDiff.at(0);
    double alpha2 = alpha + aryDiff.at(2) * aryDiff.at(2) + aryDiff.at(1) * aryDiff.at(1) + aryDiff.at(size-3) * aryDiff.at(size-3);
    double alpha3 = alpha + aryDiff.at(2) * aryDiff.at(2) + aryDiff.at(size-3) * aryDiff.at(size-3) + aryDiff.at(size-2) * aryDiff.at(size-2);

    double beta2 = beta + aryDiff.at(2) * aryDiff.at(1) + aryDiff.at(1) * aryDiff.at(0) + aryDiff.at(size-3) * aryDiff.at(size-4);
    double beta3 = beta + aryDiff.at(2) * aryDiff.at(1) + aryDiff.at(size-3) * aryDiff.at(size-4) + aryDiff.at(size-2) * aryDiff.at(size-3);
    double beta4 = beta + aryDiff.at(size-3) * aryDiff.at(size-4) + aryDiff.at(size-2) * aryDiff.at(size-3) + aryDiff.at(size-1) * aryDiff.at(size-2);

    double gamma3 = gamma + aryDiff.at(2) * aryDiff.at(0) + aryDiff.at(size-3) * aryDiff.at(size-5) + aryDiff.at(size-2) * aryDiff.at(size-4);
    double gamma4 = gamma + aryDiff.at(size-3) * aryDiff.at(size-5) + aryDiff.at(size-2) * aryDiff.at(size-4) + aryDiff.at(size-1) * aryDiff.at(size-3);

    double zeta4 = zeta + aryDiff.at(size-3) * aryDiff.at(size-6) + aryDiff.at(size-2) * aryDiff.at(size-5) + aryDiff.at(size-1) * aryDiff.at(size-4);

    double A = alpha2 * alpha1 - beta2 * beta2;
    double B = -(beta3 * alpha1 - beta2 * gamma3);
    double C = beta3 * beta2 - alpha2 * gamma3;
    double D = -(beta3 * alpha1 - gamma3 * beta2);
    double E = alpha3 * alpha1 - gamma3 * gamma3;
    double F = -(alpha3 * beta2 - beta3 * gamma3);
    double G = beta3 * beta2 - gamma3 * alpha2;
    double H = -(alpha3 * beta2 - gamma3 * beta3);
    double I = alpha3 * alpha2 - beta3 * beta3;
    double det = alpha3 * A + beta3 * B + gamma3 * C;

    double alphachange = 0;
    matrix phiArray(3, 1);

    // begin iteration
    int index = -1;
    int iterations = 1;
    do
    {
        phiArray(0, 0) = (beta4 * A + gamma4 * D + zeta4 * G) / det;
        phiArray(1, 0) = (beta4 * B + gamma4 * E + zeta4 * H) / det;
        phiArray(2, 0) = (beta4 * C + gamma4 * F + zeta4 * I) / det;
        
        matrix phi(phiArray);
        matrix yhatMatrix = vldb.combine(phi, xMatrix);

        index = vldb.repair(yhatMatrix, yMatrix);
        if (index == -1)
            break;

        double preVal = yMatrix(index, 0);
        double val = yhatMatrix(index, 0);
        yMatrix(index, 0) = val;

        // update x
        for (int j = 0; j < p; ++j)
        {
            int i = index + 1 + j;
            if (i < 0) continue;
            if (i >= rowNum) break;
            xMatrix(i, j) = val;
        }

        // update alpha
        int zPos = index + p;
        aryDiff[zPos] = val;
        alphachange = -preVal * preVal + val * val;
        if (zPos <= size - 2 && zPos >= 3)
        {
            alpha3 += alphachange;
        }
        if (zPos <= size - 3 && zPos >= 2)
        {
            alpha2 += alphachange;
        }
        if (zPos <= size - 4 && zPos >= 1)
        {
            alpha1 += alphachange;
        }

        // update beta
        double zPosPreValb2 = aryDiff.at(zPos - 1);
        double zPosPreValb3 = zPos >= 2 ? aryDiff.at(zPos - 1) : 0;
        double zPosPreValb4 = zPos >= 3 ? aryDiff.at(zPos - 1) : 0;
        double zPosNextValb2 = zPos <= size - 3 ? aryDiff.at(zPos + 1) : 0;
        double zPosNextValb3 = zPos <= size - 2 ? aryDiff.at(zPos + 1) : 0;
        double zPosNextValb4 = zPos <= size - 2 ? aryDiff.at(zPos + 1) : 0;

        beta2 = beta2 + (val - preVal) * (zPosPreValb2 + zPosNextValb2);
        beta3 = beta3 + (val - preVal) * (zPosPreValb3 + zPosNextValb3);
        beta4 = beta4 + (val - preVal) * (zPosPreValb4 + zPosNextValb4);

        double zPosPreValg3 = zPos >= 2 ? aryDiff.at(zPos - 2) : 0;
        double zPosPreValg4 = zPos >= 3 ? aryDiff.at(zPos - 2) : 0;
        double zPosNextValg3 = zPos <= size - 3 ? aryDiff.at(zPos + 2) : 0;
        double zPosNextValg4 = zPos <= size - 3 ? aryDiff.at(zPos + 2) : 0;

        gamma3 = gamma3 + (val - preVal) * (zPosPreValg3 + zPosNextValg3);
        gamma4 = gamma4 + (val - preVal) * (zPosPreValg4 + zPosNextValg4);

        double zPosPreValz4 = zPos >= 3 ? aryDiff.at(zPos - 3) : 0;
        double zPosNextValz4 = zPos <= size - 4 ? aryDiff.at(zPos + 3) : 0;

        zeta4 = zeta4 + (val - preVal) * (zPosPreValz4 + zPosNextValz4);
        A = alpha2 * alpha1 - beta2 * beta2;
        B = -(beta3 * alpha1 - beta2 * gamma3);
        C = beta3 * beta2 - alpha2 * gamma3;
        D = -(beta3 * alpha1 - gamma3 * beta2);
        E = alpha3 * alpha1 - gamma3 * gamma3;
        F = -(alpha3 * beta2 - beta3 * gamma3);
        G = beta3 * beta2 - gamma3 * alpha2;
        H = -(alpha3 * beta2 - gamma3 * beta3);
        I = alpha3 * alpha2 - beta3 * beta3;
        det = alpha3 * A + beta3 * B + gamma3 * C;

    } while (++iterations <= vldb.maxIterations);

    std::cout << "Stop after " << iterations << " iterations" << std::endl;
}

TimeSeries IMRIC_Algorithm(const TimeSeries& dirtySeries, const TimeSeries& labelSeries,
    const std::vector<bool>& labelList, int p, double delta, int maxIterations)
{
    Vldb vldb(labelList, p, delta, maxIterations);

    print_time_series(dirtySeries);
    print_time_series(labelSeries);

    int size = dirtySeries.size();
    std::vector<double> aryDiff; // form z
    for (int i = 0; i < size; ++i)
    {
        double diff = labelSeries.at(i).value - dirtySeries.at(i).value;
        aryDiff.push_back(diff);
    }

    int rows = size - p;
    matrix xMatrix(rows, p); // build x,y for params estimation
    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < p; ++j)
        {
            xMatrix(i, j) = aryDiff.at(p + i - j - 1);
        }
    }
    print_matrix(xMatrix);

    matrix yMatrix(rows, 1);
    for (int i = 0; i < rows; ++i)
    {
        yMatrix(i, 0) = (double)aryDiff.at(p + i);
    }
    print_matrix(yMatrix);

    switch (p) {
    case 1:
        incrase_1(vldb, xMatrix, yMatrix, aryDiff, size);
        break;
    case 2:
        incrase_2(vldb, xMatrix, yMatrix, aryDiff, size);
        break;
    case 3:
        incrase_3(vldb, xMatrix, yMatrix, aryDiff, size);
        break;
    default:
        Compute(vldb, xMatrix, yMatrix, aryDiff, size);
    }

    // form result series
    TimeSeries resultSeries;
    double modify = 0;

    for (int i = 0; i < size; ++i)
    {
        int64_t timestamp = labelSeries.at(i).timestamp;
        if (!labelList.at(i))
        {
            modify = dirtySeries.at(i).value + yMatrix(i - p, 0);
        }
        else
        {
            modify = labelSeries.at(i).value;
        }

        TimePoint tp(timestamp, modify);
        print_time_point(tp);
        resultSeries.push_back(tp);
    }

    print_time_series(resultSeries);
    return resultSeries;
}