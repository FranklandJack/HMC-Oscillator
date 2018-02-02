#ifndef DataArray_hpp
#define DataArray_hpp

#include <vector>
#include <cmath>
#include <random>
#include <iostream>
#include <random>
#include "Estimate.hpp"
#include "IDataFunctor.hpp"


class DataArray
{
private:
    std::vector<double> m_data;
    int m_size;

public:
    DataArray();
    DataArray(int);

    double& operator[](int);
    const double& operator[](int) const;

    void push_back(double);
    void pop_back();

    void reserve(int);

    double sum() const;
    double sum(const IDataFunctor&) const;

    double mean() const;
    double mean(const IDataFunctor&) const;

    double error() const;
    double error(const IDataFunctor&) const;

    Estimate estimate() const;
    Estimate estimate(const IDataFunctor&) const;

    
    double jackKnifeError() const;
    double jackKnifeError(const IDataFunctor&) const;

    Estimate jackKnifeEstimate() const;
    Estimate jackKnifeEstimate(const IDataFunctor&) const;

    double bootstrapError(std::default_random_engine&, int) const;
    double bootstrapError(std::default_random_engine&, int, const IDataFunctor&) const;

    Estimate bootstrapEstimate(std::default_random_engine&, int) const;
    Estimate bootstrapEstimate(std::default_random_engine&, int, const IDataFunctor&) const;

    friend std::ostream& operator<<(std::ostream&, const DataArray&);

    double autoCorrelation(int) const;

    std::vector<double> autoCorrelation(int,int) const;

    void printAutoCorrelation(int,int) const;

};

#endif /* DataArray_hpp */