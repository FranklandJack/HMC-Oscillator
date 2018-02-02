#include "DataArray.hpp"

DataArray::DataArray():m_size{0}{}

DataArray::DataArray(int size):m_size{size}
{
    m_data.reserve(size);
}

double& DataArray::operator[](int index)
{
    return m_data[index];
}

const double& DataArray::operator[](int index) const
{
    return m_data[index];
}

void DataArray::push_back(double point) 
{
    m_data.push_back(point);
    m_size++;
}

void DataArray::pop_back()
{
    m_data.pop_back();
    m_size--;
}

void DataArray::reserve(int size)
{
    m_data.reserve(size);
}

double DataArray::sum() const
{
    double sum = 0;
    for(const auto& point : m_data)
    {
        sum += point;
    }

    return sum;
}

double DataArray::sum(const IDataFunctor &fcn) const
{
    double sum = 0;
    for(const auto& point : m_data)
    {
        sum += fcn(point);
    }

    return sum;
}

double DataArray::mean() const
{
    double sum = 0;

    for(const auto& point : m_data)
    {
        sum += point;
    }

    return sum/m_size;
}

double DataArray::mean(const IDataFunctor &fcn) const
{
    double sum = 0;

    for(const auto& point : m_data)
    {
        sum += fcn(point);
    }

    return sum/m_size;
}

double DataArray::error() const
{
    double sum = 0;

    for(const auto& point : m_data)
    {
        sum += point*point;
    }

    double meanSquared = sum / m_size;
    double mean        = (*this).mean();

    return sqrt((meanSquared - mean * mean) / (m_size - 1));

}

double DataArray::error(const IDataFunctor &fcn) const
{
    double sum = 0;

    for(const auto& point : m_data)
    {
        sum += fcn(point)*fcn(point);
    }

    double meanSquared = sum / m_size;
    double mean        = (*this).mean(fcn);

    return sqrt((meanSquared - mean * mean) / (m_size - 1));

}

Estimate DataArray::estimate() const
{
    double error = (*this).error();
    double mean  = (*this).mean();

    return Estimate(mean,error);
}

Estimate DataArray::estimate(const IDataFunctor &fcn) const
{
    double error = (*this).error(fcn);
    double mean  = (*this).mean(fcn);

    return Estimate(mean,error);
}


double DataArray::jackKnifeError() const
{
    double sum = 0;
    double mean = (*this).mean();
    double total = (*this).sum();

    for(const auto& point : m_data)
    {
        double reducedMean = (total-point)/(m_size-1);
        sum += (reducedMean - mean)*(reducedMean - mean);
    }

    return sqrt(sum);
}

double DataArray::jackKnifeError(const IDataFunctor &fcn) const
{
    double sum = 0;
    double mean = (*this).mean(fcn);
    double total = (*this).sum(fcn);

    for(const auto& point : m_data)
    {
        double reducedMean = (total-fcn(point))/(m_size-1);
        sum += (reducedMean - mean)*(reducedMean - mean);
    }

    return sqrt(sum);
}


Estimate DataArray::jackKnifeEstimate() const
{
    double mean = (*this).mean();
    double error = (*this).jackKnifeError();

    return Estimate(mean,error);
}

Estimate DataArray::jackKnifeEstimate(const IDataFunctor& fcn) const
{
    double mean = (*this).mean(fcn);
    double error = (*this).jackKnifeError(fcn);

    return Estimate(mean,error);
}

double DataArray::bootstrapError(std::default_random_engine& generator, int totalIterations) const
{
    static std::uniform_int_distribution<int> distribution;

    double bsMean       = 0;
    double bsMeanSquare = 0;
    for(int iteration = 0; iteration < totalIterations; ++iteration)
    {
        double tempMean        = 0;
        double tempMeanSquared = 0;
        for(int i = 0; i < m_size; ++i)
        {
            double sample = m_data[distribution(generator)];
            tempMean += sample;
            tempMeanSquared += sample*sample;

        }

        tempMean /= m_size; 
        bsMean += tempMean;

        tempMeanSquared /= m_size;
        bsMeanSquare += tempMeanSquared;
    }

    return sqrt(bsMeanSquare-bsMean*bsMean);
}

double DataArray::bootstrapError(std::default_random_engine &generator, int totalIterations, const IDataFunctor &fcn) const
{
    static std::uniform_int_distribution<int> distribution;

    double bsMean       = 0;
    double bsMeanSquare = 0;
    for(int iteration = 0; iteration < totalIterations; ++iteration)
    {
        double tempMean        = 0;
        double tempMeanSquared = 0;
        for(int i = 0; i < m_size; ++i)
        {
            double sample = m_data[distribution(generator)];
            tempMean += fcn(sample);
            tempMeanSquared += fcn(sample)*fcn(sample);

        }

        tempMean /= m_size; 
        bsMean += tempMean;

        tempMeanSquared /= m_size;
        bsMeanSquare += tempMeanSquared;
    }

    return sqrt(bsMeanSquare-bsMean*bsMean);
}


Estimate DataArray::bootstrapEstimate(std::default_random_engine &generator, int totalIterations) const
{
    double mean = (*this).mean();
    double error = (*this).bootstrapError(generator,totalIterations);
    return Estimate(mean,error);
}

Estimate DataArray::bootstrapEstimate(std::default_random_engine &generator, int totalIterations, const IDataFunctor &fcn) const
{
    double mean = (*this).mean(fcn);
    double error = (*this).bootstrapError(generator,totalIterations,fcn);
    return Estimate(mean,error);

}

std::ostream& operator<<(std::ostream& out, const DataArray& data)
{
    int index = 0;
    for(const auto& point : data.m_data)
    {
        out << index++ << ' ' << point << '\n';
    }

    return out;
}

double DataArray::autoCorrelation(int t) const
{
    double term1          = 0;
    double mean_m         = 0;
    double mean_mSquared  = 0;

    for(int point = 0; point < m_size; ++point)
    {
        term1 += m_data[point]*m_data[(point+t)%m_size];
        mean_m += m_data[point];
        mean_mSquared += m_data[point] * m_data[point];
    }

    term1         /= m_size;
    mean_m        /= m_size;
    mean_mSquared /= m_size;

    double normalisation = mean_mSquared - mean_m * mean_m;
    double term2 = mean_m * mean_m;

    return (term1-term2)/normalisation;
}

std::vector<double> DataArray::autoCorrelation(int t1, int t2) const
{
    std::vector<double> autoCorrelationData;
    autoCorrelationData.reserve(t2-t1);
    for(int t = t1; t < t2; ++t)
    {
        autoCorrelationData.push_back(autoCorrelation(t));
    }

    return autoCorrelationData;
}