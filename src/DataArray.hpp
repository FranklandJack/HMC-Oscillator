#ifndef DataArray_hpp
#define DataArray_hpp

#include <vector>
#include <cmath>
#include <random>
#include <iostream>
#include <random>


/**
 *\file 
 *\class DataArray
 *\brief Class for statistical analysis of sample sets.
 *
 * This class essentially holds a vector but has methods that make doing statistical easy.
 */
class DataArray
{
private:
	/**
	 *\brief Member variable to hold the raw data (sample set).
	 */
    std::vector<double> m_data;
    /**
     *\brief Member variable to hold the size of the sample set.
     */
    int m_size;

public:
	/**
	 *\class IDataFunctor
	 *\brief Interface for a functor that takes a data array as its argument.
	 *
	 * This is extremely useful for calculating functions whose dependence on the data set is non-trivial.
	 * For example, a function that depends on the first and second moment can inherit from this class as 
	 * a functor where it can be customized and still act on the data in the set. This class also makes 
	 * calculating bootstrap and jackknife errors very portable with a general function for doing so.
	 */
	class IDataFunctor
	{
		public:
			/**
			 *\brief Operator() overload so the class behaves as a function (functor)
			 *
			 * Is pure virtual.
			 *
			 *\param data a const DataArray reference for the function to act on.
			 *\return a floating point value representing the result of the function.
			 */
			virtual double operator()(const DataArray &data) const = 0;
	};

	/**
	 *\brief Default constructor.
	 */
    DataArray();

    /**
     *\brief Constructor that takes reserves space.
     *\param size integer instance representing size of sample set.
     *
     * This constructor will reserve space for the vector which is more efficient. If the user
     * knows how many data points will be in the set before use this constructor should be used.
     */
    DataArray(int size);

	/** 
     *\brief Getter method for size of data.
     *\return Integer value representing the size of the DataArray.
     */
    int getSize() const;

    /**
     *\brief operator() overload to accesses samples.
     *\param index integer value representing  the index of the sample being accessed.
     *\return floating point reference to the sample so it can be used or changed.
     */
    double& operator[](int index);

    /**
     *\brief Constant operator() overload to access samples in a constant instance.
     *
     * See non constant version
     */
    const double& operator[](int index) const;

    /**
     *\brief Adds a sample to the DataArray.
     *\param sample floating point instance to be added to DataArray.
     *
     * Implemented using the vector push_back method so will be added as last element of vector
     * and will work quickest when vector has been reserved.
     */
    void push_back(double sample);

    /**
     *\brief Removes most recently added sample from DataArray.
     */
    void pop_back(); 

    /**
     *\brief Reserves memory for DataArray making it faster.
     *\param size integer value representing the number of elements to reserve space for.
     *
     * This method reserves space for the vector member variable and should be used when the number of samples
     * being recorded is known beforehand.
     */
    void reserve(int size);

    /**
     *\brief Helper function to sum the samples in the data.
	 *\param Floating point instance representing the sum of the data.
	 */
    double sum() const;

    /**
     *\brief Helper function to sum the function of the samples in the data.
	 *\param Floating point instance representing the sum of the data.
	 */
    double sum(const IDataFunctor&) const;

    /**
     *\brief Method to calculate the mean of the data.
     *\return Floating point value representing the mean of the data.
     */
    double mean() const;

    /**
     *\brief Method to calculate the square mean of the data.
     *\return Floating point value representing the square mean of the data.
     */
    double squareMean() const;

    /**
     *\brief Method to calculate the naive error of the data.
     *\return Floating point value representing the naive error of the data.
     */
    double error() const;

    /**
     *\brief operator<< overload to output the data array to a stream.
     *\param out std::ostream reference that is the stream being output to.
	 *\param dataArray constant DataArray reference to be printed.
	 *\return std::ostream reference so that the operator can be chained.
	 *
	 * The operator will print the data points in a column with their respective indices.
	 */
    friend std::ostream& operator<<(std::ostream &out , const DataArray& dataArray);

    /**
     *\brief function to calculate the autocorrelation function for a specific computer time.
     *\param t time value to compute autocorrelation of data for.
     *\return floating point value representing value of autocorrelation function.
     */
    double autoCorrelation(int t) const;

    /**
     *\brief function to calculate the autocorrelation function for a range computer times.
     *\param t1 initial time value to compute autocorrelation of data for.
     *\param t2 final time value to compute autocorrelation for.
     *\return vector of floating point values representing values of autocorrelation function indexed 
     * according to their position in the vector.
     */
    std::vector<double> autoCorrelation(int t1, int t2) const;

    /**
     *\brief Function to calculate the integrated autocorrelation time.
     *\return floating point value representing the integrated autocorrelation time.
     */
     double integratedAutocorrelationTime(int range) const;
};

#endif /* DataArray_hpp */