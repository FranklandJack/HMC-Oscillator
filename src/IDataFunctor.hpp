#ifndef IDataFunctor_hpp
#define IDataFunctor_hpp
class IDataFunctor
{
public:
	virtual double operator()(double) const;
};

#endif /* IDataFunctor */