#ifndef __UTR3FINDERFUZZY_HPP__
#define __UTR3FINDERFUZZY_HPP__

#include <string>
#include <unordered_map>
#include "seqStruct.hpp"
#include "utr3Finder.hpp"


class Utr3FinderFuzzy : public Utr3Finder {
public: 
	Utr3FinderFuzzy(const SeqStruct & sSt);
	virtual ~Utr3FinderFuzzy();

	virtual bool isMutationInMotif() const override;
	virtual std::string getSequence() const override;
	virtual std::vector<std::string> getMotifSequence() const override;
	virtual std::vector<size_t> getPolyaMotifPos() const override;
	virtual void writeInfo() const override;

protected:
	virtual void findPolyaMotif() override;

public:
	/**
	 * Class that holds values to calculate the truth value for a certain DSE location.
	 */
	class DseLocation {
	private:
		typedef std::pair<size_t, size_t> range;
		
		range positiveIntermediate; //range where TV is something between 0 and 1; positive slope of the (acute) trapezoid side (usually left side)
		range negativeIntermediate; //...; negative slope of the (acute) trapzoid side (usually right side)
		//double threshold;

	public:
		DseLocation(range p, range n);
		~DseLocation();
	
	private:
		typedef double intercept;
		typedef double slope;
		typedef std::pair<slope, intercept> straight;

		straight positiveStraight;
		straight negativeStraight;

		void calcStraights();
	};

	/**
	 * Class that holds values to calculate the truth value for a certain uracil content.
	 */
	class UracilContent {
	private:
		double upperBound;
		double lowerBound;
		double maxTruthValue;

	public:
		UracilContent(double uB, double lB, double maxTv);
		~UracilContent();

	private:
		typedef double intercept;
		typedef double slope;
		typedef std::pair<slope, intercept> straight;

		straight intermediate;

		void calcStraight();
	};

protected:
	static std::unordered_map<std::string, DseLocation> dseLocMap;
	static std::unordered_map<std::string, UracilContent> dseUracilMap;
	static std::unordered_map<std::string, UracilContent> useUracilMap;

	double getDseLocationTvalue(std::string pas, size_t pos);
	double getDseUcontentTvalue(std::string pas, double uContent);
	double getUseUcontentTvalue(std::string pas, double uContent);

};      


#endif /* __UTR3FINDERFUZZY_HPP__ */

