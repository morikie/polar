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

	double getDseLocationTvalue(std::string pas, size_t pos) const;
	double getDseUcontentTvalue(std::string pas, double uContent) const;
	double getUseUcontentTvalue(std::string pas, double uContent) const;

public:
	/**
	 * Class that holds values to calculate the truth value for a certain DSE location.
	 */
	class DseLocation {
	public:
		typedef std::pair<size_t, size_t> range;
		
		typedef double intercept;
		typedef double slope;
		typedef std::pair<slope, intercept> straight;
		
	private:
		//range where TV is something between 0 and 1; positive slope of the (acute) trapezoid side (usually left side)
		range positiveIntermediate;
		//range where TV is something between 0 and 1; negative slope of the (acute) trapzoid side (usually right side)
		range negativeIntermediate; 
		//double threshold;

		straight positiveStraight;
		straight negativeStraight;

	public:
		DseLocation(range p, range n);
		DseLocation();
		~DseLocation();
		
		range getLeftRange() const;
		range getRightRange() const;
		straight getLeftStraight() const;
		straight getRightStraight() const;

	private:
		void calcStraights();


	};

	/**
	 * Class that holds values to calculate the truth value for a certain uracil content.
	 */
	class UracilContent {
	public:
		typedef double intercept;
		typedef double slope;
		typedef std::pair<slope, intercept> straight;

	private:
		double lowerBound;
		double upperBound;
		double maxTruthValue;

		straight intermediate;

	public:
		UracilContent(double lB, double uB, double maxTv);
		UracilContent();
		~UracilContent();
		
		double getLowerBound() const;
		double getUpperBound() const;
		double getMaxTruthValue() const;
		straight getStraight() const;

	private:
		void calcStraight();

	};

protected:
	static std::unordered_map<std::string, Utr3FinderFuzzy::DseLocation> dseLocMap;
	static std::unordered_map<std::string, Utr3FinderFuzzy::UracilContent> dseUracilMap;
	static std::unordered_map<std::string, Utr3FinderFuzzy::UracilContent> useUracilMap;

};      


#endif /* __UTR3FINDERFUZZY_HPP__ */

