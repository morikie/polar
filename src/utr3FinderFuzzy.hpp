#ifndef __UTR3FINDERFUZZY_HPP__
#define __UTR3FINDERFUZZY_HPP__

#include <string>
#include <unordered_map>
#include "seqStruct.hpp"
#include "utr3Finder.hpp"


class Utr3FinderFuzzy : public Utr3Finder {
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
		//range where TV is something between 0 and 1; positive slope of the (acute) trapezoid side ("left" side)
		range positiveIntermediate;
		//range where TV is something between 0 and 1; negative slope of the (acute) trapezoid side ("right" side)
		range negativeIntermediate; 

		straight positiveStraight;
		straight negativeStraight;

	public:
		DseLocation(range p, range n);
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
		~UracilContent();
		
		double getLowerBound() const;
		double getUpperBound() const;
		double getMaxTruthValue() const;
		straight getStraight() const;

	private:
		void calcStraight();

	};

public:
	typedef std::string motifSequence;
	typedef std::unordered_map<motifSequence, DseLocation> pasToDseLocMap;
	typedef std::unordered_map<motifSequence, UracilContent> pasToUcontentMap;

public:
	static pasToDseLocMap dseLocMap;
	static pasToDseLocMap dseDinucleoLocMap;
	static pasToUcontentMap dseUracilMap;
	static pasToUcontentMap dseDinucleoMap;
	static pasToUcontentMap dseShortUracilMap;
	static pasToUcontentMap useUracilMap;
	static std::unordered_map<motifSequence, double> thresholdMap;

protected:
	std::string reversedSeq;
	const size_t windowSize = 10;
	bool searchBackward;

public: 
	Utr3FinderFuzzy(const SeqStruct & sSt, const bool & searchBackward);
	virtual ~Utr3FinderFuzzy();

	virtual bool isMutationInMotif() const override;
	virtual std::string getSequence() const override;
	std::string getRevComplementSeq() const;
	virtual std::string getMotifSequence(const Utr3FinderResult & result) const override;
	virtual std::vector<Utr3FinderResult> getPolyaMotifPos() const override;
	void setThresholdMap(std::unordered_map<std::string, double> & map);
	virtual void writeInfo() const override;

protected:
	virtual void findPolyaMotif() override;

	double locationTvalue(const std::string & pas, const size_t & pos, const pasToDseLocMap & map) const;
	double nucleotideContentTvalue(const std::string & pas, const double & uContent, const pasToUcontentMap & map) const;
	
	double calcCombinedDseTvalue(const size_t & pos, const std::string & seq);
	double calcCombinedDinucleoDseTvalue( const size_t & pos, const std::string & seq);
	double calcDseShortTvalue(const size_t & pos, const std::string & seq);
	double calcUseTvalue(const size_t & pos, const std::string & seq);
};      


#endif /* __UTR3FINDERFUZZY_HPP__ */

