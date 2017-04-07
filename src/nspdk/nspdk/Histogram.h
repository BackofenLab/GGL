/* -*- mode:c++ -*- */
#ifndef SGM_NSPDK_HISTOGRAM_H
#define SGM_NSPDK_HISTOGRAM_H


namespace nspdk {

using namespace std;

class HistogramClass{
  friend ostream& operator<<(ostream& out, const HistogramClass& aHC);
public:
  HistogramClass();
  HistogramClass(const HistogramClass& aH);
  void Insert(unsigned aValue);
  void Insert(unsigned aBin, double aValue);
  void Add(const HistogramClass& aH);
  ostream& Output(ostream& out)const;
  unsigned Size()const;
public:
  map<unsigned,double> mHistogram;
};

//---------------------------------------------------------------------------------
class SecondOrderHistogramClass{
  friend ostream& operator<<(ostream& out, const SecondOrderHistogramClass& aSOH);
public:
  SecondOrderHistogramClass();
  SecondOrderHistogramClass(const SecondOrderHistogramClass& aSOH);
  void Insert(unsigned aSecondOrderBin, unsigned aValue);
  void Insert(unsigned aSecondOrderBin, unsigned aBin, double aValue);
  void Add(const SecondOrderHistogramClass& aSOH);
  void Add(unsigned aSecondOrderBin,const HistogramClass& aH);
  ostream& Output(ostream& out)const;
  unsigned Size()const;
public:
  map<unsigned,HistogramClass> mSecondOrderHistogram;
};

//---------------------------------------------------------------------------------
class ThirdOrderHistogramClass{
  friend ostream& operator<<(ostream& out, const ThirdOrderHistogramClass& aTOH);
public:
  ThirdOrderHistogramClass();
  ThirdOrderHistogramClass(const ThirdOrderHistogramClass& aSOH);
  void Insert(unsigned aThirdOrderBin, unsigned aSecondOrderBin, unsigned aValue);
  void Insert(unsigned aThirdOrderBin, unsigned aSecondOrderBin, unsigned aBin, double aValue);
  void Add(const ThirdOrderHistogramClass& aTOH);
  void Add(unsigned aThirdOrderBin,const SecondOrderHistogramClass& aSOH);
  ostream& Output(ostream& out)const;
  unsigned Size()const;
public:
  map<unsigned,SecondOrderHistogramClass> mThirdOrderHistogram;
};

} // namespace

#endif
