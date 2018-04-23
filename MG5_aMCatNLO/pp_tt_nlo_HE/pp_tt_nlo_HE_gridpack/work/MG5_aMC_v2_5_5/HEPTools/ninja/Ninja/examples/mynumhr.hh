// This -*- C++ -*- file has been generated by the Python package
// NinjaNumGen, which is distributed with the Ninja library.

#ifndef DIAGRAM_HH
#define DIAGRAM_HH

#include <ninja/num_defs.hh>

#define NINJA_NUM Diagram
#define NINJA_NUM_NAMESPACE 0

#if NINJA_NUM_NAMESPACE
namespace NINJA_NUM_NAMESPACE {
#endif


  class NINJA_NUM : public ninja::Numerator {

  public:

    virtual ninja::Complex evaluate(const ninja::ComplexMomentum & q,
                                    const ninja::Complex & muq,
                                    int cut,
                                    const ninja::PartitionInt part[]);

    virtual void muExpansion(const ninja::ComplexMomentum v_perp[],
                             const ninja::PartitionInt part[],
                             ninja::Complex c[]);

    virtual void t3Expansion(const ninja::ComplexMomentum & a,
                             const ninja::ComplexMomentum & e3,
                             const ninja::ComplexMomentum & e4,
                             const ninja::Complex & param,
                             int mindeg,
                             int cut, const ninja::PartitionInt part[],
                             ninja::Complex c[]);

    virtual void t2Expansion(const ninja::ComplexMomentum & a0,
                             const ninja::ComplexMomentum & a1,
                             const ninja::ComplexMomentum & e3,
                             const ninja::ComplexMomentum & e4,
                             const ninja::Complex param[],
                             int mindeg,
                             int cut, const ninja::PartitionInt part[],
                             ninja::Complex c[]);
  
  public:
    ninja::ComplexMomentum v0, v1, v2, v3, v4, v5;

  private:
    // Add other private methods and data here

  };


#if NINJA_NUM_NAMESPACE
} // namespace NINJA_NUM_NAMESPACE
#endif

#undef NINJA_NUM
#undef NINJA_NUM_NAMESPACE

#endif // DIAGRAM_HH
