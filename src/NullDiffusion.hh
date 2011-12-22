#ifndef NULL_DIFFUSION_HH
#define NULL_DIFFUSION_HH

#include "Diffusion.hh"

class NullDiffusion : public Diffusion
{
 public:
   NullDiffusion(){};
   
   void calc(const std::vector<double>& Vm, std::vector<double>& dVm){};
};


#endif