#ifndef _RHS_OPERATOR
#define _RHS_OPERATOR

#include "DataStructs.h"
#include "FluxFunctions.h"

template<class T>
class RHSOperator
{
public:
  RHSOperator();
  ~RHSOperator();

  virtual void eval() = 0;
  virtual void eval(DataStruct<T> &Uin) = 0;

  virtual DataStruct<T>& ref2RHS() = 0;
};


/*
template<class T>
class Central1D : public RHSOperator<T>
{
private:

  // structure containing the RHS values
  DataStruct<T> RHS;
  
  // reference to current solution
  DataStruct<T> &U;

  // reference to mesh 
  // TODO: change to a mesh structure
  DataStruct<T> &mesh;

  // reference to flux function
  FluxFunction<T> &F;

  void evalRHS(DataStruct<T> &Uin);

public:
  Central1D(DataStruct<T> &_U, DataStruct<T> &_mesh, FluxFunction<T> &_F);
  ~Central1D();

  virtual void eval();
  virtual void eval(DataStruct<T> &Uin);

  virtual DataStruct<T>& ref2RHS();

};
*/


template<class T>
class Central1D
{
private:
    DataStruct<T> &xj;
    EulerFlux<T>  &flux_function;
    DataStruct<T> f_rho, f_rho_u, f_rho_E;
    DataStruct<T> RHS_rho, RHS_rho_u, RHS_rho_E;

public:
    Central1D(DataStruct<T> &grid, EulerFlux<T> &flux);
    ~Central1D();
    void eval(const DataStruct<T> &rho, const DataStruct<T> &rho_u, const DataStruct<T> &rho_E);
    DataStruct<T>& ref2RHS_rho() { return RHS_rho; }
    DataStruct<T>& ref2RHS_rho_u() { return RHS_rho_u; }
    DataStruct<T>& ref2RHS_rho_E() { return RHS_rho_E; }
};



#endif // _RHS_OPERATOR
