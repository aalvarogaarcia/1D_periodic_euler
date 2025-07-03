#include "FluxFunctions.h"
#include "DataStructs.h"

template<class T>
FluxFunction<T>::FluxFunction()
{

};



template<class T>
LinearFlux<T>::LinearFlux()
{
  c = 1.;
};


template<class T>
EulerFlux<T>::EulerFlux(T g) : gamma(g) {}

template<class T>
EulerFlux<T>::~EulerFlux() {}



template<class T>
void LinearFlux<T>::computeFlux(DataStruct<T> &U, DataStruct<T> &F)
{
  T *dataU = U.getData();
  T *dataF = F.getData();

  for(int n = 0; n < U.getSize(); n++)
  {
    dataF[n] = c * dataU[n];
  };
};





template<class T>
T LinearFlux<T>::computeFlux(const T &Ui)
{
  return c * Ui;
};





//EulerFlux empleando computeFlux
template<class T>
//void EulerFlux<T>::computeFlux(const DataStruct<T> &rho,
//                               const DataStruct<T> &rho_u,
//                               const DataStruct<T> &rho_E,
//                               DataStruct<T> &f_rho,
//                               DataStruct<T> &f_rho_u,
//                               DataStruct<T> &f_rho_E)
void EulerFlux<T>::computeFlux(const DataStruct<T> &U,
				DataStruct<T> &F)
{
//    const T* data_rho = rho.getData();
//    const T* data_rho_u = rho_u.getData();
//    const T* data_rho_E = rho_E.getData();
    const T*dataU = U.getData();

//    T* data_f_rho = f_rho.getData();
//    T* data_f_rho_u = f_rho_u.getData();
//    T* data_f_rho_E = f_rho_E.getData();
    T*dataF = F.getData();

//    for(int i = 0; i < rho.getSize(); i++)
    for(int i = 0; i < U.getSize()/3;i++)
    {
//        T rho_i   = data_rho[i];
//        T rho_u_i = data_rho_u[i];
//        T rho_E_i = data_rho_E[i];

        int base_idx = i * 3;
        T rho_i   = dataU[base_idx + 0];
        T rho_u_i = dataU[base_idx + 1];
        T rho_E_i = dataU[base_idx + 2];

        // Usamos un umbral pequeño para evitar división por cero
        if (rho_i <= 1e-9) {
         //   data_f_rho[i]   = 0;
         //   data_f_rho_u[i] = 0;
         //   data_f_rho_E[i] = 0;
            dataF = 0;
            continue;
        }

        T u_i = rho_u_i / rho_i;
        T p_i = (gamma - 1.0) * (rho_E_i - 0.5 * rho_i * u_i * u_i);

        // f1 = rho * u
//        data_f_rho[i] = rho_u_i;

        // f2 = rho * u^2 + p
//        data_f_rho_u[i] = rho_i * u_i * u_i + p_i;

        // f3 = u * (rhoE + p)
//        data_f_rho_E[i] = u_i * (rho_E_i + p_i);

        dataF[base_idx + 0] = rho_u_i;
        dataF[base_idx + 1] = rho_i * u_i * u_i + p_i;
        dataF[base_idx + 2] = u_i * (rho_E_i + p_i);


    }
}



template class FluxFunction<float>;
template class FluxFunction<double>;

template class LinearFlux<float>;
template class LinearFlux<double>;

template class EulerFlux<float>;
template class EulerFlux<double>;

