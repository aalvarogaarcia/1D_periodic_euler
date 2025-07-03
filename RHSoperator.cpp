#include "RHSoperator.h"

template<class T>
RHSOperator<T>::RHSOperator()
{

}


template<class T>
RHSOperator<T>::~RHSOperator()
{

}

/*
template<class T>
Central1D<T>::Central1D(DataStruct<T> &_U, 
                     DataStruct<T> &_mesh, 
                     FluxFunction<T> &_F):
U(_U), mesh(_mesh), F(_F)
{
  RHS.setSize(_U.getSize());
}

template<class T>
Central1D<T>::~Central1D()
{

}


SE DEBE AJUSTAR LA CLASE 1D

template<class T>
void Central1D<T>::evalRHS(DataStruct<T> &Uin)
{
  // the BC should be included in the mesh
  // momentarily done here by hand
  T *dataRHS = RHS.getData();
  const T *dataU = Uin.getData();
  const T *dataMesh = mesh.getData();
  const int len = U.getSize();

  for(int j = 0; j < len; j++)
  {
    T dx;
    if(j == 0)
    {
      dx = dataMesh[len-1] - dataMesh[len-2];
      dx += dataMesh[1] - dataMesh[0];
      dataRHS[0] = -(F.computeFlux(dataU[1]) - F.computeFlux(dataU[len-2]))/dx;
    }
    else
    {
      dx = dataMesh[j+1] - dataMesh[j-1];
      dataRHS[j] = -(F.computeFlux(dataU[j+1]) - F.computeFlux(dataU[j-1]))/dx;
    }
  }

  dataRHS[len-1] = dataRHS[0];
}

template<class T>
void Central1D<T>::eval()
{
  evalRHS(U);
}

template<class T>
void Central1D<T>::eval(DataStruct<T> &Uin)
{
  evalRHS(Uin);
}

template<class T>
DataStruct<T>& Central1D<T>::ref2RHS()
{
  return RHS;
}
*/




// ***IMPLEMENTACION DEL CONSTRUCTOR***
template<class T>
Central1D<T>::Central1D(DataStruct<T> &grid, EulerFlux<T> &flux) :
    xj(grid), flux_function(flux),
    F_internal(grid.getSize()*3),
    RHS(grid.getSize()*3)
//    f_rho(grid.getSize()), f_rho_u(grid.getSize()), f_rho_E(grid.getSize()),
//    RHS_rho(grid.getSize()), RHS_rho_u(grid.getSize()), RHS_rho_E(grid.getSize())
{
    // El cuerpo puede estar vacío, la inicialización se hace arriba
}

// --- IMPLEMENTACIÓN DEL DESTRUCTOR ---
template<class T>
Central1D<T>::~Central1D()
{
    // El cuerpo puede estar vacío, pero la función debe existir.
}


// *** IMPLEMENTACION DEL METODO EVAL***

template<class T>
void Central1D<T>::eval(const DataStruct<T>& U_in)
{
    int numPoints = xj.getSize();
    T dx = xj.getData()[1] - xj.getData()[0];


    flux_function.computeFlux(U_in, F_internal);

    const T* F_data = F_internal.getData();
    T* RHS_data = RHS.getData();


    for(int i = 1; i < numPoints - 1; i++)
    {
        for (int var = 0; var < 3; ++var) // Bucle para rho, rho_u, rho_E
        {
            int idx_i   = i * 3 + var;
            int idx_im1 = (i - 1) * 3 + var;
            int idx_ip1 = (i + 1) * 3 + var;
            RHS_data[idx_i] = -(F_data[idx_ip1] - F_data[idx_im1]) / (2.0 * dx);
        }
    }


    for (int var = 0; var < 3; ++var) {
        RHS_data[0 * 3 + var] = 0;
        RHS_data[(numPoints - 1) * 3 + var] = 0;
    }
}


/*
template<class T>
void Central1D<T>::eval(const DataStruct<T> &rho, const DataStruct<T> &rho_u, const DataStruct<T> &rho_E)
{
    int size = xj.getSize();
    T dx = xj.getData()[1] - xj.getData()[0];

    // 1. Calcular los 3 componentes del flujo en todos los puntos
    flux_function.computeFlux(rho, rho_u, rho_E, f_rho, f_rho_u, f_rho_E);

    const T* F1 = f_rho.getData();
    const T* F2 = f_rho_u.getData();
    const T* F3 = f_rho_E.getData();
    T* R1 = RHS_rho.getData();
    T* R2 = RHS_rho_u.getData();
    T* R3 = RHS_rho_E.getData();

    // 2. Calcular la derivada centrada para los puntos interiores
    for(int i = 1; i < size - 1; i++)
    {
        R1[i] = -(F1[i+1] - F1[i-1]) / (2.0 * dx);
        R2[i] = -(F2[i+1] - F2[i-1]) / (2.0 * dx);
        R3[i] = -(F3[i+1] - F3[i-1]) / (2.0 * dx);
    }

    // 3. Aplicar condiciones de contorno periódicas
    // Punto i = 0
    R1[0] = -(F1[1] - F1[size-1]) / (2.0 * dx);
    R2[0] = -(F2[1] - F2[size-1]) / (2.0 * dx);
    R3[0] = -(F3[1] - F3[size-1]) / (2.0 * dx);

    // Punto i = size - 1
    R1[size-1] = -(F1[0] - F1[size-2]) / (2.0 * dx);
    R2[size-1] = -(F2[0] - F2[size-2]) / (2.0 * dx);
    R3[size-1] = -(F3[0] - F3[size-2]) / (2.0 * dx);
}
*/









template class RHSOperator<float>;
template class RHSOperator<double>;

template class Central1D<float>;
template class Central1D<double>;


