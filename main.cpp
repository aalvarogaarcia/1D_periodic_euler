#include <iostream> // std::cout
#include <fstream> // ofile
#include <string>
#include <math.h>
#include <iomanip> // set precision
#include <mpi.h>

#include "DataStructs.h"
#include "rk4.h"
#include "FluxFunctions.h"
#include "RHSoperator.h"

#ifdef _DOUBLE_
#define FLOATTYPE double
#else
#define FLOATTYPE float
#endif

// declare supporting functions
void write2File(DataStruct<FLOATTYPE> &X,
		DataStruct<FLOATTYPE> &U, //Se sustituye U por 3 variables (SoA)
		//DataStruct<FLOATTYPE> &rho,
		//DataStruct<FLOATTYPE> &rho_U,
		//DataStruct<FLOATTYPE> &rho_E,
		std::string name);
FLOATTYPE calcL2norm(DataStruct<FLOATTYPE> &u, DataStruct<FLOATTYPE> &uinit);


int main(int narg, char **argv)
{
  int numPoints =  80;
  //FLOATTYPE k = 2.; // wave number (NO SE NECESITA PARA RESOLVER EULER)

  //if(narg != 3) //Ahora solo se necesita el numero de puntos
  if(narg != 2)
  {
    std::cout<< "Wrong number of arguments. You should include:" << std::endl;
    std::cout<< "    Num points" << std::endl;
    //std::cout<< "    Wave number" << std::endl;
    return 1;
  }else
  {
    numPoints = std::stoi(argv[1]);
    //k         = std::stod(argv[2]);
  }

  // solution data
  DataStruct<FLOATTYPE> U(numPoints * 3), xj(numPoints);
  // DataStruct<FLOATTYPE> rho(numPoints), rho_U(numPoints), rho_E(numPoints); //Se declaran las 3 variables solucion
  // DataStruct<FLOATTYPE> xj(numPoints); //Se declara la malla

  // DataStruct<FLOATTYPE> f_rho(numPoints), f_rho_U(numPoints), f_rho_E(numPoints); //Se declaran las 3 variables para los flujos


  // flux function
  // LinearFlux<FLOATTYPE> lf; Se modifica (No se usa resolucion lineal)
  const FLOATTYPE gamma = 1.4; //Constante para las ecuaciones de Euler
  EulerFlux <FLOATTYPE> ef(gamma);


  // time solver
  RungeKutta4<FLOATTYPE> rk(U); //(Se resuelve para las 3 variables)
  // RungeKutta4<FLOATTYPE> rk_rho(rho);
  // RungeKutta4<FLOATTYPE> rk_rho_U(rho_U);
  // RungeKutta4<FLOATTYPE> rk_rho_E(rho_E);

  // Initial Condition
  FLOATTYPE *datax = xj.getData();
  FLOATTYPE *dataU = U.getData(); //(Pulso inicial de Euler)
  // FLOATTYPE *dataRho = rho.getData();
  // FLOATTYPE *dataRho_U = rho_U.getData();
  // FLOATTYPE *dataRho_E = rho_E.getData();

//	--- CONDICIONES INICIALES ESTATICAS ---
  const FLOATTYPE p_init = 1.8; // Presion inicial
  for(int j = 0; j < numPoints; j++)
  {
    int base_idx = j * 3;
    // xj
    datax[j] = FLOATTYPE(j)/FLOATTYPE(numPoints-1);
    FLOATTYPE rho = 1.0+0.2*exp(-(datax[j]-0.5)*(datax[j]-0.5)/0.01), rho_U = 0.2, rho_E = p_init/(gamma-1.0);
    // init Uj
    // dataU[j] = sin(k*2. * M_PI * datax[j]);

    // Pulso Gaussiano para la densidad
    // dataRho[j]=1.0+0.2*exp(-(datax[j]-0.5)*(datax[j]-0.5)/0.01);
    dataU[base_idx] = rho;

    // Fluido en reposo
    // dataRho_U[j] = 0.2;
    dataU[base_idx+1] = rho_U;

    // Energia de un fluido en reposo
    // dataRho_E[j] = p_init/(gamma-1.0);
    dataU[base_idx+2] = rho_E;
}
/*

//	--- CONDICIONES INICIALES DE CHOQUE ---

for (int j = 0; j < numPoints; j++)
{
    datax[j] = static_cast<FLOATTYPE>(j) / static_cast<FLOATTYPE>(numPoints - 1);

    if (datax[j] < 0.5) { // Estado izquierdo (alta presión)
        FLOATTYPE rho_left = 1.0;
        FLOATTYPE p_left = 1.0;
        FLOATTYPE u_left = 0.0;

        dataRho[j] = rho_left;
        dataRho_U[j] = rho_left * u_left;
        dataRho_E[j] = p_left / (gamma - 1.0) + 0.5 * rho_left * u_left * u_left;
    } else { // Estado derecho (baja presión)
        FLOATTYPE rho_right = 0.125;
        FLOATTYPE p_right = 0.1;
        FLOATTYPE u_right = 0.0;

        dataRho[j] = rho_right;
        dataRho_U[j] = rho_right * u_right;
        dataRho_E[j] = p_right / (gamma - 1.0) + 0.5 * rho_right * u_right * u_right;
    }
}
*/

  DataStruct<FLOATTYPE> Uinit(numPoints*3);
  Uinit = U;  //(Manejamos Rho)
  //DataStruct <FLOATTYPE> rho_init(numPoints);
  //rho_init = rho;

  // Operator
  // Central1D<FLOATTYPE> rhs(xj,ef); (Ahora se emplean todas las variables nuevas)
  Central1D<FLOATTYPE> rhs(xj, ef);


  // FLOATTYPE CFL = 2.4; (Para Euler se recomienda CFL <1.0)
  FLOATTYPE CFL = 0.8;
  FLOATTYPE max_vel = 1.2; // Estimacion de la velocidad maxima de onda
  FLOATTYPE dt = CFL*(datax[1]-datax[0])/max_vel;

  // Output Initial Condition
  write2File(xj, U, "initialCondition.csv"); //(Output de las 3 variables)
  // write2File(xj, rho, rho_U, rho_E, "initialCondition.csv");

  //FLOATTYPE t_final = 1.;
  FLOATTYPE t_final = 0.2;
  FLOATTYPE time = 0.;
  DataStruct<FLOATTYPE> Ui(U.getSize()); // temp. data (Se incluyen todos los datos)
  //DataStruct<FLOATTYPE> rho_i(numPoints), rho_U_i(numPoints), rho_E_i(numPoints);

  // init timer
  double compTime = MPI_Wtime();

  // main loop
  while(time < t_final)
  {
    if(time+dt >= t_final) dt = t_final - time;

    // take RK step
    rk.initRK();
    // rk_rho.initRK();
    // rk_rho_U.initRK();
    // rk_rho_E.initRK();
    // for(int s = 0; s < rk_rho.getNumSteps(); s++)
    for(int s = 0;s < rk.getNumSteps(); s++)
    {
      rk.stepUi(dt);
      //rk_rho.stepUi(dt);
      //rk_rho_U.stepUi(dt);
      //rk_rho_E.stepUi(dt);


      Ui = *rk.currentU();
      //rho_i = *rk_rho.currentU();
      //rho_U_i = *rk_rho_U.currentU();
      //rho_E_i = *rk_rho_E.currentU();

      rhs.eval(Ui);
      //rhs.eval(rho_i, rho_U_i, rho_E_i);


      rk.setFi(rhs.ref2RHS());
      //rk_rho.setFi(rhs.ref2RHS_rho());
      //rk_rho_U.setFi(rhs.ref2RHS_rho_u());
      //rk_rho_E.setFi(rhs.ref2RHS_rho_E());


    }
    rk.finalizeRK(dt);
    //rk_rho.finalizeRK(dt);
    //rk_rho_U.finalizeRK(dt);
    //rk_rho_E.finalizeRK(dt);


    time += dt;
  }

  // finishe timer
  compTime = MPI_Wtime() - compTime;

  write2File(xj, U, "final.csv");
  //write2File(xj, rho, rho_U, rho_E, "final.csv");

  // L2 norm
  FLOATTYPE err = calcL2norm(Uinit, U);
  //FLOATTYPE err = calcL2norm(rho_init, rho);
  std::cout << std::setprecision(4) << "Comp. time: " << compTime;
  //std::cout << " sec. Error: " << err/k;
  std::cout << " sec.L2 Error (rho): " << err;
  //std::cout << " kdx: " << k*datax[1]*2.*M_PI; CORRESPONDE AL PROBLEMA ORIGINAL
  std::cout << std::endl;

  return 0;

}


// ==================================================================
// AUXILIARY FUNCTIONS
// ==================================================================
void write2File(DataStruct<FLOATTYPE> &X, DataStruct<FLOATTYPE> &U, std::string name)
//void write2File(DataStruct<FLOATTYPE> &X,
//		DataStruct<FLOATTYPE> &rho,
//		DataStruct<FLOATTYPE> &rho_U,
//		DataStruct<FLOATTYPE> &rho_E,
//		std::string name)
{
  std::ofstream file;
  file.open(name,std::ios_base::trunc);
  if(!file.is_open())
  {
    std::cout << "Couldn't open file for Initial Condition" << std::endl;
    exit(1);
  }

  for(int j = 0; j < U.getSize()/3; j++)
  {
    int base_idx = j * 3;
    file << X.getData()[j] << " ," << U.getData()[base_idx] << " ," << U.getData()[base_idx+1] << " ," <<  U.getData()[base_idx+2] << std::endl;
  }

  file.close();
}

FLOATTYPE calcL2norm(DataStruct<FLOATTYPE> &u, DataStruct<FLOATTYPE> &uinit)
{
  FLOATTYPE err = 0.;
  const FLOATTYPE *dataU = u.getData();
  const FLOATTYPE *dataInit = uinit.getData();

  for(int n = 0; n < u.getSize()/3; n++)
  {
    int base_idx = n * 3;
    err += (dataU[base_idx] - dataInit[base_idx])*(dataU[base_idx] - dataInit[base_idx]);
  }

  //return sqrt( err );
  return sqrt(err/u.getSize()); //CALCULO ERROR RMS
}
