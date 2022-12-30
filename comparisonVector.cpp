#include "DataStructs.h"
#include <vector>
#include <iostream>
#include <mpi.h>

#ifdef _DOUBLE_
#define FTYPE double
#else
#define FTYPE float
#endif

using namespace std;

int main(int narg, char **argv)
{
  int Narray = 1;

  if(narg < 2)
  {
    cout << "usage: compatisonVector.p Narray" << endl;
  }else
  {
    Narray = stoi( argv[1] );
  }

  // create vectors
  vector<FTYPE> stdVec(Narray);
  DataStruct<FTYPE> vecStruct(Narray);
  FTYPE* dataVec = vecStruct.getData();

  for(int n = 0; n < Narray; n++)
  {
    stdVec[n] = n;
    dataVec[n] = n;
  }

  double Wtimes, tmp = 0.;

  // measure accesses DataStruct
  tmp = 0;
  {
    FTYPE* data = vecStruct.getData();
    Wtimes = MPI_Wtime();
    for(int n = 0; n < vecStruct.getSize(); n++)
    {
      tmp += data[n];
      //data[n] = n;
    }
  }
  Wtimes = MPI_Wtime() - Wtimes;
  cout << tmp << " Accesses DataStruct " << Wtimes << " s. ";

  // measure accesses std::vector
  tmp = 0;
  Wtimes = MPI_Wtime();
  {
    for(int n = 0; n < stdVec.size(); n++)
    {
      tmp += stdVec[n];
      //stdVec[n] = n;
    }
  }
  Wtimes = MPI_Wtime() - Wtimes;
  cout << tmp << " Accesses std::vector " << Wtimes << " s." << endl;
  

  return 0;
}