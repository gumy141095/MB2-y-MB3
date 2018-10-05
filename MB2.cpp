#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <ctime>
//#include <numeric>
#include <math.h>
//#include <random>
//#include <chrono>
#include <fstream>

#include "Seb.h"
#include "Seb_debug.h"

using namespace std;

typedef double FT;
typedef Seb::Point<FT> Point;
typedef std::vector<Point> PointVector;
typedef Seb::Smallest_enclosing_ball<FT> Miniball; 

double Distancia( double *a,  double *b, int dim);
void Diferencia( double *a,  double *b,  double *dif, int dim);
//void Cap_Dim(int &dim);

int main(int argn,char **argv) {
    
    clock_t t;
    double tiempo, dist, max, radio, producto_ab, dist_ac, dist_bc, dist_12, dist_13, dist_23;
    double lambda_1, lambda_2, deter;
   	int ren=3, experimentos=10000;
   	int a, b, c, col;
   	
   	double *tiempos_1 = NULL, *tiempos_2 = NULL, *desviacion = NULL, *promedio = NULL;
   	
   	tiempos_1 = new double [experimentos+1];
   	tiempos_2 = new double [experimentos+1];
   	desviacion = new double [2];
   	promedio = new double [2];
   	
	
	PointVector S;
	col = std::atoi(argv[1]);
    col = col*50;
    //Cap_Dim(col);
    
    double **N = NULL;
    
    N = new double*[ren+2];
	for (int i = 0; i < ren+2; ++i)
    	N[i] = new double[col];
    
    vector<double> coords(col);
    srand(time(NULL));
    for(int q=0 ; q<experimentos ; ++q){
        /*
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        std::default_random_engine generator (seed);

        std::uniform_real_distribution<double> distribution(0.0,10.0);
        */
        S.clear();
        for(int i=0 ; i<ren ; ++i){
            for(int j=0 ; j<col ; ++j){
                N[i][j] = static_cast<FT>(10.0*rand()/RAND_MAX);
                coords[j] = N[i][j];
            }
            S.push_back(Point(col,coords.begin()));
        }
                            
        t = clock();
        
        max = 0;
        
        dist_12 = Distancia(N[0],N[1],col);
        dist_13 = Distancia(N[0],N[2],col);
        dist_23 = Distancia(N[1],N[2],col);
        
        if(dist_12 < dist_13){
            b = 2;
            if(dist_13 < dist_23){
                a = 1;
                c = 0;
                max = dist_23;
                dist_ac = dist_12;
                dist_bc = dist_13;
            }else{
            	a = 0;
            	c = 1;
                max = dist_13;
                dist_ac = dist_12;
                dist_bc = dist_23;
            }
        }else{
        	b = 1;
            if(dist_12 < dist_23){
            	a = 2;
            	c = 0;
	            max = dist_23;
                dist_ac = dist_13;
                dist_bc = dist_12;
            }else{
            	a = 0;
            	c = 2;
                max = dist_12;
                dist_ac = dist_13;
                dist_bc = dist_23;
            }
        }
        
        
        for(int i=0 ; i<ren-1 ; ++i)
            for (int j=ren-1 ; j>i ; --j){
                dist = Distancia(N[i],N[j],col);
                if(dist>max){
                    max=dist;
                    a=i;
                    b=j;
                    c=((j+i)*2)%3;
                }
            }
        
        dist_ac = Distancia(N[a],N[c],col);
        dist_bc = Distancia(N[b],N[c],col);
        
        // CAMBIAR EL PRODUCTO ESCALAR
        producto_ab = 0.5 * (dist_ac + dist_bc - max);
        
        if(producto_ab<=0){
            radio = sqrt(max)/2;
            for(int i=0 ; i<col ; ++i)
                N[ren][i] = (N[a][i] + N[b][i])/2;            
            
        }else{
            
            N[ren][0] = sqrt(dist_ac)/2.0;
            N[ren][1] = sqrt(dist_ac)*(dist_bc-producto_ab)/(2*sqrt(dist_ac*dist_bc-(producto_ab*producto_ab)));
            
            radio = sqrt((N[ren][0]*N[ren][0])+(N[ren][1]*N[ren][1]));
            
            deter = 2 * ( ( dist_ac * dist_bc ) - ( producto_ab*producto_ab ) );
            
            lambda_1 = ( dist_bc * ( dist_ac - producto_ab ) ) / deter;
            lambda_2 = ( dist_ac * ( dist_bc - producto_ab ) ) / deter;
            
            for(int i=0 ; i<col ; ++i)
                N[ren][i] = ( lambda_1 * N[a][i] ) + ( lambda_2 * N[b][i] ) + ( ( 1.0 - lambda_1 - lambda_2 ) * N[c][i] );
            
        }
        
        tiempo = ( clock() - t )/( double)CLOCKS_PER_SEC;
        tiempos_1[q] = tiempo;
        promedio[0]+=tiempo;
        
        t = clock();
        
        Miniball mb(col, S);
        FT rad = mb.radius();
        Miniball::Coordinate_iterator center_it = mb.center_begin();
        
        tiempo = ( clock() - t )/( double)CLOCKS_PER_SEC;
        tiempos_2[q] = tiempo;
        promedio[1]+=tiempo; 
        
    }
    
    for (int i = 0; i < ren+2; ++i)
    	delete [] N[i];
	delete [] N;
    
    for(int i=0 ; i<2 ; ++i){
    	promedio[i] /= experimentos;
	    tiempo=0;
	    for(int j=0 ; j<experimentos ; ++j){
	    	tiempo += (tiempos_1[j]-promedio[i])*(tiempos_1[j]-promedio[i]);
		}
		tiempo /= experimentos;
		desviacion[i] = sqrt(tiempo);
	}
	
	delete[] tiempos_1;
   	delete[] tiempos_2;
   	
	/*
    cout << col << '\t' << scientific << promedio[0] << " s"
				<< '\t' << scientific << desviacion[0] << " s" 
		 		<< '\t' << scientific << promedio[1] << " s"
				<< '\t' << scientific << desviacion[1] << " s"  
				<< endl;
    */
	ofstream salida;
   
   	salida.open("Tiempos.csv", std::ofstream::out | std::ofstream::app );
   
   	salida << col << ',' << scientific << promedio[0] << ','
						 << scientific << desviacion[0] << ',' 
		 			     << scientific << promedio[1] << ','
				 	     << scientific << desviacion[1] << endl;
   
    salida.close();					
   	
   	
   	delete[] desviacion;
   	delete[] promedio;
 	
    return 0;
}

 double Distancia( double *a,  double *b, int dim){
	long double dist=0;
	
	for( int i = 0 ; i < dim ; ++i)
		dist += (a[i]-b[i])*(a[i]-b[i]);
		
	return dist;
}

void Diferencia( double *a,  double *b,  double *dif, int dim){

	for( int i = 0 ; i < dim ; ++i)
		dif[i] = a[i]-b[i];

}
/*
void Cap_Dim(int &dim){
   
    ifstream entrada;
   
    entrada.open("value_1.txt");
    
    
    if(!entrada){
        cout << "Error: no se pudo abrir el archivo..." << endl;
        return;
    }
    
   
    entrada >> dim;

    entrada.close();
    dim = 3;
}
*/
