#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <cstring>
#include <fstream>

#include "Seb.h"
#include "Seb_debug.h"

using namespace std;
//Definimos las funciones a usar.
double Distancia(double a[], double b[], int dim);
double Encontrar_Max(double valor_1, double valor_2, int a, int &b);
void Diferencia(double a[], double b[], double dif[], int dim);
//void Cap_Dim(int &dim);

double MB2_Simple(double a[], double b[], double c[], double dist_ab, double dist_ac, double dist_bc, double centro[], double p[], double lambda[], double prod_esc);

//Codigo principal
int main(int argn,char **argv) {
	typedef double FT;
  	typedef Seb::Point<FT> Point;
  	typedef std::vector<Point> PointVector;
  	typedef Seb::Smallest_enclosing_ball<FT> Miniball;

   	int ren=4, col = std::atoi(argv[1]), experimentos=1000;
   	const int n=ren;
    col = col * 50;
    
    clock_t t;

    double **N = NULL;
    double *centro = NULL;
    double *tiempos_1 = NULL, *tiempos_2 = NULL, *desviacion = NULL, *promedio = NULL;
    
    tiempos_1 = new double [experimentos+1];
   	tiempos_2 = new double [experimentos+1];
   	desviacion = new double [2];
   	promedio = new double [2];
    
    promedio[0] = 0;
    promedio[1] = 0;
    
    double distancia, tiempo, prod_esc_temp;
    
    double norma_sq[2];

    //Se enlistan las posibles combinaciones de los puntos a,b,c,d.
    int Posibilidades[6][4]={	{0,1,2,3},
                                {0,2,1,3},
                                {0,3,1,2},
                                {1,2,0,3},
                                {1,3,0,2},
                                {2,3,0,1}
                            };
    srand(time(NULL));
    //Cap_Dim(col);
    for(int pruebas=0 ; pruebas<experimentos ; ++pruebas ){
        /* GENERAR PUNTOS AL AZAR EN EL RANGO 0-rango*/
        ////////////////////////////////////////////////////////////////////////
        // /*
        double rango=10;
        //inicializa una semilla random.
        
        //S.clear();

        N = new double*[ren];
        for (int i = 0; i < ren; ++i)
            N[i] = new double[col];

        PointVector S;
        vector<double> coords(col);

        for(int i=0 ; i<ren ; ++i){
            for(int j=0 ; j<col ; ++j){
                //Regresa un numero pseudo-aleatorio entre 0-RAND_MAX, se multiplica por 10 que es el rango deseado y se divide entre RAND_MAX.
                N[i][j] = static_cast<FT>(rango*rand()/RAND_MAX);
                coords[j] = N[i][j];
            }
            //Se agregan los puntos a la nube de puntos de uno en uno.
            S.push_back(Point(col,coords.begin()));
        }
        
        centro = new double [col];
        // */
        ////////////////////////////////////////////////////////////////////////

        /* LEER NUBE DE ARCHIVO */
        ////////////////////////////////////////////////////////////////////////
        /*
        ifstream entrada;
        //Se abre el archivo a leer.
        entrada.open("Point_Cloud.txt");
        /*
        if(!entrada){
            cout << "Error: no se pudo abrir el archivo..." << endl;
            return 1;
        }
        */
        /*
        entrada >> ren >> col;
        //Creamos la nube de puntos segun los valores obtenidos en el archivo de la cantidad de puntos y dimension.

        N = new double*[ren];
        for (int i = 0; i < ren; ++i)
            N[i] = new double[col];

        PointVector S;
        vector<double> coords(col);

        for(int i=0;i<ren;++i){
            for(int j=0;j<col;++j){
                entrada >> N[i][j];
                coords[j] = N[i][j];
            }
            //Se agregan los puntos a la nube de puntos de uno en uno.
            S.push_back(Point(col,coords.begin()));
        }

        entrada.close();
        */
        ////////////////////////////////////////////////////////////////////////

        /* INGRESAR NUBE MANUALMENTE */
        ////////////////////////////////////////////////////////////////////////
        /*
        col=2;

        N = new double*[ren];
        for (int i = 0; i < ren; ++i)
            N[i] = new double[col];

        PointVector S;
        vector<double> coords(col);

        N[0][0] = 7.99752;
        N[0][1] = 6.97715;

        coords[0] = N[0][0];
        coords[1] = N[0][1];

        S.push_back(Point(col,coords.begin()));

        N[1][0] = 9.5424;
        N[1][1] = 9.53253;

        coords[0] = N[1][0];
        coords[1] = N[1][1];

        S.push_back(Point(col,coords.begin()));

        N[2][0] = 5.65238;
        N[2][1] = 2.26605;

        coords[0] = N[2][0];
        coords[1] = N[2][1];

        S.push_back(Point(col,coords.begin()));

        N[3][0] = 1.97713;
        N[3][1] = 8.11956;

        coords[0] = N[3][0];
        coords[1] = N[3][1];

        S.push_back(Point(col,coords.begin()));
        */
        ////////////////////////////////////////////////////////////////////////
        t = clock();
        // Compute the miniball by inserting each value
        Miniball mb(col, S);
        FT rad = mb.radius();
            
        Miniball::Coordinate_iterator center_it = mb.center_begin();
        
        tiempo = ( clock() - t )/( double)CLOCKS_PER_SEC;
        tiempos_2[pruebas] = tiempo;
        promedio[1]+=tiempo;
        ////////////////////////////////////////////////////////////
        
        t = clock();
        //Definimos las variables para las distancias entre los puntos de nuestra nube.
        double max=0, dist_12, dist_13, dist_14, dist_23, dist_24, dist_34;

        //Calculamos las distancias entre los puntos enumerados segun se agregaron a la nube.
        dist_12 = Distancia(N[0],N[1],col);
        dist_13 = Distancia(N[0],N[2],col);
        dist_14 = Distancia(N[0],N[3],col);
        dist_23 = Distancia(N[1],N[2],col);
        dist_24 = Distancia(N[1],N[3],col);
        dist_34 = Distancia(N[2],N[3],col);

        //Creamos una matriz de distancias para facilitar el acceso a las distancias.
        double mat_dist[4][4]=	{{0,dist_12,dist_13,dist_14},
                                {dist_12,0,dist_23,dist_24},
                                {dist_13,dist_23,0,dist_34},
                                {dist_14,dist_24,dist_34,0}
                                };

        //Se definen las variables que se utilizaran durante el algoritmo.

        int orden=5, a, b, c, d;
        double radio;
        double centro_R3[5][3];
        double prod_esc[5];
        double difer_R3[3];
        double lambda[3];

        //Se usa una funcion iterativa para encontrar la combinacion adecuada segun los puntos mas alejados.
        radio = Encontrar_Max(dist_12,
                    Encontrar_Max(dist_13,
                        Encontrar_Max(dist_14,
                            Encontrar_Max(dist_23,
                                Encontrar_Max(dist_24,dist_34,4,orden),3,orden),2,orden),1,orden),0,orden);

        //Se nombran a, b, c y d con a-b la pareja de puntos mas alejados.
        a = Posibilidades[orden][0];
        b = Posibilidades[orden][1];
        c = Posibilidades[orden][2];
        d = Posibilidades[orden][3];

        /* CALCULAMOS ALGUNOS PRODUCTOS INTERIORES QUE SERAN UTILIZADOS MAS ADELANTE. NOS BASAMOS EN LA LEY DE COSENOS PARA SU CALCULO. */
        prod_esc[1] = 0.5 * (mat_dist[a][c] + mat_dist[b][c] - mat_dist[a][b]);
        prod_esc[2] = 0.5 * (mat_dist[a][d] + mat_dist[b][d] - mat_dist[a][b]);
        
        
        /////////////////////
        //Revisar caso 1.////
        /////////////////////
        if(prod_esc[1]<=0 && prod_esc[2]<=0){

            radio = sqrt(radio)/2;

            for(int i=0 ; i<col ; ++i){
                centro[i] = (N[a][i] + N[b][i])/2;
            }

            //cout << "Radio: " << radio << endl;

            //cout << "Centro: ( ";
            //for(int i=0 ; i<col ; ++i)
            //    cout << centro[i] << " ,";
            //cout << "\b)" << endl;

            //cout << "\n [a,b] es el soporte" << endl;
            // Output
            tiempo = ( clock() - t )/( double)CLOCKS_PER_SEC;
            tiempos_1[pruebas] = tiempo;
            promedio[0]+=tiempo;
            continue;
        }

        /* Para la recontruccion del tetraedo en R^3 renombramos los vertices c y d, en caso de ser necesario,
        * de tal manera que el triangulo [a,b,c] tenga mayor area que el triangulo [a,b,d].
        * Esto se realiza debido a que el valor del area del triangulo [a,b,c] aparece como denominador en algunas operaciones.
        */

        /* El area del triangulo [a,b,c] es sqrt(S_abc/16), donde S_abc es dado por la formula de Heron.
        * El area del triangulo [a,b,d] es sqrt(S_abd/16), donde S_abd es dado por la formula de Heron.
        */
        /* ES CONVENIENTE ESCRIBIR LO SIGUIENTE EN TERMINOS DE PRODUCTOS ESCALARES? */
        double S_abc = (mat_dist[a][b] + mat_dist[a][c] + mat_dist[b][c])*(mat_dist[a][b] + mat_dist[a][c] - mat_dist[b][c])*(mat_dist[a][b] - mat_dist[a][c] + mat_dist[b][c])*(mat_dist[a][c] + mat_dist[b][c] - mat_dist[a][b]);
        double S_abd = (mat_dist[a][b] + mat_dist[a][d] + mat_dist[b][d])*(mat_dist[a][b] + mat_dist[a][d] - mat_dist[b][d])*(mat_dist[a][b] - mat_dist[a][d] + mat_dist[b][d])*(mat_dist[a][d] + mat_dist[b][d] - mat_dist[a][b]);

        if(S_abd > S_abc){
            d = Posibilidades[orden][2];
            c = Posibilidades[orden][3];
            S_abc = S_abd;
            
            prod_esc_temp = prod_esc[1];
            prod_esc[1] = prod_esc[2];
            prod_esc[2] = prod_esc_temp;
        }


        /* LA CONDICION DE PARO DEBE ESTABLECERSE EN TERMINOS DE ALGUN CRITERIO DE COLINEALIDAD */
        /* REVISAR LA SIGUIENTE CONDICION DE PARO, ¿ES APROPIADA? */

        /* La precision double llega hasta alrededor de 1e-16 por lo que esa condicion no haria mucho,
        creo que seria recomendable ponerla en 1e-14 o alrededor.*/
        //if( S_abc <= 1e-14 ){
        if( S_abc <= 1e-15 ){
            radio = sqrt(radio)/2;

            //cout << "\n [a,b] es el soporte" << endl;

            for(int i=0 ; i<col ; ++i){
                centro[i] = (N[a][i] + N[b][i])/2;
            }
            
            /*
            cout << "Radio: " << radio << endl;

            cout << "Centro: ( ";
            for(int i=0 ; i<col ; ++i)
                cout << centro[i] << " ,";
            cout << "\b)" << endl;
            */
            
            tiempo = ( clock() - t )/( double)CLOCKS_PER_SEC;
            tiempos_1[pruebas] = tiempo;
            promedio[0]+=tiempo;
            continue;
        }
        
        /* CALCULAMOS ALGUNOS PRODUCTOS INTERIORES QUE SERAN UTILIZADOS MAS ADELANTE. NOS BASAMOS EN LA LEY DE COSENOS PARA SU CALCULO. */
        prod_esc[0] = 0.5 * (mat_dist[b][a] + mat_dist[c][a] - mat_dist[b][c]);
        prod_esc[3] = 0.5 * (mat_dist[a][b] + mat_dist[a][d] - mat_dist[b][d]);
        
        //Se calcula el producto escalar <b-a,c-a> en R^d, mediante la Ley de Cosenos:
        // prod_esc[0] = 0.5 * (mat_dist[b][a] + mat_dist[c][a] - mat_dist[b][c]);
        /* ¿ES APLICABLE ESTO EN EL ALGORITMO MB_2? */

        //Creamos la nube de puntos en R^3 que tendra al tetraedro.
        double NR3[4][3];
        double p[3];

        //Fijamos al punto a en el origen.
        NR3[a][0]=0;
        NR3[a][1]=0;
        NR3[a][2]=0;

        //Colocamos al punto b sobre el eje x.
        NR3[b][0]=sqrt(mat_dist[a][b]);
        NR3[b][1]=0;
        NR3[b][2]=0;

        //Colocamos al punto c en el plano x-y.
        NR3[c][0]=prod_esc[0]/sqrt(mat_dist[b][a]);
        NR3[c][1]=(sqrt(mat_dist[c][a]))*(sqrt(1.0- pow(prod_esc[0]/(sqrt(mat_dist[a][c])*sqrt(mat_dist[a][b])),2) ));
        NR3[c][2]=0;

        //Encontramos el vertice faltante del tetraedro.
        //Nota, en la coordenada y se pueden generar errores si el angulo es de 0°.
        NR3[d][0]=(mat_dist[a][b] + mat_dist[a][d] - mat_dist[b][d])/(2*sqrt(mat_dist[a][b]));

        // Cambiar a S_abc????
        //double Area_abc = sqrt( (mat_dist[a][b] * mat_dist[a][c]) - pow(prod_esc[0],2))/2;

        NR3[d][1]=((mat_dist[a][b]*(mat_dist[a][c] + mat_dist[a][d] - mat_dist[c][d])) - (prod_esc[0]*(mat_dist[a][b] + mat_dist[a][d] - mat_dist[b][d])))/(4 * sqrt(mat_dist[a][b]) * (2*sqrt(S_abc)));
        //NR3[d][2]= sqrt( mat_dist[a][d] - pow(NR3[d][0],2) - pow(NR3[d][1],2) );

        double denominador = mat_dist[b][c] * mat_dist[a][d] * (mat_dist[a][c] + mat_dist[a][b] + mat_dist[c][d] + mat_dist[b][d] - mat_dist[b][c] - mat_dist[a][d]) + mat_dist[a][c] * mat_dist[b][d] * (mat_dist[b][c] + mat_dist[a][b] + mat_dist[c][d] + mat_dist[a][d] - mat_dist[a][c] - mat_dist[b][d]) + mat_dist[a][b] * mat_dist[c][d] * (mat_dist[b][c] + mat_dist[a][c] + mat_dist[a][d] + mat_dist[b][d] - mat_dist[a][b] - mat_dist[c][d]) - mat_dist[b][c] * mat_dist[a][c] * mat_dist[a][b] - mat_dist[b][c] * mat_dist[c][d] * mat_dist[b][d] - mat_dist[a][c] * mat_dist[c][d] * mat_dist[a][d] - mat_dist[a][b] * mat_dist[a][d] * mat_dist[b][d];

        if(denominador <= 1e-15){
            NR3[d][2] = 0;
        }else{
            //NR3[d][2] = sqrt(verificar)/(6*Area_abc);
            NR3[d][2] = sqrt(denominador)/(3*sqrt(S_abc));
        }


        /////////////////////
        //Revisar caso 2.////
        //Triangulo abc  ////
        /////////////////////

        //Bandera para saber si se calculo el triangulo abc.
        bool triangulo_abc=false;
        if(prod_esc[1]>0){
            //Se señala que si se construyo el triangulo abc.
            triangulo_abc =true;

            //Se llama al algoritmo MB2_Simple para obtener el radio, el centro y los lambdas.
            radio = MB2_Simple(NR3[a], NR3[b], NR3[c], mat_dist[a][b], mat_dist[a][c], mat_dist[b][c], centro_R3[0], p, lambda, prod_esc[0]);

            //Se calcula el vector del centro de CB[a,b,c]->d.
            Diferencia(centro_R3[0],NR3[d],difer_R3,3);

            //Se obtiene la norma cuadrada del vector.
            prod_esc[4] = inner_product(difer_R3,difer_R3+3,difer_R3,double(0));

            //Se compara la norma con el radio.
            if(prod_esc[4]<=radio){
                //En caso de que el punto este en CB[a,b,c], se devuelve el punto a R^n usando los lambdas.
                for(int i=0 ; i<col ; ++i)
                    centro[i] = ( lambda[0] * N[b][i] ) + ( lambda[1] * N[c][i] ) + ( lambda[2] * N[a][i] );
                //cout << "\n [a,b,c] es soporte" << endl;
                radio = sqrt(radio);

                //cout << "Radio: " << radio << endl;

                //cout << "Centro: ( ";
                //for(int i=0 ; i<col ; ++i)
                //    cout << centro[i] << " ,";
                //cout << "\b)" << endl;

                // Output
               
                tiempo = ( clock() - t )/( double)CLOCKS_PER_SEC;
                tiempos_1[pruebas] = tiempo;
                promedio[0]+=tiempo;
                continue;

            }
        }

        /////////////////////
        //Revisar caso 2.////
        //Triangulo abd  ////
        /////////////////////
        if(prod_esc[2]>0){
            //caso 2 abd.

            //Se llama al algoritmo MB2_Simple para obtener el radio, el centro y los lambdas.
            radio = MB2_Simple(NR3[a], NR3[b], NR3[d], mat_dist[a][b], mat_dist[a][d], mat_dist[b][d], centro_R3[1], p, lambda, prod_esc[3]);

            //Se calcula el vector del centro de CB[a,b,d]->c.
            Diferencia(centro_R3[1],NR3[c],difer_R3,3);

            //Se obtiene la norma cuadrada del vector.
            prod_esc[4] = inner_product(difer_R3,difer_R3+3,difer_R3,double(0));

            //Se compara la norma con el radio.
            if(prod_esc[4]<=radio){
                //En caso de que el punto este en CB[a,b,d], se devuelve el punto a R^n usando los lambdas.
                for(int i=0 ; i<col ; ++i)
                    centro[i] = ( lambda[0] * N[b][i] ) + ( lambda[1] * N[d][i] ) + ( lambda[2] * N[a][i] );
                //cout << "\n [a,b,d] es soporte" << endl;
                radio = sqrt(radio);

                //cout << "Radio: " << radio << endl;

                //cout << "Centro: ( ";
                //for(int i=0 ; i<col ; ++i)
                //    cout << centro[i] << " ,";
                //cout << "\b)" << endl;
                
                tiempo = ( clock() - t )/( double)CLOCKS_PER_SEC;
                tiempos_1[pruebas] = tiempo;
                promedio[0]+=tiempo;
                

                continue;
            }
        }

        /////////////////////
        //Revisar caso 2.////
        //Triangulo acd  ////
        /////////////////////
        //Se encuentra al par mas alejado entre estos 3 puntos.
        int a_temp, b_temp, c_temp, d_temp;
        orden = 3;
        distancia = Encontrar_Max(mat_dist[a][c],
                            Encontrar_Max(mat_dist[a][d],mat_dist[c][d],2,orden),1,orden);
        d_temp = b;
        if(orden==1){
            a_temp = a;
            b_temp = c;
            c_temp = d;
        }
        if(orden==2){
            a_temp = a;
            b_temp = d;
            c_temp = c;
        }
        if(orden==3){
            a_temp = c;
            b_temp = d;
            c_temp = a;
        }

        //Se encuientran los vectores entre uno de los 2 puntos mas alejados y los otros 2.
        //Se obtiene el producto interior entre los vectores.
        prod_esc[3] = 0.5 * ( mat_dist[a_temp][b_temp] + mat_dist[a_temp][c_temp] - mat_dist[b_temp][c_temp] );
        // AGREGAR EL "IF"
        if(prod_esc[3]>0){
            //Se llama al algoritmo MB2_Simple para obtener el radio, el centro y los lambdas.
            radio = MB2_Simple(NR3[a_temp], NR3[b_temp], NR3[c_temp], mat_dist[a_temp][b_temp], mat_dist[a_temp][c_temp], mat_dist[b_temp][c_temp], centro_R3[2], p, lambda, prod_esc[3]);

            //Se calcula el vector del centro de CB[a,c,d]->b.
            Diferencia(centro_R3[2],NR3[d_temp],difer_R3,3);

            //Se obtiene la norma cuadrada del vector.
            prod_esc[4] = inner_product(difer_R3,difer_R3+3,difer_R3,double(0));

            //Se compara la norma con el radio.
            if(prod_esc[4]<=radio){

                //En caso de que el punto este en CB[a,c,d], se devuelve el punto a R^n usando los lambdas.
                for(int i=0 ; i<col ; ++i)
                    centro[i] = ( lambda[0] * N[b_temp][i] ) + ( lambda[1] * N[c_temp][i] ) + ( lambda[2] * N[a_temp][i] );
                //cout << "\n [a,c,d] es soporte" << endl;
                radio = sqrt(radio);

                //cout << "Radio: " << radio << endl;

                //cout << "Centro: ( ";
                //for(int i=0 ; i<col ; ++i)
                //    cout << centro[i] << " ,";
                //cout << "\b)" << endl;
                tiempo = ( clock() - t )/( double)CLOCKS_PER_SEC;
                tiempos_1[pruebas] = tiempo;
                promedio[0]+=tiempo;

                continue;
            }
        }
        /*
        //Dependiendo de cual fue el par mas alejado, se procede.
        switch(orden){
            case 1:
                //Se encuientran los vectores entre uno de los 2 puntos mas alejados y los otros 2.
                //Se obtiene el producto interior entre los vectores.
                prod_esc[3] = 0.5 * (mat_dist[a][c] + mat_dist[a][d] - mat_dist[a][d]);
                // AGREGAR EL "IF"
                if(prod_esc[3]>0){
                    //Se llama al algoritmo MB2_Simple para obtener el radio, el centro y los lambdas.
                    radio = MB2_Simple(NR3[a], NR3[c], NR3[d], mat_dist[a][c], mat_dist[a][d], mat_dist[c][d], centro_R3[2], p, lambda, prod_esc[3]);

                    //Se calcula el vector del centro de CB[a,c,d]->b.
                    Diferencia(centro_R3[2],NR3[b],difer_R3,3);

                    //Se obtiene la norma cuadrada del vector.
                    prod_esc[4] = inner_product(difer_R3,difer_R3+3,difer_R3,double(0));

                    //Se compara la norma con el radio.
                    if(prod_esc[4]<=radio){

                        //En caso de que el punto este en CB[a,c,d], se devuelve el punto a R^n usando los lambdas.
                        for(int i=0 ; i<col ; ++i)
                            centro[i] = ( lambda[0] * N[c][i] ) + ( lambda[1] * N[d][i] ) + ( lambda[2] * N[a][i] );
                        cout << "\n [a,c,d] es soporte" << endl;
                        radio = sqrt(radio);

                        cout << "Radio: " << radio << endl;

                        cout << "Centro: ( ";
                        for(int i=0 ; i<col ; ++i)
                            cout << centro[i] << " ,";
                        cout << "\b)" << endl;

                        // Output
                        FT rad = mb.radius();
                        cout << "Radius = " << rad << endl
                            << "Center:" << endl;
                        Miniball::Coordinate_iterator center_it = mb.center_begin();
                        for (int j=0; j<col; ++j)
                            cout << "  " << center_it[j] << endl;

                        for (int i = 0; i < ren; ++i)
                            delete [] N[i];
                        delete [] N;

                        delete[] centro;
                        return 0;
                    }
                }
                break;
            case 2:
                //Analogamente a la posibilidad anterior.

                prod_esc[3] = 0.5 * (mat_dist[a][c] + mat_dist[a][d] - mat_dist[a][d]);

                if(prod_esc[3]>0){
                    radio = MB2_Simple(NR3[a], NR3[d], NR3[c], mat_dist[a][d], mat_dist[a][c], mat_dist[c][d], centro_R3[2], p, lambda, prod_esc[3]);

                    //Se calcula el vector del centro de CB[a,c,d]->b.
                    Diferencia(centro_R3[2],NR3[b],difer_R3,3);

                    //Se obtiene la norma cuadrada del vector.
                    prod_esc[4] = inner_product(difer_R3,difer_R3+3,difer_R3,double(0));

                    //Se compara la norma con el radio.
                    if(prod_esc[4]<=radio){

                        //En caso de que el punto este en CB[a,c,d], se devuelve el punto a R^n usando los lambdas.
                        for(int i=0 ; i<col ; ++i)
                            centro[i] = ( lambda[0] * N[d][i] ) + ( lambda[1] * N[c][i] ) + ( lambda[2] * N[a][i] );
                        cout << "\n [a,d,c] es soporte" << endl;
                        radio = sqrt(radio);

                        cout << "Radio: " << radio << endl;

                        cout << "Centro: ( ";
                        for(int i=0 ; i<col ; ++i)
                            cout << centro[i] << " ,";
                        cout << "\b)" << endl;

                        // Output
                        FT rad = mb.radius();
                        cout << "Radius = " << rad << endl
                            << "Center:" << endl;
                        Miniball::Coordinate_iterator center_it = mb.center_begin();
                        for (int j=0; j<col; ++j)
                            cout << "  " << center_it[j] << endl;

                        for (int i = 0; i < ren; ++i)
                            delete [] N[i];
                        delete [] N;

                        delete[] centro;
                        return 0;

                    }
                }

                break;
            case 3:
                //Analogamente a la primer posibilidad.

                prod_esc[3] = 0.5 * (mat_dist[a][c] + mat_dist[d][c] - mat_dist[a][d]);

                if(prod_esc[3]>0){

                    radio = MB2_Simple(NR3[c], NR3[d], NR3[a], mat_dist[c][d], mat_dist[a][c], mat_dist[a][d], centro_R3[2], p, lambda, prod_esc[3]);


                    //Se calcula el vector del centro de CB[a,c,d]->b.
                    Diferencia(centro_R3[2],NR3[b],difer_R3,3);

                    //Se obtiene la norma cuadrada del vector.
                    prod_esc[4] = inner_product(difer_R3,difer_R3+3,difer_R3,double(0));

                    //Se compara la norma con el radio.
                    if(prod_esc[4]<=radio){

                        //En caso de que el punto este en CB[a,c,d], se devuelve el punto a R^n usando los lambdas.
                        for(int i=0 ; i<col ; ++i)
                            centro[i] = ( lambda[0] * N[d][i] ) + ( lambda[1] * N[a][i] ) + ( lambda[2] * N[c][i] );
                        cout << "\n [c,d,a] es soporte" << endl;
                        radio = sqrt(radio);

                        cout << "Radio: " << radio << endl;

                        cout << "Centro: ( ";
                        for(int i=0 ; i<col ; ++i)
                            cout << centro[i] << " ,";
                        cout << "\b)" << endl;

                        // Output
                        FT rad = mb.radius();
                        cout << "Radius = " << rad << endl
                            << "Center:" << endl;
                        Miniball::Coordinate_iterator center_it = mb.center_begin();
                        for (int j=0; j<col; ++j)
                            cout << "  " << center_it[j] << endl;

                        for (int i = 0; i < ren; ++i)
                            delete [] N[i];
                        delete [] N;

                        delete[] centro;
                        return 0;

                    }
                }

                break;
        };
        */

        /////////////////////
        //Revisar caso 2.////
        //Triangulo bcd  ////
        /////////////////////

        //El procedimiento es analogo que con el triangulo acd.
        orden = 3;
        distancia = Encontrar_Max(mat_dist[b][c],
                            Encontrar_Max(mat_dist[b][d],mat_dist[c][d],2,orden),1,orden);
        d_temp = a;
        if(orden==1){
            a_temp = b;
            b_temp = c;
            c_temp = d;
        }
        if(orden==2){
            a_temp = b;
            b_temp = d;
            c_temp = c;
        }
        if(orden==3){
            a_temp = c;
            b_temp = d;
            c_temp = b;
        }

        //Se encuientran los vectores entre uno de los 2 puntos mas alejados y los otros 2.
        //Se obtiene el producto interior entre los vectores.
        prod_esc[3] = 0.5 * (mat_dist[a_temp][b_temp] + mat_dist[a_temp][c_temp] - mat_dist[b_temp][c_temp]);
        // AGREGAR EL "IF"
        if(prod_esc[3]>0){
            //Se llama al algoritmo MB2_Simple para obtener el radio, el centro y los lambdas.
            radio = MB2_Simple(NR3[a_temp], NR3[b_temp], NR3[c_temp], mat_dist[a_temp][b_temp], mat_dist[a_temp][c_temp], mat_dist[b_temp][c_temp], centro_R3[2], p, lambda, prod_esc[3]);

            //Se calcula el vector del centro de CB[a,c,d]->b.
            Diferencia(centro_R3[2],NR3[d_temp],difer_R3,3);

            //Se obtiene la norma cuadrada del vector.
            prod_esc[4] = inner_product(difer_R3,difer_R3+3,difer_R3,double(0));

            //Se compara la norma con el radio.
            if(prod_esc[4]<=radio){

                //En caso de que el punto este en CB[a,c,d], se devuelve el punto a R^n usando los lambdas.
                for(int i=0 ; i<col ; ++i)
                    centro[i] = ( lambda[0] * N[b_temp][i] ) + ( lambda[1] * N[c_temp][i] ) + ( lambda[2] * N[a_temp][i] );
                //cout << "\n [a,c,d] es soporte" << endl;
                radio = sqrt(radio);

                //cout << "Radio: " << radio << endl;

                //cout << "Centro: ( ";
                //for(int i=0 ; i<col ; ++i)
                //    cout << centro[i] << " ,";
                //cout << "\b)" << endl;
                tiempo = ( clock() - t )/( double)CLOCKS_PER_SEC;
                tiempos_1[pruebas] = tiempo;
                promedio[0]+=tiempo;

                continue;
            }
        }
        /*
        switch(orden){
            case 1:

                prod_esc[3] = 0.5 * (mat_dist[c][b] + mat_dist[d][b] - mat_dist[c][d]);
                if(prod_esc[3]>0){

                    radio = MB2_Simple(NR3[b], NR3[c], NR3[d], mat_dist[b][c], mat_dist[b][d], mat_dist[c][d], centro_R3[3], p, lambda, prod_esc[3]);

                    Diferencia(centro_R3[3],NR3[a],difer_R3,3);

                    prod_esc[4] = inner_product(difer_R3,difer_R3+3,difer_R3,double(0));

                    if(prod_esc[4]<=radio){

                        for(int i=0 ; i<col ; ++i)
                            centro[i] = ( lambda[0] * N[c][i] ) + ( lambda[1] * N[d][i] ) + ( lambda[2] * N[b][i] );
                        cout << "\n [b,c,d] es soporte" << endl;
                        radio = sqrt(radio);

                        cout << "Radio: " << radio << endl;

                        cout << "Centro: ( ";
                        for(int i=0 ; i<col ; ++i)
                            cout << centro[i] << " ,";
                        cout << "\b)" << endl;

                        // Output
                        FT rad = mb.radius();
                        cout << "Radius = " << rad << endl
                            << "Center:" << endl;
                        Miniball::Coordinate_iterator center_it = mb.center_begin();
                        for (int j=0; j<col; ++j)
                            cout << "  " << center_it[j] << endl;

                        for (int i = 0; i < ren; ++i)
                            delete [] N[i];
                        delete [] N;

                        delete[] centro;
                        return 0;

                    }
                }

                break;
            case 2:

                prod_esc[3] = 0.5 * (mat_dist[c][b] + mat_dist[d][b] - mat_dist[c][d]);
                if(prod_esc[3]>0){

                    radio = MB2_Simple(NR3[b], NR3[d], NR3[c], mat_dist[b][d], mat_dist[b][c], mat_dist[c][d], centro_R3[3], p, lambda, prod_esc[3]);

                    Diferencia(centro_R3[3],NR3[a],difer_R3,3);

                    prod_esc[4] = inner_product(difer_R3,difer_R3+3,difer_R3,double(0));

                    if(prod_esc[4]<=radio){

                        for(int i=0 ; i<col ; ++i)
                            centro[i] = ( lambda[0] * N[d][i] ) + ( lambda[1] * N[c][i] ) + ( lambda[2] * N[b][i] );
                        cout << "\n [b,d,c] es soporte" << endl;
                        radio = sqrt(radio);

                        cout << "Radio: " << radio << endl;

                        cout << "Centro: ( ";
                        for(int i=0 ; i<col ; ++i)
                            cout << centro[i] << " ,";
                        cout << "\b)" << endl;

                        // Output
                        FT rad = mb.radius();
                        cout << "Radius = " << rad << endl
                            << "Center:" << endl;
                        Miniball::Coordinate_iterator center_it = mb.center_begin();
                        for (int j=0; j<col; ++j)
                            cout << "  " << center_it[j] << endl;

                        for (int i = 0; i < ren; ++i)
                            delete [] N[i];
                        delete [] N;

                        delete[] centro;
                        return 0;

                    }
                }

                break;
            case 3:

                prod_esc[3] = 0.5 * (mat_dist[b][c] + mat_dist[d][c] - mat_dist[b][d]);

                if(prod_esc[3]>0){

                    radio = MB2_Simple(NR3[c], NR3[d], NR3[b], mat_dist[c][d], mat_dist[b][c], mat_dist[b][d], centro_R3[3], p, lambda, prod_esc[3]);

                    Diferencia(centro_R3[3],NR3[a],difer_R3,3);

                    prod_esc[4] = inner_product(difer_R3,difer_R3+3,difer_R3,double(0));

                    if(prod_esc[4]<=radio){

                        for(int i=0 ; i<col ; ++i)
                            centro[i] = ( lambda[0] * N[d][i] ) + ( lambda[1] * N[b][i] ) + ( lambda[2] * N[c][i] );
                        cout << "\n [c,d,b] es soporte" << endl;
                        radio = sqrt(radio);

                        cout << "Radio: " << radio << endl;

                        cout << "Centro: ( ";
                        for(int i=0 ; i<col ; ++i)
                            cout << centro[i] << " ,";
                        cout << "\b)" << endl;

                        // Output
                        FT rad = mb.radius();
                        cout << "Radius = " << rad << endl
                            << "Center:" << endl;
                        Miniball::Coordinate_iterator center_it = mb.center_begin();
                        for (int j=0; j<col; ++j)
                            cout << "  " << center_it[j] << endl;

                        for (int i = 0; i < ren; ++i)
                            delete [] N[i];
                        delete [] N;

                        delete[] centro;
                        return 0;
                    }
                }
                break;
        };
        */
        //////////////////////
        /// Revisar caso 3 ///
        ///   Tetraedro    ///
        //////////////////////



        //Se verifica si se obtuvo el centro del triangulo abc, en caso de que no, se usa el algoritmo MB2_Simple para calcularlo.
        if(!triangulo_abc){
            radio = MB2_Simple(NR3[a], NR3[b], NR3[c], mat_dist[a][b], mat_dist[a][c], mat_dist[b][c], centro_R3[0], p, lambda, prod_esc[0]);
        }

        //Se obtienen las cordenadas del centro donde solo cambia z en comparacion con el centro de abc.
        centro_R3[4][0] = centro_R3[0][0];
        centro_R3[4][1] = centro_R3[0][1];
        centro_R3[4][2] = (mat_dist[a][d] - 2*( NR3[d][0]*centro_R3[4][0] - NR3[d][1]*centro_R3[4][1]))/(2*NR3[d][2]);

        //Se obtiene el radio cuadrado.
        radio = sqrt(inner_product(centro_R3[4],centro_R3[4]+3,centro_R3[4],double(0)));

        //Se calculan los lambdas con sustitucion hacia atras.
        lambda[2] = centro_R3[4][2]/NR3[d][2];
        lambda[1] = (centro_R3[4][1] - lambda[2]*NR3[d][1])/NR3[c][1];
        lambda[0] = ( centro_R3[4][0] - lambda[1]*NR3[c][1] - lambda[2]*NR3[d][0] )/NR3[b][0];

        //Se devuelve el centro a R^n.
        for(int i=0 ; i<col ; ++i)
            centro[i] = ( 1 - ( lambda[0] + lambda[1] + lambda[2] ) )*N[a][i] + ( lambda[0] * N[b][i] ) + ( lambda[1] * N[c][i] ) + ( lambda[2] * N[d][i] );

        //cout << "\n El [a,b,c,d] es el soporte" << endl;

        radio = sqrt(radio);

        //cout << "Radio: " << radio << endl;

        //cout << "Centro: ( ";
        //for(int i=0 ; i<col ; ++i)
        //    cout << centro[i] << " ,";
        //cout << "\b)" << endl;
        
        tiempo = ( clock() - t )/( double)CLOCKS_PER_SEC;
        tiempos_1[pruebas] = tiempo;
        promedio[0]+=tiempo;
    }
    
    delete [] centro;
	
    for (int i = 0; i < ren; ++i)
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
						 << scientific <<  desviacion[0] << ',' 
		 			     << scientific <<  promedio[1] << ','
				 	     << scientific <<  desviacion[1] << endl;
   
    salida.close();					
   	
   	
   	delete[] desviacion;
   	delete[] promedio;

	return 0;
}
//Funcion que encuentra el valor maximo entre 2 numeros y regresa tambien un indicador de cual valor se regreso.
double Encontrar_Max(double valor_1, double valor_2, int a, int &b){
	if(valor_1<=valor_2){
		return valor_2;
	}else{
		b=a;
		return valor_1;
	}
}
//Calcula la distancia cuadrada entre 2 puntos.
double Distancia(double a[], double b[], int dim){
	double dist=0;

	for( int i = 0 ; i < dim ; ++i)
		dist += (a[i]-b[i])*(a[i]-b[i]);

	return dist;
}

//Calcula la diferencia entre 2 puntos (el vector entre ellos).
void Diferencia(double a[], double b[], double dif[], int dim){

	for( int i = 0 ; i < dim ; ++i)
		dif[i] = a[i]-b[i];

}

// El algoritmo MB2_Simple recibe tres puntos (a,b,c) y regresa: centro y radio de la bola circunscrita en R^d, y coordenadas afines del centro.
double MB2_Simple(double a[], double b[], double c[], double dist_ab, double dist_ac, double dist_bc, double centro[], double p[], double lambda[], double prod_esc){
	double radio;

	p[0] = sqrt(dist_ab)/(double)2;
	p[1] = sqrt(dist_ab)*(dist_ac - prod_esc)/(2 * sqrt((dist_ab*dist_ac)- pow(prod_esc,2)));

	radio = (p[0]*p[0])+(p[1]*p[1]);

    double deter = 2 * ( ( dist_ab * dist_ac ) - ( prod_esc*prod_esc ) );

    lambda[0] = ( dist_ac * ( dist_ab - prod_esc ) ) / deter;
    lambda[1] = ( dist_ab * ( dist_ac - prod_esc ) ) / deter;
    lambda[2] = ( 1.0 - lambda[0] - lambda[1] );

    for(int i=0 ; i<3 ; ++i)
        centro[i] = ( lambda[0] * b[i] ) + ( lambda[1] * c[i] ) + ( lambda[2] * a[i] );

	return radio;
}

/*
void Cap_Dim(int &dim){
   
    ifstream entrada;
   
    entrada.open("value.txt");
    
    
    if(!entrada){
        cout << "Error: no se pudo abrir el archivo..." << endl;
        return;
    }
    
   
    entrada >> dim;

    entrada.close();
    
    dim = dim*50;
}

*/
