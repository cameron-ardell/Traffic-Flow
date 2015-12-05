#include "nr3.h"
#include "odeint.h"
#include "stepper.h"
#include "stepperdopr5.h"

using namespace std;


//64 cells wide (for distance)
//
struct RoadRage {

	void operator() (const Doub x, VecDoub_I &y, VecDoub_O &dydx){
		const Doub rho = y[0];
		const Doub kappa = y[1];

		//assuming dydx implies levels of derivatives
		//how will this affect integration?
		dydx[0] = kappa;
		//dydx[1] = -kappa;
	}

};

struct Traffic{
	const static Doub x_center = 0.5;
	const static Doub lambda = 0.1;
	const static Doub u_max = 1.0;
	const static Doub p_max = 1.0;
	const static Doub x_max = 1.0;



	Doub perturb;
	Doub rho_bar;

	int dim;

	MatDoub rho_matrix;
	VecDoub road;
	Doub size_road = 12;
	Doub t = 6;

	Doub delta = 1.0/ (Doub)size_road;



	Traffic(Doub perturb_in, Doub rho_bar_in) {
		perturb = perturb_in * p_max;
		rho_bar = rho_bar_in * p_max;
	}

	void init(){
		rho_matrix.assign(t,size_road, 1.0);

		//initalize values at t = 0
		int index = 0;
		for(Doub x = 0.0; x < x_max; x += delta){
			
			Doub cur_rho = rho_bar + perturb
				* exp(-( x - x_center)* ( x - x_center) / (lambda * lambda));

			road[index] = cur_rho;
			rho_matrix[index][0] = cur_rho;

			index++;
		}
	}

	~Traffic(){}

	void print_matrix(){

		for(int i = 0; i < t; i++){
			

			for(int j = 0; j < size_road; j++){

				cout << setw(12) << setprecision(3) << rho_matrix[i][j];

			}
			cout << "\n";
		}
	}

	void print_road(){
		for(int i = 0; i < size_road; i++){
			cout << setw(5) << setprecision(3) << road[i];
		}
	}

};



int main(){

	Traffic traf(10e-3, 0.5);

	//traf.print_matrix();
	//traf.print_road();

	traf.init();


	traf.print_matrix();
	//traf.print_road();

	return 0;
}