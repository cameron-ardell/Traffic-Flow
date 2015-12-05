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

	Doub delta;


	Doub perturb;
	Doub rho_bar;


	MatDoub rho_matrix;
	VecDoub road;
	Doub size_road;
	Doub t;


	Traffic(Doub perturb_in, Doub rho_bar_in, Doub size_road_in, Doub t_in) {
		perturb = perturb_in * p_max;
		size_road = size_road_in;
		rho_bar = rho_bar_in * p_max;
		t = t_in;
		delta = 1.0/ (size_road + 1.0);

		//road.resize(size_road);
		rho_matrix.assign(size_road, t, 0);
	}

	void init(){

		//initalize values at t = 0
		int index = 0;
		for(Doub x= 0.0; x < x_max; x+= delta){
			
			Doub cur_rho = rho_bar + perturb
				* exp(-( x - x_center)* ( x - x_center) / (lambda * lambda));

			//road[index] = cur_rho;
			rho_matrix[index][0] = cur_rho;

			index++;
		}
	}

	~Traffic(){}

	void print_matrix(){

		for(int i = 0; i < size_road; i++){
			
			cout << "\n";

			for(int j = 0; j < t; j++){

				cout << setw(5) << rho_matrix[size_road][j];

			}
		}
	}

	void print_road(){
		for(int i = 0; i < size_road; i++){
			cout << setw(5) << road[i];
		}
	}

};



int main(){

	Traffic traf(10e-3, 0.5, 3.0, 4.0);

	traf.print_matrix();
	traf.print_road();

	traf.init();


	traf.print_matrix();
	traf.print_road();

	return 0;
}