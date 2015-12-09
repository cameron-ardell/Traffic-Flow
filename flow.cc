#include "nr3.h"
#include "odeint.h"
#include "stepper.h"
#include "stepperdopr5.h"

using namespace std;


struct RoadRage {


    Doub u_max,p_max,eta;
    int size_road;

    //constructor (get and set rho vector)
    RoadRage(Doub u_max_in, Doub p_max_in,int size_road_in,Doub eta_in){
        u_max = u_max_in;
        p_max = p_max_in;
        size_road = size_road_in;
        eta = eta_in;
  };

  ~RoadRage() {};

  //function to integrate over
  void operator() (const Doub x, VecDoub_I &rho, VecDoub_O &dpdt){

    //Use finite differencing to find df(p)/dx
    Doub f_right,f_left,f_cent,visc;
    Doub delta = 1.0 / (Doub)size_road;

    for(int i = 0; i <= size_road; i++){

      //Boundaries (Use ghost points at boundaries)
      //Left boundary
        if(i==0) {
            f_left = rho[size_road] * u_max * (1.0 - rho[size_road]/ p_max);
            f_right = rho[i+1] * u_max * (1.0 - rho[i+1])/ p_max;
        } 
        //Right boundary
        else if (i==size_road) {
            f_right = rho[0] * u_max * (1.0 - rho[0]/ p_max);
            f_left = rho[i-1] * u_max * (1.0 - rho[i-1]/ p_max);
        }
        //Not on boundaries
        else{
            f_right = rho[i+1] * u_max * (1.0 - rho[i+1]/ p_max);
            f_left = rho[i-1] * u_max * (1.0 - rho[i-1]/ p_max);
        }

      //Calculate f at i
      f_cent = rho[i]*u_max*(1.0 - rho[i]/p_max);  

      //Find artificial viscosity using 2nd derivative
      visc = - eta * (f_right - 2*f_cent + f_left) / delta;

      //Find dpdt
      dpdt[i] = - (f_right - f_left) / (2*delta) + visc;

    }
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

    VecDoub road;
    int size_road;
    int t;
    Doub delta;
    Doub eta;


    //constructor
    Traffic(Doub perturb_in, Doub rho_bar_in, int t_in, int size_road_in, Doub eta_in) {
        //Multiply by rho_max to get real values
        perturb = perturb_in * p_max;
        rho_bar = rho_bar_in * p_max;

        //Initialize t and x and eta
        size_road = size_road_in;
        t = t_in;
        eta = eta_in;

        delta = 1.0/ (Doub)size_road;
    }

  void init(){
    //Assign 
    road.assign(size_road + 1, 1.0);
  }

  void update(){

    //sets absolute and relative tolerances
    //also sets minimum step size and first guess of step size
    const Doub atol = 1.e-5, rtol = atol, h1=0.001, hmin=0.0;
    Output out;
    //width of cells in data files
    int width = 17;


    //Initialize data file of rho values over t and x for surface plot
    ofstream outfile_p;
    outfile_p.open("Rho over x and t");
    outfile_p.setf(ios::left);
    outfile_p << setw(width) << "# t " << setw(width) << " x " << setw(width) << " rho " << endl;
    outfile_p << "#====================================================" << endl;


    //Write rho values to file at time t = 0
    ofstream outfile_t;
    ostringstream filename;
    filename << "Rho over x at t=" << 0 << ".out" << ends;
    outfile_t.open(filename.str().c_str());  
    outfile_t.setf(ios::left);
    outfile_t << setw(width) << "#Rho"<< setw(width) << "x   " << endl;


    //initialize file for making heat map
    ofstream outfile_h;
    outfile_h.open("heatmap.dat");
    outfile_h.setf(ios::left);


    //Set initial rho values as guess and write to file
    for(int i = 0; i <= size_road; i++){
        Doub x = (Doub)i*delta;

        //inital guess equation
        Doub cur_rho = rho_bar + perturb * exp(-(x-x_center)*(x-x_center) / (lambda*lambda));

        //Write to road vector (which is used for ODE)
        road[i] = cur_rho;


        //Write rho and x to file
        outfile_t << setw(width) << x << setw(width) << cur_rho << endl;
        //Write rho and x and t to file
        outfile_p << setw(width) << 0 << setw(width) << x << setw(width) << cur_rho << endl;
        //writing to heatmap file
        outfile_h << setw(width) << setprecision(5) << cur_rho;
    }
    outfile_h << endl;

     

    //Initialize odeint object
    RoadRage anger(u_max, p_max, size_road,eta);


    //Integrate at each value of x over time (not including initial guess)
    for(int i = 0; i < t; i++){
        Odeint<StepperDopr5<RoadRage> > ode(road, i, i+1, atol, rtol, h1, hmin,out,anger);
        ode.integrate();

        //Write rho values to file if t is in quarters of tmax
        if(i == t/4 || i==t/2 || i==3*t/4 || i==t-1){

            //Create file with rho values over time
            ofstream outfile_t;
            ostringstream filename;
            filename << "Rho over x at t=" << i << ".out" << ends;
            outfile_t.open(filename.str().c_str());

            outfile_t.setf(ios::left);
            outfile_t << setw(width) << "#Rho" << setw(width) << "x" << endl;


            //Write rho vs x to file
            for(int s = 0; s <= size_road; s++){
                Doub x = (Doub)s*delta;
                outfile_t << setw(width) << x << setw(width) << road[s] << endl;
            }
        }

        //write current road data to splot and heat map
        for(int r = 0; r <= size_road; r++){
            Doub x = (Doub)r*delta;

            //to heatmap
            outfile_h << setw(width) << setprecision(5) << road[r];

            //only write to splot if not first guess (already did that)
            if(i!=0){
                outfile_p << setw(width) << i << setw(width) << x << setw(width) << road[r] << endl;
            }
        }

        outfile_h << endl;
        outfile_p << endl;
    }
  }

  ~Traffic(){}


  void print_road(){
    for(int i = 0; i <= size_road; i++){
        cout << setw(5) << setprecision(3) << road[i];
    }
    cout << endl;
  }

};


int main(){

    Doub delta_rho = .5;
    Doub rho_bar = 0.5;
    int t = 25;
    int size_road = 100;
    Doub eta = 1.0;

    Traffic traf(delta_rho, rho_bar, t, size_road,eta);


    traf.init();
    traf.update();


    return 0;
}
