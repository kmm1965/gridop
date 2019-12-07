#include "average_operator.hpp"
#include "first_derivative_operator.hpp"
#include "cyclic_index.hpp"
#include "symmetric_index.hpp"

int main1(int argc, char* argv[])
{   
    std::ofstream os_enrg;
    os_enrg << std::setprecision(16);
    os_enrg << std::scientific;
    os_enrg.open("energy.dat", std::ios::out);

    std::cout << std::setprecision(16);
    std::cout << std::scientific;
    //const int TN = 10;
    const int max_iter = 100;
	const size_t N = 100;
    const double L = 1e-2; // [m]
	const double h = L/static_cast<double>(N); // [m]
    const double eps = 0.0;
    const double dt = 0.25e-7; // [s]
    const double CS = 200; // [m/s]
    const double tau = 0.5*h/CS; //[s]
    const double eta = 5e-4; // [Pa*s]
    const double lambda1 = 6e-4;
    //const double M = 1e-10;
    const double M = 1e-9;
    const double A_psi = 1e+4;
    const double beta = std::sqrt(2*A_psi/lambda1);

    const size_t CELL_N = N;
    const size_t EDGE_N = CELL_N; // number of edges
	vector_grid_function<tag_aux, double> w(EDGE_N), jm(EDGE_N), P(EDGE_N);
	vector_grid_function<tag_main, double> C(CELL_N), mu(CELL_N), u(CELL_N), rho(CELL_N, 0, 1.0),
                                           Psi0(CELL_N), E_lmd(CELL_N), Psi1(CELL_N), Qw(CELL_N),
                                           rho_new(CELL_N), rho_C_new(CELL_N), rho_u_new(CELL_N), C_new(CELL_N),
                                           totEnrgGF(CELL_N),
                                           lmdEnrgGF(CELL_N),
                                           kinEnrgGF(CELL_N),
                                           psiEnrgGF(CELL_N);

    //double x = -1;
#ifdef __CUDACC__
    host_vector<double> hC(C.size());
#else
    MATH_VECTOR_BASE_CLASS<double> &hC = C;
#endif
    for (int i = 0; i < CELL_N; ++i){
            //x = h*(i+0.5);
            //C[i] = 0.5 + 0.5*std::tanh(beta*(x-L*0.25)*0.5) + 0.5*std::tanh(-beta*(x-L*0.75)*0.5);//(i < N/2)?eps:(1-eps);
           // C[i] = 0.5 + 0.5*std::tanh(beta*(x-L*0.5)*0.5);
            hC[i] = (i > N/2)?eps:(1.0-eps);
    }
#ifdef __CUDACC__
    C.assign(hC);
#endif

    const cyclic_index index(CELL_N);
	const average_operator<tag_main> avr(index);   // operator maps from cell to edge
	const average_operator<tag_aux> avrSt(index); // operator maps from edge to cell
	const first_derivative_operator<double, tag_main> drv(index, h);
	const first_derivative_operator<double, tag_aux> drvSt(index, h);

    const int WIter = 1;

//    for (int i = 0; i < 1e+7; i++)
    for (int i = 0; i <= max_iter; i++)
    {
//        std::cout << "iter " << i << std::endl;



        Psi1  = (CS*CS)*log(rho);

        E_lmd = (lambda1*0.5)*avrSt(sqr(drv(C)));
        Psi0  = A_psi*sqr(C*(1.0-C)) + Psi1;

        lmdEnrgGF = rho*E_lmd;
        psiEnrgGF = rho*Psi0;
     //   psiEnrgGF = rho*A_psi*func::sqr(C*(1.0-C));
        kinEnrgGF = 0.5*rho*sqr(u);
        totEnrgGF = lmdEnrgGF + psiEnrgGF  + kinEnrgGF;
        const double
            totEnrg = MATH_ACCUMULATE(totEnrgGF.begin(), totEnrgGF.end(), 0.0),
            kinEnrg = MATH_ACCUMULATE(kinEnrgGF.begin(), kinEnrgGF.end(), 0.0),
            lmdEnrg = MATH_ACCUMULATE(lmdEnrgGF.begin(), lmdEnrgGF.end(), 0.0),
            psiEnrg = MATH_ACCUMULATE(psiEnrgGF.begin(), psiEnrgGF.end(), 0.0);

        std::cout << i << " Energy " << totEnrg*h << std::endl;
        os_enrg << i*dt << ' ' <<  totEnrg*h << ' ' 
                << lmdEnrg*h << ' ' << psiEnrg*h << ' ' << kinEnrg*h <<  std::endl;

////////   NEW   ///////
        mu = A_psi*2.0*C*(1.0-C)*(1.0-2.0*C) - lambda1*drvSt(avr(rho)*drv(C))/rho;

        w  = tau*(avr(u)*drv(u) + drv(Psi0) + drv(E_lmd) - avr(mu)*drv(C));
        P  = (4.0/3.0)*eta*drv(u) + avr(rho)*avr(u)*w;
        jm = avr(rho)*(avr(u) - w);

        rho_new = rho - dt*drvSt(jm);
        rho_u_new = rho*u + dt*( 
                                 - drvSt(jm*avr(u)) - avrSt(avr(rho)*drv(Psi0 + E_lmd)) 
                                 + drvSt(P) + mu*avrSt(avr(rho)*drv(C))
                               ); 
        rho_C_new = rho*C + dt*(
                                 - drvSt(jm*avr(C) - 0.25*h*h*avr(rho)*drv(C)*drv(u))
                                 + drvSt(M*drv(mu))
                               );

 // NEW 


 // OLD
 
 /*       mu = A_psi*2.0*C*(1.0-C)*(1.0-2.0*C) - lambda1*drvSt(avr(rho)*drv(C))/rho;
        Qw = lambda1*avrSt(avr(rho)*drv(C))*drvSt(avr(C));
        w  = tau*(avr(u)*drv(u) + drv(Psi1) + drv(Qw)/avr(rho));
        P  = -CS*CS*avr(rho) + (4.0/3.0)*eta*drv(u) + avr(rho)*avr(u)*w - avr(Qw);
        jm = avr(rho)*(avr(u) - w);

        rho_new = rho - dt*drvSt(jm);
        rho_u_new = rho*u + dt*( 
                                 - drvSt(jm*avr(u))
                                 + drvSt(P)
                               ); 
        rho_C_new = rho*C + dt*(
                                 - drvSt(jm*avr(C))
                                 + drvSt(M*drv(mu))
                               );*/

//        C_new = C + dt*(
//                                 - avrSt((avr(u) - w)*drv(C))
//                                 + drvSt(M*drv(mu))/rho
//                               );
 // OLD END

        rho = rho_new;
        u = rho_u_new/rho_new;
        C = rho_C_new/rho_new;

        if (rho[5] != rho[5])
            abort();

        if (i%WIter == 0)
        {
            std::ofstream os;
            os << std::setprecision(16);
            os << std::scientific;
            char buf[256];
            sprintf(buf, "results_%d.dat", i/WIter);
            os.open(buf, std::ios::out);

            for (int i = 0; i < N; ++i)
            {
                os << h*(i+0.5) << ' ' <<  C[i] << ' ' << u[i] << ' ' << rho[i]  << std::endl;
            }
            os.close();
        }

        
    }

    //	gfc = deriv1(gfv1);

    /*	cyclic_grid_function<tag_main, double> ud1(N), ud2(N);
        ud1 = u + 1.25;
        ud2 = drvSt(avr(u));

	for (int i = 0; i < N; ++i)
		std::cout << ud1[i] << ' ';

    std::cout << std::endl;

	for (int i = 0; i < N; ++i)
		std::cout << ud2[i] << ' ';

    std::cout << std::endl;*/
/*    std::cout << ' ';
	for (int i = 0; i < N; ++i)
		std::cout << Ca[i] << ' ';

    std::cout << std::endl;

    std::cout << ' ';
	for (int i = 0; i < N; ++i)
		std::cout << Ca2[i] << ' ';

    std::cout << std::endl;*/
/*	const forward_average_operator<tag_aux> av_v;
	gfv2 = (av_v & deriv1)(gfv1);

	//gfc = (av_v & deriv1)(gfv1);

	gfc = (av_f & av_v & deriv1)(gfv1);

	gfv2 = (av_v & (av_f + deriv1))(gfv1);
	gfv2 = (av_v & (av_f - deriv1))(gfv1);
	gfv2 = (av_v & (av_f * deriv1))(gfv1);
	gfv2 = (av_v & (av_f / deriv1))(gfv1);

	//gfv2 = (av_v & (av_v / deriv1))(gfv1);
    
	const auto op1 = av_v & deriv1;

	gfv2 = op1(gfv1);

	const decltype(av_v & deriv1) op2 = av_v & deriv1;

	gfv2 = op2(gfv1);

	const auto op3 = operator&(av_v, deriv1);

	gfv2 = op3(gfv1);*/
	return 0;
}
