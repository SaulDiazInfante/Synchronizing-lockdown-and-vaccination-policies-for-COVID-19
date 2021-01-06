#include "header_criterion"
{
    double a_d = constants[16];
    double a_e = constants[15];

    Tdouble D_zero = initial_state[5];
    Tdouble Y_zero = initial_state[8];

    Tdouble D_T = final_state[5];
    Tdouble Y_T = final_state[8];

    Tdouble x_zero_T = final_state[9];
    criterion = x_zero_T + a_d*(D_T - D_zero) + a_e*(Y_T - Y_zero);
}
