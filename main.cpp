#include <iostream>
#include <fstream>
#include "cmath"
#include "stdlib.h"
#include "time.h"

const float X_A = 0.0;
const float X_B = 4.0;
const float Y_A = 0.0;
const float Y_B = 4.0;

const int N_x = 600;
const int N_y = 600;
const int N_t = 100;

float tao = (N_x > 1000 || N_y > 1000) ? 1e-3 : 1e-2;

float H_x = (X_B - X_A)/(N_x - 1);
float H_y = (Y_B - Y_A)/(N_y - 1);

int src_i = 1;
int src_j = 1;


float func_src_nm(int i, int j, int n_moment){
    const float f_nul = 1.0;
    const float t_nul = 1.5;
    const float gamma = 4.0;

    float d = (2 * M_PI * f_nul * (n_moment * tao - t_nul));

    if(i == src_i && j == src_j){
        float q = (float)exp((double)((-1 * d*d) / (gamma * gamma)));
        float s = (float)sin((double)(d));
        return q * s;
    }else{
        return 0.0;
    }
}


float* fill_nul(float* arr){
    for(int i = 0; i < N_x; i++){
        for(int j = 0; j < N_y;j++){
            arr[i*N_y + j] = 0;
        }
    }
    return arr;
}

float* init_phase_arr(){
    float* phase_arr = new float[N_x* N_y];
    for(int i = 0; i < N_x - 1; i++)
        for(int j = 0; j < N_y - 1; j++){
            phase_arr[i * N_y + j] = (j >= (N_y / 2))? 0.04 : 0.01;
        }
    return phase_arr;
}


float** initialize(){
    float** storage = new float*[3];

    float* u_min = (float*)calloc(N_y * N_x,sizeof(float)); //new float[N_x * N_y];
        storage[0] = fill_nul(u_min);
    float* u_nul = (float*)calloc(N_y * N_x,sizeof(float));
        storage[1] = fill_nul(u_nul);
    float* u_plus = (float*)calloc(N_y * N_x,sizeof(float));
        storage[2] = fill_nul(u_plus);

    return storage;
}

float k_x(float** storage,int i , int j, float* phase_arr){
    float diff_u_syl1 = (storage[1][i * N_y + j + 1] - storage[1][i * N_y + j]);
    float diff_phase_syl1 = (phase_arr[(i - 1) * N_y + j] + phase_arr[i * N_y + j]);

    float diff_u_syl2 = (storage[1][i * N_y + j - 1] - storage[1][i * N_y + j]);
    float diff_phase_syl2 = (phase_arr[(i - 1) * N_y + j - 1] + phase_arr[i * N_y + j - 1]);

    return (diff_phase_syl1 * diff_u_syl1 + diff_phase_syl2 * diff_u_syl2)
            / (2 * H_x * H_x);
}

float k_y(float** storage,int i , int j, float* phase_arr){
    float diff_u_syl1 = (storage[1][(i + 1) * N_y + j] - storage[1][i * N_y + j]);
    float diff_phase_syl1 = (phase_arr[i * N_y + (j - 1)] + phase_arr[i * N_y + j]);

    float diff_u_syl2 = (storage[1][(i - 1) * N_y + j] - storage[1][i * N_y + j]);
    float diff_phase_syl2 = (phase_arr[(i - 1) * N_y + (j - 1)] + phase_arr[(i - 1) * N_y + j]);

    return (diff_phase_syl1 * diff_u_syl1 + diff_phase_syl2 * diff_u_syl2)
            / (2 * H_y * H_y);
}

float sub_func(float** storage,int i, int j, int count_step, float* phase_arr){
    float f = func_src_nm(i,j,count_step);
    float k_xx = k_x(storage, i, j, phase_arr);
    float k_yy = k_y(storage, i, j, phase_arr);
    return f + k_xx + k_yy;
}

void print_result(float* arr){
    for(int i = 0;i < N_x;i++) {
        for (int j = 0; j < N_y; j++) {
            printf("%.3f ", arr[i * N_y + j]);
        }
        printf("\n");
    }
}

void step(float** storage, float* phase_arr,int count_step){
    for(int i = 1; i < N_x - 2;i++)
        for(int j = 1; j < N_y - 2;j++) {
            storage[2][i * N_y + j] = 2 * storage[1][i * N_y + j] - storage[0][i * N_y + j] + (tao * tao * sub_func(storage, i, j, count_step, phase_arr));
        }

    float* tmp = storage[0];
        storage[0] = storage[1];
        storage[1] = storage[2];
        storage[2] = tmp;

}


float* max_elem_matrix(float* matrix){
    float* p_max_elem = (float*)malloc(sizeof(float));
    *p_max_elem = -6.0;
    for(int i = 1; i < N_x - 1; i++)
        for(int j = 1; j < N_y - 1; j++)
            if(*p_max_elem < matrix[i * N_y + j])
                *p_max_elem = matrix[i * N_y + j];
    return p_max_elem;
}


void mko(){
    float* service_check_arr = (float*)calloc(N_t, sizeof(float));
    float** storage = initialize();
    float* phase_arr = init_phase_arr();

    clock_t start_time = clock();
    for(int s_t = 2; s_t < N_t; s_t++){
        step(storage, phase_arr, s_t);
        service_check_arr[s_t] = *(max_elem_matrix(storage[2]));//!!!
    }
    clock_t end_time = clock();
    clock_t diff = end_time - start_time;
    double current_time = (double)diff / CLOCKS_PER_SEC;
    std::cout << "Time: " << current_time << std::endl;

    for(int i = 0; i < N_t; i++)
        printf("%.8f ",service_check_arr[i]);


    std::ofstream out;
    out.open("float600x600.txt");
    if (out.is_open()) {
        std::cout << "Opened";
        for (unsigned it = 0; it < N_y; ++it) {
            for (unsigned j = 0; j < N_x; ++j) {
                out << storage[1][it * N_x + j] << " ";
            }
            out << "\n";
        }
    }
    out.close();


    FILE * fdat = fopen("data600-11.dat", "wb");
    fwrite(storage[1], sizeof(float), N_x * N_y, fdat);
    fclose(fdat);
}

int main(int argc, char* argv[]) {
          mko();

    return 0;
}


