#include <iostream>
#include <fstream>
#include <set>
#include <vector>
#include <iterator>
#include "mpi.h"
#include <omp.h>
#include <stdio.h>
#include <math.h>  
#include <stdlib.h>

using namespace std;

MPI_Comm MCW = MPI_COMM_WORLD;
int Nx,Ny,K3,K4,Px,Py,subd_num,proc_num;

int comp (const void *i,const void *j)
{
    return (*(int*)i - *(int*)j);
}

void resort_local(int* mass,int dim){
    qsort(mass,dim,sizeof(mass[0]),comp);
}

void create_EN(int * &DEC_ext_IA,int* DEC_ext_JA,int Nx_local,int* &EN_IA,int* &EN_JA)
{
    int i,j,triangle,square,num;
    int halo_l,halo_t,halo_r,halo_b;
    bool halo = false,new_side = false,start = false;
    triangle = 3;
    square = 4;

    EN_IA[0] = 0;

    resort_local(DEC_ext_JA,DEC_ext_IA[2]);

    j = 1;
    for(i = 0; i < DEC_ext_IA[6]; i++) {
        if(i >= DEC_ext_IA[2]) halo = true;
        if(DEC_ext_JA[i] % (K3 + K4) < K3) {
            EN_IA[j] = EN_IA[j - 1] + triangle;
            j++;
            if(!halo) {
                EN_IA[j] = EN_IA[j - 1] + triangle;
                j++;
            }
        }
        else {
            EN_IA[j] = EN_IA[j - 1] + square;
            j++;
        }
    }
    i = 0;
    num = halo_l = halo_t = halo_r = halo_b = 0;
    halo = false;
    for(j = 0; j < DEC_ext_IA[6]; j++) {
        if((j == DEC_ext_IA[3] or j == DEC_ext_IA[4] or j == DEC_ext_IA[5] or j == DEC_ext_IA[6]) and j != DEC_ext_IA[2]) { i++; new_side = true; }
        else new_side = false;
        if((DEC_ext_JA[j - 1] % (K3 + K4)) < K3 and (DEC_ext_JA[j] % (K3 + K4)) < K3 and !new_side and halo) i++;
        if(j == DEC_ext_IA[2]) {
            halo = true;
            i = Nx_local + i;
            halo_l = 0;
            halo_t = 0;
            halo_r = Nx_local - 1;
            halo_b = i - Nx_local + 1;
        }
        if((DEC_ext_JA[j] % (K3 + K4)) < K3) {
            if(halo) {
                if(j < DEC_ext_IA[3]) {
                    if(((DEC_ext_JA[j - 1] % (K3 + K4)) >= K3) and !new_side and start) i++;
                    EN_JA[num] = halo_t;
                    EN_JA[num + 1] = i + 1;
                    EN_JA[num + 2] = halo_t + 1;
                    halo_t += 1;
                }
                if(j >= DEC_ext_IA[3] and j < DEC_ext_IA[4]){
                    EN_JA[num] = halo_l;
                    EN_JA[num + 1] = i + 1;
                    EN_JA[num + 2] = halo_l + Nx_local;
                    halo_l += Nx_local;
                }
                if(j >= DEC_ext_IA[4] and j < DEC_ext_IA[5]) {
                    EN_JA[num] = halo_r;
                    EN_JA[num + 1] = i + 1;
                    EN_JA[num + 2] = halo_r + Nx_local;
                    halo_r += Nx_local;
                }
                if(j >= DEC_ext_IA[5] and j < DEC_ext_IA[6]) {
                    if(((DEC_ext_JA[j - 1] % (K3 + K4)) >= K3) and !new_side and start) i++;
                    EN_JA[num] = halo_b;
                    EN_JA[num + 1] = i + 1;
                    EN_JA[num + 2] = halo_b + 1;
                    halo_b += 1;
                }
                num += 3;
                start = true;
            }
            else {
                if((i + 1) % Nx_local == 0) i++;
                EN_JA[num] = i;
                EN_JA[num + 1] = EN_JA[num + 3] = i + 1;
                EN_JA[num + 2] = EN_JA[num + 5] = Nx_local + i;
                EN_JA[num + 4] = Nx_local + i + 1;
                i++;
                num += 6;
            }
        }
        else {
            if(halo) {
                if(j < DEC_ext_IA[3]) {
                    if(((DEC_ext_JA[j - 1] % (K3 + K4)) < K3) and !new_side and start) i++;
                    EN_JA[num] = halo_t;
                    EN_JA[num + 3] = i + 1;
                    EN_JA[num + 1] = halo_t + 1;
                    EN_JA[num + 2] = i + 2;
                    halo_t += 1;
                }
                if(j >= DEC_ext_IA[3] and j < DEC_ext_IA[4]) {
                    EN_JA[num] = halo_l;
                    EN_JA[num + 3] = i + 1;
                    EN_JA[num + 1] = halo_l + Nx_local;
                    EN_JA[num + 2] = i + 2;
                    halo_l += Nx_local;
                }
                if(j >= DEC_ext_IA[4] and j < DEC_ext_IA[5]) {
                    if(((DEC_ext_JA[j - 1] % (K3 + K4)) < K3) and !new_side and start) i++;
                    EN_JA[num] = halo_r;
                    EN_JA[num + 3] = i + 1;
                    EN_JA[num + 1] = halo_r + Nx_local;
                    EN_JA[num + 2] = i + 2;
                    halo_r += Nx_local;
                }
                if(j >= DEC_ext_IA[5] and j < DEC_ext_IA[6]) {
                    EN_JA[num] = halo_b;
                    EN_JA[num + 3] = i + 1;
                    EN_JA[num + 1] = halo_b + 1;
                    EN_JA[num + 2] = i + 2;
                    halo_b += 1;
                }
                start = true;
            }
            else {
                if((i + 1) % Nx_local == 0) {
                    i++;
                }
                EN_JA[num] = i;
                EN_JA[num + 1] = i + 1;
                EN_JA[num + 2] = Nx_local + i + 1;
                EN_JA[num + 3] = Nx_local + i;
            }
            i++;
            num += 4;
        }
    }
}

void EN_to_NE(int* &EN_IA,int* &EN_JA,int* &NE_IA,int* &NE_JA,int EN_IA_dim,int EN_JA_dim)
{
    ofstream fout;
    fout.open("file.txt",std::ios::app);
    int i,j,k,NE_IA_dim,NE_JA_dim;
    NE_IA_dim = EN_JA[EN_JA_dim - 2] + 2;
    NE_JA_dim = EN_IA[EN_IA_dim - 1];
    for(j = 0; j < NE_IA_dim; j++) {
        NE_IA[j] = 0;
    }
    for(i = 0; i < EN_JA_dim; i++) NE_IA[EN_JA[i]+1]++;

    for(j = 1; j < NE_IA_dim; j++) NE_IA[j] = NE_IA[j] + NE_IA[j - 1];

    for(i = 0; i < NE_JA_dim; i++) {
        NE_JA[i] = -1;
    }
    k = 0;
    for(i = 0; i < EN_IA_dim - 1; i++) {
        for(j = EN_IA[i]; j < EN_IA[i + 1]; j++) {
            while(NE_JA[NE_IA[EN_JA[j]] + k] > -1) {
                k++;
            }
            NE_JA[NE_IA[EN_JA[j]] + k] = i;
            k = 0;
        }
    }
}

void EN_to_NN(int* EN_IA,int* EN_JA,int* NE_IA,int* NE_JA,int* NN_IA,int* NN_JA,int NE_IA_dim)
{
    int k,i,j,l,triangle,num;
    triangle = 3;
    EN_IA[0] = 0;
    int* vect;
    vect = new int[NE_IA_dim-1];
    for(i = 0; i < NE_IA_dim; i++) {
        NN_IA[i] = 0;
    }
    for(i = 0; i < NE_IA_dim-1; i++) {
        vect[i] = -1;
    }
    l = 0;
    for(k = 0; k < NE_IA_dim-1; k++) {
        for(i = NE_IA[k]; i < NE_IA[k + 1]; i++) {
            for(j = EN_IA[NE_JA[i]]; j < EN_IA[NE_JA[i]+1]; j++) {
                vect[EN_JA[j]] = EN_JA[j];
                if(EN_JA[j] == k) {
                    num = j;
                }
            }
            if(EN_IA[NE_JA[i] + 1] - EN_IA[NE_JA[i]] > triangle){
                if(num < EN_IA[NE_JA[i]] + 2) {
                    vect[EN_JA[num + 2]] = -1;
                }
                else {
                    vect[EN_JA[num - 2]] = -1;
                }
            }
        }
        for(i = 0; i < NE_IA_dim-1; i++) {
            if(vect[i] != -1) {
                NN_IA[k+1]++;
                NN_JA[l] = vect[i];
                l++;
            }
        }
        for(i = 0; i < NE_IA_dim-1; i++) {
            vect[i] = -1;
        }
    }
    NN_IA[0] = 0;
    for(i = 1; i < NE_IA_dim; i++) {
        NN_IA[i] = NN_IA[i - 1] + NN_IA[i];
    }
    delete[] vect;
}

void EN_to_EE(int* &EN_IA,int* &EN_JA,int* &NE_IA,int* &NE_JA,int* &EE_IA,int* &EE_JA,int EN_IA_dim)
{
    int k,i,j,l,m=0,vect_max=0;
    int* vect;
    int* vect_dim;
    vect_dim = new int[EN_IA_dim-1];
    for(k = 0; k < EN_IA_dim-1; k++) {
        vect_dim[k]=0;
        for(i = EN_IA[k]; i < EN_IA[k + 1]; i++) {
            vect_dim[k]+= (NE_IA[EN_JA[i]+1]-NE_IA[EN_JA[i]]-1);
        }
        if(vect_max<vect_dim[k]) vect_max = vect_dim[k];
    }
    vect = new int[vect_max+1];
    for(i = 0; i < EN_IA_dim; i++) {
        EE_IA[i] = 0;
    }
    for(k = 0; k < EN_IA_dim-1; k++) {
        for(i=0; i<vect_max; i++){
            vect[i]=0;
        }
        l = 0;
        vect[l]=k;
        l++;
        for(i = EN_IA[k]; i < EN_IA[k + 1]; i++) {
            for(j = NE_IA[EN_JA[i]]; j < NE_IA[EN_JA[i] + 1]; j++) {
                if(NE_JA[j] != k) { vect[l] = NE_JA[j]; l++; }
            }
        }
        resort_local(vect,vect_dim[k]+1);
        for(i=0; i<vect_dim[k]; i++){
            if(vect[i] == vect[i+1] or vect[i]==k){
                EE_IA[k+1]++;
                EE_JA[m]=vect[i];
                m++;
            }
        }
        if(vect[vect_dim[k]]==k) {
            EE_IA[k+1]++;
            EE_JA[m]=vect[i];
            m++;
        }
        l=0;
    }
    delete[] vect;
    delete[] vect_dim;
    EE_IA[0] = 0;
    for(i = 1; i < EN_IA_dim; i++) {
        EE_IA[i] = EE_IA[i - 1] + EE_IA[i];
    }
}

void print_topology(int* TOP_IA,int* TOP_JA,int TOP_IA_dim) {
    int i,j;
    ofstream fout("file.txt",std::ios::app);
    for(i = 0; i < TOP_IA_dim; i++)
    {
        fout << TOP_IA[i] << ' ';
    }
    fout << endl;
    for(j = 0; j < TOP_IA_dim-1; j++) {
        for(i = TOP_IA[j]; i < TOP_IA[j + 1]; i++)
        {
            fout << TOP_JA[i] << ' ';
        }
        fout << endl;
    }
    fout << endl;
    fout.close();
}

void matrix_fill(int* &TOP_IA,int* &TOP_JA,double* &A,int TOP_IA_dim,int Nx,int* &L2G,int* &Part,double* &b,int IA_dim) {
    int i,j,k,t,l;
    l = 0;
    for(i = 0; i < TOP_IA_dim; i++) {
        for(j = TOP_IA[i]; j < TOP_IA[i + 1]; j++)
        {
            k = L2G[TOP_JA[j]];
            t = L2G[i];
            A[l] = sin(t) + cos(k);
            l++;
        }
    }
    for(i = 0; i < IA_dim-1; i++){
        if(Part[i] == subd_num){
            b[i] = sin(L2G[i]*L2G[i]);
        }
    }
}

void matrix_vector_product(int* &IA,int* &JA,double* &A,double* &b,double* &res,int IA_dim,int* &L2G) {
    #pragma omp parallel for
        for(int i = 0; i < IA_dim; i++) {
            double result = 0;
            for(int j = IA[i]; j < IA[i + 1]; j++)
                result += A[j]*b[JA[j]];
            res[i] = result;
        }
}

void control_values(double& C_norm,double& L2_norm,double& elements_sum,double* &vect,int vect_dim) {
    int i;
    C_norm = L2_norm = elements_sum = 0;
    for(i = 0; i < vect_dim; i++) {
        if(vect[i] > C_norm) {
            C_norm = vect[i];
        }
        L2_norm += vect[i] * vect[i];
        elements_sum += vect[i];
    }
}

void decompose(int &Ne,int* (&elements)){
    int i,j,k,n,l;
    int* razm_x,*razm_y;
    razm_x = new int[Px+1];
    razm_y = new int[Py+1];
    razm_x[0] = razm_y[0] = Ne = 0;
    for(i = 1; i < Px+1; i++) {
        if((i-1) < (Nx - 1) % Px) {
            razm_x[i] = (Nx - 1) / Px + 1;
        }
        else {
            razm_x[i] = (Nx - 1) / Px;
        }
    }
    for(i = 1; i < Py+1; i++) {
        if((i-1) < (Ny - 1) % Py) {
            razm_y[i] = (Ny - 1) / Py + 1;
        }
        else {
            razm_y[i] = (Ny - 1) / Py;
        }
    }
    k = 0;
    i = subd_num/Px;
    j = (subd_num>=Px)?subd_num%Px:subd_num;
    Ne = razm_x[j+1] * razm_y[i+1];

    elements = new int[Ne];
    for(i = 1; i < Py+1; i++) {
        razm_y[i] += razm_y[i - 1];
    }
    for(i = 1; i < Px+1; i++) {
        razm_x[i] += razm_x[i - 1];
    }
    n = 0;
    i = subd_num/Px;
    j = (subd_num>=Px)?subd_num%Px:subd_num;
    for(k = razm_y[i]; k < razm_y[i + 1]; k++) {
        for(l = razm_x[j]; l < razm_x[j + 1]; l++) {
            elements[n] = k * (Nx - 1) + l;
            n++;
        }
    }
    delete[] razm_x;
    delete[] razm_y;
}

//функция формирует расширенную подобласть по массиву декомпозиции
void subdomain_ext(int Ne,int* &DEC_JA,int* &DEC_ext_IA,int *(&DEC_ext),int Nx_local) {
    int i,if_n,in_n,ha_r,ha_l,ha_t,ha_b;
    int* iface,*inner,*halo_r,*halo_l,*halo_t,*halo_b;
    iface = new int[Ne];
    inner = new int[Ne];
    halo_r = new int[(Ne+1)/(Nx_local-1)];
    halo_l = new int[(Ne+1)/(Nx_local-1)];
    halo_t = new int[Nx_local-1];
    halo_b = new int[Nx_local-1];
    if_n = in_n = ha_r = ha_l = ha_t = ha_b = 0;
    bool left,right,top,botton,iface_flg;
    for(i = 0; i < Ne; i++) {
        left = right = top = botton = true;
        iface_flg = false;
        fflush(stdout);

        //проверяем, если ли соседи у элемента в подобласти (или в принципе не может быть гало)
        if((i+1) % (Nx_local-1) == 0) right = false;
        if(i % (Nx_local-1) == 0) left = false;
        if(i < Nx_local-1) top = false;
        if(Ne-i < Nx_local) botton = false;
        if(!top) {
            if(DEC_JA[i] - Nx + 1 > 0) {
                halo_t[ha_t] = DEC_JA[i] - Nx + 1;
                ha_t++;
                if(!iface_flg) { iface[if_n] = DEC_JA[i]; if_n++; iface_flg = true; }
            }
        }
        if(!right) {
            if(DEC_JA[i] + 1 < (Nx - 1)*(Ny - 1)) {
                if((DEC_JA[i] + 1) % (Nx - 1) != 0) {
                    halo_r[ha_r] = DEC_JA[i] + 1;
                    ha_r++;
                }
                if(!iface_flg) { iface[if_n] = DEC_JA[i]; if_n++; iface_flg = true; }
            }
        }
        if(!botton) {
            if(DEC_JA[i] + Nx - 1 < (Nx - 1)*(Ny - 1)) {
                halo_b[ha_b] = DEC_JA[i] + Nx - 1;
                ha_b++;
                if(!iface_flg) { iface[if_n] = DEC_JA[i]; if_n++; iface_flg = true; }
            }
        }
        if(!left) {
            if(DEC_JA[i] - 1 >= 0) {
                if((DEC_JA[i]) % (Nx - 1) != 0) {
                    halo_l[ha_l] = DEC_JA[i] - 1;
                    ha_l++;
                }
                if(!iface_flg) { iface[if_n] = DEC_JA[i]; if_n++; iface_flg = true; }
            }
        }
        if(!iface_flg) {
            inner[in_n] = DEC_JA[i];
            in_n++;
        }
    }
    DEC_ext_IA[0] = 0;
    DEC_ext_IA[1] = in_n;
    DEC_ext_IA[2] = in_n + if_n;
    DEC_ext_IA[3] = in_n + if_n + ha_t;
    DEC_ext_IA[4] = DEC_ext_IA[3] + ha_l;
    DEC_ext_IA[5] = DEC_ext_IA[4] + ha_r;
    DEC_ext_IA[6] = DEC_ext_IA[5] + ha_b;
    DEC_ext = new int[DEC_ext_IA[6]];
    for(i = 0; i < in_n; i++) {
        DEC_ext[i] = inner[i];
    }
    for(i = 0; i < if_n; i++) {
        DEC_ext[in_n + i] = iface[i];
    }
    for(i = 0; i < ha_t; i++) {
        DEC_ext[in_n + if_n + i] = halo_t[i];
    }
    for(i = 0; i < ha_l; i++) {
        DEC_ext[DEC_ext_IA[3] + i] = halo_l[i];
    }
    for(i = 0; i < ha_r; i++) {
        DEC_ext[DEC_ext_IA[4] + i] = halo_r[i];
    }
    for(i = 0; i < ha_b; i++) {
        DEC_ext[DEC_ext_IA[5] + i] = halo_b[i];
    }
    delete[] inner;
    delete[] iface;
    delete[] halo_l;
    delete[] halo_r;
    delete[] halo_t;
    delete[] halo_b;
}

//функция возвращает массив, в котором каждой ячейки расширенной подобласти сопоставляется эта ячейка в глобальной нумерации
void get_L2G(int* &DEC_ext_IA,int* DEC_ext_JA,int dim_L2G,int* &L2G) {
    int i,j,tria_num;
    bool halo = false;
    j = 0;
    L2G = new int[dim_L2G];
    for(i = 0; i < DEC_ext_IA[6]; i++) {
        tria_num = (DEC_ext_JA[i] / (K3 + K4))*K3;
        if(DEC_ext_JA[i] % (K3 + K4) > K3) {
            tria_num += K3;
        }
        else {
            tria_num += (DEC_ext_JA[i] % (K3 + K4));
        }
        if(i >= DEC_ext_IA[2]) halo = true;
        if(DEC_ext_JA[i] % (K3 + K4) < K3) {
            if(halo) {
                if(i < DEC_ext_IA[4]) {
                    L2G[j] = DEC_ext_JA[i] + tria_num + 1;
                }
                else {
                    L2G[j] = DEC_ext_JA[i] + tria_num;
                }
            }
            else {
                L2G[j] = DEC_ext_JA[i] + tria_num;
                L2G[j + 1] = L2G[j] + 1;
                j++;
            }
        }
        else {
            L2G[j] = DEC_ext_JA[i] + tria_num;
        }
        j++;
    }
}

//функция возвращает массив, в котором каждой ячейке расширенной подобласти сопоставлена подобласть, которой она принадлежит
void get_Part(int* &DEC_ext_IA,int* DEC_ext_JA,int dim_Part,int* &Part) {
    int i,j,num;
    Part = new int[dim_Part];
    j = 0;
    num = subd_num;
    for(i = 0; i < DEC_ext_IA[6]; i++) {
        if(i >= DEC_ext_IA[2] and i < DEC_ext_IA[3]) num = subd_num - Px; //область сверху
        if(i >= DEC_ext_IA[3] and i < DEC_ext_IA[4]) num = subd_num - 1; //область слева
        if(i >= DEC_ext_IA[4] and i < DEC_ext_IA[5]) num = subd_num + 1; //область справа
        if(i >= DEC_ext_IA[5] and i < DEC_ext_IA[6]) num = subd_num + Px; //область снизу
        Part[j] = num;
        j++;
        if(DEC_ext_JA[i] % (K3 + K4) < K3 and num == subd_num) {
            Part[j] = num;
            j++;
        }
    }
}

//процедура обрабатывает подобласть
void init(int* &EN_IA,int* &EN_JA,int*(&EE_IA),int*(&EE_JA),int EN_IA_dim,int EN_JA_dim) {
    int* NE_IA,*NE_JA;
    int NE_IA_dim,NE_JA_dim;\
        NE_IA_dim = EN_JA[EN_JA_dim - 2] + 2;
    NE_JA_dim = EN_IA[EN_IA_dim - 1];
    NE_IA = new int[NE_IA_dim];
    NE_JA = new int[NE_JA_dim];

    EN_to_NE(EN_IA,EN_JA,NE_IA,NE_JA,EN_IA_dim,EN_JA_dim);

    EE_IA = new int[EN_IA_dim];
    EE_JA = new int[EN_IA_dim * 5];

    EN_to_EE(EN_IA,EN_JA,NE_IA,NE_JA,EE_IA,EE_JA,EN_IA_dim);
    delete[] NE_IA;
    delete[] NE_JA;
}

void post(int* &IA,int* &JA,int EE_IA_dim,int alien,double* &b,int* &cur_nb,int neighbours,vector<set<int> > Send_Buff,vector<set<int> > Recieve_Buff){
    int i,k,mpi_res;
    size_t res_size = 0;
    if(neighbours>0){
        MPI_Request* req;
        MPI_Status*  sts;
        req = new MPI_Request[neighbours*2];
        sts = new MPI_Status[neighbours*2];
        for(i = 0; i < neighbours; i++){
            res_size += Recieve_Buff[i].size();
        }
        double **send,**recieve;
        double* recieve_all;
        recieve_all = new double[res_size];
        send = new double*[neighbours];
        recieve = new double*[neighbours];
        for(i = 0; i < neighbours; i++){
            send[i] = new double[Send_Buff[i].size()];
            recieve[i] = new double[Recieve_Buff[i].size()];
            k = 0;
            for(std::set<int>::iterator j=Send_Buff[i].begin(); j!=Send_Buff[i].end(); ++j) {
                send[i][k] = b[*j];
                k++;
            }
            mpi_res = MPI_Isend(send[i],(int)Send_Buff[i].size(),MPI_DOUBLE,cur_nb[i],cur_nb[i]+subd_num,MCW,&req[i*2]);
            if(mpi_res!= MPI_SUCCESS)
                printf("MPI_Isend left failed (code %d)\n",mpi_res);
            mpi_res = MPI_Irecv(recieve[i],(int)Recieve_Buff[i].size(),MPI_DOUBLE,cur_nb[i],cur_nb[i]+subd_num,MCW,&req[i*2+1]);
            if(mpi_res!= MPI_SUCCESS)
                printf("MPI_Recv left failed (code %d)\n",mpi_res);
        }
        mpi_res = MPI_Waitall(neighbours*2,req,sts);
        if(mpi_res!= MPI_SUCCESS)
            printf("MPI_Waitall failed!\n");

        int m = 0;
        for(i=0; i<neighbours; i++){
            for(int ii = 0; ii < Recieve_Buff[i].size(); ii++){
                recieve_all[m] = recieve[i][ii];
                m++;
            }
            delete[] send[i];
            delete[] recieve[i];
        }
        k = 0;
        for(i = 0; i < EE_IA_dim-1; i++)
            if(i>=alien) { b[i] = recieve_all[k]; k++; }
        delete[] recieve_all;
        delete[] req;
        delete[] sts;
    }
}

void generation(int* &L2G,int* &Part,int* &EN_IA,int* &EN_JA,int &EN_IA_dim,int &Nx_local,int &EN_JA_dim){
    int i,DEC_IA_dim,DEC_JA_dim,triangle,square,dim_L2G,dim_Part,Ne=0;
    int *DEC_ext_JA;
    int* elements,*DEC_ext_IA;
    bool halo = false;
    triangle = 3;
    square = 4;
    DEC_IA_dim = Px * Py + 1;
    DEC_JA_dim = (Nx - 1)*(Ny - 1);
    decompose(Ne,elements);

    int tria_num = (DEC_JA_dim / (K3 + K4))*K3;
    if(DEC_JA_dim % (K3 + K4) > K3){
        tria_num += K3;
    }
    else {
        tria_num += ((Nx - 1)*(Ny - 1) % (K3 + K4));
    }
    DEC_ext_IA = new int[7];
    Nx_local = (Nx - 1) / Px + 1;
    if((subd_num % Px) < (Nx - 1) % Px) Nx_local++;
    subdomain_ext(Ne,elements,DEC_ext_IA,DEC_ext_JA,Nx_local);
    EN_IA_dim = 1;
    EN_JA_dim = 0;
    for(i = 0; i < DEC_ext_IA[6]; i++) {
        if(i >= DEC_ext_IA[2]) halo = true;
        if(DEC_ext_JA[i] % (K3 + K4) < K3) {
            if(halo) {
                EN_IA_dim++;
                EN_JA_dim += triangle;
            }
            else {
                EN_IA_dim += 2;
                EN_JA_dim += triangle * 2;
            }
        }
        else {
            EN_IA_dim++;
            EN_JA_dim += square;
        }
    }
    EN_IA = new int[EN_IA_dim];
    EN_JA = new int[EN_JA_dim];
    dim_L2G = EN_IA_dim - 1;
    dim_Part = dim_L2G;
    create_EN(DEC_ext_IA,DEC_ext_JA,Nx_local,EN_IA,EN_JA);
    get_L2G(DEC_ext_IA,DEC_ext_JA,dim_L2G,L2G);
    get_Part(DEC_ext_IA,DEC_ext_JA,dim_Part,Part);
    delete[] DEC_ext_IA;
    delete[] DEC_ext_JA;
    delete[] elements;
}

void create_data(int EE_IA_dim,int* &Part,int* &L2G,int* &cur_nb,int neighbours,vector<set<int> > &Send_Buff,vector<set<int> > &Recieve_Buff,int* &IA,int* &JA){
    int my_cells = 0,Neib_dem = 0;
    int* Neib_Num;
    for(int i=0; i<EE_IA_dim - 1; i++){
        if(Part[i] == subd_num)  my_cells++;
        if(Neib_dem < (Part[i]+1)) Neib_dem = Part[i]+1;
    }
    Neib_Num = new int[Neib_dem];
    for(int i = 0; i < Neib_dem; i++){
        Neib_Num[i] = -1;
    }
    int j = 0;
    for(int i=1; i<EE_IA_dim - 1; i++){
        if(Part[i] != Part[i-1]) {
            Neib_Num[Part[i]] = j;
            j++;
        }
    }
    for(int i=0; i < my_cells; i++){
        for(int j=IA[i]; j<IA[i+1]; j++){
            if(Part[JA[j]] != subd_num){
                Recieve_Buff[Neib_Num[Part[JA[j]]]].insert(JA[j]);
                Send_Buff[Neib_Num[Part[JA[j]]]].insert(i);
            }
        }
    }
    if(neighbours>0) cur_nb = new int[neighbours];
    for(int i=0; i<neighbours; i++)
        for(j = 0; j < Neib_dem; j++)
            if(Neib_Num[j] == i) cur_nb[i] = j;

    delete[] Neib_Num;
}


void SpMV_readme(int firstline = 0){
    printf("- SpMV: parallel program for matrix-vector product\n");
    if(firstline) return;
    printf("- obligatory parameters: \n <Nx> - numder of nodes in x axis \n <Ny> - number of nodes in y axis\n");
    printf("- additional parameters: \n [K3=2] - number of triangles \n [K4=2] - number of squares \n [T=1] - number of threads\n");
    printf("- Example: EN-NE.exe 1000 1000 3 4 2\n");
}


int main(int argc,char** argv)
{
    double C,L2,L1;
    int T=1;
    int mpi_res = MPI_Init(&argc,&argv);
    if(mpi_res!= MPI_SUCCESS)
        printf("MPI_Init failed (code %d)\n",mpi_res);

    mpi_res = MPI_Comm_size(MCW,&proc_num);
    if(mpi_res!= MPI_SUCCESS)
        printf("MPI_Init failed (code %d)\n",mpi_res);

    mpi_res = MPI_Comm_rank(MCW,&subd_num);
    if(mpi_res!= MPI_SUCCESS)
        printf("MPI_Init failed (code %d)\n",mpi_res);

    if(argc<3) { if(!subd_num) SpMV_readme(); MPI_Finalize(); return 0; }
    K3 = K4 = 2;
    Nx = atoi(argv[1]);
    Ny = atoi(argv[2]);
    if(argc>=4) K3 = atoi(argv[3]);
    if(argc>=5) K4 = atoi(argv[4]);
    if(argc>=6) T = atoi(argv[5]);

    if(Nx<=1 || Ny<=1 || K3==0 || K4==0 || proc_num > (Nx-1)*(Ny-1) || T==0){
        if(!subd_num) {
            cout << endl << "Wrong input parametrs!" << endl << endl;
            SpMV_readme();
        }
        MPI_Finalize(); return 0;
    }
    omp_set_num_threads(T);

    if(!subd_num) SpMV_readme(1);
    if(!subd_num) cout << "- Nx: "  << Nx << " Ny = " << Ny << " K3 = " << K3 << " K4 = " << K4 << endl;

    ofstream fout;
    fout.open("file.txt",std::ios::app);

    if(proc_num > 1){
        double a = (double(Ny-1)/double(Nx-1))*proc_num;
        Py = (int)sqrt(a);
        Py = (Py==0)?1:Py;
        Px=proc_num/Py;
        Px = (Px==0)?1:Px;
        Py=proc_num/Px;
        Py = (Py==0)?1:Py;
        if(Px>Py){
            while(proc_num%Px != 0) Px--;
            Py=proc_num/Px;
        }
        else {
            while(proc_num%Py != 0) Py--;
            Px=proc_num/Py;
        }
    }
    else Px=Py=1;
    if(subd_num==0) cout << "- Num_proc: "  << proc_num << " Px = " << Px << " Py = " << Py << " " << endl;

    int *EN_IA,*EN_JA,*L2G,*Part;
    int EN_IA_dim,Nx_local,EN_JA_dim;
    if(!subd_num) cout << "- Calculation area generation... " << endl;
    generation(L2G,Part,EN_IA,EN_JA,EN_IA_dim,Nx_local,EN_JA_dim);

    int *EE_IA,*EE_JA;
    if(!subd_num) cout << "- Matrix portrait generation... " << endl;
    init(EN_IA,EN_JA,EE_IA,EE_JA,EN_IA_dim,EN_JA_dim);

    delete[] EN_IA;
    delete[] EN_JA;

    int local_matrix_dim = 0;
    int alien = EN_IA_dim-1;
    for(int j = 0; j < EN_IA_dim-1; j++){
        if(Part[j] == subd_num) local_matrix_dim++;
        if(Part[j] != subd_num) {
            alien=j; break;
        }
    }

    double *res,*A,*b;
    b = new double[EN_IA_dim-1];
    A = new double[EE_IA[local_matrix_dim]];

    if(!subd_num) cout << "- Martix filling... " << endl;
    matrix_fill(EE_IA,EE_JA,A,local_matrix_dim,Nx_local,L2G,Part,b,EN_IA_dim);

    int neighbours=0;
    int *cur_nb;
    for(int i=0; i<EN_IA_dim - 1; i++) if(i > 0) if(Part[i] != Part[i-1]) neighbours++;
    vector<set<int> > Send_Buff(neighbours);
    vector<set<int> > Recieve_Buff(neighbours);
    if(!subd_num && proc_num>1) cout << "- Creating data for messaging... " << endl;
    create_data(EN_IA_dim,Part,L2G,cur_nb,neighbours,Send_Buff,Recieve_Buff,EE_IA,EE_JA);
    delete[] L2G;
    delete[] Part;
    res = new double[local_matrix_dim];


    int times = 50;
    if(!subd_num) cout << "- Calculation of the product... ("<< times << " times) " << endl;
    MPI_Barrier(MCW);
    double tcomp=0.0;
    double tbeg = MPI_Wtime();
    for(int i=0; i<times; i++){
        post(EE_IA,EE_JA,EN_IA_dim,alien,b,cur_nb,neighbours,Send_Buff,Recieve_Buff);
        matrix_vector_product(EE_IA,EE_JA,A,b,res,local_matrix_dim,L2G);
    }
    tcomp += MPI_Wtime() - tbeg;
    if(proc_num>1) delete[] cur_nb;
    double C_norm=0,L2_norm=0,elements_sum=0;
    if(!subd_num) cout << "- Control values calculation... " << endl;
    control_values(C_norm,L2_norm,elements_sum,res,local_matrix_dim);
    MPI_Barrier(MCW);
    MPI_Reduce(&L2_norm,&L2,1,MPI_DOUBLE,MPI_SUM,0,MCW);
    MPI_Reduce(&C_norm,&C,1,MPI_DOUBLE,MPI_MAX,0,MCW);
    MPI_Reduce(&elements_sum,&L1,1,MPI_DOUBLE,MPI_SUM,0,MCW);

    MPI_Barrier(MCW);
    if(!subd_num) cout << endl << "Results:" << endl;
    if(!subd_num) cout << "Time: "  << tcomp/double(times) << " " << endl;
    if(!subd_num) cout << "C = " << C << " L2 = " << L2 << " L1 = " << L1 << endl;
    delete[] EE_IA;
    delete[] EE_JA;
    delete[] A;
    delete[] b;
    delete[] res;
    MPI_Finalize();
    fout.close();
    return 0;
}
