#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<vector>
#include<time.h>
using namespace std;

//Units used here:  meV for energy, Angstrom for distance

//default value
const double Pi=3.1415927;
const double k=8.62e-5*1000;  // Boltzmann's constant in meV/K 
const double Latt=7.68;   // unit in Angstrom
// Latt is 7.68 A on surface. but for E/atom for step formation energy, one atom => 3.84 A. 
const int N_plank_llimit=10;  // lower limit of plank number
const int N_plank_ulimit=3000;  // upper limit of plank number

double d_sigma;  // meV/A^2
double C2;   // meV/A 
//E_elastic=-Latt*C2*ln((l/Pi/a)*cos(Pi*p/2)), A/B domains' widths are (1+p)l and (1-p)l
//C2=d_sigma^2*(1-v)/(2*Pi*miu)  where d_sigma is the stress anisotropy, v is Poisson's ratio, miu is Bulk modulus
// Latt*C2 is about 0.023 eV 和kink与SA形成能在同一数量级

double Ef_kink;//=60 meV;       // kink formation energy for each kink
double Ef_SA;  //=40 meV/atom;  // formation energy of step A per site; Ef_SA in eV/A is 0.04/7.68 
int istep,N_print;
int T; //temperature
int H[2];  // terrace width
int A[2];  // amplitude width
int N_plank_init;  // initial plank number
int repetition;  // repetition of planks.
int N_plank; // = N_plank_init*repetition

int N_event[4];  // 0 for fluc, 1 for add, 2 for delete, 3 for shift
int N_acc[4];  // number of accept
int N_rej[4];  // number of reject
vector<int> plank_list0[2];   // trail
vector<int> plank_list1[2];   // real

FILE* fec;  // output current energy
FILE* fet;  // output trial energy
FILE* fall; // output trial energy per plank for all_energy
FILE* feall; // output all energy

void   accept_or_not(int flag, int Tind, int Pind);
double get_energy(int flag,vector<int> list[],FILE* fp,int Tind,int Pind); // for event 1/2/3, Tind/Pind are not necessary
double get_Ev_i(int h1,int h2);             // Elastic energy along vertical direction,    edge direction
double get_Eh_i(vector<int> list,int Pind); // Elastic energy along horizontal direction,  step-descending direction
int rand_int(int max);
int Pleft(int Pind,int size);
int Pright(int Pind,int size);
void print(int istep,FILE* fp);

int main(int argc, char *argv[]){
    int i,steps,Pind_c;
    int ii,j,Tind,Pind;
    int rep;
    int h_min,h_max,h_c,h_flag;
    int N;
    int P_left,P_right;
    srand((unsigned)time(NULL));
    fec=fopen("E_current.res","w");
    fet=fopen("E_trail.res","w");
    fall=fopen("Event.res","w");
    feall=fopen("E.res","w");
    fprintf(fec, "#  istep  flag  E_kink  E_SA  E_v  E_h  E_all  E_all/area  # for flag=3 only output E_v\n");
    fprintf(fet, "#  istep  flag  E_kink  E_SA  E_v  E_h  E_all  E_all/area  # for flag=3 only output E_v\n");
    fprintf(fall,"#  istep flag acc/rej\n");
    fprintf(feall,"#  istep  E_kink  E_SA  E_v  E_h  E_all  E_all/area \n");
    steps=atoi(argv[1]);
    N_print=steps/100;
    T=atoi(argv[2]);
    H[0]=atoi(argv[3]);
    H[1]=atoi(argv[4]);
    N_plank_init=atoi(argv[5]);
    repetition=atoi(argv[6]);
    N_plank=N_plank_init*repetition;
    d_sigma=atof(argv[7]);   // stress anisotropy in meV/A^2
    C2=d_sigma*d_sigma*(1-0.27)/(2*Pi*98*6.24);   // unit in meV/A 
    //  to meet the result of Phys. Rev. Lett. 61 (1988) 1973;
    //C2=d_sigma^2*(1-v)/(2*pi*miu) where d_sigma is the surface stress anisotropy, v is Poisson's ratio (=0.27)
    //miu is Bulk modulus 98 GPa, while 1 GPa =6.24 meV/A^3
      
    Ef_kink=atof(argv[8]);  // kink formation energy for each kink, in meV/site
    Ef_SA=atof(argv[9]);   //formation energy of step A per site, in meV/atom
    Ef_SA=Ef_SA*2;  // meV/atom to meV/latt

    A[0]=int(H[0]*atof(argv[10])); // Amplitude
    A[1]=int(H[1]*atof(argv[11])); // Amplitude
    if(A[0]<6) A[0]=8;
    if(A[1]<6) A[1]=8;
    printf("#steps T A0 H0 A1 H1 N_plank_init repetition C2 Ef_kink(meV/site) Ef_SA(meV/atom)\n");
    printf("%d %d %d %d %d %d %d %d %f %f %f\n",steps,T,A[0],H[0],A[1],H[1],N_plank_init,repetition,C2,Ef_kink,Ef_SA);

    // ./a.out steps T H N_plank_init elasit_scale
    // initial
    for(ii=0;ii<2;ii++){ 
        for(rep=0;rep<repetition;rep++){
            for(i=0;i<N_plank_init;i++){
                //plank_list1.push_back(H/2);   // straight init
                if(i<N_plank_init/2)
                    plank_list0[ii].push_back((int)((i*(H[ii]-6)*2/N_plank_init+3)*A[ii]/H[ii]+(H[ii]-A[ii])/2));   // ZZ init
                else
                    plank_list0[ii].push_back((int)(((N_plank_init-i)*(H[ii]-6)*2/N_plank_init+3)*A[ii]/H[ii]+(H[ii]-A[ii])/2));   // ZZ init
            }
        }
        plank_list1[ii]=plank_list0[ii];
    }

    for(ii=0;ii<4;ii++){
        N_event[ii]=0;
        N_acc[ii]  =0;
        N_rej[ii]  =0;
    }

    // run kmc
    for(istep=0;istep<=steps;istep++){
        if(istep%N_print==0)
            print(istep,feall);

        // dimer row elongation/shrink
        Tind=rand_int(2);
        Pind=rand_int(plank_list0[0].size());
        if((double)rand()/RAND_MAX>0.5)
            plank_list0[Tind][Pind]++;  // checked in get_Ev_i so that h1 and h2 cannot be 0
        else
            plank_list0[Tind][Pind]--;
        accept_or_not(0,Tind,Pind);   // 0 means fluc

        if(istep<steps/10){
            if(istep%100==0) { // shift every 100 kmc steps; always let Pind=0 terrace shift
                N=plank_list0[0].size();
                if((double)rand()/RAND_MAX>0.5) { // left shift
                    for(i=0;i<N-1;i++)
                        plank_list0[0][i]=plank_list1[0][i+1];
                    plank_list0[0][N-1]=plank_list1[0][0];
                }
                else{  // right shift
                    for(i=1;i<N;i++)
                        plank_list0[0][i]=plank_list1[0][i-1];
                    plank_list0[0][0]=plank_list1[0][N-1];
                }
                accept_or_not(3,0,0);   // 3 means shift
            }
            if(istep%100==50){   // speed up the amplitude calc
                N=plank_list0[0].size();
                Tind=rand_int(2);
                // get h_min,hmax of list[Tind]
                h_max=0; h_min=H[Tind];
                for(i=0;i<N;i++){
                    if(plank_list0[Tind][i]>h_max)  h_max=plank_list0[Tind][i];
                    if(plank_list0[Tind][i]<h_min)  h_min=plank_list0[Tind][i];
                }
                if(h_max-h_min<20)  continue;  // if dh<20, not necessary to speed up
                printf("h_max, h_min %d  %d\n", h_max, h_min);    
                // get h_c
                j=0;
                while(true){
                    j++;
                    h_c=rand_int(h_max-h_min-10)+h_min+5;
                    h_flag=0;
                    for(i=0;i<N;i++)
                        if(abs(plank_list0[Tind][i]-h_c)==0){
                            h_flag=1;
                            break;
                        }
                    if(h_flag==0 || j==100)   break;
                }
                if(h_flag==1) continue; // if not easy to find possible site, then just give up the trail.
                if(h_c>H[Tind]/2){  // h_c in top half part
                    if((double)rand()/RAND_MAX>0.5) { //expand
                        for(i=0;i<N;i++)
                            if(plank_list0[Tind][i]>h_c)
                                plank_list0[Tind][i]++;
                    }
                    else{                            // shrink
                        for(i=0;i<N;i++)
                            if(plank_list0[Tind][i]>h_c)
                                plank_list0[Tind][i]--;
                    }
                }
                else{  // h_c in bottom half part
                    if((double)rand()/RAND_MAX>0.5) { //expand
                        for(i=0;i<N;i++)
                            if(plank_list0[Tind][i]<h_c)
                                plank_list0[Tind][i]--;
                    }
                    else{                            // shrink
                        for(i=0;i<N;i++)
                            if(plank_list0[Tind][i]<h_c)
                                plank_list0[Tind][i]++;
                    }
                }
                accept_or_not(4,Tind,0);
            }
        }

        //if(istep%1000==0) { // Add/Del plank every 1000 mc steps
        //    N=plank_list0[0].size();
        //    j=0;
        //    while(true){
        //        Pind_c=rand_int(N);   //Pind_c is Pind_choose
        //        P_left=Pleft(Pind_c,N);
        //        P_right=Pright(Pind_c,N);

        //        j++;
        //        if(plank_list0[0][Pind_c]==plank_list0[0][P_left] || plank_list0[0][Pind_c]==plank_list0[0][P_right])
        //            if(plank_list0[1][Pind_c]==plank_list0[1][P_left] || plank_list0[1][Pind_c]==plank_list0[1][P_right])
        //                break;

        //        if(j==100)
        //            break;
        //    }
        //    if(j<100){
        //        if((double)rand()/RAND_MAX>0.5){   //add one plank at random position
        //            if(N<N_plank_ulimit){
        //                plank_list0[0].insert(plank_list0[0].begin()+Pind_c,plank_list0[0][Pind_c]);
        //                plank_list0[1].insert(plank_list0[1].begin()+Pind_c,plank_list0[1][Pind_c]);
        //                accept_or_not(1,0,0);
        //            } 
        //        }
        //        else{  // delete one plank at random position
        //            if(N>N_plank_llimit){ 
        //                plank_list0[0].erase(plank_list0[0].begin()+Pind_c);
        //                plank_list0[1].erase(plank_list0[1].begin()+Pind_c);
        //                accept_or_not(2,0,0);
        //            }
        //        }
        //    }
        //}

    }
    printf("Fluc(all/acc/rej):   %d %d %d\n",N_event[0],N_acc[0],N_rej[0]);
    printf("Add (all/acc/rej):   %d %d %d\n",N_event[1],N_acc[1],N_rej[1]);
    printf("Del (all/acc/rej):   %d %d %d\n",N_event[2],N_acc[2],N_rej[2]);
    printf("Shift (all/acc/rej): %d %d %d\n",N_event[3],N_acc[3],N_rej[3]);
    fclose(fec); fclose(fet); fclose(fall);
    return 0;
}

void accept_or_not(int flag,int Tind,int Pind){
    double E0,E1,dE;
    double poss,poss2;
    int Tnb=(Tind+1)%2;   // the index of another terrace; nb means neighbor
    N_event[flag]++;
    E0=get_energy(flag,plank_list0,fet,Tind,Pind); // trail
    E1=get_energy(flag,plank_list1,fec,Tind,Pind);  // current

    if(flag==0||flag==3||flag==4)     //fluc or shift or amplitude
        dE=E0-E1;
    else if(flag==1)   // add    E0 is E_new, E1 is E_old
        dE=(E0-E1)-E1/plank_list1[0].size();
    else if(flag==2)   // del    E0 is E_new, E1 is E_old
        dE=(E0-E1)+E1/plank_list1[0].size();
    else
        printf("Error in accept_or_not func!\n");

    if(istep%N_print==0)
        fprintf(fall,"step: %12d flag: %1d ----- ",istep,flag);

    if(dE<0){   // accept;   E1 for current step, E0 for trail step
        if(flag==0)  // for fluc, only updata [Tind][Pind]
            plank_list1[Tind][Pind]=plank_list0[Tind][Pind];
        else if(flag==3 || flag==4)  //For shift, only Tind=0, for amplitude, it can be 1
            plank_list1[Tind]=plank_list0[Tind];
        else if(flag==1 || flag==2){   // add/del should update all terrace;
            plank_list1[0]=plank_list0[0];
            plank_list1[1]=plank_list0[1];
        }
        N_acc[flag]++;
        if(istep%N_print==0)
            fprintf(fall,"Accept!\n");
        return;
    }
    poss=(double)rand()/RAND_MAX;
    poss2=exp(-dE/(k*T));

    if(poss<poss2){   // accept
        if(flag==0)  // for fluc, only updata [Tind][Pind]
            plank_list1[Tind][Pind]=plank_list0[Tind][Pind];
        else if(flag==3 || flag==4)  //For shift, only Tind=0, for amplitude, it can be 1
            plank_list1[Tind]=plank_list0[Tind];
        else if(flag==1 || flag==2){   // add/del should update all terrace;
            plank_list1[0]=plank_list0[0];
            plank_list1[1]=plank_list0[1];
        }
        N_acc[flag]++;
        if(istep%N_print==0)
            fprintf(fall,"Accept!  poss,poss2:  %f   %f\n",poss,poss2);
    }
    else{       // reject
        if(flag==0)  // for fluc, only updata [Tind][Pind]
            plank_list0[Tind][Pind]=plank_list1[Tind][Pind];
        else if(flag==3 || flag==4)  //For shift, only Tind=0, for amplitude, it can be 1
            plank_list0[Tind]=plank_list1[Tind];
        else if(flag==1 || flag==2){   // add/del should update all terrace;
            plank_list0[0]=plank_list1[0];
            plank_list0[1]=plank_list1[1];
        }
        N_rej[flag]++;
        if(istep%N_print==0)
            fprintf(fall,"Reject!  poss,poss2:  %f   %f\n",poss,poss2);
    }
}

double get_energy(int flag,vector<int> list[],FILE* fp,int Tind,int Pind){
    int    i,dh;
    double E_all, E_stepA, E_kink;
    double E_v;  // elastic energy along step downward direction
    double E_h;  // elastic energy paralle to step 
    int N=list[Tind].size();
    int Tnb=(Tind+1)%2;   // the index of another terrace; nb means neighbor
    int P_left,P_right;

    E_stepA=0;E_kink=0;E_v=0;E_h=0;E_all=0;
    
    if(flag==0){  // fluc,   need to calc part of E_v, part of E_h, and part of E_kink/E_stepA
        //E_v   for only Pind related three cases
        E_v= get_Ev_i(list[Tind][Pind], H[Tind]-list[Tind][Pind]);
        E_v+=get_Ev_i(list[Tind][Pind],  H[Tnb]-list[Tnb][Pind]);
        E_v+=get_Ev_i(list[Tnb][Pind],  H[Tind]-list[Tind][Pind]);
        // not necessary to calc (list[Tnb], H-list[Tnb])

        //E_kink and E_stepA
        P_left= Pleft(Pind,N);
        P_right=Pright(Pind,N);
        if(list[Tind][P_left]!=list[Tind][Pind]){
                E_kink+=Ef_kink;
                E_stepA+=Ef_SA*abs(list[Tind][P_left]-list[Tind][Pind]);
        }
        if(list[Tind][P_right]!=list[Tind][Pind]){
                E_kink+=Ef_kink;
                E_stepA+=Ef_SA*abs(list[Tind][P_right]-list[Tind][Pind]);
        }

        // E_h; elastic energy along horizontal direction
        for(i=0;i<N;i++)    // not necessary to calc list[Tnb]
            E_h+=get_Eh_i(list[Tind],i);

        E_all=E_stepA+E_kink+E_v+E_h;
        if(istep%N_print==0)  // flag==3
            fprintf(fp,"%10d  %1d  %.1f  %.1f  %.2f  %.2f  %.2f None\n",istep,flag,E_kink,E_stepA,E_v,E_h,E_all);
        return E_all;
    }
    else if(flag==1 || flag==2 || flag==4 || flag==5){  // add/del plank, calc all energy 
        //  Eventhough the E_v doesn't change at this process
        //  the total energy is necessary for evaluating the accept possibility
        //  for amplitude change, not necessary to calc all energy , but don't bring large amount extra calc
        //  so let it be.
        // flag==5 for print() function to output data to E.res file
        for(i=0;i<N;i++){
            P_right=Pright(i,N);
            // for Tind=0
            dh=abs(list[0][P_right]-list[0][i]);
            if(dh>0){
                E_kink+=Ef_kink;
                E_stepA+=Ef_SA*dh;
            }   
            // for Tind=1
            dh=abs(list[1][P_right]-list[1][i]);
            if(dh>0){
                E_kink+=Ef_kink;
                E_stepA+=Ef_SA*dh;
            }   
            // E_v
            E_v+=get_Ev_i(list[0][i],H[0]-list[0][i]);
            E_v+=get_Ev_i(list[0][i],H[1]-list[1][i]);
            E_v+=get_Ev_i(list[1][i],H[0]-list[0][i]);
            E_v+=get_Ev_i(list[1][i],H[1]-list[1][i]);
            // E_h
            E_h+=get_Eh_i(list[0],i);
            E_h+=get_Eh_i(list[1],i);
        }

        E_all=E_stepA+E_kink+E_v+E_h;

        if(flag==5)
            //fprintf(feall,"#   istep E_kink  E_SA  E_v E_h  E_all E_all/area \n");
            fprintf(fp,"%10d  %.1f  %.1f  %.2f  %.2f  %.2f  %f\n",istep,E_kink,E_stepA,E_v,E_h,E_all,E_all/(Latt*Latt*2*H[0]*N_plank));
        else if(istep%N_print==0 )
            //fprintf(fec, "#   istep  flag  E_kink  E_SA  E_v  E_h  E_all  E_all/area  # for flag=3 only output E_v\n");
            fprintf(fp,"%10d  %1d  %.1f  %.1f  %.2f  %.2f  %.2f  %f\n",istep,flag,E_kink,E_stepA,E_v,E_h,E_all,E_all/(Latt*Latt*2*H[0]*N_plank));
        return E_all;
    }
    else if(flag==3){ // shift, only calc part of E_v,since E_h/E_kink/E_stepA donot change.
        for(i=0;i<N;i++){
            E_v+=get_Ev_i(list[1][i],H[0]-list[0][i]);
            E_v+=get_Ev_i(list[0][i],H[1]-list[1][i]);
        }

        return E_v;
    }
    else
        printf("Error in calc E_all!");
}

double get_Ev_i(int h1,int h2){  
    // elastic energy for one dimer chain along step downwald direction
    // elastic*log((HH/2/Pi)*cos((Pi/2)*((double)(2*h_i-HH)/HH)))/2;  
    // log is ln in C language
    // C2=3meV/A; a=7.69A; refer to PRB1993,47,13432; 
    double W_temp,p,E_gap;
    if(h1==0 || h2==0)   // don't let h1 or h2 be 0;
        return 10000000;
    else if(h1==1 || h2==1)
        E_gap=150;
    else if(h1==2 || h2==2)
        E_gap=50;
    else
        E_gap=0;
    // -a*C2*ln((l/Pi/a)*cos(Pi*p/2)), A/B domains' widths are (1+p)l and (1-p)l
    // the variable l=(h1+h2)/2 in paper Phys. Rev. Lett. 1990, 64, 2406.
    // l/Pi/a = ( (h1+h2)*Latt/2 ) /Pi/(Latt/2) = (h1+h2)/Pi, has no unit. 'a' in formula is same as Latt/2
    // p=(h1-h2)/(h1+h2) is a ratio and has no unit
    W_temp=(double)(h1+h2);
    p=(double)(h1-h2)/(h1+h2);
    return E_gap-(2*Latt*C2)*log(2*W_temp/2/Pi)*cos(p*Pi/2)/2;  //in meV
    // return -(2*latt*C2)*np.log(((Latt/a)*W_temp/2/np.pi)*np.cos(p*np.pi/2))/2 # in meV
    // a is the Lorentzian broadening of the force density, a=Latt/2
    // the final divided by 2 for double-counting-like
}

double get_Eh_i(vector<int> list, int Pind){   
    // elastic energy for one dimer unit parallel to step direction
    double v=0;
    int Ly1,Ly2,height,Pind_temp;
    int size=list.size();
    if(list[Pind]>list[Pright(Pind,size)]) {   //Pright side edge
        for(height=list[Pright(Pind,size)]+1;height<=list[Pind];height++){   
            // get Ly1; Pleft side
            Ly1=0;
            Pind_temp=Pleft(Pind,size);
            while(true){
                Ly1++;
                if(list[Pind_temp]<height)
                    break;
                Pind_temp=Pleft(Pind_temp,size);
            }
            // get Ly2; Pright side
            Ly2=-1;
            Pind_temp=Pright(Pind,size);
            while(true){
                Ly2++;
                if(list[Pind_temp]>=height)
                    break;
                Pind_temp=Pright(Pind_temp,size);
            }
            v+=get_Ev_i(Ly1,Ly2); 
        }
    }
    if(list[Pind]>list[Pleft(Pind,size)]) {   //Pleft side edge
        for(height=list[Pleft(Pind,size)]+1;height<=list[Pind];height++){   
            // get Ly1; Pleft side
            Ly1=-1;
            Pind_temp=Pleft(Pind,size);
            while(true){
                Ly1++;
                if(list[Pind_temp]>=height)
                    break;
                Pind_temp=Pleft(Pind_temp,size);
            }
            // get Ly2
            Ly2=0;
            Pind_temp=Pright(Pind,size);
            while(true){
                Ly2++;
                if(list[Pind_temp]<height)
                    break;
                Pind_temp=Pright(Pind_temp,size);
            }
            v+=get_Ev_i(Ly1,Ly2);   // double counting -like
        }
    }
    return v;
}

int Pleft(int Pind,int size){
    if(Pind>0)
        return Pind-1;
    else if(Pind==0)
        return size-1;
    else
        printf("Error in fPinding Pleft!");
}

int Pright(int Pind,int size){
    if(Pind<size-1)
        return Pind+1;
    else if(Pind==size-1)
        return 0;
    else
        printf("Error in fPinding Pright!");
}

int rand_int(int max){ // return a random integer in the range [0,max-1]
      //return int(max*(double)rand()/RAND_MAX);
      return rand()%max;
}

void print(int istep,FILE* fp){
    int i;
    FILE* f;
    char s[255];
    sprintf(s,"%d.txt",istep);
    f=fopen(s,"w");
    for(i=0;i<plank_list1[0].size();i++)
        fprintf(f,"%d   %d   %d \n",i,plank_list1[0][i],plank_list1[1][i]+H[0]);
    fclose(f);

    get_energy(5,plank_list1,fp,0,0);
}
