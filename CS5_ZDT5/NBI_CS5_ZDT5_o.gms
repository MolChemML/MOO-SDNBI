$title NBI test5 (NBI_CS5_ZDT5_o.gms)

$gdxin inputcs5_nbi.gdx

* ----------------------------------------------
*  Declare sets
*  : i, j, m will be assgiend by external interface
* ----------------------------------------------
Sets
i        number of objectivse
j        number of variables
m        number of points
k        iteration number /k1*k21/
;
$load    i j m 

Table  w_data(k,i)  weight vector
    i1      i2
k1  0        1
k2  0.05     0.95
k3  0.1      0.9
k4  0.15     0.85
k5  0.2      0.8
k6  0.25     0.75
k7  0.3      0.7
k8  0.35     0.65
k9  0.4      0.6
k10 0.45     0.55
k11 0.5      0.5
k12 0.55     0.45
k13 0.6      0.4
k14 0.65     0.35
k15 0.7      0.3
k16 0.75     0.25
k17 0.8      0.2
k18 0.85     0.15
k19 0.9      0.1
k20 0.95     0.05

;

* ----------------------------------------------
*  Declare paramters
* ----------------------------------------------
Parameters
beta(m)
f_utopia(i) utopian points in objective space /i1=1, i2=0.3225806452/
f_nadir(i)  nadir points in objective space   /i1=31, i2=10/
w(i)        weights for weighted sum method   /i1=1, i2=1/
phi(i,m)    NBI parameters for choosing referece points
normal(i)   direction vector for NBI
f_o(i)      origin to shift the objective space
yo_i
y_i(j)
*for integer cut
*yv(c,j) store x values from previous solutions
;
$load    beta normal f_o phi y_i yo_i

* ----------------------------------------------
*  Store results in the txt file
* ----------------------------------------------
file res /NBI_CS5_result.txt/;
put res;
put '---------SIMULATION START------------------'/;

Variables
t            objective function to maximize
objective    nbi objective function
objective_w  scalarized objective fucntion
obj_f(i)     each objective
obj_g        objective multiplier
yo           variables1
y(j)         variables2


v(j)
s(j)
y_1(j)
y_2(j)



nbi_01
nbi_02
;

Binary variable s;
Integer variable yo, y, y_1, y_2;

* ------------------------------------------------
*  Provide upper and lower bounds of variables
* ------------------------------------------------
* lower and upper bound of x
y.LO(j)=0;
y.UP(j)=5;
yo.UP=30;
*x.LO('j2')=0;
*x.UP('j2')=3;
t.LO=-10;
t.UP=10;



* initial gusses of x
obj_g.L=10;
*x.L(j)=x_i(j);
yo.L = yo_i;
y.L(j)=y_i(j);

Equations
OBJ       objective function of nbi
OBJ_w     weighted sum objective fuction
OBJ_G2     multiplier for 2nd objective
OBJ_F1    first objective
OBJ_F2    second objective
OBJ_F1_NBI  first objective
OBJ_F2_NBI  second objective

EQ_V
AX_EQ1
AX_EQ2
AX_EQ3


INEQ_NBI_01  constraints2
INEQ_NBI_02  constraints2

*NORMAL_SUM

;

OBJ      .. objective =e= t ;
OBJ_w    .. objective_w =e= sum(i,w(i)*obj_f(i));
OBJ_G2   .. obj_g =e=  sum(j,v(j));

EQ_V(j)  .. v(j) =e= (2+y_1(j))*s(j) + (1-s(j));
AX_EQ1(j).. y_1(j)+y_2(j) =e= y(j);
AX_EQ2(j).. y_1(j)-4*s(j) =l= 0;
AX_EQ3(j).. y_2(j)-5*(1-s(j)) =e= 0;

OBJ_F1   .. obj_f('i1') =e= ( (1+yo) -f_utopia('i1') )/(f_nadir('i1')-f_utopia('i1'));
OBJ_F2   .. obj_f('i2') =e= ( (1/(1+yo))*obj_g -f_utopia('i2') )/(f_nadir('i2')-f_utopia('i2'));

OBJ_F1_NBI .. nbi_01 =e= sum(m, phi('i1',m)*beta(m)) + t*normal('i1') - (obj_f('i1')-f_o('i1'));
OBJ_F2_NBI .. nbi_02 =e= sum(m, phi('i2',m)*beta(m)) + t*normal('i2') - (obj_f('i2')-f_o('i2'));


*-------------------------------
*  Constraints
* ------------------------------
INEQ_NBI_01 .. nbi_01 =e=0;
INEQ_NBI_02 .. nbi_02 =e=0;





model NBI_cs5 /OBJ,OBJ_G2,OBJ_F1,OBJ_F2,OBJ_F1_NBI,OBJ_F2_NBI,INEQ_NBI_01,INEQ_NBI_02,EQ_V,AX_EQ1,AX_EQ2,AX_EQ3/;

*-------------------------------
*  Optimisation paramters
* ------------------------------
OPTION minlp = DICOPT;
OPTION OPTCA=1e-7;
OPTION OPTCR=1e-7 ;

*NBI_cs2.optfile=1;




solve NBI_cs5 using minlp maximising objective;
put obj_f.L('i1'):15:10,':',obj_f.L('i2'):15:10,':',yo.L:15:10,':',y.L('j1'):15:10,':',y.L('j2'):15:10,':',y.L('j3'):15:10,':',y.L('j4'):15:10/
put y.L('j5'):15:10,':',y.L('j6'):15:10,':',y.L('j7'):15:10,':',y.L('j8'):15:10,':',y.L('j9'):15:10,':',y.L('j10'):15:10,':',t.L:15:10/
put nbi_01.L:15:10,':',nbi_02.L:15:10,':',INEQ_NBI_01.M:15:10,':',INEQ_NBI_02.M:15:10,':',NBI_cs5.etSolve:15:10,':',NBI_cs5.etSolver:15:10,':',NBI_cs5.modelStat:15:10/;

put '---------SIMULATION END------------------'/;

