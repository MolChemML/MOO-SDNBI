$title NBI test3 (NBI_CS3_WS.gms)

$gdxin inputcs3_ws.gdx

sets
i        number of objectivse
j        number of variables
m        number of points /m1*m2/
k        iteration number /k1*k21/
*c        number of integer cuts
;
$load    i j

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

Table phi(i,m)
    m1   m2
i1  0    2
i2  3    0
;

PARAMETERS
*beta(i) /i1=0.5, i2=0.5/
*w(i) /i1=0, i2=1/
*phi(i) /i1=2, i2=1/
f_utopia(i) /i1=0.0416641269000000, i2=0.0416641269000000/
f_nadir(i)  /i1= 1.03844983740000, i2=1.03844983740000 /
normal(i) /i1=-1, i2=-2/
f_o(i) /i1=0, i2=0/
beta(m)  /m1=0,m2=0/
w(i)
*/i1=0.25,    i2= 0.75/
phi(i,m)
x_i(j)
*normal(i)
*f_o(i)
*for integer cut
*yv(c,j) store x values from previous solutions
;
$load    w x_i

* store results in the txt file
file res /NBI_CS3_WS_result.txt/;
put res;
put '---------SIMULATION START------------------'/;

variables
t            maximize
objective    nbi objective function
objective_w  scalarized objective fucntion
obj_f(i)     each objective
x(j)         variables
*alpha(c)     for intercut
eq_const1    equality constraint
eq_const2    equality constraint
ineq_const   inequality constratins

nbi_01
nbi_02
;

*binary variable alpha;

* lower and upper bound of x

x.LO('j1')=1e-6;
x.UP('j1')=3.14159265358;
x.LO('j2')=1e-6;
x.UP('j2')=3.14159265358;
t.LO=-10000;
t.UP=10000;

* initial gusses of x
x.L(j)=x_i(j);


equations
OBJ       objective nbi
OBJ_w     objectives
OBJ_F1  first objective
OBJ_F2  second objective
OBJ_F1_NBI  first objective
OBJ_F2_NBI  second objective

CONST_01 constraints1
CONST_02 constraints1


INEQ_01
INEQ_02
INEQ_NBI_01  constraints2
INEQ_NBI_02  constraints2


;

OBJ   .. objective =e= t;
OBJ_w .. objective_w =e= sum(i,w(i)*obj_f(i));
OBJ_F1 .. obj_f('i1') =e= (x('j1')-f_utopia('i1'))/(f_nadir('i1')-f_utopia('i1'));
OBJ_F2 .. obj_f('i2') =e= (x('j2')-f_utopia('i2'))/(f_nadir('i2')-f_utopia('i2'));

OBJ_F1_NBI .. nbi_01 =e= sum(m, phi('i1',m)*beta(m)) + t*normal('i1') - (obj_f('i1')-f_o('i1'));
OBJ_F2_NBI .. nbi_02 =e= sum(m, phi('i2',m)*beta(m)) + t*normal('i2') - (obj_f('i2')-f_o('i2'));

CONST_01 .. eq_const1 =e= power((x('j1')-0.5),2) + power((x('j2')-0.5),2)-0.5 ;
CONST_02 .. eq_const2 =e= -power(x('j1'),2)-power(x('j2'),2)+1+0.1*cos( 16*arctan(x('j1')/x('j2')) );


*-------------------------------
*  constraints
* ------------------------------
INEQ_01 .. eq_const1 =l=0;
INEQ_02 .. eq_const2 =l=0;
INEQ_NBI_01 .. nbi_01 =g=0;
INEQ_NBI_02 .. nbi_02 =g=0;

*model NBI_cs2 /OBJ,OBJ_F1,OBJ_F2,OBJ_F1_NBI,OBJ_F2_NBI,INEQ_NBI_01,INEQ_NBI_02,CONST_01,CONST_02,EQ_01,EQ_02,IntCut1,IntCut2/;
model Weight_cs3  /OBJ_w,OBJ_F1,OBJ_F2,CONST_01,CONST_02,INEQ_01,INEQ_02/;

* Optimality conditions:
*CPLEX use dual SIMPLEX
*BARON    CONOPT or SNOPT
OPTION LP= CPLEX;
OPTION NLP=SNOPT;
OPTION OPTCA=1e-6;
OPTION OPTCR=1e-6 ;


solve Weight_cs3 using nlp minimising objective_w;
put obj_f.L('i1'):15:10,':',obj_f.L('i2'):15:10,':',x.L('j1'):15:10,':',x.L('j2'):15:10,':', objective_w.l:15:10,':',eq_const1.L:15:10,':',eq_const2.L:15:10,':',Weight_cs3.etSolve:15:10,':',Weight_cs3.etSolver:15:10/;
*);
obj_f.L('i1')=0;
obj_f.L('i2')=0;

********************************************

put '---------SIMULATION END------------------'/;
