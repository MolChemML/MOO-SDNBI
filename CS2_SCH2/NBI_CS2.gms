$title NBI test2 (NBI_CS2.gms)

$gdxin inputcs2_nbi.gdx


* ----------------------------------------------
*  Declare sets
*  : i, j, m will be assgiend by external interface
* ----------------------------------------------
Sets
i        number of objectivse
j        number of variables
m        number of points
k        iteration number /k1*k21/
c        number of integer cuts
;
$load    i j m c

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

Parameters
beta(m)
f_utopia(i) utopian points in objective space /i1=7.2509860, i2=0.5/
f_nadir(i)  nadir points in objective space   /i1=18.50000091, i2=3.04535252/
w(i)        weights for weighted sum method   /i1=1, i2=1/
phi(i,m)    NBI parameters for choosing referece points
normal(i)   direction vector for NBI
f_o(i)      origin to shift the objective space
x_i(j)      initial guesses of x
*for integer cut
yv(c,j) store x values from previous solutions
;
$load    beta normal f_o phi yv x_i


* ----------------------------------------------
*  Store results in the txt file
* ----------------------------------------------
file res /NBI_CS2_result.txt/;
put res;
put '---------SIMULATION START------------------'/;

Variables
t            objective function to maximize
objective    nbi objective function
objective_w  scalarized objective fucntion
obj_f(i)     each objective
x(j)         variables
eq_const1    equality constraint
eq_const2    equality constraint
ineq_const   inequality constratins

nbi_01
nbi_02
;


* ------------------------------------------------
*  Provide upper and lower bounds of variables
* ------------------------------------------------
* lower and upper bound of x
x.LO('j1')=0;
x.UP('j1')=5;
x.LO('j2')=0;
x.UP('j2')=3;
t.LO=-10000;
t.UP=10000;

* initial gusses of x
x.L('j1')=x_i('j1');
x.L('j2')=x_i('j2');


* ------------------------------------------------
*  Declaire equations
* ------------------------------------------------
Equations
OBJ         objective function of nbi
OBJ_w       weighted sum objective fuction
OBJ_F1      first objective
OBJ_F2      second objective
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
OBJ_F1 .. obj_f('i1') =e= ((power((x('j1')+x('j2')-7.5),2)+power((x('j2')-x('j1')+3),2)/4)-f_utopia('i1'))/(f_nadir('i1')-f_utopia('i1'));
OBJ_F2 .. obj_f('i2') =e= (power((x('j1')-1),2)/4 + power((x('j2')-4),2)/2-f_utopia('i2'))/(f_nadir('i2')-f_utopia('i2'));

OBJ_F1_NBI .. nbi_01 =e= sum(m, phi('i1',m)*beta(m)) + t*normal('i1') - (obj_f('i1')-f_o('i1'));
OBJ_F2_NBI .. nbi_02 =e= sum(m, phi('i2',m)*beta(m)) + t*normal('i2') - (obj_f('i2')-f_o('i2'));

CONST_01 .. eq_const1 =e= power((x('j1')-2),3)/2 + x('j2') -2.5 ;
CONST_02 .. eq_const2 =e= x('j1')+x('j2')-8*power((x('j2')-x('j1')+0.65),2)-3.85;




*-------------------------------
*  constraints
* ------------------------------
INEQ_01 .. eq_const1 =l=0;
INEQ_02 .. eq_const2 =l=0;
INEQ_NBI_01 .. nbi_01 =g=0;
INEQ_NBI_02 .. nbi_02 =g=0;

model NBI_cs2 /OBJ,OBJ_F1,OBJ_F2,OBJ_F1_NBI,OBJ_F2_NBI,INEQ_NBI_01,INEQ_NBI_02,CONST_01,CONST_02,INEQ_01,INEQ_02/;

* --------------------------------
* Solver and Optimality conditions
* --------------------------------
* Optimality conditions:
*CPLEX use dual SIMPLEX
*BARON   Branch and reduced for MINLP NLP
*CONOPT  GRG NLP  /  SNOPT: SQP NLP
OPTION LP= CPLEX;
OPTION NLP=CONOPT;
OPTION OPTCA=1e-7;
OPTION OPTCR=1e-7 ;

*NBI_cs2.optfile=1;



solve NBI_cs2 using nlp maximising objective;
put obj_f.L('i1'):15:10,':',obj_f.L('i2'):15:10,':',x.L('j1'):15:10,':',x.L('j2'):15:10,':', objective.l:15:10,':', nbi_01.L:15:10,':',nbi_02.L:15:10,':',eq_const1.L:15:10,':',eq_const2.L:15:10,':',INEQ_NBI_01.M:15:10,':',INEQ_NBI_02.M:15:10,':',INEQ_01.m:15:10,':',INEQ_02.m:15:10,':',NBI_cs2.etSolve:15:10,':',NBI_cs2.etSolver:15:10,':',NBI_cs2.modelStat:15:10/;
obj_f.L('i1')=0;
obj_f.L('i2')=0;

********************************************
put '---------SIMULATION END------------------'/;

*DISPLAY  OBJ_F1_NBI.M, OBJ_F2_NBI.m, INEQ_NBI_01.m, INEQ_NBI_02.m, INEQ_01.m, INEQ_02.m;
*TRACE=%gams.sysdir%auditlog\log.txt
