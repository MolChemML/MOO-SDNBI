$title NBI test3 (NBI_CS3.gms)

$gdxin inputcs3_nbi.gdx

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

* ----------------------------------------------
*  Declare paramters
* ----------------------------------------------
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
f_utopia(i) /i1=0.0416641269000000, i2=0.0416641269000000/
f_nadir(i)  /i1=1.03844983740000, i2=1.03844983740000 /
w(i) /i1=1,i2=1/
beta(m)
normal(i)
f_o(i)
f_k(i,m)
phi(i,m)
x_i(j)

;
$load    beta normal f_o f_k phi x_i

* ----------------------------------------------
*  Store results in the txt file
* ----------------------------------------------
file res /NBI_CS3_result.txt/;
put res;
put '---------SIMULATION START------------------'/;

Variables
t            objective function to maximize
objective    nbi objective function
objective_w  scalarized objective fucntion
obj_f(i)     each objective
x(j)         variables
eq_const1    equality constraint 1
eq_const2    equality constraint 2
ineq_const   inequality constratint

nbi_01
nbi_02
;



* ------------------------------------------------
*  Provide upper and lower bounds of variables
* ------------------------------------------------
* lower and upper bound of x
x.LO('j1')=0;
x.UP('j1')=3.14159265358;
x.LO('j2')=0;
x.UP('j2')=3.14159265358;


* lower and upper bound of t
t.LO=-10000;
t.UP=10000;

* ------------------------------------------------
*  Provide initial guesses
* ------------------------------------------------
x.L('j1')=x_i('j1');
x.L('j2')=x_i('j2');


*f.LO('i1')=f_k('i1','m1');
*f.UP('i1')=f_k('i2','m1');
*f.LO('i2')=f_k('i2','m2');
*f.UP('i2')=f_k('i1','m2');

* ------------------------------------------------
*  Declaire equations
* ------------------------------------------------
Equations
OBJ         objective nbi
OBJ_w       objectives
OBJ_F1      first objective
OBJ_F2      second objective
OBJ_F1_NBI  first objective for NBI
OBJ_F2_NBI  second objective for NBI

CONST_01    constraints equation 1
CONST_02    constraints equation 2


INEQ_01     ineqaulity contraint 1
INEQ_02     ineqaulity contraint 2
INEQ_NBI_01 inequality constraint 1 for NBI
INEQ_NBI_02 inequality constraint 2 for NBI
;

OBJ         .. objective =e= t;
OBJ_w       .. objective_w =e= sum(i,w(i)*obj_f(i));
OBJ_F1      .. obj_f('i1') =e= (x('j1')-f_utopia('i1'))/(f_nadir('i1')-f_utopia('i1'));
OBJ_F2      .. obj_f('i2') =e= (x('j2')-f_utopia('i2'))/(f_nadir('i2')-f_utopia('i2'));

OBJ_F1_NBI  .. nbi_01 =e= sum(m, phi('i1',m)*beta(m)) + t*normal('i1') - (obj_f('i1')-f_o('i1'));
OBJ_F2_NBI  .. nbi_02 =e= sum(m, phi('i2',m)*beta(m)) + t*normal('i2') - (obj_f('i2')-f_o('i2'));

CONST_01    .. eq_const1 =e= power((x('j1')-0.5),2) + power((x('j2')-0.5),2)-0.5 ;
CONST_02    .. eq_const2 =e= -power(x('j1'),2)-power(x('j2'),2)+1+0.1*cos( 16*arctan(x('j1')/x('j2')) );



*-------------------------------
*  constraints
* ------------------------------
INEQ_01     .. eq_const1 =l= 0;
INEQ_02     .. eq_const2 =l= 0;
INEQ_NBI_01 .. nbi_01 =g= 0;
INEQ_NBI_02 .. nbi_02 =g= 0;



model NBI_cs3 /OBJ,OBJ_F1,OBJ_F2,OBJ_F1_NBI,OBJ_F2_NBI,INEQ_NBI_01,INEQ_NBI_02,CONST_01,CONST_02,INEQ_01,INEQ_02/;

* Optimality conditions:
*CPLEX use dual SIMPLEX
*BARON   Branch and reduced for MINLP NLP
*CONOPT  GRG NLP  /  SNOPT: SQP NLP
*OPTION LP= CPLEX;
OPTION NLP=SNOPT;
OPTION OPTCA=1e-6;
OPTION OPTCR=1e-6 ;

*NBI_cs3.optfile=1;



solve NBI_cs3 using nlp maximising objective;
put obj_f.L('i1'):15:10,':',obj_f.L('i2'):15:10,':',x.L('j1'):15:10,':',x.L('j2'):15:10,':', objective.l:15:10,':', nbi_01.L:15:10,':',nbi_02.L:15:10,':',eq_const1.L:15:10,':',eq_const2.L:15:10,':',INEQ_NBI_01.M:15:10,':',INEQ_NBI_02.M:15:10,':',INEQ_01.m:15:10,':',INEQ_02.m:15:10,':',NBI_cs3.etSolve:15:10,':',NBI_cs3.etSolver:15:10,':',NBI_cs3.modelStat:15:10/;
*);
obj_f.L('i1')=0;
obj_f.L('i2')=0;

put '---------SIMULATION END------------------'/;

DISPLAY  OBJ_F1_NBI.M, OBJ_F2_NBI.m, INEQ_NBI_01.m, INEQ_NBI_02.m, INEQ_01.m, INEQ_02.m;
*TRACE=%gams.sysdir%auditlog\log.txt
