$Title Fluxo de Potencia Otimo (AC-OPF) - Minimizacao de Custos

Set
   i        'barras'            / 1*4   /
   slack(i) barra de referencia / 1     /
   GB(i)    'barras de geracao' / 1,2,3 /;

Scalar
   Sbase / 100 /;

Alias (i,j);

Table GenD(i,*) 'Caracteristicas das unidades geradoras (Exercicio 1)'
       pmax  pmin    a      b    c     Qmax  Qmin
   1   450   200     0.004  5.2  500   9999  -9999
   2   350   150     0.006  5.0  400   9999  -9999
   3   235   100     0.008  5.0  200   9999  -9999;

Table BD(i,*) 'Demanda de cada barra em MW e MVAr (Exercicio 2)'
       Pd   Qd
   1   0    0
   2   0    0
   3   220  136.34
   4   280  173.52;

Table LN(i,j,*) 'Caracteristicas das linhas de transmissao'
          r         x        b       limit
   1 .4   0.00744   0.0372   0.0775  9999
   1 .3   0.01008   0.0504   0.1025  9999
   2 .3   0.00744   0.0372   0.0775  9999
   2 .4   0.01272   0.0636   0.1275  9999;

* --- PRE-PROCESSAMENTO DOS PARAMETROS DA REDE ---
LN(i,j,'x')$(LN(i,j,'x')=0) = LN(j,i,'x');
LN(i,j,'r')$(LN(i,j,'r')=0) = LN(j,i,'r');
LN(i,j,'b')$(LN(i,j,'b')=0) = LN(j,i,'b');
LN(i,j,'Limit')$(LN(i,j,'Limit')=0) = LN(j,i,'Limit');

LN(i,j,'z')$LN(i,j,'Limit') = sqrt(sqr(LN(i,j,'x')) + sqr(LN(i,j,'r')));
LN(j,i,'z')$(LN(i,j,'z')=0) = LN(i,j,'z');

LN(i,j,'th')$(LN(i,j,'Limit') and LN(i,j,'x') and LN(i,j,'r'))   = arctan(LN(i,j,'x')/LN(i,j,'r'));
LN(i,j,'th')$(LN(i,j,'Limit') and LN(i,j,'x') and LN(i,j,'r')=0) = pi/2;
LN(i,j,'th')$(LN(i,j,'Limit') and LN(i,j,'r') and LN(i,j,'x')=0) = 0;
LN(j,i,'th')$LN(i,j,'Limit') = LN(i,j,'th');

Parameter cx(i,j);
cx(i,j)$(LN(i,j,'limit') and LN(j,i,'limit')) = 1;
cx(i,j)$(cx(j,i)) = 1;

* --- VARIAVEIS E EQUACOES ---
Variable CT, Pij(i,j), Qij(i,j), Pg(i), Qg(i), Va(i), V(i);
Equation eq1, eq2, eq_P, eq_Q, eq_cost, eq6;

* Fluxo de Potencia Ativa nas linhas
eq1(i,j)$cx(i,j)..
   Pij(i,j) =e= (sqr(V(i))*cos(LN(j,i,'th')) - V(i)*V(j)*cos(Va(i) - Va(j) + LN(j,i,'th'))) / LN(j,i,'z');

* Fluxo de Potencia Reativa nas linhas
eq2(i,j)$cx(i,j)..
   Qij(i,j) =e= (sqr(V(i))*sin(LN(j,i,'th')) - V(i)*V(j)*sin(Va(i) - Va(j) + LN(j,i,'th'))) / LN(j,i,'z') - LN(j,i,'b')*sqr(V(i))/2;

* Balanco Nodal de Potencia Ativa (Geracao - Carga = Fluxo injetado)
eq_P(i)..
    Pg(i)$GB(i) - BD(i,'Pd')/Sbase =e= sum(j$cx(j,i), Pij(i,j));
               
* Balanco Nodal de Potencia Reativa
eq_Q(i)..
    Qg(i)$GB(i) - BD(i,'Qd')/Sbase =e= sum(j$cx(j,i), Qij(i,j));
 
* Funcao Custo (Convertendo Pg de p.u. para MW para o calculo financeiro)
eq_cost..
   CT =e= sum(i$GB(i), GenD(i,'c') + GenD(i,'b')*(Pg(i)*Sbase) + GenD(i,'a')*sqr(Pg(i)*Sbase));

* Limite termico das linhas
eq6(i,j)$cx(i,j).. sqr(Pij(i,j)) + sqr(Qij(i,j)) =l= sqr(LN(i,j,'Limit')/Sbase);

Model loadflow /all/;

* --- LIMITES DAS VARIAVEIS ---
Pg.lo(i)$GB(i) = GenD(i,'Pmin')/Sbase;
Pg.up(i)$GB(i) = GenD(i,'Pmax')/Sbase;
Qg.lo(i)$GB(i) = GenD(i,'Qmin')/Sbase;
Qg.up(i)$GB(i) = GenD(i,'Qmax')/Sbase;

Va.up(i)     = pi/2;
Va.lo(i)     =-pi/2;
Va.l(i)      = 0;
Va.fx(slack) = 0;

V.lo(i)  = 0.90;
V.up(i)  = 1.10;
V.l(i)   = 1;

* Fixando a tensao da barra slack para o referencial de fase
V.fx('1') = 1.0;

* Obs: Deixei as tensoes V2 e V3 livres entre 0.9 e 1.1 para que o 
* solver otimize tambem a injecao de reativos minimizando perdas.

* --- SOLUCAO ---
solve loadflow minimizing CT using nlp;

* --- RELATORIO DE EXPORTACAO ---
Parameter report(i,*), PerdasTotais;
report(i,'V (pu)')    = V.l(i);
report(i,'Ang (grau)')= Va.l(i) * 180 / pi;
report(i,'Pg (MW)')   = Pg.l(i) * Sbase;
report(i,'Qg (MVAr)') = Qg.l(i) * Sbase;

PerdasTotais = sum(i, Pg.l(i)*Sbase) - sum(i, BD(i,'Pd'));

display CT.l, PerdasTotais;
display report;
Parameter IC(i) 'Custo Incremental da maquina ($/MWh)';
Parameter LMP(i) 'Preco Marginal Nodal / Lambda Local ($/MWh)';

* Derivada da funcao custo: dF/dP = 2*a*P + b
IC(i)$GB(i) = 2 * GenD(i,'a') * (Pg.l(i)*Sbase) + GenD(i,'b');

* Extraindo o multiplicador de Lagrange (shadow price) da equacao de balanco de potencia
* Divide-se por Sbase porque a equacao foi montada em p.u.
LMP(i)$GB(i) = eq_P.m(i) / Sbase;
Parameter FatorPenalizacao(i) 'Fator de Penalizacao equivalente (referencia na Barra Slack)';

* Fator de Penalizacao = LMP da barra / LMP da barra de referencia (Slack)
FatorPenalizacao(i)$GB(i) = LMP(i) / LMP('1');

display FatorPenalizacao;
display IC, LMP;

